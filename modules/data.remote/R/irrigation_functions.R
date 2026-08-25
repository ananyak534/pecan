# Irrigation helper functions for PEcAn.data.remote.
# Functions called directly by package scripts are exported with roxygen2.
# Internal helpers remain available within the package namespace.

# ============================================================
# CAL-ADAPT / PENMAN-MONTEITH HELPERS
# ============================================================

#' Convert specific humidity to actual vapor pressure
#'
#' @param q Specific humidity.
#' @param pressure_pa Air pressure in pascals.
#' @return Actual vapor pressure in kPa.
#' @importFrom rlang .data .env
#' @keywords internal
specific_humidity_to_ea = function(q, pressure_pa) {
  P_kpa = pressure_pa / 1000
  (q * P_kpa) / (0.622 + 0.378 * q)
}


#' Calculate net radiation
#'
#' @param sw_rad Downwelling shortwave radiation.
#' @param lw_rad Downwelling longwave radiation.
#' @param temp_c Air temperature in degrees Celsius.
#' @param albedo Surface albedo.
#' @param emissivity Surface emissivity.
#' @return Net radiation in MJ m-2 for the hourly timestep.
#' @keywords internal
calculate_net_radiation = function(sw_rad, lw_rad, temp_c, albedo = 0.23, emissivity = 0.98) {
  sigma = 5.67e-8
  Rns = (1 - albedo) * sw_rad
  Rnl = lw_rad - emissivity * sigma * (temp_c + 273.15)^4
  (Rns + Rnl) * 0.0036
}


#' Parse NetCDF time coordinates
#'
#' @param time_steps Numeric NetCDF time values.
#' @param time_units NetCDF time-unit string.
#' @param year_val Fallback year used when the time origin cannot be parsed.
#' @return A vector of POSIXct timestamps.
#' @keywords internal
parse_nc_time = function(time_steps, time_units, year_val) {
  origin_string = sub("^.*since\\s+", "", time_units)
  origin_time = lubridate::ymd_hms(origin_string, tz = "UTC", quiet = TRUE)
  
  if (is.na(origin_time)) {
    origin_time = lubridate::ymd_hms(sprintf("%d-01-01 00:00:00", year_val), tz = "UTC")
  }
  
  if (grepl("^seconds", time_units, ignore.case = TRUE)) {
    origin_time + lubridate::seconds(time_steps)
  } else if (grepl("^hours", time_units, ignore.case = TRUE)) {
    origin_time + lubridate::hours(time_steps)
  } else if (grepl("^days", time_units, ignore.case = TRUE)) {
    origin_time + lubridate::days(time_steps)
  } else {
    stop("Unknown NetCDF time units: ", time_units)
  }
}


#' Calculate daily reference evapotranspiration
#'
#' Aggregates hourly climate data to daily values and calculates FAO-56
#' Penman-Monteith reference evapotranspiration.
#'
#' @param df_hourly Data frame containing hourly climate variables.
#' @param wind_height_m Wind measurement height in meters.
#' @return A tibble containing daily climate summaries and ET0 in millimeters.
#' @keywords internal
calculate_daily_et0 = function(df_hourly, wind_height_m = 10) {
  df_daily = df_hourly |>
    dplyr::mutate(
      ea_kpa = specific_humidity_to_ea(.data$spec_hum, .data$pressure),
      Rn_MJ = calculate_net_radiation(.data$sw_rad, .data$lw_rad, .data$temp_c),
      date = as.Date(.data$datetime)
    ) |>
    dplyr::group_by(.data$date) |>
    dplyr::summarise(
      T_min = min(.data$temp_c, na.rm = TRUE),
      T_max = max(.data$temp_c, na.rm = TRUE),
      mean_temp_c = mean(.data$temp_c, na.rm = TRUE),
      R_n = sum(.data$Rn_MJ, na.rm = TRUE),
      e_a = mean(.data$ea_kpa, na.rm = TRUE),
      P_kpa = mean(.data$pressure / 1000, na.rm = TRUE),
      wind = mean(.data$wind_speed, na.rm = TRUE),
      precip_mm = sum(.data$precip_mm, na.rm = TRUE),
      .groups = "drop"
    )
  
  df_daily |>
    dplyr::rowwise() |>
    dplyr::mutate(
      ET0_mm = {
        gamma = 0.000665 * .data$P_kpa
        if (wind_height_m == 2) {
          FAO56::ETo_FPM(
            R_n = .data$R_n, G = 0, gamma = gamma, u_2 = .data$wind,
            e_a = .data$e_a, T_min = .data$T_min, T_max = .data$T_max
          )
        } else {
          FAO56::ETo_FPM(
            R_n = .data$R_n, G = 0, gamma = gamma, u_z = .data$wind,
            z = wind_height_m, e_a = .data$e_a, T_min = .data$T_min, T_max = .data$T_max
          )
        }
      }
    ) |>
    dplyr::ungroup()
}


#' Process one Cal-Adapt NetCDF file
#'
#' Reads hourly Cal-Adapt climate variables from a NetCDF file and returns
#' daily reference evapotranspiration and climate summaries.
#'
#' @param nc_path Path to a NetCDF climate file.
#' @param wind_height_m Wind measurement height in meters.
#' @return A tibble of daily climate values, or `NULL` when the file cannot be read.
#' @export
process_nc_file = function(nc_path, wind_height_m = 10) {
  nc = tryCatch(
    ncdf4::nc_open(nc_path),
    error = function(e) {
      warning("Could not open: ", nc_path)
      NULL
    }
  )
  
  if (is.null(nc)) return(NULL)
  on.exit(ncdf4::nc_close(nc))
  
  file_name = basename(nc_path)
  site_hash = basename(dirname(nc_path))
  parts = unlist(strsplit(gsub("\\.nc$", "", file_name), "\\."))
  
  if (length(parts) < 3) {
    warning("Unexpected filename: ", file_name)
    return(NULL)
  }
  
  model_name = parts[1]
  scenario_value = parts[2]
  year_val = as.integer(parts[3])
  
  time_steps = as.numeric(ncdf4::ncvar_get(nc, "time"))
  temp_k = as.numeric(ncdf4::ncvar_get(nc, "air_temperature"))
  wind_speed = as.numeric(ncdf4::ncvar_get(nc, "wind_speed"))
  precip = as.numeric(ncdf4::ncvar_get(nc, "precipitation_flux"))
  spec_hum = as.numeric(ncdf4::ncvar_get(nc, "specific_humidity"))
  lw_rad = as.numeric(ncdf4::ncvar_get(nc, "surface_downwelling_longwave_flux_in_air"))
  sw_rad = as.numeric(ncdf4::ncvar_get(nc, "surface_downwelling_shortwave_flux_in_air"))
  pressure = as.numeric(ncdf4::ncvar_get(nc, "air_pressure"))
  time_units = ncdf4::ncatt_get(nc, "time", "units")$value
  datetime = parse_nc_time(time_steps, time_units, year_val)
  
  df_hourly = tibble::tibble(
    datetime = datetime, temp_c = temp_k - 273.15, wind_speed = wind_speed,
    sw_rad = sw_rad, lw_rad = lw_rad, spec_hum = spec_hum,
    pressure = pressure, precip_mm = precip * 3600
  )
  
  calculate_daily_et0(df_hourly, wind_height_m) |>
    dplyr::mutate(
      site_hash = .env$site_hash,
      model = .env$model_name,
      scenario = .env$scenario_value,
      year = lubridate::year(.data$date),
      day_of_year = lubridate::yday(.data$date)
    ) |>
    dplyr::select(
      "date", "site_hash", "model", "scenario", "ET0_mm", "precip_mm",
      "mean_temp_c", "T_min", "T_max", dplyr::everything()
    )
}


# ============================================================
# SSURGO / AVAILABLE WATER CAPACITY HELPERS
# ============================================================

# Package-local cache so SSURGO tables are read only once per R session.
ssurgo_table_cache = new.env(parent = emptyenv())
ssurgo_table_cache$path = NULL
ssurgo_table_cache$component = NULL
ssurgo_table_cache$chorizon = NULL


#' Calculate effective available water capacity
#'
#' Clips soil horizons to the crop rooting depth and calculates the available
#' water capacity of the resulting soil profile.
#'
#' @param hzdept_r_cm Horizon top depths in centimeters.
#' @param hzdepb_r_cm Horizon bottom depths in centimeters.
#' @param awc_r Available water capacity in centimeters of water per centimeter of soil.
#' @param rooting_depth_cm Crop rooting depth in centimeters.
#' @return Effective available water capacity in millimeters.
#' @keywords internal
calc_effective_awc = function(hzdept_r_cm, hzdepb_r_cm, awc_r, rooting_depth_cm) {
  effective_top = pmin(hzdept_r_cm, rooting_depth_cm)
  effective_bottom = pmin(hzdepb_r_cm, rooting_depth_cm)
  thickness_cm = pmax(0, effective_bottom - effective_top)
  sum(awc_r * thickness_cm, na.rm = TRUE) * 10
}


#' Load SSURGO tables into the irrigation cache
#'
#' Reads the SSURGO component and horizon tables once and stores them in a
#' package-local cache for reuse across county irrigation calculations.
#'
#' @param ssurgo_gdb_path Path to the gSSURGO geodatabase.
#' @param ssurgo_weights Data frame containing parcel-to-SSURGO mapping weights.
#' @return Invisibly returns `NULL`.
#' @export
load_ssurgo_tables_once = function(ssurgo_gdb_path, ssurgo_weights) {
  current_path = normalizePath(ssurgo_gdb_path, winslash = "/", mustWork = TRUE)
  already_loaded = !is.null(ssurgo_table_cache$path) &&
    identical(ssurgo_table_cache$path, current_path)
  
  if (already_loaded) return(invisible(NULL))
  
  message("Loading SSURGO component/chorizon tables once...")
  
  component = sf::read_sf(ssurgo_gdb_path, layer = "component", as_tibble = TRUE)
  if (inherits(component, "sf")) component = sf::st_drop_geometry(component)
  
  component = component |>
    dplyr::select("mukey", "cokey", "comppct_r") |>
    dplyr::semi_join(ssurgo_weights |> dplyr::distinct(.data$mukey), by = "mukey")
  
  chorizon = sf::read_sf(ssurgo_gdb_path, layer = "chorizon", as_tibble = TRUE)
  if (inherits(chorizon, "sf")) chorizon = sf::st_drop_geometry(chorizon)
  
  chorizon = chorizon |>
    dplyr::select("cokey", "hzdept_r", "hzdepb_r", "awc_r") |>
    dplyr::semi_join(component |> dplyr::distinct(.data$cokey), by = "cokey")
  
  ssurgo_table_cache$path = current_path
  ssurgo_table_cache$component = component
  ssurgo_table_cache$chorizon = chorizon
  
  message("SSURGO tables loaded and cached.")
  invisible(NULL)
}


#' Compute parcel soil available water capacity
#'
#' @param crop_info Crop records containing parcel IDs and rooting depths.
#' @param ssurgo_weights Parcel-to-SSURGO mapping weights.
#' @param ssurgo_gdb_path Path to the gSSURGO geodatabase.
#' @return `crop_info` with calculated `whc_mm`.
#' @keywords internal
compute_soil_awc_projection = function(crop_info, ssurgo_weights, ssurgo_gdb_path) {
  load_ssurgo_tables_once(ssurgo_gdb_path, ssurgo_weights)
  
  parcel_ids = unique(crop_info[["parcel_id"]])
  weights = ssurgo_weights |> dplyr::filter(.data$parcel_id %in% parcel_ids)
  
  if (nrow(weights) == 0) {
    return(crop_info |> dplyr::mutate(whc_mm = NA_real_))
  }
  
  component = ssurgo_table_cache$component |>
    dplyr::semi_join(weights, by = "mukey")
  
  chorizon = ssurgo_table_cache$chorizon |>
    dplyr::semi_join(component, by = "cokey")
  
  crop_depths = crop_info |>
    dplyr::distinct(.data$parcel_id, .data$rooting_depth_m)
  
  awc = weights |>
    dplyr::inner_join(crop_depths, by = "parcel_id", relationship = "many-to-many") |>
    dplyr::left_join(component, by = "mukey", relationship = "many-to-many") |>
    dplyr::left_join(chorizon, by = "cokey", relationship = "many-to-many") |>
    dplyr::filter(
      !is.na(.data$awc_r), !is.na(.data$hzdept_r),
      !is.na(.data$hzdepb_r), !is.na(.data$rooting_depth_m)
    ) |>
    dplyr::mutate(rooting_depth_cm = .data$rooting_depth_m * 100) |>
    dplyr::summarise(
      whc_mm_cmp = calc_effective_awc(
        .data$hzdept_r, .data$hzdepb_r, .data$awc_r, .data$rooting_depth_cm
      ),
      .by = c("parcel_id", "rooting_depth_m", "mukey", "cokey", "area_m2", "weight", "comppct_r")
    ) |>
    dplyr::summarise(
      whc_mm_mu = sum(.data$whc_mm_cmp * .data$comppct_r / sum(.data$comppct_r)),
      .by = c("parcel_id", "rooting_depth_m", "mukey", "area_m2", "weight")
    ) |>
    dplyr::summarise(
      whc_mm = sum(.data$whc_mm_mu * .data$weight),
      .by = c("parcel_id", "rooting_depth_m")
    )
  
  duplicate_awc = awc |>
    dplyr::count(.data$parcel_id, .data$rooting_depth_m) |>
    dplyr::filter(.data$n > 1)
  
  if (nrow(duplicate_awc) > 0) {
    stop("SSURGO calculation produced duplicate parcel/root-depth AWC rows.")
  }
  
  crop_info |>
    dplyr::left_join(awc, by = c("parcel_id", "rooting_depth_m"))
}


#' Add cached soil available water capacity to crop records
#'
#' @param crop_info Crop records containing parcel IDs and rooting depths.
#' @param ssurgo_weights Parcel-to-SSURGO mapping weights.
#' @param ssurgo_gdb_path Path to the gSSURGO geodatabase.
#' @param cache_file Optional parquet cache file.
#' @return `crop_info` with `whc_mm`.
#' @keywords internal
add_soil_awc_projection = function(crop_info, ssurgo_weights, ssurgo_gdb_path, cache_file = NULL) {
  if (is.null(cache_file)) {
    return(compute_soil_awc_projection(crop_info, ssurgo_weights, ssurgo_gdb_path))
  }
  
  dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)
  needed = crop_info |> dplyr::distinct(.data$parcel_id, .data$rooting_depth_m)
  
  if (file.exists(cache_file)) {
    cache = arrow::read_parquet(cache_file) |>
      dplyr::select("parcel_id", "rooting_depth_m", "whc_mm")
  } else {
    cache = tibble::tibble(
      parcel_id = needed$parcel_id[0],
      rooting_depth_m = needed$rooting_depth_m[0],
      whc_mm = numeric(0)
    )
  }
  
  missing = needed |>
    dplyr::anti_join(cache, by = c("parcel_id", "rooting_depth_m"))
  
  if (nrow(missing) > 0) {
    message("  Calculating ", nrow(missing), " new parcel/root-depth AWC combinations.")
    
    new_awc = compute_soil_awc_projection(missing, ssurgo_weights, ssurgo_gdb_path) |>
      dplyr::select("parcel_id", "rooting_depth_m", "whc_mm")
    
    cache = dplyr::bind_rows(cache, new_awc) |>
      dplyr::distinct(.data$parcel_id, .data$rooting_depth_m, .keep_all = TRUE)
    
    arrow::write_parquet(cache, cache_file)
  }
  
  crop_info |>
    dplyr::left_join(cache, by = c("parcel_id", "rooting_depth_m"))
}


# ============================================================
# IRRIGATION WATER-BALANCE HELPERS
# ============================================================

#' Summarize an irrigation event vector
#'
#' @param irr Daily irrigation amounts.
#' @param dates Dates corresponding to `irr`.
#' @param projection_row_id Prediction-row identifier.
#' @param parcel_id Parcel identifier.
#' @param projection_year Projection year.
#' @param crop_name Crop name.
#' @param keep_event_details Whether to retain event dates and amounts.
#' @return A one-row irrigation summary tibble.
#' @keywords internal
summarize_irrigation_vector = function(irr, dates, projection_row_id, parcel_id,
                                       projection_year, crop_name, keep_event_details = TRUE) {
  event_idx = which(!is.na(irr) & irr > 0)
  
  if (length(event_idx) > 0) {
    first_date = dates[min(event_idx)]
    last_date = dates[max(event_idx)]
    max_event = max(irr[event_idx], na.rm = TRUE)
  } else {
    first_date = as.Date(NA)
    last_date = as.Date(NA)
    max_event = 0
  }
  
  if (keep_event_details && length(event_idx) > 0) {
    event_dates = paste(format(dates[event_idx], "%Y-%m-%d"), collapse = ";")
    event_amounts = paste(round(irr[event_idx], 3), collapse = ";")
  } else {
    event_dates = NA_character_
    event_amounts = NA_character_
  }
  
  tibble::tibble(
    projection_row_id = projection_row_id, parcel_id = parcel_id,
    projection_year = projection_year, irrigation_total_mm = sum(irr, na.rm = TRUE),
    irrigation_n_events = length(event_idx), irrigation_first_date = first_date,
    irrigation_last_date = last_date, irrigation_max_event_mm = max_event,
    irrigation_method = ifelse(crop_name == "Rice", "flood", "canopy"),
    irrigation_event_dates = event_dates, irrigation_event_amounts_mm = event_amounts
  )
}


#' Run a batched non-rice irrigation water balance
#'
#' @param g Crop records sharing a phenology schedule.
#' @param dates Irrigation-period dates.
#' @param etc_mm_day Daily crop evapotranspiration.
#' @param precip_mm_day Daily precipitation.
#' @param irrigation_max Maximum irrigation applied in one event, in millimeters.
#' @param keep_event_details Whether to retain event dates and amounts.
#' @return A tibble containing irrigation summaries for all parcels in `g`.
#' @keywords internal
run_nonrice_water_balance_batch = function(g, dates, etc_mm_day, precip_mm_day,
                                           irrigation_max = 150, keep_event_details = TRUE) {
  n_parcels = nrow(g)
  n_days = length(dates)
  
  if (n_parcels == 0) return(tibble::tibble())
  if (length(etc_mm_day) != n_days || length(precip_mm_day) != n_days) {
    stop("ET/precipitation length does not match dates.")
  }
  if (anyNA(etc_mm_day) || anyNA(precip_mm_day)) {
    stop("Missing ET or precipitation in irrigation period.")
  }
  
  whc = as.numeric(g$whc_mm)
  whc_min_frac = as.numeric(g$whc_min_frac)
  
  if (anyNA(whc) || anyNA(whc_min_frac)) {
    stop("Missing WHC or MAD entered batch water balance.")
  }
  
  w_min = whc * whc_min_frac
  W_prev = whc
  irrigation_total = numeric(n_parcels)
  irrigation_n_events = integer(n_parcels)
  first_event_idx = rep(NA_integer_, n_parcels)
  last_event_idx = rep(NA_integer_, n_parcels)
  max_event = numeric(n_parcels)
  
  if (keep_event_details) {
    irrigation_history = matrix(0, nrow = n_days, ncol = n_parcels)
  }
  
  for (t in seq_len(n_days)) {
    W0 = W_prev + precip_mm_day[t] - etc_mm_day[t]
    irrigate = W0 < w_min
    irr_t = numeric(n_parcels)
    
    if (any(irrigate)) {
      irr_t[irrigate] = pmin(whc[irrigate] - W0[irrigate], irrigation_max)
      W0[irrigate] = W0[irrigate] + irr_t[irrigate]
    }
    
    above_capacity = W0 > whc
    W_t = pmax(W0, w_min)
    W_t[above_capacity] = whc[above_capacity]
    W_prev = W_t
    
    irrigation_total = irrigation_total + irr_t
    event = irr_t > 0
    irrigation_n_events[event] = irrigation_n_events[event] + 1L
    first_event_idx[event & is.na(first_event_idx)] = t
    last_event_idx[event] = t
    max_event = pmax(max_event, irr_t)
    
    if (keep_event_details) irrigation_history[t, ] = irr_t
  }
  
  first_date = as.Date(rep(NA_character_, n_parcels))
  last_date = as.Date(rep(NA_character_, n_parcels))
  has_first = !is.na(first_event_idx)
  has_last = !is.na(last_event_idx)
  first_date[has_first] = dates[first_event_idx[has_first]]
  last_date[has_last] = dates[last_event_idx[has_last]]
  
  if (keep_event_details) {
    event_dates = vapply(seq_len(n_parcels), function(i) {
      idx = which(irrigation_history[, i] > 0)
      if (length(idx) == 0) return(NA_character_)
      paste(format(dates[idx], "%Y-%m-%d"), collapse = ";")
    }, character(1))
    
    event_amounts = vapply(seq_len(n_parcels), function(i) {
      idx = which(irrigation_history[, i] > 0)
      if (length(idx) == 0) return(NA_character_)
      paste(round(irrigation_history[idx, i], 3), collapse = ";")
    }, character(1))
  } else {
    event_dates = rep(NA_character_, n_parcels)
    event_amounts = rep(NA_character_, n_parcels)
  }
  
  tibble::tibble(
    projection_row_id = g$projection_row_id, parcel_id = g$parcel_id,
    projection_year = g$projection_year, irrigation_total_mm = irrigation_total,
    irrigation_n_events = irrigation_n_events, irrigation_first_date = first_date,
    irrigation_last_date = last_date, irrigation_max_event_mm = max_event,
    irrigation_method = "canopy", irrigation_event_dates = event_dates,
    irrigation_event_amounts_mm = event_amounts
  )
}


#' Run irrigation calculations for modeled crop records
#'
#' @param modeled_crops Crop records ready for irrigation modeling.
#' @param climate Daily county climate data.
#' @param rice_flood_target_mm Target rice flood depth in millimeters.
#' @param rice_flood_min_mm Minimum rice flood depth in millimeters.
#' @param rice_flood_max_mm Maximum rice flood depth in millimeters.
#' @param rice_seepage_mm_day Daily rice seepage in millimeters.
#' @param irrigation_max_mm Maximum non-rice irrigation event in millimeters.
#' @param keep_event_details Whether to retain event dates and amounts.
#' @return A tibble containing irrigation summaries.
#' @keywords internal
run_irrigation_fast = function(modeled_crops, climate, rice_flood_target_mm, rice_flood_min_mm,
                               rice_flood_max_mm, rice_seepage_mm_day, irrigation_max_mm = 150,
                               keep_event_details = TRUE) {
  if (nrow(modeled_crops) == 0) return(tibble::tibble())
  
  schedule_groups = modeled_crops |>
    dplyr::group_by(.data$crop_name, .data$planting_date, .data$peak_date, .data$harvest_date) |>
    dplyr::group_split(.keep = TRUE)
  
  message("  Unique crop/phenology schedules: ", length(schedule_groups))
  schedule_results = vector("list", length(schedule_groups))
  
  for (s in seq_along(schedule_groups)) {
    g = schedule_groups[[s]]
    crop = g$crop_name[[1]]
    plant = g$planting_date[[1]]
    peak = g$peak_date[[1]]
    harvest = g$harvest_date[[1]]
    dates = seq.Date(plant, harvest, by = "day")
    
    climate_idx = match(dates, climate$date)
    
    if (anyNA(climate_idx)) {
      missing_dates = dates[is.na(climate_idx)]
      stop(
        "Missing climate for ", crop, " from ", plant, " to ", harvest,
        ". First missing date: ", missing_dates[[1]]
      )
    }
    
    clim = climate[climate_idx, , drop = FALSE]
    days_to_peak = max(1, as.numeric(peak - plant))
    days_after_peak = max(1, as.numeric(harvest - peak))
    
    canopy_cover = ifelse(
      dates <= peak,
      0.15 + 0.85 * as.numeric(dates - plant) / days_to_peak,
      1 - 0.85 * as.numeric(dates - peak) / days_after_peak
    )
    canopy_cover = pmin(1, pmax(0, canopy_cover))
    
    etc_mm_day = PEcAn.data.land::eto_to_etc_bism(
      eto = clim$ET0_mm, crop_name = crop, canopy_cover = canopy_cover
    )
    precip_mm_day = clim$precip_mm
    
    if (crop == "Rice") {
      wb = PEcAn.data.land::calc_water_balance_rice(
        et = etc_mm_day, precip = precip_mm_day,
        flood_target = rice_flood_target_mm, flood_min = rice_flood_min_mm,
        flood_max = rice_flood_max_mm, seepage = rice_seepage_mm_day
      )
      
      one_group = vector("list", nrow(g))
      for (i in seq_len(nrow(g))) {
        one_group[[i]] = summarize_irrigation_vector(
          wb$irr, dates, g$projection_row_id[[i]], g$parcel_id[[i]],
          g$projection_year[[i]], crop, keep_event_details
        )
      }
      
      schedule_results[[s]] = dplyr::bind_rows(one_group)
      next
    }
    
    schedule_results[[s]] = run_nonrice_water_balance_batch(
      g = g, dates = dates, etc_mm_day = etc_mm_day, precip_mm_day = precip_mm_day,
      irrigation_max = irrigation_max_mm, keep_event_details = keep_event_details
    )
    
    if (s %% 25 == 0) {
      message("    Finished schedule ", s, " / ", length(schedule_groups))
    }
  }
  
  dplyr::bind_rows(schedule_results)
}


# ============================================================
# RUN IRRIGATION FOR ONE COUNTY FILE
# ============================================================

#' Run irrigation predictions for one county file
#'
#' Adds irrigation predictions to one county crop-projection file using
#' Cal-Adapt climate data, BISM crop parameters, and SSURGO soil properties.
#'
#' @param crop_file Path to a county crop-prediction CSV.
#' @param target_scenario Scenario label used for logging.
#' @param output_dir Directory for the completed prediction CSV.
#' @param county_daily_climate County-level daily climate data.
#' @param bism_crop_unique LandIQ-to-BISM crop lookup.
#' @param crop_whc_sub Crop water-holding-capacity lookup.
#' @param ssurgo_weights Parcel-to-SSURGO mapping weights.
#' @param ssurgo_gdb_path Path to the gSSURGO geodatabase.
#' @param awc_cache_dir Directory for cached parcel AWC calculations.
#' @param gcm_use Climate-model identifier written to the output.
#' @param ssp_use Climate-scenario identifier written to the output.
#' @param rice_flood_target_mm Target rice flood depth in millimeters.
#' @param rice_flood_min_mm Minimum rice flood depth in millimeters.
#' @param rice_flood_max_mm Maximum rice flood depth in millimeters.
#' @param rice_seepage_mm_day Daily rice seepage in millimeters.
#' @param irrigation_max_mm Maximum non-rice irrigation event in millimeters.
#' @param keep_event_details Whether to retain event-level irrigation details.
#' @return Invisibly returns the completed county prediction data frame.
#' @export
run_irrigation_file = function(crop_file, target_scenario, output_dir, county_daily_climate,
                               bism_crop_unique, crop_whc_sub, ssurgo_weights, ssurgo_gdb_path,
                               awc_cache_dir, gcm_use, ssp_use, rice_flood_target_mm,
                               rice_flood_min_mm, rice_flood_max_mm, rice_seepage_mm_day,
                               irrigation_max_mm = 150, keep_event_details = TRUE) {
  message("Processing: ", basename(crop_file))
  message("Scenario:   ", target_scenario)
  
  crops = read.csv(crop_file, stringsAsFactors = FALSE) |>
    dplyr::mutate(
      projection_row_id = dplyr::row_number(),
      planting_date = as.Date(.data$planting_date),
      peak_date = as.Date(.data$peak_date),
      harvest_date = as.Date(.data$harvest_date)
    )
  
  original_n_rows = nrow(crops)
  
  if (anyNA(crops$county) || dplyr::n_distinct(crops$county) != 1) {
    stop("Expected exactly one county in ", basename(crop_file))
  }
  
  county_name = unique(crops$county)
  message("County:     ", county_name)
  message("Rows:       ", format(original_n_rows, big.mark = ","))
  
  climate = county_daily_climate |>
    dplyr::filter(.data$County == .env$county_name) |>
    dplyr::select("date", "ET0_mm", "precip_mm") |>
    dplyr::arrange(.data$date)
  
  if (nrow(climate) == 0) {
    stop("No Cal-Adapt climate for ", county_name)
  }
  
  duplicate_climate = climate |>
    dplyr::count(.data$date) |>
    dplyr::filter(.data$n > 1)
  
  if (nrow(duplicate_climate) > 0) {
    stop("Climate contains duplicate dates for ", county_name, ". Need exactly one county-level climate row per date.")
  }
  
  climate_min_date = min(climate$date)
  climate_max_date = max(climate$date)
  message("Climate:    ", climate_min_date, " through ", climate_max_date)
  
  crop_info_all = crops |>
    dplyr::mutate(
      projection_year = .data$year,
      CLASS_lookup = as.character(.data$CLASS),
      SUBCLASS_lookup = suppressWarnings(as.integer(.data$SUBCLASS))
    ) |>
    dplyr::left_join(
      bism_crop_unique,
      by = c("CLASS_lookup" = "landiq_class", "SUBCLASS_lookup" = "landiq_subclass")
    ) |>
    dplyr::left_join(crop_whc_sub, by = "crop_name")
  
  if (nrow(crop_info_all) != original_n_rows) {
    stop("Crop lookup changed row count for ", basename(crop_file), ". Check duplicate BIS/WHC lookup keys.")
  }
  
  status_table = crop_info_all |>
    dplyr::transmute(
      projection_row_id = .data$projection_row_id,
      crop_name = .data$crop_name,
      irrigation_mapping_source = .data$irrigation_mapping_source,
      irrigation_status = dplyr::case_when(
        .data$CLASS_lookup %in% c("X", "I") ~ "not_applicable_non_crop",
        is.na(.data$crop_name) ~ "not_modeled_no_BIS_crop",
        .data$crop_name != "Rice" & is.na(.data$rooting_depth_m) ~ "not_modeled_missing_rooting_depth",
        .data$crop_name != "Rice" & is.na(.data$whc_min_frac) ~ "not_modeled_missing_MAD",
        is.na(.data$planting_date) | is.na(.data$peak_date) | is.na(.data$harvest_date) ~ "not_modeled_missing_dates",
        .data$peak_date < .data$planting_date | .data$harvest_date < .data$peak_date ~ "not_modeled_invalid_dates",
        .data$planting_date < .env$climate_min_date | .data$harvest_date > .env$climate_max_date ~
          "not_modeled_climate_out_of_range",
        TRUE ~ "ready"
      )
    )
  
  crop_info = crop_info_all |>
    dplyr::semi_join(
      status_table |> dplyr::filter(.data$irrigation_status == "ready"),
      by = "projection_row_id"
    )
  
  message(
    "Rows ready before soil: ", format(nrow(crop_info), big.mark = ","),
    " / ", format(original_n_rows, big.mark = ",")
  )
  
  non_rice = crop_info |> dplyr::filter(.data$crop_name != "Rice")
  rice = crop_info |>
    dplyr::filter(.data$crop_name == "Rice") |>
    dplyr::mutate(whc_mm = NA_real_)
  
  if (nrow(non_rice) > 0) {
    awc_cache_file = file.path(
      awc_cache_dir,
      paste0(gsub("[^A-Za-z0-9]+", "_", county_name), "_awc.parquet")
    )
    
    non_rice = add_soil_awc_projection(
      crop_info = non_rice, ssurgo_weights = ssurgo_weights,
      ssurgo_gdb_path = ssurgo_gdb_path, cache_file = awc_cache_file
    )
  } else {
    non_rice$whc_mm = numeric(0)
  }
  
  crops_with_soil = dplyr::bind_rows(non_rice, rice)
  
  if (nrow(crops_with_soil) != nrow(crop_info)) {
    stop("Soil join changed row count for ", basename(crop_file))
  }
  
  soil_status = crops_with_soil |>
    dplyr::select("projection_row_id", "whc_mm")
  
  status_table = status_table |>
    dplyr::left_join(soil_status, by = "projection_row_id") |>
    dplyr::mutate(
      irrigation_status = dplyr::case_when(
        .data$irrigation_status == "ready" & .data$crop_name != "Rice" & is.na(.data$whc_mm) ~
          "not_modeled_missing_SSURGO",
        .data$irrigation_status == "ready" ~ "modeled",
        TRUE ~ .data$irrigation_status
      )
    ) |>
    dplyr::select(-dplyr::all_of(c("crop_name", "whc_mm")))
  
  modeled_crops = crops_with_soil |>
    dplyr::filter(.data$crop_name == "Rice" | !is.na(.data$whc_mm))
  
  message("Rows modeled after soil: ", format(nrow(modeled_crops), big.mark = ","))
  
  if (nrow(modeled_crops) == 0) {
    irrigation_summary = tibble::tibble()
  } else {
    irrigation_summary = run_irrigation_fast(
      modeled_crops = modeled_crops, climate = climate,
      rice_flood_target_mm = rice_flood_target_mm, rice_flood_min_mm = rice_flood_min_mm,
      rice_flood_max_mm = rice_flood_max_mm, rice_seepage_mm_day = rice_seepage_mm_day,
      irrigation_max_mm = irrigation_max_mm, keep_event_details = keep_event_details
    )
  }
  
  final_predictions = crops |>
    dplyr::left_join(status_table, by = "projection_row_id")
  
  if (nrow(irrigation_summary) > 0) {
    final_predictions = final_predictions |>
      dplyr::left_join(
        irrigation_summary |>
          dplyr::select(-dplyr::all_of(c("parcel_id", "projection_year"))),
        by = "projection_row_id"
      )
  } else {
    final_predictions = final_predictions |>
      dplyr::mutate(
        irrigation_total_mm = NA_real_, irrigation_n_events = NA_integer_,
        irrigation_first_date = as.Date(NA), irrigation_last_date = as.Date(NA),
        irrigation_max_event_mm = NA_real_, irrigation_method = NA_character_,
        irrigation_event_dates = NA_character_, irrigation_event_amounts_mm = NA_character_
      )
  }
  
  final_predictions = final_predictions |>
    dplyr::mutate(irrigation_GCM = gcm_use, irrigation_SSP = ssp_use) |>
    dplyr::mutate(
      irrigation_total_mm = dplyr::if_else(
        .data$irrigation_status == "not_applicable_non_crop", 0, .data$irrigation_total_mm
      ),
      irrigation_n_events = dplyr::if_else(
        .data$irrigation_status == "not_applicable_non_crop", 0L, .data$irrigation_n_events
      )
    ) |>
    dplyr::select(
      -dplyr::any_of(c("irrigation_event_dates", "irrigation_event_amounts_mm"))
    )
  
  if (nrow(final_predictions) != original_n_rows) {
    stop("FINAL ROW COUNT ERROR for ", basename(crop_file))
  }
  if (!identical(final_predictions$parcel_id, crops$parcel_id)) {
    stop("parcel_id order changed for ", basename(crop_file))
  }
  if (!identical(final_predictions$year, crops$year)) {
    stop("year order changed for ", basename(crop_file))
  }
  
  final_predictions = final_predictions |>
    dplyr::select(-dplyr::all_of("projection_row_id"))
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  stem = tools::file_path_sans_ext(basename(crop_file))
  prediction_file = file.path(output_dir, paste0(stem, "_with_irrigation.csv"))
  write.csv(final_predictions, prediction_file, row.names = FALSE)
  
  message("Saved predictions: ", prediction_file)
  invisible(final_predictions)
}
