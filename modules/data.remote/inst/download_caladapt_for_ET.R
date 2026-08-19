#using the 198 design points that has caladapt 2024-45 data already downloaded, and collapsing them into 
#county-level daily climate to use the spatial assumption that every parcel in a county experiences the same daily climate

pacman::p_load(ncdf4, sf, dplyr, tidyr, fs, purrr, furrr, future, lubridate, tigris, FAO56)

##this part is fairly short because its mostly sourcing helpers from irrigation_function.R  
source("/projectnb/dietzelab/ananyak/irrigation_functions.R")

##the caladapt 2024-45 info for irrigation estimates is already downloaded
base_dir = "/projectnb/dietzelab/ccmmf/ensemble/CalAdapt_runs/data_raw/CalAdaptWRF"

output_dir = '/projectnb/dietzelab/ananyak'


#---- geo assignments ----
#convert the Cal-Adapt lat/lon locations into county assignments, attach the county to the climate records, and create a lookup table for the design points 
site_dirs = dir_ls(base_dir, type = "directory")

##site_index will store each site_hash and their lat/lon info 
site_index = map_dfr(site_dirs, function(dir_path) {
  nc_files = dir_ls(dir_path, glob = "*.nc")
  
  if (length(nc_files) == 0) {return(NULL)}
  nc = nc_open(nc_files[1])
  
  lat = as.numeric(ncvar_get(nc, "latitude"))
  lon = as.numeric(ncvar_get(nc, "longitude"))
  
  nc_close(nc)
  tibble(site_hash = basename(dir_path), lat = lat, lon = lon)
})

##just for a check, should say 198 
cat(sprintf("Indexed %d unique Cal-Adapt design-point sites.\n", nrow(site_index)))

##get CA county info as of 2024
ca_counties = counties(state = "CA", cb = TRUE, class = "sf") %>%
  select(County = NAME)

##sites_sf now holds each site_hash's lat/lon and point geometry 
sites_sf = st_as_sf(site_index, coords = c("lon", "lat"), crs = 4326, remove = FALSE)

#double check both layers use same CRS
ca_counties = st_transform(ca_counties, st_crs(sites_sf))

##assign each site hash its county name 
sites_with_county = st_join(sites_sf, ca_counties, join = st_intersects, left = TRUE
) %>%
  st_drop_geometry()

##can see the counties that are still not assigned to a site 
missing_counties =
  ca_counties %>%
  filter(
    !County %in% sites_with_county$County)
print(missing_counties$County)

sites_projected = st_transform(sites_sf, 3310)
missing_counties_projected = st_transform(missing_counties, 3310)

missing_county_points = st_point_on_surface(missing_counties_projected)

##assign whats leftover using nearest neighbor
nearest_site_index = st_nearest_feature(missing_county_points, sites_projected)

missing_county_lookup = tibble(County = missing_counties$County, site_hash = sites_projected$site_hash[nearest_site_index])

internal_site_lookup = sites_with_county %>%
  st_drop_geometry() %>%
  filter(
    !is.na(County)
  ) %>%
  select(County, site_hash)

county_site_lookup = bind_rows(internal_site_lookup, missing_county_lookup
  ) %>%
  distinct()

##now have a lookup table with sites and their county + geo info
write.csv(county_site_lookup, file.path(output_dir, "caladapt_county_site_lookup.csv"), row.names = FALSE)

#double check if any sites still failed to match a county (should get 0 rows)
print(sites_with_county %>%
    filter(is.na(County)))

# ---- Process all NetCDF files ----
all_nc_files = dir_ls(base_dir, recurse = TRUE, glob = "*.nc")

cat(sprintf("Found %d NetCDF files.\n", length(all_nc_files)))

plan(multisession, workers = max(1, availableCores() - 1))

daily_climate_dataset = future_map_dfr(all_nc_files, process_nc_file, wind_height_m = 10, .progress = TRUE, .options = furrr_options(seed = TRUE))

# ---- Attach counties to daily climate records ----

daily_climate_county =
  daily_climate_dataset %>%
  inner_join(county_site_lookup, by = "site_hash")

county_daily_climate =
  daily_climate_county %>%
  group_by(County, date, model, scenario
  ) %>%
  summarise(ET0_mm = mean(ET0_mm, na.rm = TRUE
      ),
    precip_mm = mean(precip_mm, na.rm = TRUE
      ),
    mean_temp_c = mean(mean_temp_c, na.rm = TRUE
      ),
    min_temp_c = mean(T_min, na.rm = TRUE
      ),
    max_temp_c = mean(T_max, na.rm = TRUE
      ),
    
    n_climate_sites = n_distinct(site_hash), .groups = "drop"
  ) %>%
  
  mutate(Year = year(date), Week = isoweek(date)
  ) %>%
  
  rename(GCM = model, SSP = scenario)

write.csv(county_daily_climate, file.path(output_dir, "caladapt_county_daily_climate.csv"), row.names = FALSE)