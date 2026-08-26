#' Normalize county GEOIDs
#'
#' Converts county GEOIDs to five-character strings with leading zeros.
#'
#' @param x A vector of county GEOIDs.
#'
#' @return A character vector of normalized five-digit GEOIDs.
#' @export
normalize_geoid = function(x) {
  x = as.character(x)
  x = gsub("\\.0$", "", x)
  x[is.na(x) | x == "" | x == "NA"] = NA_character_
  
  ok = !is.na(x)
  
  x[ok] = sprintf("%05d", as.integer(x[ok]))
  
  x
}


#' Adjust day-of-year values that span the year boundary
#'
#' Detects distributions containing observations near both the beginning and
#' end of the year and shifts early-year observations forward by one year so
#' that summary statistics can be calculated continuously.
#'
#' @param x Numeric vector of day-of-year values.
#' @param low_cutoff Lower day-of-year threshold used to detect wrapping.
#' @param high_cutoff Upper day-of-year threshold used to detect wrapping.
#'
#' @return A numeric vector of adjusted day-of-year values.
#' @export
fix_wrap_doy = function(x, low_cutoff = 45, high_cutoff = 320) {
  x = x[!is.na(x)]
  
  if (length(x) == 0) {
    return(x)
  }
  
  wraps =
    stats::quantile(x, 0.05, na.rm = TRUE) <= low_cutoff &&
    stats::quantile(x, 0.95, na.rm = TRUE) >= high_cutoff
  
  if (wraps) {
    x[x <= low_cutoff] = x[x <= low_cutoff] + 365
  }
  
  x
}


#' Wrap adjusted day-of-year values back to the calendar year
#'
#' @param x Numeric vector of adjusted day-of-year values.
#'
#' @return Integer day-of-year values between 1 and 365.
#' @export
wrap_back_doy = function(x) {
  ((round(x) - 1) %% 365) + 1
}


#' Calculate the statistical mode
#'
#' @param x A vector.
#'
#' @return The most frequently occurring non-missing value, or `NA` if no
#'   non-missing values are present.
#' @export
mode_value = function(x) {
  x = x[!is.na(x)]
  
  if (length(x) == 0) {
    return(NA)
  }
  
  ux = unique(x)
  
  ux[which.max(tabulate(match(x, ux)))]
}


#' Convert year and day-of-year to calendar dates
#'
#' @param year Integer vector of calendar years.
#' @param doy Numeric vector of day-of-year values.
#'
#' @return A `data.table::IDate` vector.
#' @export
doy_to_date = function(year, doy) {
  out = rep(as.Date(NA), length(year))
  
  ok = !is.na(year) & !is.na(doy)
  
  out[ok] =
    as.Date(paste0(year[ok], "-01-01")) +
    as.integer(round(doy[ok])) -
    1L
  
  data.table::as.IDate(out)
}