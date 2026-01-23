


library(FEISTY)
library(sf)
# Some other libraries may be required for plotting.

# Convert NetCDF Results to FEISTY Simulation Object
#
# Converts NetCDF output to a FEISTY-compatible sim object that can be used
# with standard FEISTY plotting functions like plotBiomasstime().
#
# Parameters:
#   ncFile - Path to the NetCDF file
#   groupNames - Optional vector of display names for functional types.
#     Default: c("Small pelagic", "Mesopelagic", "Large pelagic", "Bathypelagic", "Demersal")
#   timeUnit - Unit of time in the NC file. Options: "year", "day", "hour", "second".
#     Time will be converted to years for plotting. Default: "day"
#
# Returns:
#   A list (sim object) compatible with FEISTY plotting functions:
#     t - Time vector (in years)
#     p - Parameters list with groupnames, colors, etc.
#     R - Resource biomass matrix
#     totBiomass - Total fish biomass matrix
#     SSB - Spawning stock biomass (same as totBiomass)
#     yield - Yield matrix (filled with zeros)
#
# Example:
#   sim <- convertNCtoSim("result.nc")
#   plotBiomasstime(sim)
convertNCtoSim <- function(ncFile,
                           groupNames = c("Small pelagic", "Mesopelagic",
                                          "Large pelagic", "Bathypelagic", "Demersal"),
                           timeUnit = "day") {

  if (!file.exists(ncFile)) {
    stop("NetCDF file not found: ", ncFile)
  }

  # =========================================================================
  # Open NetCDF file
  # =========================================================================
  nc <- ncdf4::nc_open(ncFile)
  on.exit(ncdf4::nc_close(nc), add = TRUE)

  # =========================================================================
  # Read time and convert to calendar years
  # =========================================================================
  time_raw <- nc$dim$time$vals

  # Try to get time units attribute (e.g., "seconds since 1963-01-01 00:00:00")
  time_units <- tryCatch(
    ncdf4::ncatt_get(nc, "time", "units")$value,
    error = function(e) NULL
  )

  # Parse time units to extract reference date and unit
  if (!is.null(time_units) && grepl("since", time_units)) {
    parts <- strsplit(time_units, " since ")[[1]]
    unit <- tolower(trimws(parts[1]))
    ref_date_str <- trimws(parts[2])

    ref_date <- tryCatch(
      as.POSIXct(ref_date_str, tz = "UTC"),
      error = function(e) as.POSIXct("1963-01-01", tz = "UTC")
    )

    time_seconds <- switch(unit,
      "seconds" = time_raw, "second" = time_raw,
      "days" = time_raw * 86400, "day" = time_raw * 86400,
      "hours" = time_raw * 3600, "hour" = time_raw * 3600,
      time_raw
    )

    time_posix <- ref_date + time_seconds
    time <- as.numeric(format(time_posix, "%Y")) +
            (as.numeric(format(time_posix, "%j")) - 1) / 365

    cat("  Time units from NC:", time_units, "\n")
  } else {
    time <- switch(timeUnit,
      "year" = time_raw, "day" = time_raw / 365,
      "hour" = time_raw / (365 * 24), "second" = time_raw / (365 * 24 * 3600),
      time_raw
    )
    cat("  Time units: manual conversion (", timeUnit, " -> year)\n", sep = "")
  }

  # =========================================================================
  # Find and read total biomass variables
  # =========================================================================
  totB_vars <- grep("^fish_fft_[0-9]+_totB$", names(nc$var), value = TRUE)
  totB_vars <- totB_vars[order(as.numeric(gsub("fish_fft_([0-9]+)_totB", "\\1", totB_vars)))]

  if (length(totB_vars) == 0) {
    stop("No fish_fft_X_totB variables found in NetCDF file")
  }

  nGroups <- length(totB_vars)
  nTime <- length(time)

  # Read biomass data into matrix (time x groups)
  totBiomass <- matrix(NA, nrow = nTime, ncol = nGroups)
  for (i in seq_along(totB_vars)) {
    totBiomass[, i] <- as.vector(ncdf4::ncvar_get(nc, totB_vars[i]))
  }

  # =========================================================================
  # Create group names and identifiers
  # =========================================================================
  fft_ids <- paste0("fft_", 1:nGroups)

  # Adjust groupNames length if needed
  if (length(groupNames) < nGroups) {
    groupNames <- c(groupNames, paste0("FFT_", (length(groupNames)+1):nGroups))
  } else if (length(groupNames) > nGroups) {
    groupNames <- groupNames[1:nGroups]
  }

  # =========================================================================
  # Create color palette (matching FEISTY style)
  # =========================================================================
  default_colors <- c(
    "fft_1" = "#33BBEE",
    "fft_2" = "#009988",
    "fft_3" = "#EE7733",
    "fft_4" = "#CC3311",
    "fft_5" = "#EE3377"
  )
  my_palette <- default_colors[1:nGroups]
  names(my_palette) <- fft_ids

  my_names <- groupNames
  names(my_names) <- fft_ids

  # =========================================================================
  # Create dummy resource data (required by getTimeseries)
  # =========================================================================
  nResources <- 1
  R <- matrix(0, nrow = nTime, ncol = nResources)

  # =========================================================================
  # Create u0 (initial values)
  # =========================================================================
  u0_resources <- 0
  names(u0_resources) <- "dummy_resource"

  u0_fish <- rep(1, nGroups)
  names(u0_fish) <- paste0(fft_ids, "_1")

  u0 <- c(u0_resources, u0_fish)

  # =========================================================================
  # Build parameters list (p)
  # =========================================================================
  all_groupnames <- c("dummy_resource", fft_ids)

  full_palette <- c("dummy_resource" = "#CCCCCC", my_palette)
  full_names <- c("dummy_resource" = "Resource", my_names)

  p <- list(
    nResources = nResources,
    nGroups = nGroups,
    groupnames = all_groupnames,
    my_palette = full_palette,
    my_names = full_names,
    ixR = 1,
    u0 = u0
  )

  # =========================================================================
  # Build simulation object (sim)
  # =========================================================================
  sim <- list(
    t = time,
    p = p,
    R = R,
    totBiomass = totBiomass,
    SSB = totBiomass,
    yield = matrix(0, nrow = nTime, ncol = nGroups),
    nTime = nTime
  )

  # Add class for S3 method dispatch
  class(sim) <- c("FEISTY", "list")

  cat("Converted NC file to FEISTY sim object:\n")
  cat("  Time steps:", nTime, "\n")
  cat("  Time range:", round(min(time), 1), "-", round(max(time), 1), "(years)\n")
  cat("  Functional types:", nGroups, "\n")
  cat("  Groups:", paste(groupNames, collapse = ", "), "\n")

  return(sim)
}


# Plot Global Fish Biomass Map from NetCDF
#
# Creates a global map of total fish biomass using data from a NetCDF file.
# Uses Mollweide projection with viridis color scale.
#
# Parameters:
#   ncFile - Path to the NetCDF file
#   groupNames - Optional vector of group names to include. If NULL, uses all
#     fish_fft_X_totB variables found. Default: NULL (use all groups)
#   logScale - Whether to use log10 scale for biomass. Default: TRUE
#   cellSize - Grid cell size in degrees. Default: NULL (auto-detect from data)
#   title - Optional plot title. Default: NULL
#
# Returns:
#   A ggplot object of the global fish biomass map
#
# Example:
#   library(FEISTY)
#   plotGlobalBiomassMap("output.nc")
plotGlobalBiomassMap <- function(ncFile,
                                  groupNames = NULL,
                                  logScale = TRUE,
                                  cellSize = NULL,
                                  title = NULL) {

  if (!file.exists(ncFile)) {
    stop("NetCDF file not found: ", ncFile)
  }

  # =========================================================================
  # Open NetCDF file
  # =========================================================================
  nc <- ncdf4::nc_open(ncFile)
  on.exit(ncdf4::nc_close(nc), add = TRUE)

  # =========================================================================
  # Read coordinates and convert to real degrees
  # =========================================================================
  # Get grid dimensions from x and y
  nx <- nc$dim$x$len
  ny <- nc$dim$y$len

  # Try to read lon/lat variables
  lon_raw <- NULL
  lat_raw <- NULL

  if ("lon" %in% names(nc$var)) {
    lon_raw <- ncdf4::ncvar_get(nc, "lon")
  } else if ("longitude" %in% names(nc$var)) {
    lon_raw <- ncdf4::ncvar_get(nc, "longitude")
  }

  if ("lat" %in% names(nc$var)) {
    lat_raw <- ncdf4::ncvar_get(nc, "lat")
  } else if ("latitude" %in% names(nc$var)) {
    lat_raw <- ncdf4::ncvar_get(nc, "latitude")
  }

  # Convert to actual degree coordinates
  # Calculate grid resolution
  dlon <- 360.0 / nx  # degrees per grid cell in longitude
  dlat <- 180.0 / ny  # degrees per grid cell in latitude

  # Create 1D coordinate vectors (cell centers)
  lon <- seq(from = dlon/2, by = dlon, length.out = nx)        # 0 to 360

  lat <- seq(from = -90 + dlat/2, by = dlat, length.out = ny)  # -90 to 90

  cat("  Grid resolution:", round(dlon, 2), "x", round(dlat, 2), "degrees\n")
  cat("  Grid size:", nx, "x", ny, "\n")

  # =========================================================================
  # Find and read total biomass variables
  # =========================================================================
  totB_vars <- grep("^fish_fft_[0-9]+_totB$", names(nc$var), value = TRUE)
  totB_vars <- totB_vars[order(as.numeric(gsub("fish_fft_([0-9]+)_totB", "\\1", totB_vars)))]

  if (length(totB_vars) == 0) {
    stop("No fish_fft_X_totB variables found in NetCDF file")
  }

  # =========================================================================
  # Read biomass data and calculate total
  # =========================================================================
  # Get dimensions - assuming [x, y] or [x, y, time]
  var_dims <- nc$var[[totB_vars[1]]]$dim
  ndims <- length(var_dims)

  # Get fill value
  fill_val <- ncdf4::ncatt_get(nc, totB_vars[1], "_FillValue")
  fill_value <- if (fill_val$hasatt) fill_val$value else -2e+20

  # Read all groups and sum
  total_biomass <- NULL
  for (var_name in totB_vars) {
    biomass_data <- ncdf4::ncvar_get(nc, var_name)

    # Replace fill values with NA
    biomass_data[biomass_data < -1e+10] <- NA

    # If 3D (x, y, time), take the last time step or mean of last portion
    if (ndims == 3) {
      ntime <- dim(biomass_data)[3]
      # Use mean of last 40% of time steps (like Global_run_illustration.R)
      start_idx <- max(1, round(0.6 * ntime))
      biomass_data <- apply(biomass_data[, , start_idx:ntime, drop = FALSE], c(1, 2), mean, na.rm = TRUE)
    }

    # Replace NaN with NA
    biomass_data[is.nan(biomass_data)] <- NA

    if (is.null(total_biomass)) {
      total_biomass <- biomass_data
    } else {
      # Sum, treating NA as 0 for addition (NA + value = NA otherwise)
      total_biomass <- ifelse(is.na(total_biomass) & is.na(biomass_data), NA,
                              ifelse(is.na(total_biomass), 0, total_biomass) +
                              ifelse(is.na(biomass_data), 0, biomass_data))
    }
  }

  # =========================================================================
  # Create data frame with lon, lat, biomass
  # =========================================================================
  # Use expand.grid to create coordinate pairs (matches column-major order of biomass)
  out <- expand.grid(lon = lon, lat = lat)
  out$tot <- as.vector(total_biomass)

  # Remove NA and zero/negative values
  out <- out[!is.na(out$tot) & out$tot > 0, ]

  # Adjust longitude for plotting (shift to -180 to 180 range)
  out$lon <- ifelse(out$lon > 179.5, out$lon - 360, out$lon)

  # =========================================================================
  # Convert to sf and create grid (matching Figure 4 approach)
  # =========================================================================
  out_sf <- sf::st_as_sf(out, coords = c("lon", "lat"), crs = 4326)

  # Auto-detect cell size from data resolution if not specified
  if (is.null(cellSize)) {
    lon_sorted <- sort(unique(out$lon))
    lat_sorted <- sort(unique(out$lat))
    if (length(lon_sorted) > 1) {
      cellSize <- median(diff(lon_sorted))
    } else if (length(lat_sorted) > 1) {
      cellSize <- median(diff(lat_sorted))
    } else {
      cellSize <- 2.8125  # default for typical global model
    }
    cat("  Auto-detected cell size:", round(cellSize, 2), "degrees\n")
  }

  # Create grid
  grid <- sf::st_as_sf(
    out_sf %>% sf::st_make_grid(cellsize = cellSize, what = "polygons")
  ) %>%
    sf::st_join(out_sf, join = sf::st_intersects, left = TRUE)

  # Load world map
  world <- rnaturalearth::ne_countries(scale = 50, returnclass = "sf")

  # =========================================================================
  # Create plot with Mollweide projection (like Figure 4)
  # =========================================================================
  if (logScale) {
    p <- ggplot2::ggplot(grid) +
      ggplot2::geom_sf(ggplot2::aes(fill = log10(tot + 0.01)), colour = NA) +
      viridis::scale_fill_viridis(
        name = bquote(atop("Fish biomass g m"^-2, "")),
        labels = c("<0.01", "1", "100"),
        breaks = c(-2, 0, 2)
      )
  } else {
    p <- ggplot2::ggplot(grid) +
      ggplot2::geom_sf(ggplot2::aes(fill = tot), colour = NA) +
      viridis::scale_fill_viridis(
        name = bquote("Fish biomass (g m"^-2*")")
      )
  }

  p <- p +
    ggplot2::geom_sf(data = world, col = "grey", fill = "grey") +
    ggplot2::coord_sf(crs = "+proj=moll") +
    ggplot2::theme(
      plot.background = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(size = 8)
    )

  if (!is.null(title)) {
    p <- p + ggplot2::ggtitle(title)
  }

  cat("Plotted global fish biomass map:\n")
  cat("  Grid points:", nrow(out), "\n")
  cat("  Lon range:", round(min(out$lon), 1), "-", round(max(out$lon), 1), "\n")
  cat("  Lat range:", round(min(out$lat), 1), "-", round(max(out$lat), 1), "\n")
  cat("  Biomass range:", round(min(out$tot, na.rm = TRUE), 4), "-",
      round(max(out$tot, na.rm = TRUE), 2), "g/m2\n")
  cat("  Functional types summed:", length(totB_vars), "\n")

  return(p)
}
