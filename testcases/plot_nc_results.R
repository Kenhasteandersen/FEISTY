


library(FEISTY)
library(sf)
# Some other libraries may be required for plotting.

# Convert NetCDF to FEISTY Simulation Object
#
# Reads NetCDF output and converts to a FEISTY-compatible sim object for plotting.
# Uses setupVertical2() to get full parameters including theta matrix, enabling
# use with both plotBiomasstime() and plotNetwork().
#
# For single-location files, lon/lat are ignored. For multi-location files,
# extracts data at the grid cell closest to the specified coordinates.
#
# Parameters:
#   ncFile - Path to the NetCDF file. Default: "output.nc"
#   lon - Target longitude in degrees (-180 to 180 or 0 to 360). Optional for single-location files.
#   lat - Target latitude in degrees (-90 to 90). Optional for single-location files.
#
# Returns:
#   A list (sim object) compatible with FEISTY plotting functions:
#     t - Time vector (in years)
#     p - Parameters list from setupVertical2 (includes theta, groupnames, colors, etc.)
#     u - Full state variable matrix (resources + fish by size class)
#     R - Resource biomass matrix
#     totBiomass - Total fish biomass matrix
#     SSB - Spawning stock biomass (same as totBiomass)
#     yield - Yield matrix (filled with zeros)
#     nTime - Number of time steps
#     location - List with actual lon/lat of extracted grid cell
#
# Example:
#   # Single-location file (lon/lat optional)
#   sim <- ncToSim("result.nc")
#
#   # Multi-location file (lon/lat required)
#   sim <- ncToSim("global_output.nc", lon = -120, lat = 35)
#
#   # Use with plotting functions
#   plotBiomasstime(sim)
#   plotNetwork(sim)
ncToSim <- function(ncFile = "output.nc",
                    lon = NULL,
                    lat = NULL) {

  if (!file.exists(ncFile)) {
    stop("NetCDF file not found: ", ncFile)
  }

  # =========================================================================
  # Open NetCDF file
  # =========================================================================
  nc <- ncdf4::nc_open(ncFile)
  on.exit(ncdf4::nc_close(nc), add = TRUE)

  # =========================================================================
  # Read grid coordinates - handle different dimension naming conventions
  # =========================================================================
  dim_names <- names(nc$dim)

  # Determine dimension names for longitude and latitude
  if ("x" %in% dim_names && "y" %in% dim_names) {
    # Format 1: x/y dimensions (global grid, need to calculate lon/lat)
    nx <- nc$dim$x$len
    ny <- nc$dim$y$len
    dlon <- 360.0 / nx
    dlat <- 180.0 / ny
    lon_grid <- seq(from = dlon/2, by = dlon, length.out = nx)
    lat_grid <- seq(from = -90 + dlat/2, by = dlat, length.out = ny)
  } else if ("lon" %in% dim_names && "lat" %in% dim_names) {
    # Format 2: lon/lat dimensions (coordinates stored directly)
    lon_grid <- nc$dim$lon$vals
    lat_grid <- nc$dim$lat$vals
    nx <- length(lon_grid)
    ny <- length(lat_grid)
  } else if ("longitude" %in% dim_names && "latitude" %in% dim_names) {
    # Format 3: longitude/latitude dimensions
    lon_grid <- nc$dim$longitude$vals
    lat_grid <- nc$dim$latitude$vals
    nx <- length(lon_grid)
    ny <- length(lat_grid)
  } else {
    stop("Cannot find spatial dimensions. Expected 'x/y', 'lon/lat', or 'longitude/latitude'")
  }

  # Check if this is a single-point file
  if (nx == 1 && ny == 1) {
    # Single-location file: ignore lon/lat input
    cat("Single-location file detected.\n")
    cat("Location: lon =", round(lon_grid[1], 2), ", lat =", round(lat_grid[1], 2), "\n")
    ix <- 1
    iy <- 1
    actual_lon <- lon_grid[1]
    actual_lat <- lat_grid[1]
    display_lon <- actual_lon
  } else {
    # Multi-location file: lon/lat required
    if (is.null(lon) || is.null(lat)) {
      stop("This is a multi-location file. Please provide 'lon' and 'lat' arguments.")
    }

    # =========================================================================
    # Convert input longitude to match grid range if needed
    # =========================================================================
    target_lon <- lon
    # If grid is 0-360 and input is negative, convert
    if (min(lon_grid) >= 0 && target_lon < 0) {
      target_lon <- target_lon + 360
    }
    # If grid is -180 to 180 and input is > 180, convert
    if (min(lon_grid) < 0 && target_lon > 180) {
      target_lon <- target_lon - 360
    }

    # =========================================================================
    # Find closest grid cell
    # =========================================================================
    ix <- which.min(abs(lon_grid - target_lon))
    iy <- which.min(abs(lat_grid - lat))

    actual_lon <- lon_grid[ix]
    actual_lat <- lat_grid[iy]

    # Convert back to -180 to 180 for display if needed
    display_lon <- ifelse(actual_lon > 180, actual_lon - 360, actual_lon)

    cat("Target location: lon =", lon, ", lat =", lat, "\n")
    cat("Closest grid cell: lon =", round(display_lon, 2), ", lat =", round(actual_lat, 2), "\n")
    cat("Grid indices: ix =", ix, ", iy =", iy, "\n")
  }

  # =========================================================================
  # Read depth from NC file
  # =========================================================================
  depth <- 500  # default

  # Try zc variable first (layer center depths, 3D: x, y, z)
  if ("zc" %in% names(nc$var)) {
    zc <- ncdf4::ncvar_get(nc, "zc")
    if (length(dim(zc)) == 3) {
      zc_local <- zc[ix, iy, ]
      # Get deepest non-NA layer
      valid_zc <- zc_local[!is.na(zc_local)]
      if (length(valid_zc) > 0) {
        depth <- abs(min(valid_zc))  # deepest layer (most negative value)
      }
    }
  } else if ("z" %in% dim_names) {
    # Fallback to z dimension if values are actual depths (negative)
    z_vals <- nc$dim$z$vals
    if (min(z_vals) < 0) {
      depth <- abs(min(z_vals))
    }
  }
  cat("Depth:", round(depth, 1), "m\n")

  # =========================================================================
  # Initialize FEISTY setupVertical2 with depth from NC file
  # =========================================================================
  p <- setupVertical2(depth = depth)

  # =========================================================================
  # Read time and convert to calendar years
  # =========================================================================
  time_raw <- nc$dim$time$vals

  time_units <- tryCatch(
    ncdf4::ncatt_get(nc, "time", "units")$value,
    error = function(e) NULL
  )

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

    cat("Time units from NC:", time_units, "\n")
  } else {
    # Fallback: assume time is in days if no units attribute found
    time <- time_raw / 365
    cat("Time units: no attribute found, assuming days -> converting to years\n")
  }

  nTime <- length(time)
  nResources <- p$nResources
  nGroups <- p$nGroups

  # =========================================================================
  # Read resource data
  # =========================================================================
  R <- matrix(0, nrow = nTime, ncol = nResources)

  # Try to read resource variables (small zoo, large zoo, benthos)
  resource_vars <- c("fish_small_zooplankton_target_integrator_zoop_integrator_result",
                     "fish_large_zooplankton_target_integrator_zoop_integrator_result",
                     "fish_benthos")

  for (i in seq_along(resource_vars)) {
    if (i <= nResources && resource_vars[i] %in% names(nc$var)) {
      res_data <- ncdf4::ncvar_get(nc, resource_vars[i])
      if (length(dim(res_data)) == 3) {
        R[, i] <- res_data[ix, iy, ]
      } else if (length(dim(res_data)) == 1) {
        R[, i] <- res_data
      }
    }
  }

  # =========================================================================
  # Read total biomass variables at the specified location
  # =========================================================================
  totB_vars <- grep("^fish_fft_[0-9]+_totB$", names(nc$var), value = TRUE)
  totB_vars <- totB_vars[order(as.numeric(gsub("fish_fft_([0-9]+)_totB", "\\1", totB_vars)))]

  if (length(totB_vars) == 0) {
    stop("No fish_fft_X_totB variables found in NetCDF file")
  }

  # Read biomass data at the specific location
  totBiomass <- matrix(NA, nrow = nTime, ncol = nGroups)

  for (i in seq_along(totB_vars)) {
    if (i > nGroups) break
    # Read 3D data: [x, y, time]
    biomass_3d <- ncdf4::ncvar_get(nc, totB_vars[i])

    # Extract time series at the specific grid cell
    if (length(dim(biomass_3d)) == 3) {
      biomass_ts <- biomass_3d[ix, iy, ]
    } else if (length(dim(biomass_3d)) == 2) {
      # If 2D, assume it's [location, time] or similar
      biomass_ts <- biomass_3d[ix, ]
    } else {
      biomass_ts <- as.vector(biomass_3d)
    }

    # Replace fill values with NA
    biomass_ts[biomass_ts < -1e+10] <- NA

    totBiomass[, i] <- biomass_ts
  }

  # =========================================================================
  # Read fish biomass by size class
  # =========================================================================
  fish_vars <- grep("^fish_fft_[0-9]+_size_[0-9]+$", names(nc$var), value = TRUE)

  # Determine size classes per group from NC file
  size_classes_nc <- list()
  for (g in 1:nGroups) {
    pattern <- paste0("^fish_fft_", g, "_size_[0-9]+$")
    vars_g <- grep(pattern, fish_vars, value = TRUE)
    size_classes_nc[[g]] <- length(vars_g)
  }

  # Total fish state variables
  nFishStates <- sum(unlist(size_classes_nc))
  fish_matrix <- matrix(0, nrow = nTime, ncol = nFishStates)

  col_idx <- 1
  for (g in 1:nGroups) {
    for (s in 1:size_classes_nc[[g]]) {
      var_name <- paste0("fish_fft_", g, "_size_", s)
      if (var_name %in% names(nc$var)) {
        fish_data <- ncdf4::ncvar_get(nc, var_name)
        if (length(dim(fish_data)) == 3) {
          fish_matrix[, col_idx] <- fish_data[ix, iy, ]
        } else if (length(dim(fish_data)) == 1) {
          fish_matrix[, col_idx] <- fish_data
        }
      }
      col_idx <- col_idx + 1
    }
  }

  # Replace NA/negative with 0
  fish_matrix[is.na(fish_matrix) | fish_matrix < 0] <- 0
  R[is.na(R) | R < 0] <- 0

  # =========================================================================
  # Build state variable matrix u (resources + fish)
  # =========================================================================
  u <- cbind(R, fish_matrix)

  # =========================================================================
  # Adjust p$ix to match NC file structure if needed
  # =========================================================================
  # Rebuild ix based on actual size classes in NC file
  p$ix <- list()
  start_idx <- nResources + 1
  for (g in 1:nGroups) {
    n_sizes <- size_classes_nc[[g]]
    p$ix[[g]] <- start_idx:(start_idx + n_sizes - 1)
    start_idx <- start_idx + n_sizes
  }

  # =========================================================================
  # Build simulation object (sim)
  # =========================================================================
  sim <- list(
    t = time,
    p = p,
    u = u,
    R = R,
    totBiomass = totBiomass,
    SSB = totBiomass,
    yield = matrix(0, nrow = nTime, ncol = nGroups),
    nTime = nTime,
    location = list(
      lon = display_lon,
      lat = actual_lat,
      ix = ix,
      iy = iy
    )
  )

  class(sim) <- c("FEISTY", "list")

  cat("\nExtracted time series from NC file:\n")
  cat("  Location: (", round(display_lon, 2), ", ", round(actual_lat, 2), ")\n", sep = "")
  cat("  Time steps:", nTime, "\n")
  cat("  Time range:", round(min(time), 1), "-", round(max(time), 1), "(years)\n")
  cat("  Resources:", nResources, "\n")
  cat("  Fish groups:", nGroups, "\n")
  cat("  Size classes per group:", paste(unlist(size_classes_nc), collapse = ", "), "\n")
  cat("  Total state variables:", ncol(u), "\n")

  # Report biomass statistics
  valid_biomass <- totBiomass[!is.na(totBiomass) & totBiomass > 0]
  if (length(valid_biomass) > 0) {
    cat("  Biomass range:", round(min(valid_biomass), 4), "-",
        round(max(valid_biomass), 2), "g/m2\n")
  } else {
    cat("  Warning: No valid biomass data at this location\n")
  }

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
