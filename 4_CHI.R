# Processing Cumulative Human Impact Data
	# Written by David Schoeman (david.schoeman@gmail.com)
		# May 2020 - March 2021
	
# Data last downloaded data on 9 March 2021 from https://knb.ecoinformatics.org/view/doi:10.5063/F19Z92TW
	# I saved these geoTIFFs in a subfolder called CHI_Tiff

# NOTE that this process takes significant time and disk space. Derived products are output at the end for use in subsequent scripts.

# Source the functions ---------------------------------------------------------

	source("0_Seabird_Helpers.R")


# Projections ---------------------------------------------------------------------------------

	lonlat <- projection(raster()) # A standard lon-lat projection
	mollCRS <- CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs") # Mollweide projection


# Reproject desired CHI geoTIFFs and write to cache ----------------------------
	
	chiFold <- "path/folder_name" # The folder where you placed the downloaded CHI geoTIFFs
	f <- dir(chiFold, pattern = ".tif", full.names = TRUE)
	registerDoMC(cores = 4)
		plyr::ldply(f, .fun = getCHI, .parallel = TRUE)


# Aggregate, reproject and save ------------------------------------------------
	
	# Give the two CHI rasters you want to use easier names
		CHI2003 <- cumulative_impact_2003
		RTrend <- cumulative_impact_trend_2003_2013

	# Aggregate 10x to reduce the data density a bit
		R2003_10 <- raster::aggregate(CHI2003, 10)
		RTrend_10 <- raster::aggregate(RTrend, 10)
	
	# Next, create a resampling grid at close to the new resolution to create a uniform grid, then aggregate to 0.5-deg		
			fn <- raster(resolution = 0.1)
			CHI2003_05deg <- resample(R2003_10, fn, method = "ngb") %>%  # Resample to regular 
				raster::aggregate(., 5) # Aggregate to 0.5-deg
			CHItrend_05deg <- resample(RTrend_10, fn, method = "ngb") %>%  # Resample to regular 
				raster::aggregate(., 5) # Aggregate to 0.5-deg
			writeRaster(CHI2003_05deg, "Data/CHI2003.grd")
			  saveRDS(CHI2003_05deg, "Data/CHI2003.rda")
			writeRaster(CHItrend_05deg, "Data/CHItrend.grd")
			  saveRDS(CHItrend_05deg, "Data/CHItrend.rda")
			

# Goto 5_HexGrid
			