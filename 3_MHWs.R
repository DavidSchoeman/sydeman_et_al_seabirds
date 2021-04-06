# Computation of Marine Heatwave (MHW) statistics
	# Written by David Schoeman (david.schoeman@gmail.com)
			# May 2020 - March 2021

# Heavily based on code from https://robwschlegel.github.io/heatwaveR/
# Requires the installation of Climate Data Operators (cdo) on host machine https://code.mpimet.mpg.de/projects/cdo/wiki
# OISST data would need to be downloaded by the user from https://www.ncdc.noaa.gov/oisst

# NOTE that this process takes significant time and disk space. Derived products are output at the end for use in subsequent scripts.

# Source the functions ---------------------------------------------------------

	source("0_Seabird_Helpers.R")


# Folders where files can be initially be read/written--------------------------

	oisst_d <- "/Volumes/Data/OISST_download"	
	oisst <- "/Volumes/Data/OISST"	
	oisst_y <- "/Volumes/Data/OISST_year"	
	
	
# Merge daily OISST into annual bins -------------------------------------------
	# First, extract just sst to save space, then clean up
	ff <- dir(oisst_d, pattern = "oisst-avhrr-v02r01.")
	registerDoMC(cores = detectCores()-2)
		plyr::ldply(ff, .fun = dOISST, .parallel = TRUE)
	
	# Next, combine daily OISST data into years - otherwise my machine falls over
		yrs <- 1982:2017
		registerDoMC(cores = detectCores()-2)
			plyr::ldply(yrs, .fun = mOISST, .parallel = TRUE)
	# Next, merge all annual files into a single netCDF
		system(paste0("cdo mergetime ", oisst_y, "/OISST* ", oisst_y, "/OISST_1982_2017.nc")) # Merge years
		file.remove(dir(oisst_y, full.names = TRUE)[-1]) # Clean up individual years
	# Next, reproject to same 0.5º grid and clean up
		system(paste0("cdo remapbil,global_0.5 ", oisst_y, "/OISST_1982_2017.nc ", oisst_y, "/OISST_1982_2017_0_5_deg.nc")) # Bilinear interpolation


# Break global grid into blocks to process separately --------------------------

		oFold <- "mhwTemp" # A temporary output path...this is going to generate a LOT of files!
		
	# Break world into manageable blocks for processing
		f <- paste0(oisst_y, "/OISST_1982_2017_0_5_deg.nc") # The combined netCDF you created above
			plyr::ldply(f, .fun = bWorld1, .parallel = FALSE) # Break the world into blocks
		

# For each pixel in each box, extract the time series of SST --------------

		fff <- dir(oFold, pattern = "BLOCK_", full.names = TRUE, recursive = TRUE)
		dO <- paste0("cdo timmean ", f, " Mask.nc") # Get the overall climatology to use as a land mask
			system(dO)
		r <- stack("Mask.nc") # Load the mask
		pts <- as.data.frame(rasterToPoints(r)) # A data frame from just the cells we want
		names(pts)[3] <- "sst" # Rename ssts
		pts <- na.omit(pts)[,1:2] # This gives us coordinates where there are ssts
		registerDoMC(cores = detectCores()-2)	
			plyr::ldply(fff, .fun = pPixel, .parallel = TRUE) # Collect time series for each pixel in the box
		file.remove("Mask.nc") # Remove the mask

				
# For each time series, detect MHWs --------------------------------------------
		
		SST_files <- dir(oFold, pattern = ".csv", full.names = TRUE, recursive = TRUE) # List the csvs for each pixel
		d <- as.Date("1982/01/01"):as.Date("2017/12/31") # Make a dat string****change the dates here, if needed
		d <- as.Date(d, origin = "1970/01/01") # Write them as dates
		nd <- length(d) # How many dates are there? 
		Dates <- data.frame(Date = d) # Need it as a data frame - used in the function below
		registerDoMC(cores = detectCores()-2)	
			out <- plyr::ldply(SST_files, .fun = do_stats, .parallel = TRUE)
		saveRDS(out, file = "Data/MHW_0_5.Rda") # Output the resultant data
		file.remove(SST_files) # Remove the mask
		rm(out, SST_files) # Clean up
	
			
# Summarise MHW stats -----------------------------------------------------------------------------------
		
		MHW_result <- readRDS("Data/MHW_0_5.Rda")
		event_sum <- MHW_result %>%
			mutate(year = lubridate::year(date_start)) %>%
			group_by(x, y, year) %>%
			dplyr::summarise(n = n_distinct(event_no),
											 Peak_Int = mean(intensity_max),
											 Cum_Int = sum(intensity_cumulative),
											 Duration = mean(as.numeric(date_end-date_start)),
											 Cum_Dur = sum(as.numeric(date_end-date_start)))
		saveRDS(event_sum, file = "Data/MHW_stats_0_5.Rda")


# Compute HMW trends and write rasters ---------------------------------

	event_sum	<- readRDS("Data/MHW_stats_0_5.Rda")
	registerDoMC(cores = detectCores()-2)
	nTrend <- plyr::ddply(event_sum, c("x", "y"), lin_n, .parallel = TRUE) %>%
		dplyr::rename(nslope = slope, np = p)
	iTrend <- plyr::ddply(event_sum, c("x", "y"), lin_int, .parallel = TRUE) %>%
		dplyr::rename(islope = slope, ip = p)
	ciTrend <- plyr::ddply(event_sum, c("x", "y"), lin_cumi, .parallel = TRUE) %>%
		dplyr::rename(cislope = slope, cip = p)
	dTrend <- plyr::ddply(event_sum, c("x", "y"), lin_dur, .parallel = TRUE) %>%
		dplyr::rename(dslope = slope, dp = p)
	cdTrend <- plyr::ddply(event_sum, c("x", "y"), lin_cdur, .parallel = TRUE) %>%
		dplyr::rename(cdslope = slope, cdp = p)
	# Join, then rasterise
		out <- left_join(nTrend, iTrend) %>%
			left_join(., ciTrend) %>%
			left_join(., dTrend) %>%
			left_join(., cdTrend)
		r <- rasterFromXYZ(out, crs = projection(raster())) # Convert to a raster
		writeRaster(r, file = "Data/MHW_trends.grd", overwrite = TRUE)	# Write the raster
		  saveRDS(r, "Data/MHW_trends.rda")

		
# Clean up to save space in the Data folder on GitHub ***Skip if y --------
	file.remove(c("Data/MHW_0_5.Rda", "Data/MHW_stats_0_5.Rda"))
								
# Goto 4_CHI.R
