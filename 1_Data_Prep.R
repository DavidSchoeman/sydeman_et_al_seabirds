# Preparing data for the analysis of global seabird breeding success
	# Written by David Schoeman (david.schoeman@gmail.com) and Brian Hoover (bhoover@faralloninstitute.org)
		# Input from Bill Sydeman (wsydeman@comcast.net) and Sarah-Ann Thompson (sathompson@faralloninstitute.org)
			# June 2020 - March 2021


# Source the functions ---------------------------------------------------------
	source("0_Seabird_Helpers.R")


# Load the data, combine and prep ----------------------------------------------

	dat1 <- read.csv("Data/global_seabird_data.csv", stringsAsFactors = TRUE) %>% 
		rename(ID = idnum, # Rename variables so that they're slightly easier to deal with
					 ts = timeseries_num,
					 year = year,
					 yearno = yearnum,
					 site = site,
					 siteno = site_num,
					 lat = latitude,
					 lon = longitude ,
					 group = site_group,
					 spp = species,
					 sppcode = species_code,
					 sppno = species_num,
					 bs = breeding_success,
					 Hemisphere = hemisphere) %>% 
		mutate(sppsite = as.factor(paste(spp, site, sep = "_"))) # A new factor variable per time series
		sum(duplicated(dat1[, -c(1:2)])) # Check to see whether there are any duplicates: none found
		
	# Merge in the trophic levels
		dat2 <- left_join(dat1, read.csv("Data/Lookup_table.csv", stringsAsFactors = TRUE),
											by = c("spp" = "Species")) %>%  # Merge Biome and TL into dat
			rename(Depth = Foraging.depth, spp.name = Scientific.name, TL = Trophic.level)
			dat2$group <- reFactor(dat2$group, c(3, 1, 2, 5, 4, 7, 8, 9, 11, 10, 6)) # Rearrange levels to more coherent geographic distribution
			dat2$TL <- reFactor(dat2$TL, c(3, 1, 2)) # Rearrange levels to more coherent tropic distribution
			dat2$Depth <- reFactor(dat2$Depth, c(2, 1)) # Rearrange Depth levels
				levels(dat2$Depth) <- c("Surface", "Deep") # Capitalise the factor levels
				levels(dat2$Hemisphere) <- c("North", "South")

# Create data required for modelling -------------------------------------------
				
		out <- list() 
		for(i in levels(dat2$sppsite)) {
			d <- filter(dat2, sppsite == i) %>% # For each time series
				drop_na(bs) %>%  # Drop NAs from the response
				mutate(pr_failure = ifelse(bs > mean(bs)*.1, 0, 1), # Code breeding failure
							 stbs = scale(bs), # Compute scaled breeding success
							 nyear = n_distinct(year)) # Compute number of data points
			out[[i]] <- d
			}
		dat <- bind_rows(out) # Bind the elements of the list back into a data frame


# Check against selection criteria and select just the series we need ----------

	# Which series are shorter than the 20-year cut-off?
		dat %>% 
			group_by(sppsite) %>% 
			summarise(minYr = min(year), maxYr = max(year)) %>% 
			mutate(lengthYr = maxYr - minYr + 1) %>% 
			data.frame() %>% 
			arrange(lengthYr)
	
	# We need to drop the short time series for least auklets at St George - see SOM
		dat <- dat %>%
			filter(sppsite != "least auklet_St. George")
	
	# Save the data	
		saveRDS(dat, file = "Data/dat.Rda")

# Make lookup table of sppsite and coordinates ---------------------------------

		xy <- dat %>%
			group_by(sppsite, Hemisphere, TL) %>% # For each time series
			dplyr::summarise(spp = spp[1], # Get the species
								site = site[1], # Get the site
								syr = min(year), # Get the start year
								eyr = max(year), # Get the end year
								x = mean(lon), # Get the longitude
								y = mean(lat), # Get the laitude
								depth = Depth[1], # Get the depth of foraging for that species
								nyear = mean(nyear)) %>% # Get the number of data points
			as.data.frame() # Make it a data frame
		saveRDS(xy, file = "Data/xy.Rda")

			
# Explore data to identify best 50 years to use for RoW/VoCC -------------------
		
		hist(dat$year) # A histogram of observation years
		# Look at cumulative proportions of observations by year
			x <- dat$year
			cbind(Freq = table(x), Cumul = cumsum(table(x)), relative = prop.table(table(x)), relprop = round(100*cumsum(prop.table(table(x))),1))
		# 1968:2017 looks optimal
		filter(dat, year >= 1968 & year <= 2017) %>% 
			summarise(prop = length(year)/nrow(dat))
		# Contains 99.7% of all observations, so we'll go with that
		
		
# Compute Rate of Warming (RoW) and VoCC ---------------------------------------
	
	# NOTE: Although code and instructions are provided here to work with raw data, derived products for rate of warming and velocity of ocean warming are output at the end and can be used "as is" in subsequent scripts
	
		# Using HadISST
		# These data can be retrieved from https://www.metoffice.gov.uk/hadobs/hadisst/
		# We downloaded a netCDF of the data and loaded them as a raster stack
			# SST <- stack("HadISSTnetCDF.nc") # Read the downloaded data
			# NAvalue(SST) <- -1000.0 # Replace -10000 with NA
			# Of course, we cannot host these data here, so this needs to be done by each user
		
		yrs <- as.numeric(substr(names(SST), 2, 5)) # For which years do we have data?
		yr50 <- 1968:2017 # 50 comprising the time series we're after
		yrs50 <- (1:length(yrs))[yrs %in% yr50] # The indices of months within that time series
		SST <- subset(SST, yrs50) # Select only the data we need
		rYr <- sumSeries(SST, # Compute annual means - uses a function from package VoCC
										 p = "1968-01/20017-12",
										 yr0 = "1968-01-16",
										 l = nlayers(SST),
										 fun = function(x) colMeans(x, na.rm = TRUE),
										 freqin = "months", freqout = "years") # Mean annual sst
		RoW <- tempTrend(rYr, th = 10) # Another function from VoCC: provides the long-term local climate trends - annual
		sGrad <- spatGrad(rYr, th = 0.0001, projected = FALSE) # Another function from VoCC:provides the spatial gradient - for use in computing VoCC
		vocc <- gVoCC(RoW, sGrad) # Compute VoCC
		writeRaster(RoW, file = "Data/RoW.grd", overwrite = TRUE)
		saveRDS(RoW, "Data/RoW.rda")
		writeRaster(vocc, file = "Data/vocc.grd", overwrite = TRUE)
		saveRDS(RoW, "Data/VoCC.rda")
		
# Goto 2_Modelling_Seabird_Breeding

		
			
			