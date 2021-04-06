# Helper functions for analysis of global seabird breeding success
	# Written by David Schoeman (david.schoeman@gmail.com) and Brian Hoover (bhoover@faralloninstitute.org)
		# Input from Bill Sydeman (wsydeman@comcast.net) and Sarah-Ann Thompson (sathompson@faralloninstitute.org)
			# June 2020 - March 2021

# Load all of the packages required for the analyses ---------------------------

	library(nlme)
	library(MASS)
	library(car)
	library(mblm)

	library(sf)
	library(raster)
	library(fasterize)
	library(rnaturalearth)
	library(rnaturalearthdata)

	library(gridExtra)
	library(ggeffects)
	library(ggthemes)
	library(tmap)  

	library(SOAR)
	library(tidyverse)
	library(doMC)

	library(VoCC) # devtools::install_github("JorGarMol/VoCC", dependencies = TRUE, build_vignettes = TRUE)
	library(heatwaveR)
	library(ncdf4)


# Function to rearrange factor levels ------------------------------------------

	reFactor <- function(x, y = 1:length(levels(x))) { # x is a factor, y is the desired sequence of factor levels
		l <- levels(x) # What levels are there?
		n <- length(l) # How many fators are there?
		z <- c(y, (1:n)[-which(l %in% l[y])]) # Place the factors you want at the front at the front 
		f <- factor(x, levels = l[z]) # Rearrange the factor levels 
		return(f) # Output the factor
	}
	
	
# A function to build a plotting data frame for mixed-effects mode -------------
	# Modified from https://bbolker.github.io/mixedmodels-misc/ecostats_chap.html
	
	pltmm <- function(mod, d, alpha = 0.05) { # mod is the model name, d is the data frame
		# Build prediction data frame
		m <- formula(mod,fixed.only = TRUE)[-2] # Fixed effects from mod
		fs <- all.vars(m) # Get just the names for the fixed effects
		# Figure out which effects are factor and which are continuous, and fill, as necessary with values fromt he data
		out <- list() # A list to collect output
		for(i in 1:length(fs)) { # For each fixed effect
			if(eval(parse(text = paste0("with(d,is.factor(", fs[i], "))")))) {
				out[[fs[i]]] <- eval(parse(text = paste0("with(d, ", fs[i], " <- levels(", fs[i], "))"))) # Fill in the factor levels for factors
			} else {
				out[[fs[i]]] <- eval(parse(text = paste0("with(d, ", fs[i], " <- seq(min(", fs[i], "), max(", fs[i], "), length.out = 250))")))	# If it is continuous make a sequence of values from min to max, lenght 250
			}
		}
		p <- expand.grid(out) # Make predictor data frame from the fixed effects
		mm <- model.matrix(m, p) # Model matrix for predictors
		beta <- fixef(mod) # Fixed-effects coefficients
		y <- mm %*% beta # Predicted values 
		V <- vcov(mod) # Variance-covariance matrix of beta
		pred.se <- sqrt(diag(mm %*% V %*% t(mm))) # Std errors of predictions
		if(!is.null(mod$family)) {linv <- mod$family$linkinv} # If there is a link function, get the inverse-link function
		# Construct 95% Normal CIs on the link scale and transform back to the response (probability) scale
		crit <- -qnorm(alpha/2) # Critical value for 95% CIs
		fits <- data.frame(y = y, # Predicted value from model
											 se.hi =  y + pred.se, # Upper error limit
											 se.lw =  y - pred.se, # Lower error limit
											 upr = y + crit*pred.se, # Approx upper 95% conf limit for fit
											 lwr = y - crit*pred.se # Approx lower 95% conf limit for fit
		)
		Fits <- fits
		if(!is.null(mod$family)) {Fits <- as.data.frame(lapply(fits, linv))} # Where there is an error family, pass fits through the inverse-link function
		return(cbind(p, Fits))
	}
			
		
# A function to extract sst from daily OISST netCDFs for MHW -------------------
	# Requires the installation of CDO
	
	dOISST <- function(x) {
		system(paste0("cdo select,name=sst ", oisst_d, "/", x, " ", oisst, "/", x)) # Select just sst and make a new copy
		if(file.exists(paste0(oisst, "/", x))) {file.remove(paste0(oisst_d, "/", x))} # Delete old copy to save space
	}
		
		
# A function to merge daily OISST into annual bins for MHW ---------------------
	# Requires the installation of CDO
		
	mOISST <- function(x) {
		system(paste0("cdo mergetime ", oisst, "/oisst-avhrr-v02r01.", x, "* ", oisst_y, "/OISST", x, ".nc"))
	}

							
# A function to break world into bits for MHW ----------------------------------
		
	bWorld1 <- function(i) { # For a spliced netCDF, and start and end years
		if(!isTRUE(file.info(oFold)$isdir)) dir.create(oFold, recursive = TRUE) # If the folder doesn't exist, create it			
		lons <- seq(-180, 180, 20) # Longitudinal blocks
		lats <- seq(-90, 90, 20) # Latitudinal blocks
		for(j in 1:(length(lons)-1)) {
			xmin <- lons[j]
			xmax <- lons[j + 1]
			for(k in 1:(length(lats)-1)) {
				ymin <- lats[k]
				ymax <- lats[k +1]
				oName <- paste0(oFold, "/BBLOCK_", xmin, "_", xmax, "_", ymin, "_", ymax, ".nc") # Output file name
				system(paste0("cdo -f nc -sellonlatbox,", xmin, ",", xmax, ",", ymin, ",", ymax, " ", i, " ", oName)) # Crop OISST data to box and write to netCDF
			}
		}
	}
	
			
# Function to extract single-coordinate time series for MHW---------------------
		
	pPixel <- function(i) {
		# # # First get just the cells we want
		pth <- oFold
		blk <- gsub(".nc", "", i) %>% 
			gsub("mhwTemp/BBLOCK_", "", .)
		bits <- unlist(strsplit(blk, "_"))
		xmin <- as.numeric(bits[1])
		xmax <- as.numeric(bits[2])
		ymin <- as.numeric(bits[3])
		ymax <- as.numeric(bits[4])
		toDo <- pts[which(pts$x >= xmin & pts$x <= xmax & pts$y >= ymin & pts$y <= ymax),]
		if(nrow(toDo) > 0) { # If there are some cells in this block that are in the mask, write the block, then extract sst for each pixel
			nc <- nc_open(i)
			for(j in 1:nrow(toDo)) {
				out <- ncvar_get(nc, varid = "sst",
												 start= c(which(nc$dim$lon$vals == toDo[j,1]), # look for long
												 				 which(nc$dim$lat$vals == toDo[j,2]),  # look for lat
												 				 1,1),
												 count = c(1,1,1,-1)) #count '-1' means 'all values along that dimension'that dimension'
				pName <- paste0(toDo[j,1], "_", toDo[j,2])
				pName1 <- pName # Take a copy so that file names can include a "-"
				if(length(grep("-" , pName)) > 0) {pName <- gsub("-", "neg", pName)}
				if(grepl("\\d", substr(pName, 1, 1))) {pName <- paste0("pos", pName)} # Avoid starting an object name with a number
				eval(parse(text = paste0(pName, " <- out")))
				eval(parse(text = paste0("write.csv(", pName, ", '", pth, "/", pName1, ".csv', row.names = FALSE)")))
			}
		}
		file.remove(i)
	}	
		

# Detect MHWs per time series --------------------------------------------
	# Based on https://cran.r-project.org/web/packages/heatwaveR/vignettes/gridded_event_detection.html
	# First, the wrapper function to compute baseline climatology and then detect MHWs
	event_only <- function(df){
		# First calculate the climatologies
		clim <- ts2clm(data = df, x = Date, y = SST, climatologyPeriod = c("1983-01-01", "2012-12-31")) # As per Oliver et al. and Smale et al.
		# Then the events
		event <- detect_event(data = clim, x = Date, y = SST)
		# Last, we return only the event dataframe of results
		return(event$event)
	}

				
# Compute MHW stats per time series --------------------------------------------
		
	do_stats <- function(file_name){
	sst <- read.csv(file_name)
	if(length(na.omit(sst$x)) == nd) { # If ALL values are available for the cell
		bits <- file_name %>% 
			gsub(paste0(oFold, "/"), "", .) %>% 
			gsub(".csv", "", .)
		bits <- unlist(strsplit(bits, "_"))
		sstDat <- data.frame(Date = Dates$Date,
												 x = as.numeric(bits[1]),
												 y = as.numeric(bits[2]),
												 SST = round(sst$x, 2))
		event <- event_only(sstDat) %>% 
			mutate(x = sstDat$x[1], y = sstDat$y[1]) %>% 
			dplyr::select(x, y, everything())
		}
	}
		

# Functions to compute MHW trends in parallel ----------------------------------
		# Using mblm instead of lm because it is more robust
	
	# Number of MHWs
		lin_n <- function(ev) {
			mod1 <- mblm(n ~ year, data = filter(ev, Cum_Dur <= 365))
			# extract slope coefficient and its p-value
			tr <- data.frame(slope = summary(mod1)$coefficients[2,1],
											 p = summary(mod1)$coefficients[2,4])
			return(tr)
		}
	# Intensity of MHWs
		lin_int <- function(ev) {
			mod1 <- mblm(Peak_Int ~ year, data = filter(ev, Cum_Dur <= 365))
			# extract slope coefficient and its p-value
			tr <- data.frame(slope = summary(mod1)$coefficients[2,1],
											 p = summary(mod1)$coefficients[2,4])
			return(tr)
		}
	# Cumulative intensity of MHWs
		lin_cumi <- function(ev) {
			mod1 <- mblm(Cum_Int ~ year, data = filter(ev, Cum_Dur <= 365))
			# extract slope coefficient and its p-value
			tr <- data.frame(slope = summary(mod1)$coefficients[2,1],
											 p = summary(mod1)$coefficients[2,4])
			return(tr)
		}			
	# Mean duration of MHWs
		lin_dur <- function(ev) {
			mod1 <- mblm(Duration ~ year, data = filter(ev, Cum_Dur <= 365))
			# extract slope coefficient and its p-value
			tr <- data.frame(slope = summary(mod1)$coefficients[2,1],
											 p = summary(mod1)$coefficients[2,4])
			return(tr)
		}
	# Cumulative duration of MHWs
		lin_cdur <- function(ev) {
			mod1 <- mblm(Cum_Dur ~ year, data = filter(ev, Cum_Dur <= 365))
			# extract slope coefficient and its p-value
			tr <- data.frame(slope = summary(mod1)$coefficients[2,1],
											 p = summary(mod1)$coefficients[2,4])
			return(tr)
		}		
		
# Reproject CHI geoTIFFs -------------------------------------------------------

	getCHI <- function(x) {
		m <- "ngb" # Nearest-neighbour method
		nm <- x %>% 
			gsub(".tif", "", .) %>% 
			strsplit(., "CHI_Tiff/") %>% 
			unlist() %>% 
			.[2] %>% 
			gsub("-", "_", .) # Get the name of the TIFF
		r <- stack(x) # Load the TIFF into raster
		dO <- paste0(nm, " <- projectRaster(r, crs = lonlat, method = ", "'ngb'", ", over = FALSE)") # Reproject to WGS84 lonlat
		eval(parse(text = dO))
		dO <- paste0("Store(", nm, ")") # Write the output to cache
		eval(parse(text = dO))
	}

				
# Hexgrid function -------------------------------------------------------------
	# Function to create square or heaxagonal planning units for your area of interest.
	# Inputs needed are:
	#   Bndry: An sf polygon object which outlines the limits of the study area.
	#       The code takes the bbbox so the limits are the most important.
	#       The output inherits the crs from this sf object so ensure it is in the correct projection for your needs
	#   LandMass: An sf multipolygon object which contains all the areas (ie land) that you wish to remove from the grid.
	#       The code assumes that any Planning Unit whose centroids is over land will be removed. This approximates > 50% of the PU is landward.
	#   CellArea: The area in km you wish your resultant Planning Units to be.
	#   Shape: Hexagon or Square
	#
	# Written by Isaac Brito Morales (i.britomorales@uq.edu.au) and Jason Eevertt (jason.everett@uq.edu.au) (UQ/UNSW/CSIRO)
	# Written: December 2020
	# Edited/updated: by Dave Schoeman March 2021
	
	fCreate_PlanningUnits <- function(Bndry, LandMass, CellArea, Shape){
		if(Shape %in% c("hexagon", "Hexagon")){
			sq <- FALSE
			diameter <- 2 * sqrt((CellArea*1e6)/((3*sqrt(3)/2))) * sqrt(3)/2 # Diameter in m's
		}
		if(Shape %in% c("square", "Square")){
			sq <- TRUE
			diameter <- sqrt(CellArea*1e6) # Diameter in m's
		}
		# First create planning units for the whole region
		PUs <- st_make_grid(Bndry,
												square = sq,
												cellsize = c(diameter, diameter),
												what = "polygons") %>%
			st_sf()
		# Check cell size worked ok.
		print(paste0("Range of cellsize are ",
								 round(as.numeric(range(units::set_units(st_area(PUs), "km^2")))[1])," km2 to ",
								 round(as.numeric(range(units::set_units(st_area(PUs), "km^2")))[2])," km2")) # Check area
		# First get all the PUs partially/wholly within the planning region
		logi_Reg <- st_centroid(PUs) %>%
			st_intersects(Bndry) %>%
			lengths > 0 # Get logical vector instead of sparse geometry binary
		PUs <- PUs[logi_Reg == TRUE, ]
		# Second, get all the pu's with < 50 % area on land (approximated from the centroid) ### THIS is where things slow down!!!
		logi_Ocean <- st_centroid(PUs) %>%
			st_within(st_union(LandMass)) %>%
			lengths > 0 # Get logical vector instead of sparse geometry binary
		PUs <- PUs[logi_Ocean==FALSE, ]
		return(PUs)
	}		
	
	
# A function to truncate values in rasters for plotting, where needed	 ---------
	
	tRast <- function(x, z = .975, y = stats::quantile(abs(x[]), z, na.rm = TRUE)) { # Function to truncate the values in a raster, based on percentiles
		xx <- which(abs(x[]) > y)
		x[xx] <- y * sign(x[xx])
		X <- x
		return(X)
	}
		
		
# A function for constructing masking polygons on maps -------------------------
		
	mkSF <- function(ymin, ymax, 
									 xmin = -180, xmax = 180, 
									 ny = abs(xmax-xmin)*20, nx = abs(ymax-ymin)*20) {
		x <- c(seq(xmin, xmax, length = nx),
					 rep(xmax, ny),
					 seq(xmax, xmin, length = nx),
					 rep(xmin, ny))
		y <- c(rep(ymin, ny),
					 seq(ymin, ymax, length = nx),
					 rep(ymax, ny),
					 seq(ymax, ymin, length = nx))
		return(st_polygon(list(cbind(x, y))))
	}			



# A function to plot maps of hemispheric asymmetry ----------------------------

	figPlt <-function(hex = hex, x, msk = m,
									pFill = "lightgrey",
									bgFill = "transparent",
									dCol = "white",
									pal = "-RdBu") { 
	dO <- paste0("p <- tm_shape(hex) + ",
							 "tm_fill(x, lwd = NA, n = 15, palette = pal, legend.is.portrait = FALSE) + ",
							 "tm_shape(world, proj = rob) + ",
							 "tm_polygons() + ",
							 "tm_shape(b0) + ",
							 "tm_polygons(fill = pFill, alpha = .35) + ",
							 "tm_shape(b1) + ",
							 "tm_polygons(fill = pFill, alpha = .35) + ",
							 "tm_shape(b2) + ",
							 "tm_polygons(fill = pFill, alpha = .35) + ",
							 "tm_shape(dLine1) + ",
							 "tm_lines() + ",
							 "tm_shape(dLine2) + ",
							 "tm_lines() + ",
							 "tm_shape(dLine3) + ",
							 "tm_lines() + ",
							 "tm_shape(dLine4) + ",
							 "tm_lines() + ",
							 "tm_shape(XY) + ",
							 "tm_dots(shape = 21, size = 0.2, col = dCol, alpha = .7)")
	eval(parse(text = dO))
	return(p)
	}	
		

# A function to find coordinates for centroid of nearest non-NA cell -----------

	nearestsstcell <- function(r, crds, type = "land") { # Where r is a raster layer, xy are SpatialPoints (or just a matrix), all projected lonlat, and type is either "land" (default) or "sea"
		if(type == "land") {
			coast <- boundaries(r, type = "inner")
		} else {
			coast <- boundaries(r, type = "outer")			
		}
		coast[coast[] == 0] <- NA # Replace all 0s in coast with NAs
		cmat <- data.frame(xyFromCell(r, 1:ncell(r)), Temp = r[])
		cmat <- na.omit(cmat)
		tsp <- SpatialPoints(cmat[,1:2])
		if(class(crds) == "SpatialPoints") {
			XY <- crds@coords
		} else {
			XY <- crds
		}
		ccells <- apply(XY, 1, FUN = function(x) {which.min(pointDistance(tsp, x, lonlat = TRUE))})
		return(coordinates(tsp[ccells]))
	}		
		
		
# A function to extract density data from hex grid -----------------------------

	gDens <- function(V = "RoW") {
		dsxy <- data.frame(dat = sh_dens[[V]])
		dnxy <- data.frame(dat = nh_dens[[V]])
		pdsxy <- round(quantile(dsxy, c(.5, .1, .9), na.rm = TRUE), 2)
		pdnxy <- round(quantile(dnxy, c(.5, .1, .9), na.rm = TRUE), 2)
		psth <- round(quantile(filter(sites, Hem == "South")[[V]], c(.5, .1, .9), na.rm = TRUE), 2)
		pnth <- round(quantile(filter(sites, Hem == "North")[[V]], c(.5, .1, .9), na.rm = TRUE), 2)
		return(list(nth = dnxy, sth = dsxy, all_north = pdnxy, sites_north = pnth, all_south = pdsxy, sites_south = psth))
	}
		


# A function to plot density data ----------------------------------------------

	pltDens <- function(x, ll = 0.001, ul = 0.999, xLab = "Insert x-axis label here",
											nCol = "#440154FF",
											nFill = "#440154FF",
											sCol = "#35B779FF",
											sFill = "#35B779FF",
											pBG = "transparent") {	
		d <- gDens(x)
		dO <- paste0("p <- ggplot() +
		geom_density(data = d$sth,
			aes(x = dat),
			color = sCol, fill = sFill,
			alpha = .45) +
		geom_density(data = d$nth,
			aes(x = dat),
				color = nCol, fill = nFill,
				alpha = .45) +
		labs(x = xLab) +
		xlim(quantile(c(unlist(d$sth), unlist(d$nth)), c(ll, ul), na.rm = TRUE)) +
		theme_classic() +
		theme(axis.title.y = element_blank(),
			axis.text.y = element_blank(),
			axis.ticks.y = element_blank(),
			plot.background = element_rect(fill = pBG, colour = NA),
			panel.background = element_rect(fill = pBG, colour = NA))")
		eval(parse(text = dO))
		return(p)
	}

				
# Custom function to get p-value for linaer model ------------------------------
		
	getP <- function(x) {
	  f <- summary(x)$fstatistic
	  p <- pf(f[1], f[2], f[3], lower.tail = FALSE)
	  attributes(p) <- NULL
	  return(p)
	}