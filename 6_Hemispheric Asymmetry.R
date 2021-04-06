# Plots of Hemispheric Assymetry
	# Written by David Schoeman (david.schoeman@gmail.com) and Brian Hoover (bhoover@faralloninstitute.org)
		# Input from Bill Sydeman (wsydeman@comcast.net) and Sarah-Ann Thompson (sathompson@faralloninstitute.org)
		# June 2020 - March 2021

# Depends on data output by 3_MHWs.R, 4_CHI.R and 5_HexGrid.R


# Source the functions ---------------------------------------------------------
	
	source("0_Seabird_Helpers.R")


# Make the saved data available ------------------------------------------------

	RoW <- stack("Data/RoW.grd")
	vocc <- stack("Data/vocc.grd")
	mhw <- stack("Data/MHW_trends.grd") 
	CHI2003 <- raster("Data/CHI2003.grd")
	CHItrend <- raster("Data/CHItrend.grd")
	xy <- readRDS("Data/xy.rda")


# Load world shapefile and prepare projection ----------------------------------

	rob <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs +towgs84=0,0,0" # Robinson projection
	lonlat <- projection(raster()) # Standard WGS84 lonlat
	world <- st_read("Shapes/Land_Rob/Land_Rob.shp") # Created in 5_HexGrid.R
		st_crs(world) <- rob # Assign the correct projection


# Set up masked areas and outline boundaries -----------------------------------

	# Masking boxes
		b0 <- st_sfc(mkSF(-15, 15),
								 crs = lonlat) %>%
			st_transform(crs = rob) # Robinson-projected equatorial mask polygon
		b1 <- st_sfc(mkSF(75, 90),
								 crs = lonlat) %>%
			st_transform(crs = rob) # Robinson-projected north-pole mask polygon
		b2 <- st_sfc(mkSF(-90, -75),
								 crs = lonlat) %>%
			st_transform(crs = rob)	# Robinson-projected south-pole mask polygon

	# Global bounding lines
		dLine1 <- data.frame(x = -179.999, y = seq(-90, 90, .1), attr = 1, id = "a") %>% 
			st_as_sf(coords = c("x", "y"), crs = tmaptools::get_proj4("longlat", output = "character")) %>%
			dplyr::summarize(m = mean(attr),do_union=FALSE) %>% st_cast("LINESTRING")
		dLine2 <- data.frame(x = 179.999, y = seq(-90, 90, .1), attr = 1, id = "a") %>% 
			st_as_sf(coords = c("x", "y"), crs = tmaptools::get_proj4("longlat", output = "character")) %>%
			dplyr::summarize(m = mean(attr),do_union=FALSE) %>% st_cast("LINESTRING")
		dLine3 <- data.frame(x = c(-179.999, 179.999), y = c(90, 90), attr = 1, id = "a") %>% 
			st_as_sf(coords = c("x", "y"), crs = tmaptools::get_proj4("longlat", output = "character")) %>%
			dplyr::summarize(m = mean(attr),do_union=FALSE) %>% st_cast("LINESTRING")
		dLine4 <- data.frame(x = c(-179.999, 179.999), y = c(-90, -90), attr = 1, id = "a") %>% 
			st_as_sf(coords = c("x", "y"), crs = tmaptools::get_proj4("longlat", output = "character")) %>%
			dplyr::summarize(m = mean(attr),do_union=FALSE) %>% st_cast("LINESTRING")
		
	# Hemispheric boxes for density plots	
		n_hem <- st_sfc(mkSF(15, 75),
										crs = lonlat) %>%
			st_transform(crs = rob) # Northern hemisphere, excluding masked areas
		s_hem <- st_sfc(mkSF(-75, -15),
										crs = lonlat) %>%
			st_transform(crs = rob) # Southern hemisphere, excluding masked areas
			
	# A global mask to crop repeated data after Robinson projection
		df <- cbind(c(-180, rep(180, 360), -180, rep(-180, 360)),
								c(90, seq(90, -90, length.out = 360), -90, seq(-90, 90, length.out = 360))) # Coordinates with high-res meridional outline
		m <- st_sf(data.frame(outline = "Outline",
													st_sfc(st_polygon(list(df)))),
							 crs = 4326) %>% 
			st_transform(crs = rob) %>%  # Convert matrix into a reprojected spatial polygon
			as_Spatial()


# Prepare the data for plotting-------------------------------------------------

	R1 <- RoW[[1]]*10 # RoW per decade
			R1 <- projectRaster(R1, crs = rob, method = "ngb", over = FALSE) # Project to Robinson
			r1 <- tRast(R1) # Truncate extremes for plotting
			
	R2 <- vocc[[1]]*10 # VoCC per decade
			R2 <- projectRaster(R2, crs = rob, method = "ngb", over = FALSE) # Project to Robinson
			r2 <- tRast(R2, z = .9) # Truncate extremes for plotting
			
	R3 <- mhw[[9]]*10 # Cumulative MHW duration (i.e., number of MHW days) per decade (i.e., number of MHW days)
			R3 <- projectRaster(R3, crs = rob, method = "ngb", over = FALSE) # Project to Robinson
			r3 <- tRast(R3) # Truncate extremes for plotting
			
	R4 <- mhw[[5]]*10 # Cumulative MHW intensity per decade
			R4 <- projectRaster(R4, crs = rob, method = "ngb", over = FALSE) # Project to Robinson
			r4 <- tRast(R4) # Truncate extremes for plotting
			
	R5 <- projectRaster(CHI2003, crs = rob, method = "ngb", over = FALSE) # Project CHI to Robinson
		r5 <- tRast(R5) # Truncate extremes for plotting
	
	R6 <- CHItrend *10 # CHI trend per decade
			R6 <- projectRaster(R6, crs = rob, method = "ngb", over = FALSE) # Project to Robinson		
			r6 <- tRast(R6) # Truncate extremes for plotting
			
	XY <- st_as_sf(xy, coords = c("x", "y"), crs = lonlat) %>% # Sample locations
		st_transform(crs = rob) # The locations of the sites
			
	# For each variable find centre of nearest non-NA cell to avoid NA values for each Site, where possible, then project Robinson
		rowXY <- nearestsstcell(RoW[[1]], dplyr::select(xy, x, y)) %>% 
			data.frame() %>% 
			bind_cols(dplyr::select(xy, Hemisphere)) %>% 
			st_as_sf(., coords = c("x", "y"), crs = lonlat) %>%
			st_transform(crs = rob)
		voccXY <- nearestsstcell(vocc[[1]], dplyr::select(xy, x, y)) %>% 
			data.frame() %>% 
			bind_cols(dplyr::select(xy, Hemisphere)) %>% 
			st_as_sf(., coords = c("x", "y"), crs = lonlat) %>%
			st_transform(crs = rob)
		cdurXY <- nearestsstcell(mhw[[9]], dplyr::select(xy, x, y)) %>% 
			data.frame() %>% 
			bind_cols(dplyr::select(xy, Hemisphere)) %>% 
			st_as_sf(., coords = c("x", "y"), crs = lonlat) %>%
			st_transform(crs = rob)
		cintXY <- nearestsstcell(mhw[[5]], dplyr::select(xy, x, y)) %>% 
			data.frame() %>% 
			bind_cols(dplyr::select(xy, Hemisphere)) %>% 
			st_as_sf(., coords = c("x", "y"), crs = lonlat) %>%
			st_transform(crs = rob)
		chiXY <- nearestsstcell(CHI2003, dplyr::select(xy, x, y)) %>% 
			data.frame() %>% 
			bind_cols(dplyr::select(xy, Hemisphere)) %>% 
			st_as_sf(., coords = c("x", "y"), crs = lonlat) %>%
			st_transform(crs = rob)
		chitXY <- nearestsstcell(CHItrend, dplyr::select(xy, x, y)) %>% 
			data.frame() %>% 
			bind_cols(dplyr::select(xy, Hemisphere)) %>% 
			st_as_sf(., coords = c("x", "y"), crs = lonlat) %>%
			st_transform(crs = rob)
	
	
# Regrid to centroids from hex grid for analysis and plotting ------------------

	grd <- st_read("Shapes/Hex_half_Deg/Hex_half_Deg.shp") # The half-degree hex-grid shapefile
		st_crs(grd) <- rob # Assign the correct CRS
	ctrds <- st_centroid(grd) # Get centroids of hex-grid cells
	hex <- grd %>% 
		mutate(RoW = raster::extract(r1, ctrds),
					 vocc = raster::extract(r2, ctrds),
					 mdays = raster::extract(r3, ctrds),
					 mint = raster::extract(r4, ctrds),
					 cimp = raster::extract(r5, ctrds),
					 cimptrnd = raster::extract(r6, ctrds)) %>% 
		na.omit # Extract the data for each hex-grid point and remover NAs

	# Repeat for northern hemisphere - on non-truncated data
		nh_dens <- st_intersection(hex, n_hem) %>% 
			mutate(RoW = raster::extract(R1, st_centroid(.)),
						 vocc = raster::extract(R2, st_centroid(.)),
						 mdays = raster::extract(R3, st_centroid(.)),
						 mint = raster::extract(R4, st_centroid(.)),
						 cimp = raster::extract(R5, st_centroid(.)),
						 cimptrnd = raster::extract(R6, st_centroid(.)))
		
	# Repeat for southern hemisphere - on non-truncated data
		sh_dens <- st_intersection(hex, s_hem) %>% 
			mutate(RoW = raster::extract(R1, st_centroid(.)),
						 vocc = raster::extract(R2, st_centroid(.)),
						 mdays = raster::extract(R3, st_centroid(.)),
						 mint = raster::extract(R4, st_centroid(.)),
						 cimp = raster::extract(R5, st_centroid(.)),
						 cimptrnd = raster::extract(R6, st_centroid(.)))

	# Repeat for positions of sample sites
		sites <- data.frame(Hem = XY$Hemisphere) %>%  
			mutate(RoW = raster::extract(R1, rowXY),
						 vocc = raster::extract(R2, voccXY),
						 mdays = raster::extract(R3, cdurXY),
						 mint = raster::extract(R4, cintXY),
						 cimp = raster::extract(R5, chiXY),
						 cimptrnd = raster::extract(R6, chitXY)) 
		# Each variable has one NA remaining, but these would have little effect on estimates, and likely less effect from projecting back to lonlat and trying again
	
	
# Plot the maps for each variable-----------------------------------------------
	
	# RoW
		p1 <- figPlt(hex = hex, x = "RoW") + 
			tm_layout(legend.width = 2, frame = FALSE, 
								legend.outside = TRUE, 
								legend.outside.position = "bottom",
								legend.format = list(scientific = FALSE, digits = 2, text.separator = ""),
								panel.label.bg.color = NA, bg.color = "transparent")
		tmap_save(p1, "Figs/FigRoW.pdf", dpi = 600, bg = "transparent", colormodel = "cmyk")

	# VoCC
		p2 <- figPlt(hex = hex, x = "vocc") +
			tm_layout(legend.width = 2, frame = FALSE, 
								legend.outside = TRUE, 
								legend.outside.position = "bottom",
								legend.format = list(scientific = FALSE, digits = 2, text.separator = ""),
								panel.label.bg.color = NA, bg.color = "transparent")
		tmap_save(p2, "Figs/FigVoCC.pdf", dpi = 600, bg = "transparent", colormodel = "cmyk")
		
	# CumHHWDurSlope
		p3 <- figPlt(hex = hex, x = "mdays") + 
			tm_layout(legend.width = 2, frame = FALSE, 
								legend.outside = TRUE, 
								legend.outside.position = "bottom",
								legend.format = list(scientific = FALSE, digits = 2, text.separator = ""),
								panel.label.bg.color = NA, bg.color = "transparent")
		tmap_save(p3, "Figs/FigMHWCDur.pdf", dpi = 600, bg = "transparent", colormodel = "cmyk")

	# CumHHWItensitySlope
		p4 <- figPlt(hex = hex, x = "mint") + 
			tm_layout(legend.width = 2, frame = FALSE, 
								legend.outside = TRUE, 
								legend.outside.position = "bottom",
								legend.format = list(scientific = FALSE, digits = 2, text.separator = ""),
								panel.label.bg.color = NA, bg.color = "transparent")
		tmap_save(p4, "Figs/FigMHWCIint.pdf", dpi = 600, bg = "transparent", colormodel = "cmyk")
		
	# HumanImpact2003
		p5 <- figPlt(hex = hex, x = "cimp") + 
			tm_layout(legend.width = 2, frame = FALSE, 
								legend.outside = TRUE, 
								legend.outside.position = "bottom",
								legend.format = list(scientific = FALSE, digits = 2, text.separator = ""),
								panel.label.bg.color = NA, bg.color = "transparent")
		tmap_save(p5, "Figs/FigCImp2003.pdf", dpi = 600, bg = "transparent", colormodel = "cmyk")
		

	# RecentHumanImpactSlope
		p6 <- figPlt(hex = hex, x = "cimptrnd") + 
			tm_layout(legend.width = 2, frame = FALSE, 
								legend.outside = TRUE, 
								legend.outside.position = "bottom",
								legend.format = list(scientific = FALSE, digits = 2, text.separator = ""),
								panel.label.bg.color = NA, bg.color = "transparent")
		tmap_save(p6, "Figs/FigCImpSlope.pdf", dpi = 600, bg = "transparent", colormodel = "cmyk")
		

# Density plots and stats -----------------------------------------------------------

	# RoW
		pltDens("RoW", xLab = "Rate of ocean warming (°C per decade)") + 
			theme(panel.grid.minor = element_blank(),
						panel.grid.major = element_blank(),
						plot.background = element_rect(fill = "transparent", colour = NA))
		ggsave("Figs/RoWDens.pdf", width = 4, height = 2, dpi = 600, colormodel = "cmyk")
		dDat <- gDens("RoW")
		sink("Figs/RoW.txt")
			dDat$all_north
			dDat$sites_north
			dDat$all_south
			dDat$sites_south
		sink()
	
	# VoCC
		pltDens("vocc", xLab = "Velocity of ocean warming (km per decade)", ll = .01, ul = .99) + 
			theme(panel.grid.minor = element_blank(),
						panel.grid.major = element_blank(),
						plot.background = element_rect(fill = "transparent", colour = NA))
		ggsave("Figs/VoCCDens.pdf", width = 4, height = 2, dpi = 600, colormodel = "cmyk")
		dDat <- gDens("vocc")
		sink("Figs/VoCC.txt")
		dDat$all_north
		dDat$sites_north
		dDat$all_south
		dDat$sites_south
		sink()
	
	# CumHHWDurSlope
		pltDens("mdays", xLab = "Rate of cumulative MHW extension (days per decade)") + 
			theme(panel.grid.minor = element_blank(),
						panel.grid.major = element_blank(),
						plot.background = element_rect(fill = "transparent", colour = NA))
		ggsave("Figs/CumHHWDur.pdf", width = 4, height = 2, dpi = 600, colormodel = "cmyk")
		dDat <- gDens("mdays")
		sink("Figs/CumHHWDur.txt")
		dDat$all_north
		dDat$sites_north
		dDat$all_south
		dDat$sites_south
		sink()
	
	# CumHHWIIntensitySlope
		pltDens("mint", xLab = "Rate of cumulative MHW intensification (°C.days per decade)", ll = .005, ul = .995) + 
			theme(panel.grid.minor = element_blank(),
						panel.grid.major = element_blank(),
						plot.background = element_rect(fill = "transparent", colour = NA))
		ggsave("Figs/CumHHWIntens.pdf", width = 4, height = 2, dpi = 600, colormodel = "cmyk")
		dDat <- gDens("mint")
		sink("Figs/CumHHWIntens.txt")
		dDat$all_north
		dDat$sites_north
		dDat$all_south
		dDat$sites_south
		sink()
	
	# HumanImpact2003
		pltDens("cimp", xLab = "Cumulative human impact Index 2003") + 
			theme(panel.grid.minor = element_blank(),
						panel.grid.major = element_blank(),
						plot.background = element_rect(fill = "transparent", colour = NA))
		ggsave("Figs/CumHumImp2003.pdf", width = 4, height = 2, dpi = 600, colormodel = "cmyk")
		dDat <- gDens("cimp")
		sink("Figs/CumHumImp2003.txt")
		dDat$all_north
		dDat$sites_north
		dDat$all_south
		dDat$sites_south
		sink()
	
	# RecentHumanImpactSlope
		pltDens("cimptrnd", xLab = "Change in cumulative human impact Index (per decade))") + 
			theme(panel.grid.minor = element_blank(),
						panel.grid.major = element_blank(),
						plot.background = element_rect(fill = "transparent", colour = NA))
		ggsave("Figs/CumHumImpSlope.pdf", width = 4, height = 2, dpi = 600, colormodel = "cmyk")
		dDat <- gDens("cimptrnd")
		sink("Figs/CumHumImpSlope.txt")
		dDat$all_north
		dDat$sites_north
		dDat$all_south
		dDat$sites_south
		sink()
		
# All done
