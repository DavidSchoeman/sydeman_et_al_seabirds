# Making a 0.5º equal area global grid
	# Written by David Schoeman (david.schoeman@gmail.com)
		# March 2021

# Heavily based on code by Isaac Brito Morales (i.britomorales@uq.edu.au) and Jason Eevertt (jason.everett@uq.edu.au)

# NOTE that this process takes significant time and disk space. Derived products are output at the end for use in subsequent scripts.

# Source the functions ---------------------------------------------------------

	source("0_Seabird_Helpers.R")


# Set up elements we need ------------------------------------------------------

	# Shapes
		robCRS <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
		world <- ne_countries(scale = "medium", returnclass = "sf") %>% 
			st_transform(robCRS) # Robinson-projected world
		Bndry <- tibble(x = seq(-180, 180, by = 1), y = -90) %>%
			bind_rows(tibble(x = 180, y = seq(-90, 90, by = 1))) %>%
			bind_rows(tibble(x = seq(180, -180, by = -1), y = 90)) %>%
			bind_rows(tibble(x = -180, y = seq(90, -90, by = -1))) %>%
			as.matrix() %>%
			list() %>%
			st_polygon() %>%
			st_sfc(crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") %>%
			st_transform(crs = robCRS) # Robinson-projected global boundary


# Make a landmass object --------------------------------------------------

	r <- raster::raster(res = .1) # A 0.1 deg raster
	r[] <- 1:raster::ncell(r) # Fill it
	land_rs <- fasterize(ne_countries(scale = "medium", returnclass = "sf"), r) # Rasterise the world shapefile
	land_rs[] <- ifelse(is.na(land_rs[]), NA, 1) # only land cells!
	land_rs <- raster::setValues(raster(land_rs), land_rs[])
	land <- as(land_rs,  "SpatialPolygonsDataFrame") # Convert to a spatial polygon
	land$layer <- seq(1, length(land)) # A label for land

	# Now convert to an sf object and create ONE BIG polygon that we can use to populate with PUs
		landmass <- st_as_sf(land) %>% 
			select(layer) %>% 
			summarise(total_layer = sum(layer, do_union = TRUE)) %>% 
			st_transform(landmass, crs = robCRS)
		# We can plot the object to see if it is correct
		ggplot() +
			geom_sf(data = landmass) # Seems about right
		Store(landmass) # Write to cache so we don't lose it
		

# Create the hexagons and check ------------------------------------------------

	PUs <- fCreate_PlanningUnits(Bndry = Bndry, LandMass = landmass, CellArea = (111.324/2)^2, Shape = "Hexagon") # Roughly equivalent to 0.5-deg squares
	
	# Plot the hexagons to check
	gg_pus <- ggplot() +
		geom_sf(data = Bndry, aes(colour = "Planning Boundary"), fill = NA, size = 0.1, show.legend = "line") +
		geom_sf(data = PUs, aes(colour = "Planning Units"), fill = NA, size = 0.5, show.legend = "line") +
		geom_sf(data = world, aes(colour = "Continent"), alpha = 0.5, size = 0.1, show.legend = "line") +
		scale_colour_manual(values = c("Planning Boundary" = "red",
																	 "Planning Units" = "lightblue",
																	 "Continent" = "grey20")) +
		theme_bw() +
		theme(legend.key = element_rect(fill = NA, colour = NA, size = 2),
					legend.position = c(0.15, 0.11),
					legend.background = element_blank(),
					legend.text = element_text(colour = "black", size = 12),
					legend.title = element_blank(),
					axis.title = element_blank()) +
		guides(colour = guide_legend(override.aes = list(size = 1)))
	gg_pus

	# Save the shapefiles for later use
		st_write(PUs, dsn = "Shapes/Hex_half_Deg", driver = "ESRI Shapefile")
		st_write(landmass, dsn = "Shapes/Land_Rob", driver = "ESRI Shapefile")

# Goto 6_Hemispheric_Asymmetry