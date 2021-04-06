# Data for Tables in the Supplementary Materials
	# Written by David Schoeman (david.schoeman@gmail.com) and Brian Hoover (bhoover@faralloninstitute.org)
		# Input from Bill Sydeman (wsydeman@comcast.net) and Sarah-Ann Thompson (sathompson@faralloninstitute.org)
			# June 2020 - March 2021

# Depends on data output by 1_Data_Prep.R


# Source the functions ---------------------------------------------------------
	
	source("0_Seabird_Helpers.R")


# Table S1 ---------------------------------------------------------------------
  
  S1 <- readRDS("Data/dat.rda") %>% 
    group_by(site, spp) %>% 
    summarise(site = unique(site),
              spp = unique(spp),
              slope = round(coefficients(lm(bs ~ year))[2], 4),
              p = round(getP(lm(bs ~ year)), 4),
              duration = paste(min(year), max(year), sep = "-"),
              nYrs = n()) %>% 
    mutate(sig = case_when(
      p <= 0.001 ~ "***",
      p <= 0.01 & p >= 0.001 ~ "**",
      p <= 0.05 & p >= 0.01 ~ "*",
      p >= 0.05 ~ "N.S."
      )
      ) %>% 
    arrange(site, spp) %>% 
    select(site, spp, duration, nYrs, slope, sig) %>% 
    data.frame()
  write.csv(S1, "Tables/S1.csv", row.names = FALSE)

  
# Table S2 ---------------------------------------------------------------------
  
	S2 <- readRDS("Data/dat.rda") %>% 
		group_by(Family, spp, spp.name, Depth, TL) %>% 
		summarise(Family = unique(Family),
							spp = unique(spp),
							spp.name = unique(spp.name),
							Depth = unique(Depth),
							TL = unique(TL),
							NoSeries = n_distinct(sppsite),
							NoYrs = n()) %>% 
		data.frame()
  write.csv(S2, "Tables/S2.csv", row.names = FALSE)
 

# Table S3 ---------------------------------------------------------------------

  S3 <- readRDS("Data/dat.rda") %>% 
    group_by(site) %>%
    summarise(Lat = min(lat),
              Lon = min(lon),
              NoSpp = length(unique(spp)))
  write.csv(S3, "Tables/S3.csv", row.names = FALSE)
  
  
   
