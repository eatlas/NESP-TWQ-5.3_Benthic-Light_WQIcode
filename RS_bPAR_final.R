library(ncdf4) 
library(ggplot2)
library(dplyr)
library(sp) # Spatial Polygons
library(lubridate)
library(stringi)
library(stringr)
library(gridExtra)
library(readxl)

# Load Murray's handy functions to calculate seasons, waterYears and grades
source("WQI_functions.R")

# Load the geographic polygons defining the regions
load("Polys.RData")
spatial = read.csv('spatial.csv', strip.white=TRUE)

# Location of river load data
#load_file <- 'C:\\Users\\brobson\\OneDrive - Australian Institute of Marine Science\\Documents\\NESP benthic light\\river loads from steve.xlsx'
#load_file <- '/Users/maritescanto/OneDrive/buhayphd/mydocs/thesis/Chapter5/newwqi/16mol/river loads from steve.xlsx'
load_file <- '2020-Loads-for-loading-maps.xlsx'

# Where will we find Tess' netcdf files?
#data_dir <- "C:\\Users\\brobson\\OneDrive - Australian Institute of Marine Science\\Documents\\NESP benthic light\\data\\"
data_dir <- "https://maps.eatlas.org.au/thredds/dodsC/NESP-TWQ-5-3_Benthic-light/xr_parb/orig/xr_parb_daily/"

#output/input directory
iodir <- '.'

# set working directoty
setwd(iodir)


# Which steps do we want to recalculate? (If you set any level to TRUE, subsequent levels are automatically also TRUE)
calculate_regional_maps = TRUE
calculate_reference_maps = FALSE + calculate_regional_maps
calculate_masks = FALSE + calculate_reference_maps
calculate_monthly_values = FALSE + calculate_masks
calculate_aggregates = FALSE + calculate_monthly_values
calculate_grades_and_generate_figures = FALSE + calculate_aggregates
calculate_regressions = FALSE + calculate_grades_and_generate_figures

# If the following parameters are changed, you will need to set calculate_masks and subsequent steps to TRUE
lit_threshold <- 16
seagrass_threshold = 10
coral_threshold = 6
KI = 1

if (calculate_regional_maps) {
  print('calculating regional maps')

  # Set up necessary variables and label pixels by region
  region_names <- names(Polys)
  for (i in 1:24) {
    region_names[i] <- as.character(spatial[which(spatial$GBRMPA_Zone==region_names[i]), "Zone"])
  }
  inputfile <- paste0(data_dir, 'xr_parb_daily_2002.nc')
  nc <- nc_open(inputfile)
  lon <- ncvar_get(nc, "lon")
  lat <- ncvar_get(nc, "lat")
  ll <- expand.grid(lon = lon, lat = lat)
  coordinates(ll) <- ~lon+lat
  pts <- over(ll, Polys)
  ll$region <- as.factor(region_names[pts])
  ll$region <- ordered(ll$region, levels=region_names)
  nc_close(nc)
  # Pixels that haven't been labelled with a region name are outside the GBRMP. Set up a mask,
  # "in_domain" which we can use later to avoid processing pixels outside the marine park.
  in_domain <- which(!is.na(ll$region))
  save(file='region_map.rda', list=c('ll', 'in_domain'))
} else {
  print('loading regional maps from previously saved file')
  load('region_map.rda')
}
if (calculate_reference_maps) {
  print('calculating reference maps')
  parb_95th_monthly <- NULL
  for (m in 1:12) {
    print(paste('month', m))
    parb_m <- NULL
    pb <- txtProgressBar(min = 0, max = 10, style = 1, width=as.integer(getOption("width")*2/3))
    # One decade is a nice round period of time to use when calculating the monthly 95th percentiles.
    # We could use the full time-series, but then we'd need to decide whether to recalculate this when
    # we add more years of data in future.
    for (y in 2003:2012) { 
      inputfile <- paste0(data_dir, 'xr_parb_daily_', y, '.nc')
      nc <- nc_open(inputfile)
      ds <- ncvar_get(nc, "time")
      ds_units <- ncatt_get(nc, 'time')
      orig <- stri_datetime_parse(ncatt_get(nc ,'time')$units, "'days since 'yyyy-MM-dd HH:mm:ss")[1]
      ds <- ds*86400 + orig
      dm <- which(month(ds) == m)
      parb_d <- ncvar_get(nc, "parb", start=c(1, 1, dm[1]), count=c(-1, -1, length(dm)))
      parb_d<- array(parb_d, dim = c(prod(dim(parb_d)[1:2]), dim(parb_d)[3]))[in_domain, ]
      parb_m <- cbind(parb_m, parb_d)
      nc_close(nc)
	    setTxtProgressBar(pb,y-2002)
    }
    cat('\n')
    parb_95th_m <- data.frame(reference_parb = apply(parb_m, 1, quantile, probs=0.95, na.rm=TRUE), month = m)
    #parb_95th_m <- data.frame(reference_parb = apply(parb_m, 1, quantile, probs=0.9, na.rm=TRUE), month = m)
    parb_95th_monthly <- rbind(parb_95th_monthly, parb_95th_m)
  }
  save(file='reference_maps.rda', list = 'parb_95th_monthly')
} else {
  print('loading previously calculated reference maps from file')
  load('reference_maps.rda')
}
if (calculate_masks) {
  # use 2006 as the reference year since that is the driest year in the record and prevents a shifting goalpost when the record gets longer
  print('calculating regional masks')
  # Find 95th percentile for 2006 (a very dry year)
  inputfile <- paste0(data_dir, 'xr_parb_daily_2006.nc')
  nc <- nc_open(inputfile)
  parb_2016 <- ncvar_get(nc, "parb")
  nc_close(nc)
  parb_2016[parb_2016>29] <- NA
  #parb_2016[parb_2016>lit_threshold] <- lit_threshold
  parb_2016 <- array(parb_2016, dim = c(prod(dim(parb_2016)[1:2]), dim(parb_2016)[3]))[in_domain, ]
  parb_95th <- apply(parb_2016, 1, quantile, probs=0.95, na.rm=TRUE)
  parb_50th <- apply(parb_2016, 1, quantile, probs=0.5, na.rm=TRUE)
  # These masks need to be applied AFTER in_domain has been applied to the data
  mask <- which(parb_95th >= lit_threshold)
  coral_mask <- which(parb_50th >= coral_threshold)
  seagrass_mask <- which(parb_50th >= seagrass_threshold)

  save(file='masks.rda', list=c('parb_95th', 'mask', 'll', 'seagrass_mask', 'coral_mask', 'in_domain'))
} else {
  print('loading masks from previously saved file')
  load('masks.rda')
}

if (calculate_monthly_values) {
  print('Applying mask')
  ll <- ll[in_domain, ]

  reference_parb <- parb_95th
  reference_parb[reference_parb > lit_threshold] <- lit_threshold

  print('calculating monthly values')
  parb_m <- NULL
  parb_m_seagrass <- NULL
  parb_m_coral <- NULL

  for (y in 2002:2019) {
    print(y)
    inputfile <- paste0(data_dir, 'xr_parb_daily_', y, '.nc')
    nc <- nc_open(inputfile)
    ds <- ncvar_get(nc, 'time')
    ds_units <- ncatt_get(nc, 'time')
    orig <- stri_datetime_parse(ncatt_get(nc ,'time')$units, "'days since 'yyyy-MM-dd HH:mm:ss")[1]
    ds <- ds*86400 + orig

    pb <- txtProgressBar(min = 0, max = length(ds), style = 1, width=as.integer(getOption("width")*2/3))

    for (i in 1:length(ds)) {
      d <- ds[i]
      parb_d = c(ncvar_get(nc, "parb", start=c(1,1,i), count=c(-1,-1,1)))[in_domain]

      # filter out bad values
      parb_d[parb_d>29] <- NA
      valid_pixels <- !is.na(parb_d)

      # Set NAs to zero so they don't add to the cumulative sum
      parb_d[is.na(parb_d)] <- 0
      parb_capped_d <- parb_d
      parb_capped_d[parb_capped_d > lit_threshold] <- lit_threshold

      if ((mday(d) == 1)|((y==2002)&(i==1))) { # It's the first day of the month
        valid_days <- as.integer(valid_pixels)
        valid_pixels <- which(valid_pixels)
        parb_cumulative <- parb_d
        parb_capped_cumulative <- parb_capped_d
      } else {
        valid_days <- valid_days + as.integer(valid_pixels)
        if (length(valid_days)!=(length(valid_pixels))) browser
        valid_pixels <- which(valid_pixels)
        parb_cumulative[valid_pixels] <- parb_cumulative[valid_pixels] + parb_d[valid_pixels]
        parb_capped_cumulative[valid_pixels] <- parb_capped_cumulative[valid_pixels] + parb_capped_d[valid_pixels]

        daysinmonth <- days_in_month(d)
        if (daysinmonth == mday(d)) { # It's the last day of the month
          season <- WQ_season(d)
          waterYear <- WQ_waterYear(d)
          m <- month(d)

          # Calculate monthly means (considering only valid pixels) and total missing photons
          # We read one day at a time to avoid using too much memory.
          parb <- parb_cumulative/valid_days
          parb_capped <- parb_capped_cumulative/valid_days
          reference_parb <- parb_95th_monthly$reference_parb[parb_95th_monthly$month == m]

          reference_parb[reference_parb > lit_threshold] <- lit_threshold
          missing_photons <- daysinmonth * (reference_parb - parb_capped)
          missing_photons[missing_photons<0] <- 0
          worst_case <- reference_parb * daysinmonth

          # Calculate final values for the current month.
          parb_i <- data.frame(region = ll$region[mask],
                               parb = parb[mask],
                               parb_capped = parb_capped[mask],
                               missing_photons = missing_photons[mask],
                               stress1 = missing_photons[mask] / worst_case[mask],
                               lat = ll$lat[mask],
                               lon = ll$lon[mask],
                               mon = m,
                               season = season,
                               waterYear = waterYear,
                               daysinmonth = daysinmonth)

          parb_i_seagrass <- data.frame(region = ll$region[seagrass_mask],
                               parb = parb[seagrass_mask],
                               parb_capped = parb_capped[seagrass_mask],
                               missing_photons = missing_photons[seagrass_mask],
                               stress1 = missing_photons[seagrass_mask] / worst_case[seagrass_mask],
                               lat = ll$lat[seagrass_mask],
                               lon = ll$lon[seagrass_mask],
                               mon = m,
                               season = season,
                               waterYear = waterYear,
                               daysinmonth = daysinmonth)

          parb_i_coral <- data.frame(region = ll$region[coral_mask],
                               parb = parb[coral_mask],
                               parb_capped = parb_capped[coral_mask],
                               missing_photons = missing_photons[coral_mask],
                               stress1 = missing_photons[coral_mask] / worst_case[coral_mask],
                               lat = ll$lat[coral_mask],
                               lon = ll$lon[coral_mask],
                               mon = m,
                               season = season,
                               waterYear = waterYear,
                               daysinmonth = daysinmonth)

          # Add current month to record of all months
          # Note that growing variables inside a loop is slow. It would be faster if the variables were initialised at their full size before starting,
          # and values added with inidices instead of using rbind. This would also make it possible to parallelise this loop.
          parb_m <- rbind(parb_m, parb_i)
          parb_m_seagrass <- rbind(parb_m_seagrass, parb_i_seagrass)
          parb_m_coral <- rbind(parb_m_coral, parb_i_coral)
        }
      }
	    setTxtProgressBar(pb,i)
    }
    cat('\n')
    nc_close(nc)
  }
  close(pb)
  save(file='monthly_parb.rda', list=c('parb_m', 'parb_m_seagrass', 'parb_m_coral', 'worst_case'))
} else {
  print('Loading monthly mean parb and parb_capped from file')
  load('monthly_parb.rda')
}

if (calculate_aggregates) {
  print('Calculating annual and seasonal values')

  # dplyr does all the work for us here. If we had unlimited memory, we could have used the same approach to calculate monthly
  # means and 95th percentiles from daily values.
  annual <- parb_m %>% group_by(region, waterYear) %>%
                       summarise(parb = mean(parb, na.rm=TRUE), 
                                 parb_capped = mean(parb_capped, na.rm=TRUE), 
                                 missing_photons = sum(missing_photons, na.rm = TRUE),
                                 stress1 = mean(stress1, na.rm = TRUE),
                                 .groups="drop")
  seasonal <- parb_m %>% group_by(region, waterYear, season) %>%
                         summarise(parb = mean(parb, na.rm=TRUE), 
                                   parb_capped = mean(parb_capped, na.rm=TRUE),
                                   missing_photons = sum(missing_photons, na.rm = TRUE),
                                   stress1 = mean(stress1, na.rm = TRUE),
                                   .groups="drop")
  
  annual_seagrass <- parb_m_seagrass %>% group_by(region, waterYear) %>%
                                 summarise(parb = mean(parb, na.rm=TRUE), 
                                           parb_capped = mean(parb_capped, na.rm=TRUE), 
                                           missing_photons = sum(missing_photons, na.rm = TRUE),
                                           stress1 = mean(stress1, na.rm = TRUE),
                                           .groups="drop")
  seasonal_seagrass <- parb_m_seagrass %>% group_by(region, waterYear, season) %>%
                                   summarise(parb = mean(parb, na.rm=TRUE), 
                                             parb_capped = mean(parb_capped, na.rm=TRUE),
                                             missing_photons = sum(missing_photons, na.rm = TRUE),
                                             stress1 = mean(stress1, na.rm = TRUE),
                                             .groups="drop") 
  
  annual_coral <- parb_m_coral %>% group_by(region, waterYear) %>%
                              summarise(parb = mean(parb, na.rm=TRUE), 
                                        parb_capped = mean(parb_capped, na.rm=TRUE), 
                                        missing_photons = sum(missing_photons, na.rm = TRUE),
                                        stress1 = mean(stress1, na.rm = TRUE),
                                        .groups="drop")
  seasonal_coral <- parb_m_coral %>% group_by(region, waterYear, season) %>%
                                summarise(parb = mean(parb, na.rm=TRUE), 
                                          parb_capped = mean(parb_capped, na.rm=TRUE),
                                          missing_photons = sum(missing_photons, na.rm = TRUE),
                                          stress1 = mean(stress1, na.rm = TRUE),
                                          .groups="drop") 

  save(file='annual_parb.rda', list=ls())
} else {
  print('loading annual and seasonal values and WQIs from file.')
  load('annual_parb.rda')
  #load('annual_parb (old).rda')
}
if (calculate_grades_and_generate_figures) {
  print('Calculating grades and generating time-series figures')
  # Cycle through to calculate and plot several versions of the WQI.
  # Save the one we want to use for the regressions vs river loads for last -- this is the only one that will be saved at the end.

  max_stress = 1
  annual <- annual %>% group_by(region) %>% mutate(min_stress = min(stress1, na.rm=TRUE)) %>% ungroup() %>% mutate(max_stress = max_stress)
  annual_coral <- annual_coral %>% group_by(region) %>% mutate(min_stress = min(stress1, na.rm=TRUE)) %>% ungroup() %>% mutate(max_stress = max_stress)
  annual_seagrass <- annual_seagrass %>% group_by(region) %>% mutate(min_stress = min(stress1, na.rm=TRUE)) %>% ungroup() %>% mutate(max_stress = max_stress)

  seasonal <- seasonal %>% group_by(region) %>% mutate(min_stress = min(stress1, na.rm=TRUE)) %>% ungroup() %>% mutate(max_stress = max_stress)
  seasonal_coral <- seasonal_coral %>% group_by(region) %>% mutate(min_stress = min(stress1, na.rm=TRUE)) %>% ungroup() %>% mutate(max_stress = max_stress)
  seasonal_seagrass <- seasonal_seagrass %>% group_by(region) %>% mutate(min_stress = min(stress1, na.rm=TRUE)) %>% ungroup() %>% mutate(max_stress = max_stress)

  for (i in c(1, 3:9, 2)) {
    if (i==1) {
      mydf <- annual %>% filter((waterYear>2002)&(waterYear<2020))
      WQI_label_a <- "Annual WQI calculated over 95th percentile mask"
      WQI_label <- expression("Annual WQI"[bPAR]*" calculated over 95"^{th}*" percentile mask")
    } else if (i==2) {
      mydf <- annual_coral %>% filter((waterYear>2002)&(waterYear<2020))
      WQI_label_a <- "Annual WQI calculated over coral habitat mask"
      WQI_label <- expression("Annual WQI"[bPAR]*" calculated over coral habitat mask")
    } else if (i==3) {
      mydf <- annual_seagrass %>% filter((waterYear>2002)&(waterYear<2020))
      WQI_label_a <- "Annual WQI calculated over seagrass habitat mask"
      WQI_label <- expression("Annual WQI"[bPAR]*" calculated over seagrass habitat mask")
    } else if (i==4) {
      mydf <- seasonal %>% filter((waterYear>2002)&(waterYear<2020)&(season=="Wet"))
      WQI_label_a <- "Wet season WQI calculated over 95th percentile mask"
      WQI_label <- expression("Wet season WQI"[bPAR]*" calculated over 95"^{th}*" percentile mask")
    } else if (i==5) {
      mydf <- seasonal %>% filter((waterYear>2002)&(waterYear<2020)&(season=="Dry"))
      WQI_label_a <- "Dry season WQI calculated over 95th percentile mask"
      WQI_label <- expression("Dry season WQI"[bPAR]*" calculated over 95"^{th}*" percentile mask")
    } else if (i==6) {
      mydf <- seasonal_coral %>% filter((waterYear>2002)&(waterYear<2020)&(season=="Wet"))
      WQI_label_a <- "Wet season WQI calculated over coral habitat mask"
      WQI_label <- expression("Wet season WQI"[bPAR]*" calculated over coral habitat mask")
    } else if (i==7) {
      mydf <- seasonal_coral %>% filter((waterYear>2002)&(waterYear<2020)&(season=="Dry"))
      WQI_label_a <- "Dry season WQI calculated over coral habitat mask"
      WQI_label <- expression("Dry season WQI"[bPAR]*" calculated over coral habitat mask")
    } else if (i==8) {
      mydf <- seasonal_seagrass %>% filter((waterYear>2002)&(waterYear<2020)&(season=="Wet"))
      WQI_label_a <- "Wet season WQI calculated over seagrass habitat mask"
      WQI_label <- expression("Wet season WQI"[bPAR]*" calculated over seagrass habitat mask")
    } else if (i==9) {
      mydf <- seasonal_seagrass %>% filter((waterYear>2002)&(waterYear<2020)&(season=="Dry"))
      WQI_label_a <- "Dry season WQI calculated over seagrass habitat mask"
      WQI_label <- expression("Dry season WQI"[bPAR]*" calculated over seagrass habitat mask")
    }
                           
    mydf$stress <- mydf$stress1
    # If we wanted a stricter score, we'd set the minimum score to zero, but instead we use the minimum light stress for the region as the
    # benchmark for a perfect score. This means that even naturally turbid coastal areas can get an A if they achieve their best possible WQ, 
    # and so can deep offshore regions where even a little material in the water will reduce light at the bottom.
    #mydf$min_stress <- 0

    {
    mydf <- mydf %>% mutate(WaterBody=case_when(stringr:::str_detect(region, 'Enc') ~ 'Enclosed Coastal',
                                    stringr:::str_detect(region, 'Open') ~ 'Open Coastal',
                                    stringr:::str_detect(region, 'Midshelf') ~ 'Midshelf',
                                    stringr:::str_detect(region, 'Offshore') ~ 'Offshore'),
                     Region=case_when(stringr:::str_detect(region, 'Cape') ~ 'Cape York',
                                 stringr:::str_detect(region, 'Wet') ~ 'Wet Tropics',
                                 stringr:::str_detect(region, 'Dry') ~ 'Dry Tropics',
                                 stringr:::str_detect(region, 'Mackay') ~ 'Mackay Whitsunday',
                                 stringr:::str_detect(region, 'Fitzroy') ~ 'Fitzroy',
                                 stringr:::str_detect(region, 'Burnett') ~ 'Burnett Mary')) %>%
               mutate(Region=factor(Region, levels=unique(spatial$Region)),
                      WaterBody=factor(WaterBody, levels=unique(spatial$WaterBody))) %>%
               group_by(region) %>%
               mutate(WQI = scales::squish(scales::rescale(stress, from=c(min_stress[1], max_stress[1]), to = c(1, 0)))) %>%
               ungroup() %>%
               mutate(grade = WQ_generateGrades(WQI))

    GradeType='Uniform'
    gradeMids = WQ_gradeMids(type=GradeType)
    gradeBoundaries = WQ_gradeBoundaries(type=GradeType)

    g=ggplot(mydf, aes(y=WQI, x=waterYear)) +
        annotate(geom='rect', xmin=-Inf,xmax=Inf,ymin=gradeBoundaries[2],ymax=gradeBoundaries[1], fill='#00734D10') +
        annotate(geom='rect', xmin=-Inf,xmax=Inf,ymin=gradeBoundaries[3],ymax=gradeBoundaries[2], fill='#B0D23510') + 
        annotate(geom='rect', xmin=-Inf,xmax=Inf,ymin=gradeBoundaries[4],ymax=gradeBoundaries[3], fill='#F0C91810') +
        annotate(geom='rect', xmin=-Inf,xmax=Inf,ymin=gradeBoundaries[5],ymax=gradeBoundaries[4], fill='#F4772110') +
        annotate(geom='rect', xmin=-Inf,xmax=Inf,ymin=gradeBoundaries[6],ymax=gradeBoundaries[5], fill='#ED1C2410') +
        geom_line() +
        geom_point(aes(fill=grade), shape=21, color='black', size=3) +
        scale_fill_manual('Grade\n', breaks=LETTERS[1:5], labels=LETTERS[1:5],values=c('#00734D','#B0D235','#F0C918','#F47721','#ED1C24'), limits=LETTERS[1:5]) + 
        scale_x_continuous('Water year') +
        scale_y_continuous(expression("WQI"[bPAR])) +
        facet_grid(Region~WaterBody, scales='free_y') + theme_bw()+
        theme(strip.background=element_blank(), legend.key=element_blank(),
                                        #strip.text.y=element_text(angle=0, hjust=0)
              axis.title.y=element_text(margin = margin(0,0.5,0,0,'cm'))) +
        ggtitle(WQI_label)
    print(g)
    }
    ggsave(file=paste0(WQI_label_a, '.jpg'), g, width=10, height=10, units='in',dpi=300) 
  }
  save(file='calculated_grades.rda', mydf)
} else {
  print('loading calculated grades from a file')
  load('calculated_grades.rda')
}
if (calculate_regressions) {
  river_names <- as.character(t(read_excel(load_file, range="A2:A36", col_names=FALSE)))
  river_names <- stri_replace_all(river_names, '_', regex='  *')
  river_names <- stri_replace_all(river_names, '_', fixed='-')

  annual_discharge <- as.data.frame(t(read_excel(load_file, sheet="Discharge", col_names = FALSE, col_type = c("skip", rep("numeric", 17)))))
  #names(annual_discharge) <- c('waterYear', paste(river_names, "Discharge", sep="_"))
  names(annual_discharge) <- c('waterYear', river_names)

  annual_DIN <- as.data.frame(t(read_excel(load_file, sheet="DIN", col_names = FALSE, col_type = c("skip", rep("numeric", 17)))))
  names(annual_DIN) <- c('waterYear', river_names)

  annual_TSS <- as.data.frame(t(read_excel(load_file, sheet="TSS", col_names = FALSE, col_type = c("skip", rep("numeric", 17)))))
  names(annual_TSS) <- c('waterYear', river_names)

  annual_PN <- as.data.frame(t(read_excel(load_file, sheet="PN", col_names = FALSE, col_type = c("skip", rep("numeric", 17)))))
  names(annual_PN) <- c('waterYear', river_names)

  # These loads are as calculated in the 2016 Fabricius et al. paper (except we don't separate North and South Wet Tropics and handle the Fitzroy slightly differently):
  annual_loads <- data.frame(waterYear = annual_discharge$waterYear,
                             Cape_York_discharge = annual_discharge$Normanby_R + annual_discharge$Endeavour_River + annual_discharge$Stewart_River,
                             Cape_York_DIN = annual_DIN$Normanby_R + annual_DIN$Endeavour_River + annual_DIN$Stewart_River,
                             Cape_York_TSS = annual_TSS$Normanby_R + annual_TSS$Endeavour_River + annual_TSS$Stewart_River,
                             Cape_York_PN = annual_PN$Normanby_R + annual_PN$Endeavour_River + annual_PN$Stewart_River,
                             NWT_discharge = annual_discharge$Daintree + annual_discharge$Barron + annual_discharge$Russell_Mulgrave + annual_discharge$Johnstone + 0.3*annual_discharge$Burdekin,
                             NWT_DIN = annual_DIN$Daintree + annual_DIN$Barron + annual_DIN$Russell_Mulgrave + annual_DIN$Johnstone + 0.3*annual_DIN$Burdekin,
                             NWT_TSS = annual_TSS$Daintree + annual_TSS$Barron + annual_TSS$Russell_Mulgrave + annual_TSS$Johnstone + 0.3*annual_TSS$Burdekin,
                             NWT_PN = annual_PN$Daintree + annual_PN$Barron + annual_PN$Russell_Mulgrave + annual_PN$Johnstone + 0.3*annual_PN$Burdekin,
                             SWT_discharge = annual_discharge$Russell_Mulgrave + annual_discharge$Johnstone + 0.5*annual_discharge$Burdekin + annual_discharge$Tully + annual_discharge$Herbert,
                             SWT_DIN = annual_DIN$Russell_Mulgrave + annual_DIN$Johnstone + 0.5*annual_DIN$Burdekin + annual_DIN$Tully + annual_DIN$Herbert,
                             SWT_TSS = annual_TSS$Russell_Mulgrave + annual_TSS$Johnstone + 0.5*annual_TSS$Burdekin + annual_TSS$Tully + annual_TSS$Herbert,
                             SWT_PN = annual_PN$Barron + annual_PN$Russell_Mulgrave + annual_PN$Johnstone + 0.5*annual_PN$Burdekin + annual_PN$Tully + annual_PN$Herbert,
                             Burdekin_discharge = annual_discharge$Burdekin,
                             Burdekin_DIN = annual_DIN$Burdekin,
                             Burdekin_TSS = annual_TSS$Burdekin,
                             Burdekin_PN = annual_PN$Burdekin,
                             Whitsundays_discharge = annual_discharge$Proserpine + annual_discharge$"O'Connell" + annual_discharge$Pioneer,
                             Whitsundays_DIN = annual_DIN$Proserpine + annual_DIN$"O'Connell" + annual_DIN$Pioneer,
                             Whitsundays_TSS = annual_TSS$Proserpine + annual_TSS$"O'Connell" + annual_TSS$Pioneer,
                             Whitsundays_PN = annual_PN$Proserpine + annual_PN$"O'Connell" + annual_PN$Pioneer,
                             Fitzroy_discharge = annual_discharge$Fitzroy,
                             Fitzroy_DIN = annual_DIN$Fitzroy,
                             Fitzroy_TSS = annual_TSS$Fitzroy,
                             Fitzroy_PN = annual_PN$Fitzroy,
                             FitzBurnett_discharge = annual_discharge$Fitzroy + annual_discharge$Burnett_River,
                             FitzBurnett_DIN = annual_DIN$Fitzroy + annual_DIN$Burnett_River,
                             FitzBurnett_TSS = annual_TSS$Fitzroy + annual_DIN$Burnett_River,
                             FitzBurnett_PN = annual_PN$Fitzroy + annual_DIN$Burnett_River)

  mydf <- left_join(mydf, annual_loads, by="waterYear")

  mydf <- mydf %>% mutate(local_TSS = case_when(Region=="Fitzroy" ~ scales::rescale(Fitzroy_TSS),
                                            Region=="Dry Tropics" ~ scales::rescale(Burdekin_TSS),
                                            Region=="Mackay Whitsunday" ~ scales::rescale(Whitsundays_TSS),
                                            Region=="Wet Tropics" ~ scales::rescale(NWT_TSS + SWT_TSS),
                                            Region=="Burnett Mary" ~ scales::rescale(FitzBurnett_TSS),
                                            Region=="Cape York" ~ scales::rescale(Cape_York_TSS)),
                          local_discharge = case_when(Region=="Fitzroy" ~ scales::rescale(Fitzroy_discharge),
                                            Region=="Dry Tropics" ~ scales::rescale(Burdekin_discharge),
                                            Region=="Mackay Whitsunday" ~ scales::rescale(Whitsundays_discharge),
                                            Region=="Wet Tropics" ~ scales::rescale(NWT_discharge + SWT_discharge),
                                            Region=="Burnett Mary" ~ scales::rescale(FitzBurnett_discharge),
                                            Region=="Cape York" ~ scales::rescale(Cape_York_discharge)),
                          local_DIN = case_when(Region=="Fitzroy" ~ scales::rescale(Fitzroy_DIN),
                                            Region=="Dry Tropics" ~ scales::rescale(Burdekin_DIN),
                                            Region=="Mackay Whitsunday" ~ scales::rescale(Whitsundays_DIN),
                                            Region=="Wet Tropics" ~ scales::rescale(NWT_DIN + SWT_DIN),
                                            Region=="Burnett Mary" ~ scales::rescale(FitzBurnett_DIN),
                                            Region=="Cape York" ~ scales::rescale(Cape_York_DIN)),
                          local_PN = case_when(Region=="Fitzroy" ~ scales::rescale(Fitzroy_PN),
                                            Region=="Dry Tropics" ~ scales::rescale(Burdekin_PN),
                                            Region=="Mackay Whitsunday" ~ scales::rescale(Whitsundays_PN),
                                            Region=="Wet Tropics" ~ scales::rescale(NWT_PN + SWT_PN),
                                            Region=="Burnett Mary" ~ scales::rescale(FitzBurnett_PN),
                                            Region=="Cape York" ~ scales::rescale(Cape_York_PN))) %>%
                   select(-Fitzroy_TSS, -Fitzroy_DIN, -Fitzroy_discharge, -Fitzroy_PN,
                          -FitzBurnett_TSS, -FitzBurnett_DIN, -FitzBurnett_discharge, -FitzBurnett_PN,
                          -Burdekin_TSS, -Burdekin_DIN, -Burdekin_discharge, -Burdekin_PN,
                          -NWT_TSS, -NWT_DIN, -NWT_discharge, -NWT_PN,
                          -SWT_TSS, -SWT_DIN, -SWT_discharge, -SWT_PN,
                          -Cape_York_TSS, -Cape_York_DIN, -Cape_York_discharge, -Cape_York_PN,
                          -Whitsundays_TSS,-Whitsundays_DIN, -Whitsundays_discharge, -Whitsundays_PN)
  
  #g2 <- ggplot(mydf %>% filter(Region == "Fitzroy" | Region =="Dry Tropics" | Region == "Mackay Whitsunday" | Region=="Wet Tropics" | Region=="Burnett Mary"), 
  g2 <- ggplot(mydf,
               aes(y = WQI)) + 
               geom_point(aes(x=local_TSS), colour='red3') +
               geom_smooth(aes(x=local_TSS), method='lm', colour='red3', linetype='dashed', fill="red3", alpha=0.2) +
               geom_point(aes(x=local_discharge), colour="steelblue") +
               geom_smooth(aes(x=local_discharge), colour="steelblue", method='lm', fill="steelblue") +
               annotate(geom='rect', xmin=-Inf,xmax=Inf,ymin=gradeBoundaries[2],ymax=gradeBoundaries[1], fill='#00734D10') +
               annotate(geom='rect', xmin=-Inf,xmax=Inf,ymin=gradeBoundaries[3],ymax=gradeBoundaries[2], fill='#B0D23510') + 
               annotate(geom='rect', xmin=-Inf,xmax=Inf,ymin=gradeBoundaries[4],ymax=gradeBoundaries[3], fill='#F0C91810') +
               annotate(geom='rect', xmin=-Inf,xmax=Inf,ymin=gradeBoundaries[5],ymax=gradeBoundaries[4], fill='#F4772110') +
               annotate(geom='rect', xmin=-Inf,xmax=Inf,ymin=gradeBoundaries[6],ymax=gradeBoundaries[5], fill='#ED1C2410') +
               xlab('River load (proportion of maximum load)') +
               ylab(expression("WQI"[bPAR])) +
               facet_grid(Region~WaterBody, scales="free") + theme_bw() +
               ylim(0,1)
  plot(g2)
  
  #ggsave(file="WQI_vs_loads.png", g2)
  ggsave(file="WQI_vs_loads_edit_v3.png", g2, width=10, height=12, units='in',dpi=300) 
  #save(file='WQI_vs_loads_edit.rda', mydf)
  

  #Some color-blind friendly colors
  #cbPalette <- c("#abd9e9", "#2c7bb6", "#ffffbf","#fdae61", "#d7191c") #make burdekin dark blue
  g3 <- ggplot(mydf,
              aes(y = WQI)) + 
              geom_point(aes(x=local_TSS), colour='#d7191c') + #red
              geom_smooth(aes(x=local_TSS), method='lm', colour='#d7191c', linetype='dashed', fill="#d7191c",lwd=0.5) +
              geom_point(aes(x=local_DIN), colour='#4dac26') + #green
              geom_smooth(aes(x=local_DIN), method='lm', colour='#4dac26', linetype='dashed', fill="#4dac26",lwd=0.5) +
              geom_point(aes(x=local_discharge), colour="#2c7bb6") + #blue
              geom_smooth(aes(x=local_discharge), colour="#2c7bb6", method='lm', fill="#2c7bb6", alpha=0.1,lwd=0.5) +
              annotate(geom='rect', xmin=-Inf,xmax=Inf,ymin=gradeBoundaries[2],ymax=gradeBoundaries[1], fill='#00734D10') +
              annotate(geom='rect', xmin=-Inf,xmax=Inf,ymin=gradeBoundaries[3],ymax=gradeBoundaries[2], fill='#B0D23510') + 
              annotate(geom='rect', xmin=-Inf,xmax=Inf,ymin=gradeBoundaries[4],ymax=gradeBoundaries[3], fill='#F0C91810') +
              annotate(geom='rect', xmin=-Inf,xmax=Inf,ymin=gradeBoundaries[5],ymax=gradeBoundaries[4], fill='#F4772110') +
              annotate(geom='rect', xmin=-Inf,xmax=Inf,ymin=gradeBoundaries[6],ymax=gradeBoundaries[5], fill='#ED1C2410') +
              xlab('River load (proportion of maximum load)') +
              ylab(expression("WQI"[bPAR])) +
              facet_grid(Region~WaterBody, scales="free") + theme_bw() +
              ylim(0,1)
            plot(g3)
  
  #ggsave(file="WQI_vs_loads.png", g2)
  ggsave(file="WQI_vs_TSSDINdischarge_edit_v2.png", g3, width=10, height=12, units='in',dpi=300) 
  
  
  #TO get the linear regression coefficients
  fitted_models_TSS = mydf %>% group_by(Region,WaterBody) %>% do(model = lm(WQI ~ local_TSS, data = .)) 
  fitted_models_DIN = mydf %>% group_by(Region,WaterBody) %>% do(model = lm(WQI ~ local_DIN, data = .))  
  fitted_models_Discharge = mydf %>% group_by(Region,WaterBody) %>% do(model = lm(WQI ~ local_discharge, data = .))  
  
  #Tidy up the results to get coefficients and r-squared values
  library(broom)
  TSS_fit<-fitted_models_TSS %>% mutate(rsquared = summary(model)$r.squared) %>% 
                                 mutate(intercept = summary(model)$coefficients[1]) %>% 
                                 mutate(slope = summary(model)$coefficients[2])  
    
  DIN_fit<-fitted_models_DIN %>% mutate(rsquared = summary(model)$r.squared) %>% 
                                 mutate(intercept = summary(model)$coefficients[1]) %>% 
                                 mutate(slope = summary(model)$coefficients[2])  
  
  Discharge_fit<-fitted_models_Discharge %>% mutate(rsquared = summary(model)$r.squared) %>% 
                                             mutate(intercept = summary(model)$coefficients[1]) %>% 
                                             mutate(slope = summary(model)$coefficients[2])  
  
  write.csv(apply(TSS_fit,2,as.character),"tss_fit_TullyHerbert.csv")  
  write.csv(apply(DIN_fit,2,as.character),"din_fit_TullyHerbert.csv")  
  write.csv(apply(Discharge_fit,2,as.character),"discharge_fit_TullyHerbert.csv")  
}
