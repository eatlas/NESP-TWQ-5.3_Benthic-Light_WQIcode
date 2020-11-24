# NESP-TWQ-5.3_Benthic-Light_WQIcode

This dataset contains code used to generate the daily benthic light (bPAR) data product provided at https://eatlas.org.au/data/uuid/356e7b3c-1508-432e-9d85
-263ec8a67cef and the bPAR index for water quality in the Great Barrier Reef (GBR). It can also be used to calculate photosynthetically active radiation (PAR) at any specified depth in GBR waters.

More information on the project can be found at: https://eatlas.org.au/nesp-twq-5/benthic-light-5-3 and metadata record https://eatlas.org.au/data/uuid/61a4bac5-79d1-4c1f-9358-a7bb587e07df

This work is available under Creative Commons licence 3.0 Attribution (CC-BY) https://creativecommons.org/licenses/by/3.0/ Robson B, Magno-Canto M, McKinna L, Logan M, Fabricius K, Collier C, Garcia R (2020) Code used to generate benthic photosynthetically active radiation (bPAR) and the bPAR index as reported in (NESP 5.3, AIMS, JCU and Go2Q)

The amount of light available for photosynthesis (photosynthetically active radiation, or PAR) is an important determinant of ecosystem health. PAR reaching the bottom of the water column is known as benthic PAR (bPAR). Where there is sufficient light reaching the bottom, seagrasses and corals may thrive. bPAR varies seasonally as a function of surface PAR, but also varies function of both water depth and water quality.

The functions needed to calculate benthic PAR or PAR at any specified depth is provided as c code in the file “get_bpar.c”. This is designed to be executed as part of NASA’s SeaDAS software.

The script to calculate the bPAR index, “RS_bPAR_final.R” can be executed in the R programming language (v 3.4 or later). This script calls on functions defined in a second R script, “WQI_functions.R” and spatial data and region labels provided in the R data files “Polys.rda” and “spatial.csv”, which are also provided.
