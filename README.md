# NESP-TWQ-5.3_Benthic-Light_WQIcode

This dataset contains code used to generate the daily benthic light (bPAR) data product provided at https://eatlas.org.au/data/uuid/356e7b3c-1508-432e-9d85-263ec8a67cef and the bPAR index for water quality in the Great Barrier Reef (GBR). It can also be used to calculate photosynthetically active radiation (PAR) at any specified depth in GBR waters.

The amount of light available for photosynthesis (photosynthetically active radiation, or PAR) is an important determinant of ecosystem health. PAR reaching the bottom of the water column is known as benthic PAR (bPAR). Where there is sufficient light reaching the bottom, seagrasses and corals may thrive. bPAR varies seasonally as a function of surface PAR, but also varies function of both water depth and water quality.

More information on the project can be found at: https://eatlas.org.au/nesp-twq-5/benthic-light-5-3 and metadata record https://eatlas.org.au/data/uuid/61a4bac5-79d1-4c1f-9358-a7bb587e07df

This work is available under Creative Commons licence 3.0 Attribution (CC-BY) https://creativecommons.org/licenses/by/3.0/ Robson B, Magno-Canto M, McKinna L, Logan M, Fabricius K, Collier C, Garcia R (2020) Code used to generate benthic photosynthetically active radiation (bPAR) and the bPAR index as reported in (NESP 5.3, AIMS, JCU and Go2Q)


## Files
|Filename|Description
|---|---
|[get_bpar.c](get_bpar.c)|The functions, written in C, for calculating benthic PAR or PAR at any specified depth. This code is designed to be executed as part of [NASA’s SeaDAS software](https://seadas.gsfc.nasa.gov/).
|[RS_bPAR_final.R](RS_bPAR_final.R)|R script to calculate the bPAR index. This script can executed in R v3.4 or later.
|[WQI_functions.R](WQI_functions.R)|R helper functions written by Murray Logan and called from the [RS_bPAR_final.R](RS_bPAR_final.R) script.
|[2020-Loads-for-loading-maps.xlsx](data/2020-Loads-for-loading-maps.xlsx)|River discharge and load data provided by Dr Stephen Lewis (JCU), included in this repository with permission.
|[spatial.csv](data/spatial.csv)|Region/zone information used for labelling.
|[Polys.RData](data/)|Geographic polygons defining the regions/zones.

## Running the code

This code can be run within a [Docker container](https://www.docker.com/resources/what-container), eliminating the need to have an existing R development environment. The following instructions will create an R environment in a Docker container.

### Install Docker
This step is only required if Docker is not already installed on your computer. There are too many options to discuss here how to install Docker, so instead refer to the [Docker website](https://hub.docker.com/search?q=&type=edition&offering=community).

### Extract this repository from Github.
Utilise your preferred method to extract this repository from Github. This will vary depending on your development tools. As an example, this is done from Linux as follows:

```
$ cd <parent directory of preferred location>
$ git clone https://github.com/eatlas/NESP-TWQ-5.3_Benthic-Light_WQIcode.git
```

### Start the Docker container.

```
$ cd <local directory of repository, contains RB_bPAR_final.R file>
$ docker run -ti --rm -v `pwd`:/workdir -w /workdir r-base bash
```

### Prepare the Docker container

```
$ apt-get update
$ apt-get upgrade -y
$ apt-get install libnetcdf-dev -y
```

### Prepare the R environment

```
$ R
> install.packages('ncdf4')
> install.packages('ggplot2')
> install.packages('dplyr')
> install.packages('sp')
> install.packages('lubridate')
> install.packages('stringi')
> install.packages('stringr')
> install.packages('gridExtra')
> install.packages('readxl')
> install.packages('zoo')
> q()
Save workspace image? [y/n/c]: n
```

### Execute the script

```
Rscript RS_bPAR_final.R
```
