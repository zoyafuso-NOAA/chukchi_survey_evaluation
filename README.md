This repository contains preliminary operating models for Arctic (Chukchi Sea) 
groundfishes. Currently, vector autoregressive spatiotemporal models using the
[VAST package](https://github.com/James-Thorson-NOAA/VAST).

Contributors: Zack Oyafuso, Lewis Barnett, and Stan Kotwicki

## Versions
R version 4.0.2 (2020-06-22) was used for this analysis. Here are some relevant
package versions used:

| Package Name    | Version        |
|-----------------|----------------|
| VAST            | 3.6.1          | 
| TMB             | 1.7.19         |   
| FishStatsUtils  | 2.8.0          | 
| Matrix          | 1.2-18         | 
| INLA            | 21.02.23       | 
| dplyr           | 1.0.5          | 
| readxl          | 1.3.1          | 
| reshape         | 0.8.8          | 
| raster          | 3.4.5          | 
| sp              | 1.4.5          | 
| rgdal           | 1.5.23         | 
|                 |            | 

## Species List

## Model Fitting Settings
Four versions of random field configurations were conducted for each 
single-species run.

| Spatial 1st Pred (Omega_1)| Spatial 2nd Pred (Omega_2)| Spatiotemporal 1st Pred (Epsilon_1)| Spatiotemporal 2nd Pred (Epsilon_2)|
|---------------------------|---------------------------|------------------------------------|------------------------------------|
| X                         | X                         |                                    |                                    | 
| X                         | X                         | X                                  |                                    | 
| X                         | X                         |                                    | X                                  | 
| X                         | X                         | X                                  | X                                  | 

200 knots were used

