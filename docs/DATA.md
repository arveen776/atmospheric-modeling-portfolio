# Data

Raw data is not stored in this repository because of its size. Figures and summary results are included, so the work can be reviewed without it.

| Project | Data | Approx. size | Source | Expected location |
|---------|------|--------------|------|-------------------|
| 01 | WRF `wrfout_d02_*`, YSU and MYNN2.5 runs | 5.4 GB | Course simulations | `data/YSU/`, `data/MYNN/` |
| 01 | `geo_em.d01.nc`, `geo_em.d02.nc` | 5 MB | WPS geogrid | `data/` |
| 01 | `era5_100m_winds.grib` (100 m u and v, hourly) | 3 MB | [ERA5 single levels](https://cds.climate.copernicus.eu/) | `data/` |
| 02 | WRF `wrfout_d01_*` with turbines, 48 hourly files | See note | Course simulations | `data/wfp/` |
| 02 | WRF `wrfout_d01_*` without turbines (MYNN), 48 hourly files | See note | Provided by the instructor | `data/nwf/` |
| 03 | `azmth_d01.nc`, `azmth_d02.nc` for five storm cases | 1 GB | Provided in the course | `notebooks/data/IdealizedSimulations/` |
| 04 | WRF `wrfout_d01_*`, `wrfout_d02_*`, 49 hourly files each | 2 GB | Course simulations | `data/` |

**Note:** the Project 02 folders were about 13 GB and 17 GB in the original coursework, including other files and archives.

Small inputs that the code needs are included: turbine locations and power tables are in `projects/02_wind_farm_wakes_awaken/inputs/`.

Projects 01 and 04 can also read data from elsewhere through the `DUNKELFLAUTE_DATA` and `KATRINA_DATA` environment variables. Project 02 scripts accept `--wfp-dir` and `--nwf-dir`.

ERA5 citation: Hersbach et al. (2020), *The ERA5 global reanalysis*, Quarterly Journal of the Royal Meteorological Society, 146, 1999–2049.
