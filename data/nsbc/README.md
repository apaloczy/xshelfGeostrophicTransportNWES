NSBC
----

The North Sea Biogeochemical Climatology (NSBC) files should be placed in this directory. The climatology can be downloaded from [https://www.cen.uni-hamburg.de/en/icdc/data/ocean/nsbc.html](https://www.cen.uni-hamburg.de/en/icdc/data/ocean/nsbc.html), archived at [https://doi.org/10.1594/WDCC/NSBClim_v1.1](https://doi.org/10.1594/WDCC/NSBClim_v1.1)

The subdirectory `level2/climatological_monthly_mean` contains the profile data, while the subdirectory `level3/all_data_mean/` contains the objectively-mapped (gridded) fields.

After downloading the files, this directory should contain four files:

* `level2/climatological_monthly_mean/NSBC_Level2_nitrate__UHAM_ICDC__v1.1__0.25x0.25deg__1960_2014.nc`
* `level2/climatological_monthly_mean/NSBC_Level2_salinity__UHAM_ICDC__v1.1__0.25x0.25deg__1960_2014.nc`
* `level3/all_data_mean/NSBC_Level3_nitrate__UHAM_ICDC__v1.1__0.25x0.25deg__OAN__all_data_mean__1960_2014.nc`
* `NSBC_Level2_land_sea_mask__UHAM_ICDC__v1.1__0.25x0.25deg.nc`

The python script `../../proc/03-calc_Nshelf-NSBC.py` processes the NSBC level2 profiles in `level2/climatological_monthly_mean/` to derive the shelf area-averaged nitrate concentrations.
