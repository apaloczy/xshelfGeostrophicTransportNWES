NSBC
----

The North Sea Biogeochemical Climatology (NSBC) files should be placed in this directory. The climatology can be downloaded from [https://www.cen
.uni-hamburg.de/en/icdc/data/ocean/nsbc.html](https://www.cen
.uni-hamburg.de/en/icdc/data/ocean/nsbc.html), archived at [https://doi.org/10.1594/WDCC/NSBClim_v1.1](https://doi.org/10.1594/WDCC/NSBClim_v1.1)

The subdirectory `level2/` contains the profile data, while the subdirectory `level3/` contains the objectively-mapped (gridded) fields.

The python script `../../proc/03-calc_Nshelf-NSBC.py` processes the NSBC level2 profiles in directory `level2/` to derive the shelf area-averaged nitrate concentrations.
