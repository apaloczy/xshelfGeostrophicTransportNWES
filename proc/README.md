proc
----

The python scripts in this directory process the different files to calculate derived fields.

* `00-interpADT_shfbrk.py`: Extracts the shelf edge isobath (200 m) from the bathymetry file, downloads and interpolates the Absolute Dynamic Topography to the isobath.
* `01-interpWOA_shfbrk.py`: Interpolates the WOA23 T/S fields to the 200 m isobath
* `02-binavgEN4_alongisobath.py`: Bin-averages the EN4 profiles along the 200 m isobath.

The file `segment_latitudes_wei_etal2024.txt` contains the (longitude, latitude) coordinates used to define the segments of the 200 m isobath (AS, CSS, PSBS, MS, HS, WSS, NENSS).
