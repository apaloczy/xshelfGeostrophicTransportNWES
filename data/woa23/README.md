WOA23
-----

The World Ocean Atlas 2023 (WOA23) climatology files should be placed in this directory.

The python scripts in this directory download and subset the relevant WOA23 fields:

* `00-subsetWOA23-TS.py`: Subsets the annual/winter/spring/summer/autumn temperature and salinity fields.
* `01-subsetWOA23-nitrate.py`: Subsets the annual/winter/spring/summer/autumn nitrate fields.
* `02-derive_densityWOA23.py`: Derives density fields from the T/S fields.
