# Download WOA23 subsets for the NWES area.
from xarray import open_dataset

#---
xmin, xmax = -16.5, 10.5
ymin, ymax = 44, 63

dll = 1
xmin -= dll
xmax += dll
ymin -= dll
ymax += dll

# 1 degree, annual 1955-2020ish averages.
urlN = "https://www.ncei.noaa.gov/thredds-ocean/dodsC/woa23/DATA/nitrate/netcdf/all/1.00/woa23_all_n00_01.nc"

# 1 degree, sesonal averages.
urlNwinter = "https://www.ncei.noaa.gov/thredds-ocean/dodsC/woa23/DATA/nitrate/netcdf/all/1.00/woa23_all_n13_01.nc"
urlNspring = "https://www.ncei.noaa.gov/thredds-ocean/dodsC/woa23/DATA/nitrate/netcdf/all/1.00/woa23_all_n14_01.nc"
urlNsummer = "https://www.ncei.noaa.gov/thredds-ocean/dodsC/woa23/DATA/nitrate/netcdf/all/1.00/woa23_all_n15_01.nc"
urlNautumn = "https://www.ncei.noaa.gov/thredds-ocean/dodsC/woa23/DATA/nitrate/netcdf/all/1.00/woa23_all_n16_01.nc"

bbox = dict(lon=slice(xmin, xmax), lat=slice(ymin, ymax))

print("Downloading nitrate seasonal climatologies.")
open_dataset(urlN, decode_times=False).squeeze().sel(bbox).to_netcdf("woa23N-annual.nc")
open_dataset(urlNwinter, decode_times=False).squeeze().sel(bbox).to_netcdf("woa23N-winter.nc")
open_dataset(urlNspring, decode_times=False).squeeze().sel(bbox).to_netcdf("woa23N-spring.nc")
open_dataset(urlNsummer, decode_times=False).squeeze().sel(bbox).to_netcdf("woa23N-summer.nc")
open_dataset(urlNautumn, decode_times=False).squeeze().sel(bbox).to_netcdf("woa23N-autumn.nc")
