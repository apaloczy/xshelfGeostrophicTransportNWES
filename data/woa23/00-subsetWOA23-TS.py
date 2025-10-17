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

# 1/4 degree, annual 1955-2020ish averages.
urlT = "https://www.ncei.noaa.gov/thredds-ocean/dodsC/woa23/DATA/temperature/netcdf/decav/0.25/woa23_decav_t00_04.nc"
urlS = "https://www.ncei.noaa.gov/thredds-ocean/dodsC/woa23/DATA/salinity/netcdf/decav/0.25/woa23_decav_s00_04.nc"

# 1/4 degree, sesonal averages.
urlTwinter = "https://www.ncei.noaa.gov/thredds-ocean/dodsC/woa23/DATA/temperature/netcdf/decav/0.25/woa23_decav_t13_04.nc"
urlTspring = "https://www.ncei.noaa.gov/thredds-ocean/dodsC/woa23/DATA/temperature/netcdf/decav/0.25/woa23_decav_t14_04.nc"
urlTsummer = "https://www.ncei.noaa.gov/thredds-ocean/dodsC/woa23/DATA/temperature/netcdf/decav/0.25/woa23_decav_t15_04.nc"
urlTautumn = "https://www.ncei.noaa.gov/thredds-ocean/dodsC/woa23/DATA/temperature/netcdf/decav/0.25/woa23_decav_t16_04.nc"

urlSwinter = "https://www.ncei.noaa.gov/thredds-ocean/dodsC/woa23/DATA/salinity/netcdf/decav/0.25/woa23_decav_s13_04.nc"
urlSspring = "https://www.ncei.noaa.gov/thredds-ocean/dodsC/woa23/DATA/salinity/netcdf/decav/0.25/woa23_decav_s14_04.nc"
urlSsummer = "https://www.ncei.noaa.gov/thredds-ocean/dodsC/woa23/DATA/salinity/netcdf/decav/0.25/woa23_decav_s15_04.nc"
urlSautumn = "https://www.ncei.noaa.gov/thredds-ocean/dodsC/woa23/DATA/salinity/netcdf/decav/0.25/woa23_decav_s16_04.nc"

bbox = dict(lon=slice(xmin, xmax), lat=slice(ymin, ymax))

print("Downloading T seasonal climatologies.")
open_dataset(urlT, decode_times=False).squeeze().sel(bbox).to_netcdf("woa23T-annual.nc")
open_dataset(urlTwinter, decode_times=False).squeeze().sel(bbox).to_netcdf("woa23T-winter.nc")
open_dataset(urlTspring, decode_times=False).squeeze().sel(bbox).to_netcdf("woa23T-spring.nc")
open_dataset(urlTsummer, decode_times=False).squeeze().sel(bbox).to_netcdf("woa23T-summer.nc")
open_dataset(urlTautumn, decode_times=False).squeeze().sel(bbox).to_netcdf("woa23T-autumn.nc")

print("Downloading S seasonal climatologies.")
open_dataset(urlS, decode_times=False).squeeze().sel(bbox).to_netcdf("woa23S-annual.nc")
open_dataset(urlSwinter, decode_times=False).squeeze().sel(bbox).to_netcdf("woa23S-winter.nc")
open_dataset(urlSspring, decode_times=False).squeeze().sel(bbox).to_netcdf("woa23S-spring.nc")
open_dataset(urlSsummer, decode_times=False).squeeze().sel(bbox).to_netcdf("woa23S-summer.nc")
open_dataset(urlSautumn, decode_times=False).squeeze().sel(bbox).to_netcdf("woa23S-autumn.nc")

# 1/4 degree, decadal averages (75-84, 85-94, 95-04, 05-14, 15-22).
decs = ["7584", "8594", "95A4", "A5B4", "B5C2"]
decs_fout = {"7584":"1975-1984", "8594":"1985-1994", "95A4":"1995-2004", "A5B4":"2005-2014", "B5C2":"2015-2022"}
print("Downloading T/S decadal climatologies.")
for dec in decs:
    urlT_dec = "https://www.ncei.noaa.gov/thredds-ocean/dodsC/woa23/DATA/temperature/netcdf/%s/0.25/woa23_%s_t00_04.nc"%(dec, dec)
    urlS_dec = "https://www.ncei.noaa.gov/thredds-ocean/dodsC/woa23/DATA/salinity/netcdf/%s/0.25/woa23_%s_s00_04.nc"%(dec, dec)
    decm = decs_fout[dec]
    open_dataset(urlT_dec, decode_times=False).squeeze().sel(bbox).to_netcdf("woa23T-%s.nc"%decm)
    open_dataset(urlS_dec, decode_times=False).squeeze().sel(bbox).to_netcdf("woa23S-%s.nc"%decm)
