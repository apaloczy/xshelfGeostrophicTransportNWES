import numpy as np
from xarray import open_dataset
from glob import glob
from os import system

#---
lonmin, lonmax = -16.5, 10.5
latmin, latmax = 44, 63

fnames = glob("EN.4.2.2.f.profiles.g10.??????.nc")
fnames.sort()

for f in fnames:
    ds = open_dataset(f)
    lats = ds["LATITUDE"].values
    lons = ds["LONGITUDE"].values
    fgud_prof = np.bool_(ds["PROFILE_POTM_QC"].values.astype("int")) # Only profiles that passed QC.

    fbb = np.logical_and(np.logical_and(lons>=lonmin, lons<=lonmax), np.logical_and(lats>=latmin, lats<=latmax))
    fbbgud = np.logical_and(fbb, fgud_prof)
    yrmo = f.split(".")[-2]
    yr, mo = yrmo[:4], yrmo[4:]
    nprofs = fbb.sum()
    nprofsgud = fbbgud.sum()
    nfrac = 100*nprofs/ds.sizes["N_PROF"]
    nbadprofs = nprofs - nprofsgud
    print("%s/%s: Number of profiles in bbox: %d (%1.1f%%). Number of NWES profiles that failed QC: %d"%(yr, mo, nprofs, nfrac, nbadprofs))

    dsout = ds.isel(N_PROF=np.where(fbbgud)[0])
    fout = f.replace(".nc", "-NWES.nc")
    dsout.to_netcdf(fout)

    ds.close()
    cmd = "rm " + f
    # _ = system(cmd) # Delete netCDF file.
    # _ = system(cmd.replace(".nc", ".zip")) # Delete .zip file.
