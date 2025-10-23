import numpy as np
from xarray import open_dataset, Variable, Coordinates, DataArray, Dataset, concat
from gsw import SA_from_SP, CT_from_pt, p_from_z
from glob import glob
from os import system

#---
fout_all = "en4profiles-NWES_all.nc"

fnames = glob("EN.4.2.2.f.profiles.g10.??????-NWES.nc")
fnames.sort()

Tattrs = dict(units="degrees Celsius", long_name="Conservative Temperature")
Sattrs = dict(units="g/kg", long_name="Absolute Salinity")

dzi = 5 # [m]
iso = 200 # [m]

zi = np.arange(0, iso + dzi, dzi, dtype=np.float32)
Nprof = 1
for f in fnames:
    print(f)
    ds = open_dataset(f)
    nx = ds.sizes["N_PROF"]

    yrmo = f.split(".")[-2]
    yr, mo = yrmo[:4], yrmo[4:]
    print("========= %s/%s ========="%(yr, mo))

    for n in range(nx):
        print("Profile %d / %d"%(n+1, nx))
        dsi = ds.isel(N_PROF=n)

        Spqc = dsi["PSAL_CORRECTED_QC"].astype(int)
        thetaqc = dsi["POTM_CORRECTED_QC"].astype(int)

        Sp = dsi["PSAL_CORRECTED"].where(Spqc==1).values
        theta = dsi["POTM_CORRECTED"].where(thetaqc==1).values
        z = dsi["DEPH_CORRECTED"].values

        if np.isnan(Sp).all() or np.isnan(theta).all():
            print("All-NaN profile in either T or S. Skipping.........")
            continue

        lat = dsi["LATITUDE"].values
        lon = dsi["LONGITUDE"].values
        t = dsi["JULD"].values

        p = p_from_z(-z, lat)
        SA = SA_from_SP(Sp, p, lon, lat)
        CT = CT_from_pt(SA, theta)

        nprof = np.array(Nprof, ndmin=1)
        lat = np.array(lat, ndmin=1)
        lon = np.array(lon, ndmin=1)
        t = np.array(t, ndmin=1)
        CT = np.array(CT, ndmin=2).T
        SA = np.array(SA, ndmin=2).T

        nprof = Variable("n", nprof, attrs=dict(long_name="Profile number"))
        t = Variable("n", t)
        lon = Variable("n", lon, attrs=dict(long_name="Longitude", units="Degrees east"))
        lat = Variable("n", lat, attrs=dict(long_name="Latitude", units="Degrees north"))
        z = Variable("z", z, attrs=dict(long_name="Depth", units="m"))
        coords = Coordinates(dict(z=z, n=nprof, t=t, lon=lon, lat=lat))

        CT = DataArray(CT, coords=coords, name="CT", attrs=Tattrs)
        SA = DataArray(SA, coords=coords, name="SA", attrs=Sattrs)

        dsii = Dataset(data_vars=dict(CT=CT, SA=SA), coords=coords).dropna(dim="z", how="all")

        if dsii.sizes["z"]==1:
            method = "nearest"
        else:
            method="linear"

        dsio = dsii.interp(z=zi, method=method)

        if Nprof==1:
            dsall = dsio
        else:
            dsall = concat((dsall, dsio), dim="n")
        Nprof += 1

    ds.close()
    cmd = "rm " + f
    _ = system(cmd) # Delete parent netCDF file.

dsall.to_netcdf(fout_all)
