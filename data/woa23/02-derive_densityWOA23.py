import numpy as np
import matplotlib.pyplot as plt
from gsw import p_from_z, SA_from_SP, CT_from_t, sigma0
from xarray import open_dataset, Dataset
from glob import glob

#---
fnamesT = glob("woa23T*.nc")
fnamesT.sort()

for fT in fnamesT:
    ff = fT.split("woa23T-")[-1].strip(".nc")
    print(ff)
    dsT = open_dataset(fT, decode_times=False)
    dsS = open_dataset(fT.replace("T", "S"), decode_times=False)

    T, Sp = dsT["t_an"], dsS["s_an"]

    z = -T["depth"].values
    lon, lat = T["lon"].values, T["lat"].values
    lon, lat = np.meshgrid(lon, lat)
    p = p_from_z(z[:, np.newaxis, np.newaxis], lat)

    # Derive CT, SA and potential density.
    SA = SA_from_SP(Sp, p, lon, lat)
    CT = CT_from_t(SA, T, p)
    rho0_an = sigma0(SA, CT) + 1000

    T, Sp = dsT["t_mn"], dsS["s_mn"]
    SA = SA_from_SP(Sp, p, lon, lat)
    CT = CT_from_t(SA, T, p)
    rho0_mn = sigma0(SA, CT) + 1000
    rho0_dd = np.minimum(dsT["t_dd"], dsS["s_dd"])
    rho0_dd.attrs = dict()

    rho0_an.name = "rho0_an"
    rho0_an.attrs = dict(units="kg/m3", long_name="Objectively-mapped potential density referenced to 0 dbar", comment="Calculated using TEOS-10 from objectively-mapped WOA23 temperature and salinity fields")
    rho0_mn.name = "rho0_mn"
    rho0_mn.attrs = dict(units="kg/m3", long_name="Bin-averaged potential density referenced to 0 dbar", comment="Calculated using TEOS-10 from bin-averaged WOA23 temperature and salinity fields")
    rho0_dd.name = "rho0_dd"
    rho0_dd.attrs = dict(units="", long_name="Number of density observations in each bin", comment="Calculated as the minimum number of temperature or salinity profiles (whichever is less) from bin-averaged WOA23 temperature and salinity fields")

    # Save potential density as a separate file.
    fout = fT.replace("woa23T", "woa23D")
    dsout = Dataset(data_vars=dict(rho0_an=rho0_an, rho0_mn=rho0_mn, rho0_dd=rho0_dd), coords=T.coords)
    dsout.to_netcdf(fout)
