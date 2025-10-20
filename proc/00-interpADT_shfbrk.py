# Interpolate CMEMS ADT to the shelfbreak isobath for different decades.
import numpy as np
import matplotlib.pyplot as plt
from gsw import distance, grav
from gsw import f as fcor
import xarray as xr
from os.path import isfile
from copernicusmarine import open_dataset
from scipy.interpolate import LinearNDInterpolator, CloughTocher2DInterpolator
import cartopy.crs as ccrs


def get_isobath(lon, lat, topo, iso):
    """
    USAGE
    -----
    lon_isob, lat_isob = get_isobath(lon, lat, topo, iso)

    Retrieves the 'lon_isob','lat_isob' coordinates of a wanted 'iso'
    isobath from a topography array 'topo', with 'lon_topo','lat_topo'
    coordinates.
    """
    lon, lat, topo = map(np.array, (lon, lat, topo))

    fig, ax = plt.subplots()
    cs = ax.contour(lon, lat, topo, [iso])
    coll = cs.allsegs[0]

    ## Test all lines to find thel longest one.
    ## This is assumed to be the wanted isobath.
    ncoll = len(coll)
    siz = np.array([])
    for n in range(ncoll):
        siz = np.append(siz, coll[n].shape[0])

    f = siz.argmax()
    coll = coll[f]
    xiso = coll[:, 0]
    yiso = coll[:, 1]
    plt.close()

    return xiso, yiso


def near2(x, y, x0, y0):
    dr2 = (x - x0)**2 + (y - y0)**2
    return int(np.nanargmin(dr2))


#---
plt.close("all")

id = "cmems_obs-sl_eur_phy-ssh_my_allsat-l4-duacs-0.0625deg_P1D" # Delayed mode, https://doi.org/10.48670/moi-00141

xmin, xmax = -16.5, 10.5
ymin, ymax = 44, 63

iso = 200 # [m]
binlen = 75 # [km]

years = [1995, 2005, 2015, 2025]
decavg_bins = ["%d-01-01"%year for year in years]
yearavg_bins = ["%d-01-01"%year for year in np.arange(1993, 2026)]

method = "cubic"
tmin, tmax = decavg_bins[0], decavg_bins[-1]
tmin_monthly = "1993-01-01"
decavg_bins = [np.datetime64(ti) for ti in decavg_bins]
yearavg_bins = [np.datetime64(ti) for ti in yearavg_bins]

fbathymetry = "../data/srtm15p/SRTM15_V2.7.nc"

fadt_decavg = "altimetry_adt_decavg.nc"
fadt_yavg = "altimetry_adt_yearavg.nc"
fadt_ssnavg = "altimetry_adt_ssnavg.nc"

# Get isobath contour.
dstopo0 = xr.open_dataset(fbathymetry)
dstopo = dstopo0.sel(lon=slice(xmin, xmax), lat=slice(ymin, ymax))
xt, yt = np.meshgrid(dstopo["lon"].values, dstopo["lat"].values)
ht = -dstopo["z"].values
xi, yi = get_isobath(xt, yt, ht, iso)

xmax_ext = xmax + 1
ymin_ext = ymin - 1
dstopo_ext = dstopo0.sel(lon=slice(xmin, xmax_ext), lat=slice(ymin_ext, ymax))
xt_ext, yt_ext = np.meshgrid(dstopo_ext["lon"].values, dstopo_ext["lat"].values)
ht_ext = -dstopo_ext["z"].values
xi_ext, yi_ext = get_isobath(xt_ext, yt_ext, ht_ext, iso)

# Replace isobath with extended one for avoiding boundary problems in the smoothing, clip later.
xi, yi = xi_ext.copy(), yi_ext.copy()

dd = np.append(0, np.cumsum(distance(xi, yi)))*1e-3 # [km]
dxi = np.median(np.diff(dd))
di = np.arange(0, dd[-1] + dxi, dxi)
xi = np.interp(di, dd, xi)
yi = np.interp(di, dd, yi)

nbins = int(di[-1]/binlen)

dix = di[1] - di[0]
nwin = int(binlen/dix)
xis = xr.DataArray(xi, coords=dict(x=di)).rolling(x=nwin, center=True).mean().dropna(dim="x").values
yis = xr.DataArray(yi, coords=dict(x=di)).rolling(x=nwin, center=True).mean().dropna(dim="x").values

##################
figm = plt.figure()
axm = plt.axes(projection=ccrs.PlateCarree())
axm.coastlines()
axm.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)
axm.plot(xi, yi, "k", linewidth=3, alpha=0.4)
axm.plot(xis, yis, "r", linewidth=1)
axm.set_extent((xmin, xmax_ext, ymin_ext, ymax))
axm.set_title("%d m isobath"%iso, fontsize=16)
# figm.savefig("alongisob_map_%dm.png"%iso, bbox_inches="tight", dpi=200)
plt.show(block=False)

##################

xi, yi = xis.copy(), yis.copy()
coordsi = dict(longitude=xi, latitude=yi)
di = np.append(0, np.cumsum(distance(xi, yi)))*1e-3 # [km]
xi0, yi0, di0 = xi.copy(), yi.copy(), di.copy()

#---

# Subsample isobath.
binlen = 20 # [km]

ns = int(di0[-1]/binlen) # Bin length increments.
di_aux = np.linspace(di0[0], di0[-1], num=ns)
xi_aux = np.interp(di_aux, di0, xi0)
yi_aux = np.interp(di_aux, di0, yi0)

################## Clip start/end to fit segments.
di_ext, xi_ext, yi_ext = di_aux.copy(), xi_aux.copy(), yi_aux.copy()

xseg0, yseg0 = np.loadtxt("segment_latitudes_wei_etal2024.txt", unpack=True)
fseg = [near2(xi_aux, yi_aux, xseg0[i], yseg0[i]) for i in range(len(yseg0))]

fl, fr = fseg[0], fseg[-1] + 1
xi, yi, di = xi_aux[fl:fr], yi_aux[fl:fr], di_aux[fl:fr]
xip, yip, dip = xi_aux[fl:fr+1], yi_aux[fl:fr+1], di_aux[fl:fr+1]
di -= di[0]
dip -= dip[0]

##################
figm = plt.figure()
axm = plt.axes(projection=ccrs.PlateCarree())
axm.coastlines()
axm.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)
axm.plot(xi_ext, yi_ext, "k", linewidth=3, alpha=0.4)
axm.plot(xi, yi, "r", linewidth=1)
axm.set_extent((xmin, xmax_ext, ymin_ext, ymax))
axm.set_title("Smoothed and clipped %d m isobath"%iso, fontsize=16)
# figm.savefig("alongisob_map_%dm_clipped.png"%iso, bbox_inches="tight", dpi=200)
plt.show(block=False)
##################

# Save isobath to file.
fisobath = "isobath_%dm.npz"%iso
np.savez(fisobath, xi=xi, yi=yi, di=di, xip=xip, yip=yip, dip=dip, xi_ext=xi_ext, yi_ext=yi_ext, di_ext=di_ext, xi0=xi0, yi0=yi0, di0=di0)

#---

# Interpolate at midpoints to compare with EN4 climatology at xim, yim points.
xi = 0.5*(xip[1:] + xip[:-1])
yi = 0.5*(yip[1:] + yip[:-1])
di = np.append(0, np.cumsum(distance(xi, yi)))*1e-3 # [km]

# Recalculate fseg after clipping isobath.
fseg = [near2(xi, yi, xseg0[i], yseg0[i]) for i in range(len(yseg0))]

# Clip start/end to fit segments.
xseg, yseg = xi[fseg], yi[fseg]
dibound = di[fseg]
txt_segs = ["AS", "CSS", "PSBS", "MS", "HS", "WSS", "NENSS"]

def annotate_segments(axs, dibound=dibound, txt_segs=txt_segs, ytxt=-190, lineplot=False, last_axs_txt_only=False):
    for axi in axs:
        _ = [axi.axvline(x=dd, color="k", linestyle="--", linewidth=1.5, zorder=2) for dd in dibound]
        for i in range(len(txt_segs)):
            xtxt = 0.5*(dibound[i] + dibound[i+1])
            if not last_axs_txt_only or (last_axs_txt_only and axi==axs[-1]):
                axi.text(xtxt, ytxt, txt_segs[i], ha="center", va="center", fontsize=12, color="k", bbox=dict(boxstyle="round,pad=0.3", fc="w", ec="k", lw=1))

#---

dyi = np.gradient(di*1e3) # [m]
g = grav(yi, 0)
f = fcor(yi)

dll = 1
xmin -= dll
xmax += dll
ymin -= dll
ymax += dll

kw = dict(dataset_id=id,
          minimum_longitude=xmin, maximum_longitude=xmax,
          minimum_latitude=ymin, maximum_latitude=ymax,
          start_datetime=tmin, end_datetime=tmax)

kw_monthly = dict(dataset_id=id,
                  minimum_longitude=xmin, maximum_longitude=xmax,
                  minimum_latitude=ymin, maximum_latitude=ymax,
                  start_datetime=tmin_monthly, end_datetime=tmax)

# Time-averaging ADT in years.
if not isfile(fadt_yavg):
    print("Downloading and time-averaging monthly ADT by year...")
    ds = open_dataset(**kw_monthly)
    adt = ds["adt"]
    adt_yavg = adt.groupby_bins("time", bins=yearavg_bins).mean()
    tbins = adt_yavg["time_bins"].values
    tl = np.array([ti.left.to_datetime64() for ti in tbins])

    adt_yavg = adt_yavg.rename(time_bins="tleft")
    adt_yavg["tleft"] = tl
    adt_yavg.to_netcdf(fadt_yavg)
else:
    adt_yavg = xr.open_dataset(fadt_yavg)["adt"]


# Time-averaging ADT in decades.
if not isfile(fadt_decavg):
    print("Downloading and time-averaging ADT...")
    ds = open_dataset(**kw)
    adt = ds["adt"]
    adt_tavg = adt.groupby_bins("time", bins=decavg_bins).mean()
    tbins = adt_tavg["time_bins"].values
    tl = np.array([ti.left.to_datetime64() for ti in tbins])

    adt_tavg = adt_tavg.rename(time_bins="tleft")
    adt_tavg["tleft"] = tl
    adt_tavg.to_netcdf(fadt_decavg)
else:
    adt_tavg = xr.open_dataset(fadt_decavg)["adt"]
xadt, yadt = np.meshgrid(adt_tavg["longitude"].values, adt_tavg["latitude"].values)


# Interpolate daily ADT maps to 200 m isobath.
print("Interpolating ADT to %d m isobath..."%iso)
if method=="linear":
    Interpolator = LinearNDInterpolator
elif method=="cubic":
    Interpolator = CloughTocher2DInterpolator

xyint = (xi, yi)
tls = []
ADTi = None
fig, ax = plt.subplots(nrows=2, sharex=True)
for n in range(adt_tavg.sizes["tleft"]):
    adtn = adt_tavg.isel(tleft=n)
    tl = str(adtn["tleft"].values).split("-")[0]
    tl = tl + "-" + str(int(tl) + 9)
    tls.append(tl)
    coordsi = dict(lon=("x", xi), lat=("x", yi), x=("x", di))
    adtn = adtn.values
    fg = np.isfinite(adtn)
    xy = (xadt[fg], yadt[fg])
    adtn = Interpolator(xy, adtn[fg], fill_value=np.nan)(xyint)
    dadtdy = np.gradient(adtn)/dyi
    uxi = -dadtdy*g/f
    ax[0].plot(di, adtn, label=tl)
    ax[1].plot(di, uxi, label=tl)
    adti = xr.DataArray(adtn, coords=coordsi, dims=("x"))
    if ADTi is not None:
        ADTi = xr.concat((ADTi, adti), dim="tl")
    else:
        ADTi = adti

ADTi["tl"] = np.array(tls)
ADTi.to_dataset(name="adt").to_netcdf("decadal_altimetry_adt_%dmisob.nc"%iso)

annotate_segments([ax[0]], lineplot=True, ytxt=0.12)
annotate_segments([ax[1]], lineplot=True, ytxt=0.05)

ax[0].set_ylim(-0.15, 0.15)
ax[1].set_ylim(-0.05, 0.06)

ax[0].set_xlim(di[0], di[-1])
ax[1].axhline(color="gray")
ax[0].legend()
ax[1].set_xlabel("Distance [km]", fontsize=14)
ax[0].set_ylabel("ADT [m]", fontsize=14)
ax[1].set_ylabel(r"Surface $u_g^\perp$ [m/s]", fontsize=14)
ax[0].set_title("Altimetry - Along-isobath ADT - Decadal variability")
# fig.savefig("decadal_altimetry_adt_%dmisob.png"%iso, bbox_inches="tight", dpi=200)

print("Done interpolating decadal ADT to isobath.")


tls = []
ADTi = None
fig, ax = plt.subplots(nrows=2, sharex=True)
nt = adt_yavg.sizes["tleft"]
for n in range(nt):
    print("%d / %d"%(n+1, nt))
    adtn = adt_yavg.isel(tleft=n)
    tl = str(adtn["tleft"].values).split("-")[0]
    tl = tl + "-" + str(int(tl) + 1)
    coordsi = dict(lon=("x", xi), lat=("x", yi), x=("x", di))
    adtn = adtn.values
    fg = np.isfinite(adtn)
    if fg.any():
        tls.append(tl)
        xy = (xadt[fg], yadt[fg])
        adtn = Interpolator(xy, adtn[fg], fill_value=np.nan)(xyint)
        dadtdy = np.gradient(adtn)/dyi
        uxi = -dadtdy*g/f
        ax[0].plot(di, adtn, label=tl)
        ax[1].plot(di, uxi, label=tl)
        adti = xr.DataArray(adtn, coords=coordsi, dims=("x"))
        if ADTi is not None:
            ADTi = xr.concat((ADTi, adti), dim="tl")
        else:
            ADTi = adti
    else:
        print("Skipping %s, all-NaN timestamp."%tl)

ADTi["tl"] = np.array(tls)
ADTi.to_dataset(name="adt").to_netcdf("yearly_altimetry_adt_%dmisob.nc"%iso)

print("Done interpolating yearly ADT to isobath.")



# Time-averaging ADT by seasons.
if not isfile(fadt_ssnavg):
    print("Downloading and seasonally-averaging ADT...")
    ds = open_dataset(**kw)
    adt = ds["adt"]

    from xarray import full_like
    months = adt.time.dt.month
    seasons = full_like(months, fill_value="none", dtype="U4")
    seasons.name = "season"

    # set values
    seasons[months.isin([1, 2, 3])] = "JFM"
    seasons[months.isin([4, 5, 6])] = "AMJ"
    seasons[months.isin([7, 8, 9])] = "JAS"
    seasons[months.isin([10, 11, 12])] = "OND"

    # seasonal mean, then reindex to get nicely ordered seasons
    adt_ssnavg = adt.groupby(seasons).mean().reindex(season=["JFM", "AMJ", "JAS", "OND"])

    adt_ssnavg = adt_ssnavg.to_dataset(name="adt")
    adt_ssnavg.to_netcdf(fadt_ssnavg)
    adt_ssnavg = adt_ssnavg["adt"]
else:
    adt_ssnavg = xr.open_dataset(fadt_ssnavg)["adt"]


tl0s = []
ADTi = None
fig, ax = plt.subplots(nrows=2, sharex=True)
for n in range(adt_ssnavg.sizes["season"]):
    adtn = adt_ssnavg.isel(season=n)
    coordsi = dict(lon=("x", xi), lat=("x", yi), x=("x", di))
    tl = adtn.season.values.flatten()[0]
    tl0s.append(tl)
    adtn = adtn.values
    fg = np.isfinite(adtn)
    xy = (xadt[fg], yadt[fg])
    adtn = Interpolator(xy, adtn[fg], fill_value=np.nan)(xyint)
    dadtdy = np.gradient(adtn)/dyi
    uxi = -dadtdy*g/f
    ax[0].plot(di, adtn, label=tl)
    ax[1].plot(di, uxi, label=tl)
    adti = xr.DataArray(adtn, coords=coordsi, dims=("x"))
    if ADTi is not None:
        ADTi = xr.concat((ADTi, adti), dim="season")
    else:
        ADTi = adti

ADTi["season"] = np.array(tl0s)
ADTi.to_dataset(name="adt").to_netcdf("seasonal_altimetry_adt_%dmisob.nc"%iso)

ax[0].set_xlim(di[0], di[-1])
ax[1].axhline(color="gray")
ax[0].legend(ncols=4)
ax[1].set_xlabel("Distance [km]")
ax[0].set_ylabel("ADT [m]")
ax[1].set_ylabel(r"$u^\perp$ [m/s]")
ax[0].set_title("Altimetry - Along-isobath ADT - Seasonal variability (1993-2023)")
# fig.savefig("seasonal_altimetry_adt_%dmisob.png"%iso, bbox_inches="tight", dpi=200)
