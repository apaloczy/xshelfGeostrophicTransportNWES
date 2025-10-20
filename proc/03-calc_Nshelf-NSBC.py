import numpy as np
import matplotlib.pyplot as plt
from xarray import open_dataset
from pandas import DataFrame, read_json
from pathlib import Path
from glob import glob
import cartopy.crs as ccrs
from cmocean.cm import matter
from cartopy.feature import LAND
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER


def point_in_poly(x, y, poly):
    """
    USAGE
    -----
    isinside = point_in_poly(x, y, poly)

    Determine if a point is inside a given polygon or not
    Polygon is a list of (x, y) pairs. This fuction
    returns True or False. The algorithm is called
    'Ray Casting Method'.

    Source: https://pseentertainmentcorp.com/smf/index.php?topic=545.0
    """
    n = len(poly)
    inside = False

    p1x, p1y = poly[0]
    for i in range(n + 1):
        p2x, p2y = poly[i % n]
        if y > min(p1y, p2y):
            if y <= max(p1y, p2y):
                if x <= max(p1x, p2x):
                    if p1y != p2y:
                        xinters = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                    if p1x == p2x or x <= xinters:
                        inside = not inside
        p1x, p1y = p2x, p2y

    return inside


def near2(x, y, x0, y0):
    dr2 = (x - x0) ** 2 + (y - y0) ** 2
    return int(np.nanargmin(dr2))


def poly_avg_plot(Nmean, Nstd, ssns, lons, lats, xis, yis, xy, axs, col="k", min_Nobs=20):
    xx, yy = np.flipud(xy[:, 0]), np.flipud(xy[:, 1])
    xs, ys = np.concatenate((xis, xx)), np.concatenate((yis, yy))
    xs, ys = np.append(xs, xs[0]), np.append(ys, ys[0])
    ps = np.array([(x0, y0) for x0, y0 in zip(xs, ys)])

    fps = []
    for xii, yii in zip(lons, lats):
        fps.append(point_in_poly(xii, yii, ps))
    lonsp, latsp = lons[fps], lats[fps]

    Nmeans_ssn = dict()
    Nstds_ssn = dict()
    Nobs_ssn = dict()
    for ssn in ssns:
        Nmeans_in = Nmean[ssn].values.ravel()[fps]
        Nmeans = np.nanmean(Nmeans_in)
        Nstds_in = Nstd[ssn].values.ravel()[fps]
        Nstds = np.nanstd(Nstds_in)
        Nnobs_in = Nnobs[ssn].values.ravel()[fps]
        Nobs_in = np.nansum(Nnobs_in)
        Nobs_ssn[ssn] = Nobs_in
        if Nobs_in < min_Nobs:
            print(
                "Not enough profiles (%d, the minimum required is %d)."
                % (Nobs_in, min_Nobs)
            )
            Nmeans = np.nan
            Nstds = np.nan
        Nmeans_ssn[ssn] = Nmeans
        Nstds_ssn[ssn] = Nstds

    _ = [axi.plot(xs, ys, color=col, linestyle="dashed", linewidth=1.5, zorder=20) for axi in axs]
    _ = [axi.plot(lonsp, latsp, linestyle="none", marker=".", ms=1, mfc=col, zorder=21) for axi in axs[4:]]

    xt, yt = lonsp.mean(), latsp.mean()
    for ssn, axi in zip(ssns, axs[:4]):
        Ntxt = "%1.1f" % Nmeans_ssn[ssn]
        axi.text(
            xt,
            yt,
            Ntxt,
            color="k",
            fontsize=10,
            fontweight="black",
            ha="center",
            va="center",
            zorder=999,
            bbox=dict(facecolor="w", edgecolor="k", boxstyle="circle"),
        )

    for ssn, axi in zip(ssns, axs[4:]):
        Nobstxt = "%d" % Nobs_ssn[ssn]
        axi.text(
            xt,
            yt,
            Nobstxt,
            color="k",
            fontsize=10,
            fontweight="black",
            ha="center",
            va="center",
            zorder=999,
            bbox=dict(facecolor="w", edgecolor="k", boxstyle="circle"),
        )

    return Nmeans_ssn, Nstds_ssn, Nobs_ssn


proj = ccrs.PlateCarree()


def bmap(fig, ax, bb, proj=proj, xticks=None, yticks=None, dlon=2, dlat=2, land=None, coastlines=True):
    ax.set_extent(bb, proj)

    if land is not None:
        ax.add_feature(land, zorder=9999999998, color="gray")

    if coastlines:
        ax.coastlines(zorder=9999999999)

    if not xticks:
        xticks = np.arange(bb[0], bb[1], dlon)
    if not yticks:
        yticks = np.arange(bb[2], bb[3], dlat)

    gl = ax.gridlines(crs=proj, linewidth=0, color="gray", alpha=0.5, linestyle="--")
    gl.top_labels = True
    gl.left_labels = True
    gl.bottom_labels = False
    gl.right_labels = False
    gl.xlines = True
    gl.ylines = True
    gl.xlocator = mticker.FixedLocator(xticks)
    gl.ylocator = mticker.FixedLocator(yticks)
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {"size": 8, "color": "gray"}
    gl.ylabel_style = {"size": 8, "color": "gray"}  # , "rotation":45}

    return ax


# ---
plt.close("all")

xmin, xmax = -16.5, 10.5
ymin, ymax = 44, 63
bb = (xmin, xmax, ymin, ymax)
ssns = ["JFM", "AMJ", "JAS", "OND"]

head = "../data/nsbc/level2/climatological_monthly_mean/"
fN = "NSBC_Level2_nitrate__UHAM_ICDC__v1.1__0.25x0.25deg__1960_2014.nc"
fS = "NSBC_Level2_salinity__UHAM_ICDC__v1.1__0.25x0.25deg__1960_2014.nc"

ds = open_dataset(head + fN).sel(depth=slice(0, 200))
ds = ds.sel(lon=slice(xmin, xmax), lat=slice(ymin, ymax))

dsS = open_dataset(head + fS).sel(depth=slice(0, 200))
dsS = dsS.sel(lon=slice(xmin, xmax), lat=slice(ymin, ymax))

msk200mshelf = ~np.isfinite(ds["nitrate_mean"].sel(depth=182).values)
msk200mshelf = np.expand_dims(msk200mshelf, 1)
Nmean0, Nnobs0 = (ds["nitrate_mean"].where(msk200mshelf), ds["nitrate_noobs"].where(msk200mshelf))

Nstd0 = ds["nitrate_stddev"].where(msk200mshelf)
Smean0 = dsS["salinity_mean"].where(
    msk200mshelf
)  # For depth-mean salinity. Anywhere too fresh is masked out for the nitrate values.
fmsk = "../data/nsbc/NSBC_Level2_land_sea_mask__UHAM_ICDC__v1.1__0.25x0.25deg.nc"
landmask0 = open_dataset(fmsk)["landmask"].sel(lon=slice(xmin, xmax), lat=slice(ymin, ymax), depth=slice(0, 200))

z = ds["depth"].values
dz = np.diff(z)[np.newaxis, :, np.newaxis, np.newaxis]

print("Calculating local depth at each grid point from the land mask.")
landmask = landmask0.values
_, Ny, Nx = landmask.shape
H = np.empty((Ny, Nx))
for ny in range(Ny):
    for nx in range(Nx):
        mskij = landmask[:, ny, nx]
        if mskij[0] == 0:
            H[ny, nx] = np.nan  # Land grid point.
        elif np.all(mskij == 1):  # Deepest data point above bottom.
            H[ny, nx] = z[-1]
        else:  # First bin above the bottom.
            fbot = np.where(mskij == 0)[0][0]
            H[ny, nx] = z[fbot]

print("Done.")

Nmean0 = Nmean0.mean("depth")
Nstd0 = Nstd0.mean("depth")
Smean0 = Smean0.mean("depth")

Nnobs0 = Nnobs0.max("depth")

# Mask low salinity grid points.
Smin = 34.5  # [psu]
mskSmean_zmean = Smean0.values > Smin
Nmean0, Nnobs0 = Nmean0.where(mskSmean_zmean), Nnobs0.where(mskSmean_zmean)

x, y = Nmean0["lon"].values, Nmean0["lat"].values
lons, lats = np.meshgrid(x, y)
lons, lats = lons.ravel(), lats.ravel()

# Seasonal averages, standard deviations, and total number of observations.
Nmean = dict()
Nmean["JFM"] = Nmean0.isel(time=slice(0, 3)).mean("time")
Nmean["AMJ"] = Nmean0.isel(time=slice(3, 6)).mean("time")
Nmean["JAS"] = Nmean0.isel(time=slice(6, 9)).mean("time")
Nmean["OND"] = Nmean0.isel(time=slice(9, 12)).mean("time")

Nstd = dict()
Nstd["JFM"] = Nstd0.isel(time=slice(0, 3)).mean("time")
Nstd["AMJ"] = Nstd0.isel(time=slice(3, 6)).mean("time")
Nstd["JAS"] = Nstd0.isel(time=slice(6, 9)).mean("time")
Nstd["OND"] = Nstd0.isel(time=slice(9, 12)).mean("time")

Nnobs = dict()
Nnobs["JFM"] = Nnobs0.isel(time=slice(0, 3)).sum("time")
Nnobs["AMJ"] = Nnobs0.isel(time=slice(3, 6)).sum("time")
Nnobs["JAS"] = Nnobs0.isel(time=slice(6, 9)).sum("time")
Nnobs["OND"] = Nnobs0.isel(time=slice(9, 12)).sum("time")

# Plot bin-averaged nitrate and number of observations per bin.
fac = 1.2
figsize = (21 * fac, 8 * fac)

fig = plt.figure(figsize=figsize)
ax11 = fig.add_subplot(241, projection=proj)
ax12 = fig.add_subplot(242, projection=proj, sharex=ax11, sharey=ax11)
ax13 = fig.add_subplot(243, projection=proj, sharex=ax11, sharey=ax11)
ax14 = fig.add_subplot(244, projection=proj, sharex=ax11, sharey=ax11)
ax21 = fig.add_subplot(245, projection=proj, sharex=ax11, sharey=ax11)
ax22 = fig.add_subplot(246, projection=proj, sharex=ax11, sharey=ax11)
ax23 = fig.add_subplot(247, projection=proj, sharex=ax11, sharey=ax11)
ax24 = fig.add_subplot(248, projection=proj, sharex=ax11, sharey=ax11)
ax11 = bmap(fig, ax11, bb, proj=proj, land=LAND)
ax12 = bmap(fig, ax12, bb, proj=proj, land=LAND)
ax13 = bmap(fig, ax13, bb, proj=proj, land=LAND)
ax14 = bmap(fig, ax14, bb, proj=proj, land=LAND)
ax21 = bmap(fig, ax21, bb, proj=proj, land=LAND)
ax22 = bmap(fig, ax22, bb, proj=proj, land=LAND)
ax23 = bmap(fig, ax23, bb, proj=proj, land=LAND)
ax24 = bmap(fig, ax24, bb, proj=proj, land=LAND)
axs = [ax11, ax12, ax13, ax14, ax21, ax22, ax23, ax24]

# fig.subplots_adjust(hspace=0, wspace=0)

head = "../data/woa23"
fbathymetry = "../data/srtm15p/SRTM15_V2.7.nc"
isob = 200

nobsmin = 0.1
nobsmax = 10
vmi, vma = 1, 13
clevs = np.arange(vmi, vma + 1, 1)
ccfmt = "%d"
cmap = matter
cmap2 = plt.cm.nipy_spectral
cmap2.set_extremes(bad="w", under="w")
cbax = [0.05, 0.2, 0.2, 0.025]
cbname1 = r"N [mmol/m$^3$]"
cbname2 = r"# obs."

# Bin-averaged fields, number of observations per bin.
ssn = "JFM"
cs11 = ax11.pcolormesh(x, y, Nmean[ssn], vmin=vmi, vmax=vma, cmap=cmap)
cs21 = ax21.pcolormesh(x, y, Nnobs[ssn], vmin=nobsmin, vmax=nobsmax, cmap=cmap2)
ax11.set_title(ssn, fontsize=14, fontweight="black")

ssn = "AMJ"
cs12 = ax12.pcolormesh(x, y, Nmean[ssn], vmin=vmi, vmax=vma, cmap=cmap)
cs22 = ax22.pcolormesh(x, y, Nnobs[ssn], vmin=nobsmin, vmax=nobsmax, cmap=cmap2)
ax12.set_title(ssn, fontsize=14, fontweight="black")

ssn = "JAS"
cs13 = ax13.pcolormesh(x, y, Nmean[ssn], vmin=vmi, vmax=vma, cmap=cmap)
cs23 = ax23.pcolormesh(x, y, Nnobs[ssn], vmin=nobsmin, vmax=nobsmax, cmap=cmap2)
ax13.set_title(ssn, fontsize=14, fontweight="black")

ssn = "OND"
cs14 = ax14.pcolormesh(x, y, Nmean[ssn], vmin=vmi, vmax=vma, cmap=cmap)
cs24 = ax24.pcolormesh(x, y, Nnobs[ssn], vmin=nobsmin, vmax=nobsmax, cmap=cmap2)
ax14.set_title(ssn, fontsize=14, fontweight="black")

# Isobaths and colorbars.
dstopo = open_dataset(fbathymetry).sel(lon=slice(xmin, xmax), lat=slice(ymin, ymax))
xt, yt = np.meshgrid(dstopo["lon"].values, dstopo["lat"].values)
ht = -dstopo["z"].values

# Load smoothed isobaths.
d = np.load("isobath_100m.npz")
xi100, yi100 = d["xi_ext"], d["yi_ext"]

d = np.load("isobath_125m.npz")
xi125, yi125 = d["xi_ext"], d["yi_ext"]

d = np.load("isobath_200m.npz")
xi, yi = d["xi_ext"], d["yi_ext"]

# Load Wei et al. (2024) sections.
xseg0, yseg0 = np.loadtxt("segment_latitudes_wei_etal2024.txt", unpack=True)
fseg = [near2(xi, yi, xseg0[i], yseg0[i]) for i in range(len(yseg0))]
xseg, yseg = xi[fseg], yi[fseg]

_ = [axi.plot(xi100, yi100, color="k", linewidth=1, zorder=20) for axi in axs]
_ = [axi.plot(xi125, yi125, color="k", linewidth=1, zorder=20) for axi in axs]
_ = [axi.plot(xi, yi, color="r", linewidth=1.5, zorder=20) for axi in axs]
_ = [axi.plot(xseg, yseg, linestyle="none", marker="o", ms=6, mfc="w", mec="r", zorder=21) for axi in axs]

cb11 = fig.colorbar(cs11, cax=ax11.inset_axes(cbax), orientation="horizontal", extend="both")
cb21 = fig.colorbar(cs21, cax=ax21.inset_axes(cbax), orientation="horizontal", extend="both")
cb11.set_label(cbname1, fontsize=12)
cb21.set_label(cbname2, fontsize=12)

fig.tight_layout()
# fig.savefig("binavg_nitrate_NSBC.png", bbox_inches="tight")

#####################################
### Create polygons for each segment.

Nmeans_allsegs = dict()
Nstds_allsegs = dict()
Nobs_allsegs = dict()

## CSS.
seg = "CSS"
xsegl, xsegr = xseg[1], xseg[2]
il, ir = np.where(xi == xsegl)[0][0], np.where(xi == xsegr)[0][0] + 1
xis, yis = xi[il:ir], yi[il:ir]

# pts = plt.ginput(n=-1, timeout=-1)
xy = np.array(
    [
        (np.float64(-4.864641986253229), np.float64(48.75028951080567)),
        (np.float64(-5.027830437645855), np.float64(49.416384885121715)),
        (np.float64(-6.758231903940402), np.float64(49.78089058700602)),
        (np.float64(-6.787513136397685), np.float64(50.29278551413694)),
        (np.float64(-6.597992641428364), np.float64(50.63222311053623)),
        (np.float64(-6.784212853602739), np.float64(50.94343275622131)),
        (np.float64(-7.944297365417542), np.float64(51.005014628799344)),
        (np.float64(-8.189220480072226), np.float64(51.207314942251216)),
        (np.float64(-9.136050548307248), np.float64(51.32331637155441)),
    ]
)

Nmeans_ssn, Nstds_ssn, Nobs_ssn = poly_avg_plot(Nmean, Nstd, ssns, lons, lats, xis, yis, xy, axs, col="c")
Nmeans_allsegs.update({seg: Nmeans_ssn})
Nstds_allsegs.update({seg: Nstds_ssn})
Nobs_allsegs.update({seg: Nobs_ssn})


## PSBS.
seg = "PSBS"
xsegl, xsegr = xseg[2], xseg[3]
il, ir = np.where(xi == xsegl)[0][0], np.where(xi == xsegr)[0][0] + 1
xis, yis = xi[il:ir], yi[il:ir]

# pts = plt.ginput(n=-1, timeout=-1)
xy = np.array(
    [
        tuple(xy[-1]),
        (np.float64(-9.27435793314636), np.float64(51.328738414782066)),
        (np.float64(-9.568526989683543), np.float64(51.30122465202925)),
        (np.float64(-9.945840139922904), np.float64(51.36042875810925)),
        (np.float64(-10.12069677613448), np.float64(51.42449836527075)),
        (np.float64(-10.368062295719682), np.float64(51.57205652638798)),
        (np.float64(-10.48268144509058), np.float64(51.704372321285334)),
        (np.float64(-10.63239248279254), np.float64(51.8911558937763)),
        (np.float64(-10.63239248279254), np.float64(52.14136321045335)),
        (np.float64(-10.56298374170092), np.float64(52.28065432548522)),
        (np.float64(-10.56298374170092), np.float64(52.46920325677567)),
        (np.float64(-10.365134383564452), np.float64(52.724233017002504)),
        (np.float64(-10.093484416103522), np.float64(52.89310465630854)),
    ]
)

Nmeans_ssn, Nstds_ssn, Nobs_ssn = poly_avg_plot(Nmean, Nstd, ssns, lons, lats, xis, yis, xy, axs, col="m")
Nmeans_allsegs.update({seg: Nmeans_ssn})
Nstds_allsegs.update({seg: Nstds_ssn})
Nobs_allsegs.update({seg: Nobs_ssn})


## MS.
seg = "MS"
xsegl, xsegr = xseg[3], xseg[4]
il, ir = np.where(xi == xsegl)[0][0], np.where(xi == xsegr)[0][0] + 1
xis, yis = xi[il:ir], yi[il:ir]

# pts = plt.ginput(n=-1, timeout=-1)
xy = np.array(
    [
        tuple(xy[-1]),
        (np.float64(-10.081244082594491), np.float64(53.04204014795244)),
        (np.float64(-10.3670516316629), np.float64(53.14169631128712)),
        (np.float64(-10.26606671948376), np.float64(53.21281601134311)),
        (np.float64(-10.372809543322237), np.float64(53.6017597802436)),
        (np.float64(-10.454369412870427), np.float64(53.905283980571504)),
        (np.float64(-10.399700888983975), np.float64(54.09042298623326)),
        (np.float64(-10.130787432366589), np.float64(54.398249801867046)),
        (np.float64(-10.013414617772412), np.float64(54.57503667193504)),
        (np.float64(-9.863898735893144), np.float64(54.782764407952655)),
        (np.float64(-9.573598750255364), np.float64(54.935697072574825)),
        (np.float64(-9.245650880691151), np.float64(55.082302361747175)),
        (np.float64(-9.080696202714321), np.float64(55.25605209159926)),
        (np.float64(-9.038176239691527), np.float64(55.42549920614546)),
        (np.float64(-8.884357742506381), np.float64(55.587100375133886)),
        (np.float64(-8.635944982346416), np.float64(55.690742785001945)),
        (np.float64(-8.177463357252396), np.float64(55.712762051567324)),
        (np.float64(-7.632423236004829), np.float64(55.84696568485803)),
        (np.float64(-7.368381858483804), np.float64(55.94725458573768)),
    ]
)

Nmeans_ssn, Nstds_ssn, Nobs_ssn = poly_avg_plot(Nmean, Nstd, ssns, lons, lats, xis, yis, xy, axs, col="c")
Nmeans_allsegs.update({seg: Nmeans_ssn})
Nstds_allsegs.update({seg: Nstds_ssn})
Nobs_allsegs.update({seg: Nobs_ssn})


## HS.
seg = "HS"
xsegl, xsegr = xseg[4], xseg[5]
il, ir = np.where(xi == xsegl)[0][0], np.where(xi == xsegr)[0][0] + 1
xis, yis = xi[il:ir], yi[il:ir]

# pts = plt.ginput(n=-1, timeout=-1)
xy = np.array(
    [
        tuple(xy[-1]),
        (np.float64(-7.247067121305481), np.float64(56.04129008911888)),
        (np.float64(-7.393708772105212), np.float64(56.19390582452306)),
        (np.float64(-7.372550555797819), np.float64(56.38687428880629)),
        (np.float64(-7.217058409287805), np.float64(56.532879810226845)),
        (np.float64(-6.7655448390679584), np.float64(56.78445515079161)),
        (np.float64(-6.767287280410921), np.float64(56.9596673080561)),
        (np.float64(-6.649299109473226), np.float64(57.136704880060826)),
        (np.float64(-6.8607429835516776), np.float64(57.243768220356145)),
        (np.float64(-6.882454355840963), np.float64(57.39552656398838)),
        (np.float64(-6.726630315741813), np.float64(57.54725724982153)),
        (np.float64(-6.686139297867273), np.float64(57.72708825953482)),
        (np.float64(-6.938710319199442), np.float64(57.67191095034103)),
        (np.float64(-7.0567261479362315), np.float64(57.49183102043589)),
        (np.float64(-7.181103270464785), np.float64(57.107470586417804)),
        (np.float64(-7.328574655237357), np.float64(56.95264222708606)),
        (np.float64(-7.453587907145089), np.float64(56.829454389918574)),
        (np.float64(-7.645754295254582), np.float64(56.67777901968362)),
        (np.float64(-7.912569083120486), np.float64(56.632115993378385)),
        (np.float64(-7.974439579695174), np.float64(56.79853297053078)),
        (np.float64(-7.976071389841758), np.float64(56.97731298387847)),
        (np.float64(-8.005388656882065), np.float64(57.157061020194476)),
        (np.float64(-8.176839353469681), np.float64(57.30854278583577)),
        (np.float64(-8.138090776938103), np.float64(57.46580503148784)),
        (np.float64(-8.006743889037702), np.float64(57.61543372458979)),
        (np.float64(-7.806418450195292), np.float64(57.77291723263462)),
        (np.float64(-7.882919922491041), np.float64(57.80812561088209)),
        (np.float64(-8.170173823887875), np.float64(57.76763459300754)),
        (np.float64(-8.128963703236877), np.float64(57.93056668747401)),
        (np.float64(-7.8553174389946), np.float64(58.19846012950461)),
        (np.float64(-7.5691698495615505), np.float64(58.393973111304554)),
        (np.float64(-7.220100767188216), np.float64(58.67215525459835)),
        (np.float64(-6.688711473183075), np.float64(58.64186996458973)),
        (np.float64(-6.193387949197332), np.float64(58.58893293712261)),
        (np.float64(-6.006089333728493), np.float64(58.44151686794822)),
        (np.float64(-6.074625359884987), np.float64(58.27089590533345)),
        (np.float64(-6.234155545062817), np.float64(58.11543141662253)),
        (np.float64(-6.19654093829412), np.float64(58.0099722286747)),
        (np.float64(-5.991596647002906), np.float64(58.125332908698404)),
        (np.float64(-5.807561651827234), np.float64(58.288458607758535)),
        (np.float64(-5.638434210363565), np.float64(58.15415233535501)),
        (np.float64(-5.530983660880926), np.float64(58.261077386654854)),
        (np.float64(-5.514416639223242), np.float64(58.43410457779086)),
        (np.float64(-5.577863630346325), np.float64(58.61393558750415)),
        (np.float64(-5.645763527123643), np.float64(58.79622814133686)),
        (np.float64(-5.86213048944094), np.float64(58.91949895190163)),
        (np.float64(-6.224088106192372), np.float64(58.931806672498745)),
        (np.float64(-6.3104634127649035), np.float64(59.06879575141445)),
    ]
)

Nmeans_ssn, Nstds_ssn, Nobs_ssn = poly_avg_plot(Nmean, Nstd, ssns, lons, lats, xis, yis, xy, axs, col="m")
Nmeans_allsegs.update({seg: Nmeans_ssn})
Nstds_allsegs.update({seg: Nstds_ssn})
Nobs_allsegs.update({seg: Nobs_ssn})


## WSS.
seg = "WSS"
xsegl, xsegr = xseg[5], xseg[6]
il, ir = np.where(xi == xsegl)[0][0], np.where(xi == xsegr)[0][0] + 1
xis, yis = xi[il:ir], yi[il:ir]

# pts = plt.ginput(n=-1, timeout=-1)
xy = np.array(
    [
        tuple(xy[-1]),
        (np.float64(-6.197577208575989), np.float64(59.13471643189948)),
        (np.float64(-6.077773844444643), np.float64(59.17838435001397)),
        (np.float64(-5.749108142997262), np.float64(59.275016056661926)),
        (np.float64(-5.24758459746635), np.float64(59.11499104832751)),
        (np.float64(-4.971655956005627), np.float64(59.2111239521725)),
        (np.float64(-4.631427098671228), np.float64(59.222369688093984)),
        (np.float64(-4.5031894325987665), np.float64(59.40298164847822)),
        (np.float64(-4.515523465544916), np.float64(59.50845576845154)),
        (np.float64(-4.589573008931172), np.float64(59.609440663198136)),
        (np.float64(-4.321443829480511), np.float64(59.716320500161636)),
        (np.float64(-4.19579086884162), np.float64(59.552939909334526)),
        (np.float64(-4.2451723463355755), np.float64(59.38221331359095)),
        (np.float64(-4.09081555167127), np.float64(59.20921943237934)),
        (np.float64(-3.8674879330691248), np.float64(59.11802921085469)),
        (np.float64(-3.5895640803965243), np.float64(59.21085187791633)),
        (np.float64(-3.4316703204020733), np.float64(59.35681971634888)),
        (np.float64(-3.3313656112958903), np.float64(59.533531945728086)),
        (np.float64(-3.49247891665496), np.float64(59.65882214069209)),
        (np.float64(-3.3569859370847652), np.float64(59.807283993139485)),
        (np.float64(-2.8846196826728736), np.float64(60.04326306465331)),
        (np.float64(-2.557223661087594), np.float64(60.11409306267494)),
        (np.float64(-2.35280520328892), np.float64(60.01650909613041)),
        (np.float64(-1.674025279866477), np.float64(60.618155967930065)),
        (np.float64(-1.084576403884892), np.float64(60.84461244047811)),
        (np.float64(-0.7265720284810442), np.float64(60.861798464325865)),
        (np.float64(-0.6599138357206087), np.float64(60.64989796448265)),
        (np.float64(-0.7603092362455115), np.float64(60.29674557998042)),
        (np.float64(-1.1721389686606791), np.float64(59.805515510474414)),
        (np.float64(-1.3329801997635826), np.float64(59.41513429858693)),
        (np.float64(-1.634438475594461), np.float64(58.9957318327085)),
        (np.float64(-1.6446866059100103), np.float64(58.693548025527846)),
        (np.float64(-1.7655782670660152), np.float64(58.41172444185021)),
        (np.float64(-1.2015683340358656), np.float64(57.85129681985957)),
        (np.float64(-0.8806114031798238), np.float64(57.911017119087795)),
        (np.float64(-0.8291893687646308), np.float64(57.819826897563146)),
        (np.float64(-0.9737968359163567), np.float64(57.67648911027352)),
        (np.float64(-0.6531119793164812), np.float64(57.56380502251183)),
        (np.float64(-0.0631643005319269), np.float64(57.53895557378209)),
        (np.float64(0.5088718230550953), np.float64(57.66619563424861)),
        (np.float64(0.8753105004000581), np.float64(57.65204777292803)),
        (np.float64(1.309858433205008), np.float64(57.864900532667825)),
    ]
)

Nmeans_ssn, Nstds_ssn, Nobs_ssn = poly_avg_plot(Nmean, Nstd, ssns, lons, lats, xis, yis, xy, axs, col="c")
Nmeans_allsegs.update({seg: Nmeans_ssn})
Nstds_allsegs.update({seg: Nstds_ssn})
Nobs_allsegs.update({seg: Nobs_ssn})


## NENSS.
seg = "NENSS"
xsegl, xsegr = xseg[6], xseg[7]
il, ir = np.where(xi == xsegl)[0][0], np.where(xi == xsegr)[0][0] + 1
xis, yis = xi[il:ir], yi[il:ir]

# pts = plt.ginput(n=-1, timeout=-1)
xy = np.array(
    [
        tuple(xy[-1]),
        (np.float64(1.5387830481325828), np.float64(58.051498389314915)),
        (np.float64(1.603917973582913), np.float64(58.238850256538925)),
        (np.float64(1.7389446591771867), np.float64(58.3716296186559)),
        (np.float64(2.1499426677835567), np.float64(58.693558706778724)),
        (np.float64(2.4688378692120296), np.float64(58.63363008071747)),
        (np.float64(3.3290009301475334), np.float64(58.286942978953086)),
        (np.float64(4.033649206454071), np.float64(58.10419812485753)),
        (np.float64(4.826654750809677), np.float64(57.941267172753484)),
        (np.float64(5.253346568365839), np.float64(57.681963498864675)),
        (np.float64(6.092984074875362), np.float64(57.492626495902385)),
        (np.float64(6.863366562892843), np.float64(57.37112120656319)),
        (np.float64(7.714240686788852), np.float64(57.401797172028296)),
        (np.float64(8.12093132539707), np.float64(57.445619979835584)),
        (np.float64(8.195055544756588), np.float64(57.46906705478206)),
    ]
)


Nmeans_ssn, Nstds_ssn, Nobs_ssn = poly_avg_plot(Nmean, Nstd, ssns, lons, lats, xis, yis, xy, axs, col="m")
Nmeans_allsegs.update({seg: Nmeans_ssn})
Nstds_allsegs.update({seg: Nstds_ssn})
Nobs_allsegs.update({seg: Nobs_ssn})

# fig.savefig("binavg_nitrate_NSBC_polygons.png", bbox_inches="tight")

# Done with all segments.
############################################################################

# Put segment-average nitrate concentrations into pandas DataFrames and save as json files.
# Columns: Segments.
# Rows: Seasons.

dfNmeans = DataFrame(Nmeans_allsegs)
dfNstds = DataFrame(Nstds_allsegs)
dfNobs = DataFrame(Nobs_allsegs, dtype=int)

dfNsterr = dfNstds / np.sqrt(dfNobs)

N_fileout = "Nmeans_shelf-NSBC.json"
Nstd_fileout = "Nstds_shelf-NSBC.json"
Nobs_fileout = "Nobs_shelf-NSBC.json"
dfNmeans.to_json(N_fileout)
dfNstds.to_json(Nstd_fileout)
dfNobs.to_json(Nobs_fileout)

# Read the json files into DataFrames.
dfNmeans_in = read_json(N_fileout)
dfNstdss_in = read_json(Nstd_fileout)
dfNobs_in = read_json(Nobs_fileout)
