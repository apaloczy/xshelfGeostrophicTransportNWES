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


proj = ccrs.PlateCarree()

def bmap(fig, ax, bb, proj=proj, TOP_LABELS=False, xticks=None, yticks=None, dlon=2, dlat=2, land=None, coastlines=True):
    ax.set_extent(bb, proj)

    if land is not None:
        ax.add_feature(land, zorder=8, color="gray")

    if coastlines:
        ax.coastlines(zorder=9)

    if not xticks:
        xticks = np.arange(bb[0], bb[1], dlon)
    if not yticks:
        yticks = np.arange(bb[2], bb[3], dlat)

    gl = ax.gridlines(crs=proj, linewidth=0, color="gray", alpha=0.5, linestyle="--")
    gl.top_labels = TOP_LABELS
    gl.left_labels = True
    gl.bottom_labels = False
    gl.right_labels = False
    gl.xlines = True
    gl.ylines = True
    gl.xlocator = mticker.FixedLocator(xticks)
    gl.ylocator = mticker.FixedLocator(yticks)
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {"size": 12, "color": "gray"}
    gl.ylabel_style = {"size": 12, "color": "gray"}  # , "rotation":45}

    return ax


def near2(x, y, x0, y0):
    dr2 = (x - x0) ** 2 + (y - y0) ** 2
    return int(np.nanargmin(dr2))


# ---
plt.close("all")

xmin, xmax = -16.5, 10.5
ymin, ymax = 44, 63
bb = (xmin, xmax, ymin, ymax)
ssns = ["JFM", "AMJ", "JAS", "OND"]
head = "../data/woa23/"
nobsmin = 1
nobsmax = 20

fJFM = head + "woa23N-winter.nc"
fAMJ = head + "woa23N-spring.nc"
fJAS = head + "woa23N-summer.nc"
fOND = head + "woa23N-autumn.nc"

dsJFM = open_dataset(fJFM, decode_times=False)
fupper200m = dsJFM["depth"]<=200

dsJFM = dsJFM.where(fupper200m)
dsAMJ = open_dataset(fAMJ, decode_times=False).where(fupper200m)
dsJAS = open_dataset(fJAS, decode_times=False).where(fupper200m)
dsOND = open_dataset(fOND, decode_times=False).where(fupper200m)

ddJFM = dsJFM["n_dd"].max("depth")
ddAMJ = dsAMJ["n_dd"].max("depth")
ddJAS = dsJAS["n_dd"].max("depth")
ddOND = dsOND["n_dd"].max("depth")

x0, y0 = ddJFM["lon"].values, ddJFM["lat"].values
x, y = np.meshgrid(x0, y0)
dxy = y0[1] - y0[0]
x0 = np.append(x0, x0[-1] + dxy)
y0 = np.append(y0, y0[-1] + dxy)

ddJFM = ddJFM.values
ddAMJ = ddAMJ.values
ddJAS = ddJAS.values
ddOND = ddOND.values

# Load smoothed isobaths.
d = np.load("isobath_100m.npz")
xi100, yi100 = d["xi_ext"], d["yi_ext"]

d = np.load("isobath_200m.npz")
xi, yi = d["xi_ext"], d["yi_ext"]

# Load Wei et al. (2024) sections.
xseg0, yseg0 = np.loadtxt("segment_latitudes_wei_etal2024.txt", unpack=True)
fseg = [near2(xi, yi, xseg0[i], yseg0[i]) for i in range(len(yseg0))]
xseg, yseg = xi[fseg], yi[fseg]

figsize = (12, 8)
cmap = plt.cm.viridis
cmap.set_extremes(bad="w", under="w")

fig = plt.figure(figsize=figsize)
ax11 = fig.add_subplot(221, projection=proj)
ax12 = fig.add_subplot(222, projection=proj, sharex=ax11, sharey=ax11)
ax13 = fig.add_subplot(223, projection=proj, sharex=ax11, sharey=ax11)
ax14 = fig.add_subplot(224, projection=proj, sharex=ax11, sharey=ax11)
ax11 = bmap(fig, ax11, bb, proj=proj, land=LAND, TOP_LABELS=True)
ax12 = bmap(fig, ax12, bb, proj=proj, land=LAND, TOP_LABELS=True)
ax13 = bmap(fig, ax13, bb, proj=proj, land=LAND, TOP_LABELS=True)
ax14 = bmap(fig, ax14, bb, proj=proj, land=LAND, TOP_LABELS=True)
axs = [ax11, ax12, ax13, ax14]

cmap = plt.cm.viridis
cmap.set_extremes(bad="w", under="w")
cbax = [0.625, 0.15, 0.35, 0.03]
cbname = r"# obs."

xtxt, ytxt = 0.7, 0.035
fig.subplots_adjust(hspace=0, wspace=0)

ssn = "JFM"
cs11 = ax11.pcolormesh(x, y, ddJFM, vmin=nobsmin, vmax=nobsmax, cmap=cmap)
ax11.text(xtxt, ytxt, ssn, fontsize=18, fontweight="black", zorder=9, transform=ax11.transAxes)

ssn = "AMJ"
cs12 = ax12.pcolormesh(x, y, ddAMJ, vmin=nobsmin, vmax=nobsmax, cmap=cmap)
ax12.text(xtxt, ytxt, ssn, fontsize=18, fontweight="black", zorder=9, transform=ax12.transAxes)

ssn = "JAS"
cs13 = ax13.pcolormesh(x, y, ddJAS, vmin=nobsmin, vmax=nobsmax, cmap=cmap)
ax13.text(xtxt, ytxt, ssn, fontsize=18, fontweight="black", zorder=9, transform=ax13.transAxes)

ssn = "OND"
cs14 = ax14.pcolormesh(x, y, ddOND, vmin=nobsmin, vmax=nobsmax, cmap=cmap)
ax14.text(xtxt, ytxt, ssn, fontsize=18, fontweight="black", zorder=9, transform=ax14.transAxes)

_ = [axi.plot(xi, yi, color="r", linewidth=3, zorder=20) for axi in axs]
_ = [axi.plot(xseg, yseg, linestyle="none", marker="o", ms=6, mfc="w", mec="r", zorder=21) for axi in axs]

for csi, axi in zip([cs11, cs12, cs13, cs14], [ax11, ax12, ax13, ax14]):
    cbi = fig.colorbar(csi, cax=axi.inset_axes(cbax, zorder=10), orientation="horizontal", extend="max")
    cbi.ax.xaxis.set_ticks_position("top"); cbi.ax.xaxis.set_label_position("top"); cbi.ax.xaxis.set_ticks(np.int32(np.linspace(nobsmin, nobsmax, 4)))
    cbi.set_label(cbname, fontsize=15)

fig.tight_layout()

fig.savefig("../plot_figs/figS1.png", bbox_inches="tight")
