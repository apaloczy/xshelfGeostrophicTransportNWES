# Interpolate seasonal/decadal WOA23 fields to shelfbreak isobath.
import numpy as np
import matplotlib.pyplot as plt
from gsw import SA_from_SP, CT_from_t, sigma0
from cmocean.cm import thermal, haline, dense, balance, curl
from xarray import open_dataset, DataArray, Variable, Coordinates, Dataset
from gsw import distance, grav, p_from_z, geo_strf_dyn_height
from gsw import f as fcor
from glob import glob

from scipy.spatial import Delaunay
from scipy.interpolate import LinearNDInterpolator, CloughTocher2DInterpolator

from scipy.interpolate import interp1d


def compute_weights(x, y, F, verbose=True):
    nk = F.sizes["depth"]
    F = F.values
    wk = []
    for k in range(nk):
        if verbose:
            print("Computing interpolation weights for level %d / %d"%(k + 1, nk))
        try:
            Fk = F[k, :, :]
            fgud = np.isfinite(Fk)
            xk, yk, Fkn = x[fgud], y[fgud], Fk[fgud]
            w = Delaunay(np.array([yk, xk]).T)
        except:
            if verbose:
                print("Warning: Skipped level %d/%d"%(k+1, nk))
            w = None

        wk.append(w)

    return wk


def interp2trk(F, xi, yi, kweights, zbot, method="cubic", verbose=True):
    if method=="linear":
        Interpolator = LinearNDInterpolator
    elif method=="cubic":
        Interpolator = CloughTocher2DInterpolator

    nk = F.sizes["depth"]
    z = F["depth"].values
    assert nk==len(kweights), "Different number of depth levels and weights."
    F = F.values
    nx = xi.size
    Fi = np.empty((nk, nx))*np.nan
    for k in range(nk):
        if verbose:
            print("Interpolating level %d / %d."%(k + 1, nk))
        Fk = F[k, :, :]
        fgud = np.isfinite(Fk)
        if fgud.any():
            try:
                Fi[k, :] = Interpolator(kweights[k], Fk[fgud])(yi, xi)
            except ValueError:
                if verbose:
                    print("ValueError - different number of points. Skipping.")
                pass
        else:
            if verbose:
                print("No data at this depth. Skipping.")

    # Mask depths below bottom at each stations.
    for i in range(nx):
        fmsk = z>zbot[i]
        Fi[fmsk, i] = np.nan

    return Fi


def near(x, x0):
    return np.nanargmin(np.abs(x - x0))


def near2(x, y, x0, y0):
    dr2 = (x - x0)**2 + (y - y0)**2
    return int(np.nanargmin(dr2))


#---
plt.close("all")
iso = 200 # [m]

rho0 = 1027.5 # [kg/m3] - Reference density at 0 dbar.

method = "cubic"
# method = "linear"

xmin, xmax = -16.5, 10.5
ymin, ymax = 44, 63

Tmi, Tma = 6, 17
Smi, Sma = 27.6, 35.9
Dmi, Dma = 26.0, 27.5
Nmi, Nma = 2, 15
uxim = 2.5
Nmi, Nma = 2, 15
uxNim = 0.3 # [(mm/kg) * (m/s)]

head_woa23 = "../data/woa23/"
fname_adt_decavg = "decadal_altimetry_adt_%dmisob.nc"%iso

fnames = glob(head_woa23 + "woa23T-????-????.nc")
fnames.sort()

dsT = open_dataset(fnames[0], decode_times=False)
x, y = np.meshgrid(dsT["lon"].values, dsT["lat"].values)
x0, y0 = x[0, :], y[:, 0]

#---

# Read isobath contour's coordinates.
fisobath = "isobath_%dm.npz"%iso
disob = np.load(fisobath)

# Interpolate at midpoints to compare with EN4 climatology at xim, yim points.
xip, yip = disob["xip"], disob["yip"] # [km]
latm = yip.mean()
xi = 0.5*(xip[1:] + xip[:-1])
yi = 0.5*(yip[1:] + yip[:-1])
di = np.append(0, np.cumsum(distance(xi, yi)))*1e-3 # [km]

#---

xseg0, yseg0 = np.loadtxt("segment_latitudes_wei_etal2024.txt", unpack=True)
fseg = [near2(xi, yi, xseg0[i], yseg0[i]) for i in range(len(yseg0))]
xseg, yseg = xi[fseg], yi[fseg]
dibound = di[fseg]
txt_segs = ["AS", "CSS", "PSBS", "MS", "HS", "WSS", "NENSS"]

def annotate_segments(axs, dibound=dibound, txt_segs=txt_segs, ytxt=-190, lineplot=False, last_axs_txt_only=False):
    if lineplot:
        ytxt = -10
    for axi in axs:
        _ = [axi.axvline(x=dd, color="k", linestyle="--", linewidth=1.5, zorder=9999) for dd in dibound]
        for i in range(len(txt_segs)):
            xtxt = 0.5*(dibound[i] + dibound[i+1])
            if not last_axs_txt_only or (last_axs_txt_only and axi==axs[-1]):
                axi.text(xtxt, ytxt, txt_segs[i], ha="center", va="center", fontsize=12, color="k", bbox=dict(boxstyle="round,pad=0.3", fc="w", ec="k", lw=1))

nlevs = 75
nlevsc = 10
nlevsc2 = 30

Tlevs = np.linspace(Tmi, Tma, num=nlevs)
Slevs = np.linspace(Smi, Sma, num=nlevs)
Dlevs = np.linspace(Dmi, Dma, num=nlevs)
Nlevs = np.linspace(Nmi, Nma, num=nlevs)
Tlevsc = np.linspace(Tmi, Tma, num=nlevsc)
Slevsc = np.linspace(Smi, Sma, num=nlevsc2)
Dlevsc = np.linspace(Dmi, Dma, num=nlevsc2)
Nlevsc = np.linspace(Nmi, Nma, num=nlevsc)
uxilevs = np.linspace(-uxim, uxim, num=nlevs)
uxNlevs = np.linspace(-uxNim, uxNim, num=nlevs)

dyi = np.gradient(di*1e3) # [m]

adt_decavg = open_dataset(fname_adt_decavg)["adt"]

# WOA's vertical grid changes dz from 5 to 15 to 25 abruptly at around 100 m. Interpolate to regular levels to avoid an artificial jump in the z-integral of rho.
dz = 5.0 # [m]

zclip = 3000
figsize = (12, 8)
count = 1
zboti = np.array([iso]*xi.size)
Tis, Sis, Dis = [], [], []
tstart, tend = [], []
for nf, f in enumerate(fnames):
    fstr = f.split("woa23T-")[-1].strip(".nc")[:4]
    fstr = fstr + "-" + str(int(fstr) + 9)
    print("")
    print(fstr)
    dsT = open_dataset(f, decode_times=False)
    dsS = open_dataset(f.replace("woa23T", "woa23S"), decode_times=False)
    if nf==0:
        z = -dsT["depth"].values
        zint = np.arange(0, -iso - dz, -dz)
        p = p_from_z(z, latm)[:, np.newaxis] # [dbar]

    # Interpolate WOA TS fields and OM standard errors.
    T, Tse = dsT["t_an"], dsT["t_sea"]
    S, Sse = dsS["s_an"], dsS["s_sea"]

    kweightsT = compute_weights(x, y, T, verbose=False)
    kweightsS = compute_weights(x, y, S, verbose=False)
    kweightsTse = compute_weights(x, y, Tse, verbose=False)
    kweightsSse = compute_weights(x, y, Sse, verbose=False)

    Ti = interp2trk(T, xi, yi, kweightsT, zboti, method=method, verbose=False)
    Si = interp2trk(S, xi, yi, kweightsS, zboti, method=method, verbose=False)

    # Derive CT, SA and potential density.
    SAi = SA_from_SP(Si, p, xi, yi)
    CTi = CT_from_t(SAi, Ti, p)

    Ti, Si = CTi.copy(), SAi.copy()
    Di = sigma0(SAi, CTi)

    nzi, nxi = zint.size, di.size
    Di_zint = np.empty((nzi, nxi))*np.nan
    Ti_zint = Di_zint.copy()
    Si_zint = Di_zint.copy()
    for n in range(nxi):
        Tin = Ti[:, n]
        fg = np.isfinite(Tin)
        Ti_zint_aux = interp1d(z[fg], Tin[fg], kind=method, bounds_error=False, fill_value=np.nan)(zint)
        Sin = Si[:, n]
        fg = np.isfinite(Sin)
        Si_zint_aux = interp1d(z[fg], Sin[fg], kind=method, bounds_error=False, fill_value=np.nan)(zint)
        Din = Di[:, n]
        fg = np.isfinite(Din)
        Di_zint_aux = interp1d(z[fg], Din[fg], kind=method, bounds_error=False, fill_value=np.nan)(zint)

        Ti_zint[:, n] = Ti_zint_aux
        Si_zint[:, n] = Si_zint_aux
        Di_zint[:, n] = Di_zint_aux

    Di = Di_zint.copy()
    Ti = Ti_zint.copy()
    Si = Si_zint.copy()

    fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, sharex=True, sharey=True, figsize=(8, 10))

    cs1 = ax1.pcolormesh(di, zint, Ti, vmin=Tmi, vmax=Tma, cmap=thermal)
    cs2 = ax2.pcolormesh(di, zint, Si, vmin=Smi, vmax=Sma, cmap=haline)
    cs3 = ax3.pcolormesh(di, zint, Di, vmin=Dmi, vmax=Dma, cmap=dense)
    cc1 = ax1.contour(di, zint, Ti, levels=Tlevsc, colors="k")
    cc2 = ax2.contour(di, zint, Si, levels=Slevsc, colors="k")
    cc3 = ax3.contour(di, zint, Di, levels=Dlevsc, colors="k")
    ax1.clabel(cc1, fmt="%1.1f")
    ax2.clabel(cc2, fmt="%1.1f", inline_spacing=0)
    ax3.clabel(cc3, fmt="%1.1f")

    fig.colorbar(cs1, label=r"$\Theta$ [$^o$C]")
    fig.colorbar(cs2, label=r"$S_A$ [g/kg]")
    fig.colorbar(cs3, label=r"$\sigma_0$ [kg/m$^3$]")

    # Add segment locations.
    axs = (ax1, ax2, ax3)
    annotate_segments(axs)

    ax3.set_xlabel("Distance [km]", fontsize=15)
    ax2.set_ylabel("Depth [m]", fontsize=15, y=0)
    ax1.set_ylim(-iso, 0)
    ax1u = ax1.twiny()
    ax1u.xaxis.set_tick_params(top=True, labeltop=True)
    ax1u.plot(yi, yi*0)
    ax1u.set_xlabel("Latitude [degrees]", fontsize=15)
    fig.subplots_adjust(hspace=0.03)
    fig.suptitle(fstr, fontsize=15, x=0.45, y=0.95, ha="center")
    # fig.savefig("alongisobTS_%s_%dm.png"%(fstr, iso), bbox_inches="tight", dpi=150)

    Tis.append(Ti)
    Sis.append(Si)

    tstart.append(int(f.split("-")[-2]))
    tend.append(int(f.split("-")[-1].strip(".nc")))

# Put variables in an xarray dataset.
xin = Variable("x", di, attrs=dict(long_name="Along-isobath distance", units="km"))
zin = Variable("z", zint, attrs=dict(long_name="Depth", units="m"))
lonin, latin = Variable("x", xi, attrs=dict(long_name="Longitude", units="Degrees east")), Variable("x", yi, attrs=dict(long_name="Latitude", units="Degrees north"))
t = Variable("t", tstart, attrs=dict(long_name="Start time for averaging"))
tend = Variable("t", tend, attrs=dict(long_name="End time for averaging"))

Dis = sigma0(np.array(Sis), np.array(Tis)) # Recalculate density from interpolated SA, CT instead of using interpolated density (density in loop just for plotting)

dims = ("t", "z", "x")
coords = Coordinates(dict(x=xin, lon=lonin, lat=latin, z=zin, t=tstart, tend=tend))
Tis = DataArray(np.array(Tis), dims=dims, coords=coords, attrs=dict(long_name="Conservative Temperature", units="Degrees Celsius"))
Sis = DataArray(np.array(Sis), dims=dims, coords=coords, attrs=dict(long_name="Absolute Salinity", units="g/kg"))
Dis = DataArray(np.array(Dis), dims=dims, coords=coords, attrs=dict(long_name="Potential density referenced to 0 dbar", units="kg/m3"))

dsout = Dataset(dict(CT=Tis, SA=Sis, rho=Dis), coords=coords).sel(z=slice(0, -iso))

# Calculate ug from geopotential anomaly.
xx, zz, lat = dsout.x.values, dsout.z.values, dsout.lat.values
g = grav(lat, 0)  # [m/s^2]
f = fcor(lat)  # [1/s]

PLOT_Phi_ug = True
tls = coords["t"].values
uxis = []
for tl in tls:
    SA, CT = dsout["SA"].sel(t=tl).values, dsout["CT"].sel(t=tl).values

    p = p_from_z(zz[:, np.newaxis], lat)
    Phi = geo_strf_dyn_height(SA, CT, p, p_ref=0)
    ug = np.gradient(Phi, xx*1e3, axis=1)/f  # [m/s] geostrophic velocity *relative to sea surface*.
    uxis.append(ug)

    if PLOT_Phi_ug:
        fig, (ax1, ax2) = plt.subplots(nrows=2, sharex=True, sharey=True, figsize=(12, 8))
        cs1 = ax1.pcolormesh(xx, zz, Phi, cmap=balance, shading="auto")
        cs2 = ax2.pcolormesh(xx, zz, ug*1e2, cmap=balance, shading="auto")
        fig.colorbar(cs1, ax=ax1, label=r"$\Phi(p, 0)$ [m2/s2]")
        fig.colorbar(cs2, ax=ax2, label=r"$u_g^\perp(p, 0)$ [cm/s]")
        ax1.set_title(r"Geopotential anomaly, $\Phi + \Phi_\eta$ (%s)"%tl)
        ax1.set_title(r"Cross-isobath geostrophic velocity relative to surface, $u_g^\perp - u_{g0}^\perp$ (%s)"%tl)
        ax2.set_xlabel(r"Along-isobath distance [km]")
        ax1.set_ylabel(r"$z$ [m]", y=0)
        plt.show()

uxis = DataArray(np.array(uxis)*1e2, dims=dims, coords=coords, attrs=dict(long_name="Geostrophic cross-isobath velocity relative to the surface", units="cm/s"))
dsout["ug"] = uxis

kwencode = dict(zlib=True, complevel=5)
dsout.to_netcdf("alongisob_decadal_WOA23_%dm.nc"%iso, encoding=dict(CT=kwencode, SA=kwencode, rho=kwencode, ug=kwencode))

# Plot rho(y, z) and ux(y, z) side by side.
xtxt, ytxt = 1.08, 0.5
years = [1975, 1985, 1995, 2005, 2015, 2023]#2025]
nrows = 5
fig, ax = plt.subplots(nrows=nrows, ncols=2, sharex=True, sharey=True, figsize=(8, 12))
for i in range(nrows):
    axl, axr = ax[i]
    Dii, uxii = Dis[i], uxis[i]
    # csl = axl.contourf(di, zint, Dii, levels=Dlevs, cmap=dense)
    csl = axl.pcolormesh(di, zint, Dii, vmin=Dlevs[0], vmax=Dlevs[-1], cmap=dense)
    cc = axl.contour(di, zint, Dii, levels=Dlevsc, colors="k")
    axl.clabel(cc, fmt="%1.1f")
    # csr = axr.contourf(di, zint, uxii, levels=uxilevs, cmap=balance)
    csr = axr.pcolormesh(di, zint, uxii, vmin=uxilevs[0], vmax=uxilevs[-1], cmap=balance)
    axr.contour(di, zint, uxii, levels=[0], colors="k")
    # csr = axr.pcolormesh(di, zint, uxi, vmin=uxilevs[0], vmax=uxilevs[-1], cmap=balance)

    yl, yr = years[i+1] - 1, years[i]
    tstr = str(yr) + " - " + str(yl)
    axl.text(xtxt, ytxt, tstr, fontsize=15, transform=axl.transAxes, ha="center", va="center", rotation=-90)
    if i==2:
        axl.set_ylabel("Depth [m]", fontsize=15)

    del Dii, uxii
# del Tis, Sis, Dis, uxis

# Add segment locations.
axs = (axl, axr)
annotate_segments(ax.ravel())

numtks_rho = 5
numtks_u = 5
axl.set_ylim(-iso, 0)
axl.set_xlabel("Distance [km]", fontsize=15, x=1)
axl, axr = ax[0]
cax = axl.inset_axes([0, 1.05, 1, 0.05])
cb = fig.colorbar(csl, cax=cax, orientation="horizontal", label=r"$\sigma_0$ [kg/m$^3$]")
cb.ax.xaxis.set_ticks_position("top"); cb.ax.xaxis.set_label_position("top")
cb.ax.xaxis.set_ticks(np.linspace(Dlevs[0], Dlevs[-1], num=numtks_rho))
cax = axr.inset_axes([0, 1.05, 1, 0.05])
cb = fig.colorbar(csr, cax=cax, orientation="horizontal", label=r"u$_{g0}^\perp$ [cm/s]")
cb.ax.xaxis.set_ticks_position("top"); cb.ax.xaxis.set_label_position("top")
cb.ax.xaxis.set_ticks(np.linspace(uxilevs[0], uxilevs[-1], num=numtks_u))
# fig.savefig("alongisobTS_rho_ux_decadal_%dm.png"%iso, bbox_inches="tight", dpi=150)








# Plot decade-to-decade differences in rho(y, z) and ux(y, z).

dDmax = 0.18
duximax = 1
dDlevs = np.linspace(-dDmax, dDmax, num=nlevs)
duxilevs = np.linspace(-duximax, duximax, num=nlevs)
nrows = 4
fig, ax = plt.subplots(nrows=nrows, ncols=2, sharex=True, sharey=True, figsize=(8, 12))
for i in range(nrows):
    axl, axr = ax[i]
    Disp, Dii = Dis[i+1], Dis[i]
    uxisp, uxii = uxis[i+1], uxis[i]
    dDii = Disp - Dii
    duxii = uxisp - uxii
    # csl = axl.contourf(di, zint, dDii, levels=dDlevs, cmap=balance)
    csl = axl.pcolormesh(di, zint, dDii, vmin=dDlevs[0], vmax=dDlevs[-1], cmap=balance)
    cc = axl.contour(di, zint, Disp, levels=Dlevsc, colors="k")
    axl.clabel(cc, fmt="%1.1f")
    # csr = axr.contourf(di, zint, duxii, levels=duxilevs, cmap=curl)
    csr = axr.pcolormesh(di, zint, duxii, vmin=duxilevs[0], vmax=duxilevs[-1], cmap=curl)
    axr.contour(di, zint, uxii, levels=[0], colors="k")
    # csr = axr.pcolormesh(di, zint, uxi, vmin=duxilevs[0], vmax=duxilevs[-1], cmap=balance)

    decr = str(years[i+2]-1)[2:] + "-" + str(years[i+1])[2:]
    decl = str(years[i+1]-1)[2:] + "-" + str(years[i])[2:]
    tstr = "(" + str(decr) + ") - " + "(" + str(decl) + ")"
    axl.text(xtxt, ytxt, tstr, fontsize=13, transform=axl.transAxes, ha="center", va="center", rotation=-90)
    if i==2:
        axl.set_ylabel("Depth [m]", fontsize=15, y=1)
    del Disp, Dii, dDii, uxisp, uxii, duxii

# Add segment locations.
# axs = (axl, axr)
annotate_segments(ax.ravel())

axl.set_ylim(-iso, 0)
axl.set_xlabel("Distance [km]", fontsize=15, x=1)
axl, axr = ax[0]
cax = axl.inset_axes([0, 1.05, 1, 0.05])
cb = fig.colorbar(csl, cax=cax, orientation="horizontal", label=r"$\sigma_0$ change [kg/m$^3$]")
cb.ax.xaxis.set_ticks_position("top"); cb.ax.xaxis.set_label_position("top")
cb.ax.xaxis.set_ticks(np.linspace(dDlevs[0], dDlevs[-1], num=numtks_rho))
cax = axr.inset_axes([0, 1.05, 1, 0.05])
cb = fig.colorbar(csr, cax=cax, orientation="horizontal", label=r"u$_{g0}^\perp$ change [cm/s]")
cb.ax.xaxis.set_ticks_position("top"); cb.ax.xaxis.set_label_position("top")
cb.ax.xaxis.set_ticks(np.linspace(duxilevs[0], duxilevs[-1], num=numtks_u))
# fig.savefig("alongisobTS_rho_ux_decadaldiff_%dm.png"%iso, bbox_inches="tight", dpi=150)


##### Seasonal plots with nitrate transports.
names_ssn = ["annual", "winter", "spring", "summer", "autumn"]
fnames_seasonal = [head_woa23 + "woa23T-%s.nc"%n for n in names_ssn]
fname_adt_ssnavg = "seasonal_altimetry_adt_200misob.nc"

adt_ssnavg = open_dataset(fname_adt_ssnavg)["adt"]

Nbbox = dict(lon=slice(xmin, xmax), lat=slice(ymin, ymax))

di = np.append(0, np.cumsum(distance(xi, yi)))*1e-3 # [km]
dsN = open_dataset(fnames_seasonal[0].replace("T", "N"), decode_times=False).sel(Nbbox).squeeze()
xN, yN = np.meshgrid(dsN["lon"].values, dsN["lat"].values)
ssndict = dict(winter="JFM", spring="AMJ", summer="JAS", autumn="OND")
Dis, uxis, Nis, uxNis = [], [], [], []
Tis, Sis, Dis, Nis = [], [], [], []
Nseis = []
ssns = []
for nf, f in enumerate(fnames_seasonal):
    fstr = f.split("woa23T-")[-1].split(".nc")[0]
    print("")
    print(fstr)
    dsT = open_dataset(f, decode_times=False)
    dsS = open_dataset(f.replace("woa23T", "woa23S"), decode_times=False)
    dsN = open_dataset(f.replace("woa23T", "woa23N"), decode_times=False).sel(Nbbox).squeeze()

    # Interpolate WOA TS fields and OM standard errors.
    T = dsT["t_an"]
    S = dsS["s_an"]
    N = dsN["n_an"]
    Nse = dsN["n_sea"]

    kweightsT = compute_weights(x, y, T, verbose=False)
    kweightsS = compute_weights(x, y, S, verbose=False)
    kweightsN = compute_weights(xN, yN, N, verbose=False)
    kweightsNse = compute_weights(xN, yN, Nse, verbose=False)

    Ti = interp2trk(T, xi, yi, kweightsT, zboti, method=method, verbose=False)
    Si = interp2trk(S, xi, yi, kweightsS, zboti, method=method, verbose=False)
    Ni = interp2trk(N, xi, yi, kweightsN, zboti, method=method, verbose=False)
    Nsei = interp2trk(Nse, xi, yi, kweightsNse, zboti, method=method, verbose=False)

    zN = -N["depth"].values # 102 z levels for annual and 43 for seasonal.
    if nf==0:
        z = -T["depth"].values
        p = p_from_z(z, latm)[:, np.newaxis] # [dbar]

        print("dz = %1.1f m"%dz)
        zint = np.arange(0, -iso - dz, -dz)
        nzi, nxi = zint.size, di.size

    # Derive CT, SA and potential density.
    SAi = SA_from_SP(Si, p, xi, yi)
    CTi = CT_from_t(SAi, Ti, p)
    Di = sigma0(SAi, CTi)
    Ti, Si = CTi.copy(), SAi.copy()

    Di_zint = np.empty((nzi, nxi))*np.nan
    Ni_zint = Di_zint.copy()
    Nsein_zint = Di_zint.copy()
    for n in range(nxi):
        Tin = Ti[:, n]
        Sin = Si[:, n]
        Din = Di[:, n]
        Nin = Ni[:, n]
        Nsein = Nsei[:, n]

        fg = np.isfinite(Tin)
        if not fg.any():
            Ti_zint_aux = zint*np.nan
        else:
            Ti_zint_aux = interp1d(z[fg], Tin[fg], kind=method, bounds_error=False, fill_value=np.nan)(zint)

        fg = np.isfinite(Sin)
        if not fg.any():
            Si_zint_aux = zint*np.nan
        else:
            Si_zint_aux = interp1d(z[fg], Sin[fg], kind=method, bounds_error=False, fill_value=np.nan)(zint)

        fg = np.isfinite(Din)
        if not fg.any():
            Di_zint_aux = zint*np.nan
        else:
            Di_zint_aux = interp1d(z[fg], Din[fg], kind=method, bounds_error=False, fill_value=np.nan)(zint)

        fg = np.isfinite(Nin)
        if not fg.any():
            Ni_zint_aux = zint*np.nan
        else:
            Ni_zint_aux = interp1d(zN[fg], Nin[fg], kind=method, bounds_error=False, fill_value=np.nan)(zint)

        fg = np.isfinite(Nsein)
        if not fg.any():
            Nise_zint_aux = zint*np.nan
        else:
            Nise_zint_aux = interp1d(zN[fg], Nsein[fg], kind=method, bounds_error=False, fill_value=np.nan)(zint)

        Ti_zint[:, n] = Ti_zint_aux
        Si_zint[:, n] = Si_zint_aux
        Di_zint[:, n] = Di_zint_aux
        Ni_zint[:, n] = Ni_zint_aux
        Nsein_zint[:, n] = Nise_zint_aux

    Di = Di_zint.copy()
    Ti = Ti_zint.copy()
    Si = Si_zint.copy()
    Ni = Ni_zint.copy()
    Nsein = Nsein_zint.copy()

    Tis.append(Ti)
    Sis.append(Si)
    Dis.append(Di)
    Nis.append(Ni*rho0*1e-3) #From umol/kg to mmol/m3. Annual nitrate has 102 z levels, the seasonal have only 43 z levels.
    Nseis.append(Nsein*rho0*1e-3)

    if fstr=="annual":
        ssnstr = "annual"
    else:
        ssnstr = ssndict[fstr]
    ssns.append(ssnstr)


# Put variables in an xarray dataset.
xin = Variable("x", di, attrs=dict(long_name="Along-isobath distance", units="km"))
zin = Variable("z", zint, attrs=dict(long_name="Depth", units="m"))
lonin, latin = Variable("x", xi, attrs=dict(long_name="Longitude", units="Degrees east")), Variable("x", yi, attrs=dict(long_name="Latitude", units="Degrees north"))
t = Variable("t", ssns, attrs=dict(long_name="Season"))

dims = ("t", "z", "x")
coords = Coordinates(dict(x=xin, lon=lonin, lat=latin, z=zin, t=t))
Tis = DataArray(np.array(Tis), dims=dims, coords=coords, attrs=dict(long_name="Conservative Temperature", units="Degrees Celsius"))
Sis = DataArray(np.array(Sis), dims=dims, coords=coords, attrs=dict(long_name="Absolute Salinity", units="g/kg"))
Dis = DataArray(np.array(Dis), dims=dims, coords=coords, attrs=dict(long_name="Potential density referenced to 0 dbar", units="kg/m3"))
Nis = DataArray(np.array(Nis), dims=dims, coords=coords, attrs=dict(long_name="Nitrate concentration", units="mmol/m3"))
Nseis = DataArray(np.array(Nseis), dims=dims, coords=coords, attrs=dict(long_name="Nitrate concentration standard error", units="mmol/m3"))

dsout = Dataset(dict(CT=Tis, SA=Sis, rho=Dis, N=Nis, Nse=Nseis), coords=coords).sel(z=slice(0, -iso))

# Calculate surface-referenced ug using TEOS-10's geopotential anomaly equation.
assert dsout.sizes["x"] == adt_ssnavg.sizes["x"], "Different along-isobath grid points"
xx, zz, lat = dsout.x.values, dsout.z.values, dsout.lat.values
g = grav(lat, 0)  # [m/s^2]
f = fcor(lat)     # [1/s]

# plt.close("all")

ssns = coords["t"].values
uxis = []
for ssn in ssns:
    print("Season: ", ssn)
    if ssn == "annual":
        adt = adt_ssnavg.mean("season").values
    else:
        adt = adt_ssnavg.sel(season=ssn).values
    SA, CT = dsout["SA"].sel(t=ssn).values, dsout["CT"].sel(t=ssn).values

    p = p_from_z(zz[:, np.newaxis], lat)

    Phi = geo_strf_dyn_height(SA, CT, p, p_ref=0)
    ugbc = np.gradient(Phi, xx*1e3, axis=1)/f  # [m/s]
    ug0 =  - g*np.gradient(adt, xx*1e3)/f # Surface geostrophic velocity [m/s].
    ug = ugbc + ug0 # Absolute geostrophic velocity [m/s].
    uxis.append(ug)

    if PLOT_Phi_ug:
        fig, (ax1, ax2) = plt.subplots(nrows=2, sharex=True, sharey=True, figsize=(12, 8))
        cs1 = ax1.pcolormesh(xx, zz, Phi, cmap=balance, shading="auto")
        cs2 = ax2.pcolormesh(xx, zz, ug*1e2, cmap=balance, shading="auto")
        fig.colorbar(cs1, ax=ax1, label=r"$\Phi(p, 0)$ [m2/s2]")
        fig.colorbar(cs2, ax=ax2, label=r"$u_g^\perp(p, 0)$ [cm/s]")
        ax1.set_title(r"Geopotential anomaly, $\Phi + \Phi_\eta$ (%s)"%ssn)
        ax1.set_title(r"Cross-isobath absolute geostrophic velocity, $u_g^\perp$ (%s)"%ssn)
        ax2.set_xlabel(r"Along-isobath distance [km]")
        ax1.set_ylabel(r"$z$ [m]", y=0)
        plt.show()

uxNis = np.array(uxis)*dsout["N"].values # m/s * mmol/m3 = [mmol/m2/s]

uxis = DataArray(np.array(uxis)*1e2, dims=dims, coords=coords, attrs=dict(long_name="Absolute geostrophic cross-isobath velocity", units="cm/s"))
uxNis = DataArray(uxNis, dims=dims, coords=coords, attrs=dict(long_name="Geostrophic cross-isobath nitrate transport per unit along-isobath distance", units="mmol/(m2*s)"))
dsout["ug"] = uxis
dsout["ugN"] = uxNis

kwencode = dict(zlib=True, complevel=5)
dsout.to_netcdf("alongisob_seasonal_WOA23_%dm.nc"%iso, encoding=dict(CT=kwencode, SA=kwencode, rho=kwencode, ug=kwencode, ugN=kwencode))

uxis, Nis, uxNis = dsout["ug"].values, dsout["N"].values, dsout["ugN"].values
di = dsout["x"].values
xi, yi = dsout["lon"].values, dsout["lat"].values

fseg = [near2(xi, yi, xseg0[i], yseg0[i]) for i in range(len(yseg0))]
xseg, yseg = xi[fseg], yi[fseg]
dibound = di[fseg]


# Integrate along the isobath to estimate transports.

uxii_zavgs, Nii_zavgs, uxNii_zints = [], [], []
for n in range(uxis.shape[0]):
    uxii, Nii, uxNii = uxis[n], Nis[n], uxNis[n]

    Niim = 0.5*(Nii[1:, :] + Nii[:-1, :])

    uxii_zavg = np.nansum(uxii*dz, axis=0)/iso # Depth-averaged cross-isobath velocity [m/s].
    Nii_zavg = np.nansum(Niim*dz, axis=0)/iso # Depth-averaged nitrate [mmol/m3].
    uxNii_zint = np.nansum(uxNii*dz, axis=0) # [mmol/m/s] Depth-integrated geostrophic nitrate transport.

    uxii_zavgs.append(uxii_zavg)
    Nii_zavgs.append(Nii_zavg)
    uxNii_zints.append(uxNii_zint)


fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, sharex=True, figsize=(12, 9))
for i in range(5):
    ax1.plot(di, uxii_zavgs[i], label=ssns[i])
    ax2.plot(di, Nii_zavgs[i], label=ssns[i])
    ax3.plot(di, uxNii_zints[i], label=ssns[i])

axs = (ax1, ax2, ax3)
# annotate_segments(axs, dibound=dibound, lineplot=False, ytxt=5, last_axs_txt_only=True)

for axi in axs:
    _ = [axi.axvline(x=dd, color="k", linestyle="--", linewidth=1.5, zorder=9999) for dd in dibound]
    for i in range(len(txt_segs)):
        xtxt = 0.5*(dibound[i] + dibound[i+1])
        if axi==axs[1]:
            ytxt = 13
        elif axi==axs[2]:
            ytxt = 90
        else:
            ytxt=7
        axi.text(xtxt, ytxt, txt_segs[i], ha="center", va="center", fontsize=12, color="k", bbox=dict(boxstyle="round,pad=0.3", fc="w", ec="k", lw=1))

ax1.set_xlim(0, di[-1])
ax1.axhline(color="gray"); ax3.axhline(color="gray")
ax2.legend(ncols=5)
ax3.set_xlabel(r"Distance [km]")
ax1.set_ylabel(r"Depth-avg. $u_g^\perp$ [cm/s]", fontsize=11)
ax2.set_ylabel(r"Depth-avg. nitrate [mmol m$^{-3}$]", fontsize=11)
ax3.set_ylabel(r"Depth-int. $u_g^\perp$ [mmol m$^{-1}$ s$^{-1}$]", fontsize=11)
# fig.savefig("seasonal_nitrate_flux_depthavg.png", bbox_inches="tight", dpi=200)
# plt.close("all")


for i in range(len(uxNii_zints)):
    uxNii = uxNii_zints[i]
    print("============")
    print(ssns[i])
    print("============")
    for n in range(len(fseg) - 1):
        seg = txt_segs[n]
        fl, fr = fseg[n], fseg[n+1]
        dx = np.diff(di[fl:fr])*1e3 # [m]
        L = dx.sum() # [m]
        uxNii_aux = uxNii[fl:fr]
        uxNii_aux = 0.5*(uxNii_aux[1:] + uxNii_aux[:-1])
        uxNii_segavg = np.nansum(uxNii_aux*dx)/L # mmol/m/s
        print("Season: " + ssns[i] + " - Segment: " + seg + " --> ugN = " + "%1.1f"%uxNii_segavg + " mmol/m/s")
