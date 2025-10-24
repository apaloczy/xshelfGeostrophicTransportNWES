import numpy as np
import matplotlib.pyplot as plt
from xarray import open_dataset, Variable, DataArray, concat, merge
from pandas import Timestamp
from pygeodesy.sphericalNvector import LatLon
from gsw import distance, grav, p_from_z, geo_strf_dyn_height, sigma0
from gsw import f as fcor
from cmocean.cm import balance

import warnings
warnings.filterwarnings("ignore", message="Degrees of freedom <= 0 for slice.")

import cartopy.crs as ccrs
from cartopy.feature import LAND


def fmt_isobath(cs, fontsize=8, fmt='%g', inline=True, inline_spacing=7, manual=False, **kw):
	isobstrH = plt.clabel(cs, fontsize=fontsize, fmt=fmt, inline=inline, inline_spacing=inline_spacing, manual=manual, **kw)
	for ih in range(0, len(isobstrH)): # Appends 'm' for meters at the end of the label.
		isobstrh = isobstrH[ih]
		isobstr = isobstrh.get_text()
		isobstr = isobstr.replace('-','') + ' m'
		isobstrh.set_text(isobstr)


proj = ccrs.PlateCarree()
def bmap(fig, ax, bb, proj=proj, xticks=None, yticks=None, dlon=5, dlat=5, land=None, coastlines=True):
    ax.set_extent(bb, proj)

    if land is not None:
        ax.add_feature(land, color="gray", zorder=999)

    if coastlines:
        ax.coastlines()

    if not xticks:
        xticks = np.arange(bb[0], bb[1], dlon)
    if not yticks:
        yticks = np.arange(bb[2], bb[3], dlat)

    ax.gridlines(draw_labels=True)

    return ax


def angle_isobath(xiso, yiso):
    R = 6371000.0 # Mean radius of the earth in meters (6371 km), from gsw.constants.earth_radius.
    deg2rad = np.pi/180 # [rad/deg]

    # From the coordinates of the isobath, find the angle it forms with the
    # zonal axis, using points k+1 and k.
    shth = yiso.size - 1
    theta = np.zeros(shth)
    for k in range(shth):
        dyk = R*(yiso[k+1] - yiso[k])
        dxk = R*(xiso[k+1] - xiso[k])*np.cos(yiso[k]*deg2rad)
        theta[k] = np.arctan2(dyk, dxk)

    xisom = 0.5*(xiso[1:] + xiso[:-1])
    yisom = 0.5*(yiso[1:] + yiso[:-1])

    return xisom, yisom, theta/deg2rad


def get_xtrackline_from_angle(lon0, lat0, ang, L=(75, 25)):
    pm = LatLon(lat0, lon0)
    ang = (90 - ang) + 90
    angb = ang + 180

    km2m = 1e3
    Ll, Lr = L
    Ll, Lr = Ll*km2m, Lr*km2m
    dLl, dLr = Ll/2, Lr/2

    # Create perpendicular half-lines starting from the midpoint.
    p = []
    if dLl>0:
        nhl = int(Ll/dLl)
        Nl = range(1, nhl + 1)
        _ = [p.append(pm.destination(dLl*n, ang)) for n in Nl]

    p.reverse()
    p.append(pm)

    if dLr>0:
        nhr = int(Lr/dLr)
        Nr = range(1, nhr + 1)
        _ = [p.append(pm.destination(dLr*n, angb)) for n in Nr]

    lon = np.array([p.lon for p in p])
    lat = np.array([p.lat for p in p])
    lon = np.array([lon[0], lon[-1]])
    lat = np.array([lat[0], lat[-1]])

    return lon, lat


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
	for i in range(n+1):
		p2x,p2y = poly[i % n]
		if y > min(p1y,p2y):
			if y <= max(p1y,p2y):
				if x <= max(p1x,p2x):
					if p1y != p2y:
						xinters = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
					if p1x == p2x or x <= xinters:
						inside = not inside
		p1x, p1y = p2x, p2y

	return inside


def near(x, x0):
    return int(np.nanargmin(np.abs(x - x0)))


def near2(x, y, x0, y0):
    dr2 = (x - x0)**2 + (y - y0)**2
    return int(np.nanargmin(dr2))


def get_profiles_polygon(xim, yim, angs, k, maxdy, nk, maxdx, lonk0, latk0, dxi, asymmetrical=True):
    lonspl, lonspr, latspl, latspr = [], [], [], []
    lonsp0, latsp0 = [], []

    if asymmetrical:
        assert len(maxdx)==2, "maxdx needs to be a 2-element tuple/list/array if asymmetrical=True."
    else:
        assert type(maxdx)==int or type(maxdx)==float, "maxdx needs to be a number if asymmetrical=False."
        maxdx = (maxdx, maxdx)

    Lh_alongisob = maxdy/2 # [km]
    nhpts_alongisob = int(np.ceil(Lh_alongisob/dxi))

    # Add backward half of strip.
    kback = np.maximum(0, k - nhpts_alongisob)
    for kb in range(kback, k):
        xixlr, yixlr = get_xtrackline_from_angle(xim[kb], yim[kb], angs[kb], L=maxdx)
        lonpli, lonpri = xixlr[0], xixlr[-1]
        latpli, latpri = yixlr[0], yixlr[-1]
        lonspl.append(lonpli); lonspr.append(lonpri); latspl.append(latpli); latspr.append(latpri)
        lonsp0.append(xim[kb]); latsp0.append(yim[kb])

    # Add central point.
    xixlr, yixlr = get_xtrackline_from_angle(xik, yik, angi, L=maxdx)
    lonplc, lonprc = xixlr[0], xixlr[-1]
    latplc, latprc = yixlr[0], yixlr[-1]
    lonspl.append(lonplc); lonspr.append(lonprc); latspl.append(latplc); latspr.append(latprc)
    lonsp0.append(xik); latsp0.append(yik)

    # Add forward half of strip.
    kfwrd = np.minimum(k + nhpts_alongisob + 1, nk)
    for kf in range(k, kfwrd):
        xixlr, yixlr = get_xtrackline_from_angle(xim[kf], yim[kf], angs[kf], L=maxdx)
        lonpli, lonpri = xixlr[0], xixlr[-1]
        latpli, latpri = yixlr[0], yixlr[-1]
        lonspl.append(lonpli); lonspr.append(lonpri); latspl.append(latpli); latspr.append(latpri)
        lonsp0.append(xim[kf]); latsp0.append(yim[kf])

    lonspl, lonspr, latspl, latspr, lonsp0, latsp0 = map(np.array, (lonspl, lonspr, latspl, latspr, lonsp0, latsp0))
    lonsp = np.concatenate((lonspl, np.flipud(lonspr)))
    latsp = np.concatenate((latspl, np.flipud(latspr)))
    lonsp = np.append(lonsp, lonsp[0])
    latsp = np.append(latsp, latsp[0])

    poly = np.array([(x0, y0) for x0, y0 in zip(lonsp, latsp)])

    fnear = []
    for xii, yii in zip(lonk0, latk0):
        fnear.append(point_in_poly(xii, yii, poly))

    fnear = np.array(fnear)
    nnear = fnear.sum()

    return poly, fnear, nnear, lonsp, latsp, lonsp0, latsp0


def gaussr(r, rmax, a=0.01):
    Lr = rmax/np.sqrt(-2*np.log(a)) # Find a length scale that produces a weight of 'a' (e.g., a = 0.01) at the edge of the strip.

    return np.exp(-(r/Lr)**2/2)


# Visualize different weighting functions.
Lmax = [100, 200]
aa = [0.1, 0.05, 0.01]
fig, ax = plt.subplots()
for lmax in Lmax:
    dr = np.linspace(0, lmax, num=100)
    for a in aa:
        ax.plot(dr, gaussr(dr, lmax, a=a), label="$l$, $a$ = %d km, %.2f"%(lmax, a))
        ax.set_ylim(0, 1); ax.set_xlim(0, lmax)
        ax.set_xlabel("Distance [km]")
        ax.set_ylabel("Weight")
    ax.plot(dr, np.exp(-(dr/50)**2/2), "--", label="$l$ = %d km"%50)
ax.grid()
ax.legend()
plt.show(block=False)


lonmin, lonmax = -16.5, 10.5
latmin, latmax = 44, 63

#---
plt.close("all")

GLOBAL_DENSINV_TEST = True # Whether to apply the full-profile rejection criterion based on a large density inversion threshold.

density_inversion_threshold = -0.03 # [kg/m3]
dzi_density_inversion_test = 20 # [dbar]

SORT_DENS_PROFILES = True # Whether to apply density sorting in the vertical.
PLOT_RANDOM_PROFILES_SORT_DENS = False

EXTRAPNN_SURFACE = True # Whether to extrapolate shallowest valid data point to the surface.
ANISOTROPIC_BIN = True
DISTANCE_WEIGHTING = True

znearsurf = 10 # [m] Deepest depth from which to extrapolate to the surface if all depths above are NaN.
gaussr_frac = 0.1 # Fraction of maximum weight at the edge of the strip.
maxdy_initial = 100 # Initial along-isobath search distance [km]

# maxdx = 50 # Cross-isobath search distance [km]
maxdx = [2, 8] # Left-right cross-isobath search distance [km]
asymmetrical = True

PLOT_BINAVGMAPS = False
SAVEFIG_BINAVGMAPS = False

nprofsin_min = 5 # Minimum number of profiles in each polygon to calculate a spatial average.
nprofsin_max = 1e30 # Maximum number of profiles in each polygon to calculate a spatial average.

nrand_plot_densinv = 500 # On average, every nrand_plot_densinv-th corrected profile is plotted.
max_rand_plots = 10

if asymmetrical:
    dy_min = np.minimum(maxdx[0], maxdx[1])*2
else:
    dy_min = maxdx*2

dy_max = maxdy_initial*3 # [km]
maxdy_increment = 50     # Incremental along-isobath search distance [km]. Should be larger than 'binlen'.

ming = 3                   # Minimum number of data points across all profiles in average to consider that depth in the spatial average
dllbin0 = 2*dy_max/111.120 # Meridional width of bbox for initial subset along each isobath point [degrees]
nstd_thresh = 5            # Maximum number of standard deviations (to filter out outliers)

f = "../data/en4/en4profiles-NWES_all.nc"
fbathymetry = "../data/srtm15p/SRTM15_V2.7.nc"
iso = 200

isobs = [iso]
ds = open_dataset(f)
rho = sigma0(ds["SA"], ds["CT"])
rho.name = "rho"
rho.attrs["units"] = "kg/m3"
rho.attrs["long_name"] = "Potential density referenced to 0 dbar"
lons, lats = ds["lon"].values, ds["lat"].values

dll = 0
xmin, xmax = np.floor(lons.min() - dll), np.ceil(lons.max() + dll)
ymin, ymax = np.floor(lats.min() - dll), np.ceil(lats.max() + dll)
bb = [xmin, xmax, ymin, ymax]

# Get bathymetry just for plotting.
dstopo = open_dataset(fbathymetry).sel(lon=slice(xmin-1, xmax+1), lat=slice(ymin-1, ymax+1))
xt, yt = np.meshgrid(dstopo["lon"].values, dstopo["lat"].values)
ht = -dstopo["z"].values

#---

# Get isobath contour.
fisobath = "isobath_%dm.npz"%iso
disob = np.load(fisobath)

xi, yi = disob["xi_ext"], disob["yi_ext"]
xip_clip, yip_clip = disob["xip"], disob["yip"]
xipm_clip = 0.5*(xip_clip[1:] + xip_clip[:-1])
yipm_clip = 0.5*(yip_clip[1:] + yip_clip[:-1])

#---
xim, yim, angs = angle_isobath(xi, yi)
di = np.append(0, np.cumsum(distance(xim, yim)))*1e-3 # [km]
dxi = np.median(np.diff(di))

# Find start and end indices in the extended isobath.
kstart, kend = near2(xim, yim, xipm_clip[0], yipm_clip[0]), near2(xim, yim, xipm_clip[-1], yipm_clip[-1])
kstart = kstart - 1 # Avoid missing the fist point.
nk = di.size

di = di - di[kstart]

ssns = ["annual", "JFM", "AMJ", "JAS", "OND"]
dmonths = dict(JFM=[1, 2, 3], AMJ=[4, 5, 6], JAS=[7, 8, 9], OND=[10, 11, 12], annual=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12])
nt = len(ssns)
fout = "alongisob_seasonal_EN4_%dm.nc"%iso
headout_figs_binavgmap = "figs_alongisobclim_bins_seasonal/"

xseg0, yseg0 = np.loadtxt("segment_latitudes_wei_etal2024.txt", unpack=True)
fseg = [near2(xim, yim, xseg0[i], yseg0[i]) for i in range(len(yseg0))]
xseg, yseg = xim[fseg], yim[fseg]
dibound = di[fseg]

figsize_map = (14, 11)
isobs_cc = [100, 200, 1000]

t = ds.t.values
months = np.array([Timestamp(ti).month for ti in t])
ntt = nt

zin_density_inversion_test = np.arange(0, iso + dzi_density_inversion_test, dzi_density_inversion_test)

Nnear_all = 0
total_badprofs_densinv = 0
dsiso = None
count_fig = 1
rand_count = 0
z0 = ds["z"].values
for i in range(ntt):
    tspan = ssns[i]
    fm = dmonths[tspan]
    ft = np.array([mo in fm for mo in months])
    if tspan=="annual":
        assert ft.sum()==t.size, "Error: Annual average does not match number of months in the year."

    loni, lati = lons[ft], lats[ft]
    nprofs = ft.sum()
    print("--------------------------------------------")
    print("%s: %d profiles"%(tspan, nprofs))
    print("--------------------------------------------")

    dsi = ds.isel(n=np.where(ft)[0])

    # Save location of all profiles.
    np.savez("en4_xyNWES-" + tspan + "_all.npz", lon=dsi["lon"].values, lat=dsi["lat"].values)
    lonin_in_strip, latin_in_strip = np.array([]), np.array([])

    if PLOT_BINAVGMAPS:
        figm, axm = plt.subplots(figsize=figsize_map, subplot_kw=dict(projection=proj))
        axm.set_extent((lonmin, lonmax, latmin, latmax))
        axm.coastlines()
        axm.add_feature(LAND)
        axm.contour(xt, yt, ht, levels=isobs_cc, colors="brown", linewidths=0.5)
        axm.gridlines(alpha=0.1, draw_labels=True)
        axm.plot(xim, yim, "m")
        axm.plot(loni, lati, linestyle="none", marker="x", ms=0.5, mfc="k", mec="k")
        axm.plot(xseg, yseg, linestyle="none", marker="o", ms=8, mfc="w", mec="b")
        figm.suptitle(tspan, fontsize=15, fontweight="black", y=0.92)

    repeated_profiles = False
    for k in range(kstart, kend):
        xik, yik, angi = xim[k], yim[k], angs[k]
        dllbin0x = dllbin0/np.cos(yik*np.pi/180)
        xikl, xikr = xik - dllbin0x, xik + dllbin0x
        yikl, yikr = yik - dllbin0, yik + dllbin0
        lonkk, latkk = dsi["lon"].values, dsi["lat"].values
        print("")

        NT_longitude_west1 = 2.1 # [degrees east] Longitude where the Norwegian Trench starts to narrow < 75 km.
        NT_longitude_west2 = 3.8 # [degrees east] Longitude where the Norwegian Trench starts to narrow < 50 km.

        # [km] Limit offshore cross-isobath search distance in the Norwegian Trench to avoid getting profiles from the Norwegian shelf.
        maxdxi = maxdx.copy()
        if xik>NT_longitude_west1:
            max_NT_dx = 60 # [km]

            if xik>NT_longitude_west2:
                max_NT_dx = 30 # [km]

            if asymmetrical:
                maxdxi[1] = np.minimum(maxdx[1], max_NT_dx)
            elif not asymmetrical and maxdx>max_NT_dx:
                maxdxi = np.minimum(maxdx, max_NT_dx)

        # Preliminary subset (bbox around point along isobath).
        fkk = np.logical_and(np.logical_and(lonkk>=xikl, lonkk<=xikr), np.logical_and(latkk>=yikl, latkk<=yikr))
        dsikk = dsi.isel(n=np.where(fkk)[0])
        lonk0, latk0 = dsikk["lon"].values, dsikk["lat"].values

        if ANISOTROPIC_BIN:
            maxdy, nnear = maxdy_initial, 0

            while np.logical_or(nnear<nprofsin_min, nnear>nprofsin_max):
                fnear = np.bool_(np.zeros(lonk0.size))
                fnear_previous = fnear.copy()
                poly, fnear, nnear, lonsp, latsp, lonsp0, latsp0 = get_profiles_polygon(xim, yim, angs, k, maxdy, nk, maxdxi, lonk0, latk0, dxi, asymmetrical=asymmetrical)

                #----- This block is just to subtract the number of rejected profiles from nnear to decide if the strip should be expanded.
                if GLOBAL_DENSINV_TEST and nnear>0:
                    dsin_aux = dsikk.isel(n=np.where(fnear)[0])
                    dsin_aux = dsin_aux.interp(z=zin_density_inversion_test, method="linear")
                    rhoaux = sigma0(dsin_aux["SA"], dsin_aux["CT"]).drop_attrs()
                    ddzrhoaux = rhoaux.diff("z")
                    badprofs_densinv = 0
                    iprof_bad_densinv = []
                    for iprof in range(nnear):
                        ddzrhoi = ddzrhoaux.isel(n=iprof)<density_inversion_threshold
                        ninversions = ddzrhoi.values.sum()
                        if ninversions>0:
                            badprofs_densinv += 1

                    nnear_test = nnear - badprofs_densinv
                    # print("----------------------------------------------")
                    # print("nnear was %d"%nnear)
                    # print("nnear after bad density inversion profiles: %d"%nnear_test)
                    # print("----------------------------------------------")
                else:
                    nnear_test = nnear
                #-----

                maxdxdystr = "[%d km left, %d km right] x %d km"%(maxdxi[1], maxdxi[0], maxdy)
                if nnear_test<nprofsin_min:
                    ttl = "%d valid profiles inside (%s) along-isobath strip centered around lon, lat (%1.1f, %1.1f), less than minimum (%d). Increasing strip size."%(nnear_test, maxdxdystr, xik, yik, nprofsin_min)
                    maxdy += maxdy_increment
                elif nnear_test>nprofsin_max:
                    ttl = "%d valid profiles inside (%s) along-isobath strip centered around lon, lat (%1.1f, %1.1f), more than maximum (%d). Decreasing strip size."%(nnear_test, maxdxdystr, xik, yik, nprofsin_max)
                    maxdy -= maxdy_increment
                else:
                    ttl = "[%s] Profiles inside (%s) along-isobath strip centered around lon, lat (%1.1f, %1.1f): -|%d|-"%(tspan, maxdxdystr, xik, yik, nnear_test)
                print(ttl)
                print("")

                # Break if strip too long or too short.
                if maxdy>dy_max:
                    print("%d km search length exceeds maximum (%d km). Stopping**********************************"%(maxdy, dy_max)); nnear = 0; break
                elif maxdy<dy_min:
                    print("%d km search length less than minimum (%d km). Stopping**********************************"%(maxdy, dy_min)); nnear = 0; break

                if PLOT_BINAVGMAPS:
                    pp = axm.plot(lonsp, latsp, "r-")
                    pc = axm.plot(lonsp0, latsp0, "k", marker="o", ms=3, mfc="w", mec="k")
                    ppin = axm.plot(lonk0[fnear], latk0[fnear], linestyle="none", marker="x", ms=1.5, mfc="r", mec="r")
                    pc2 = axm.plot(xik, yik, marker="o", ms=5, mfc="m", mec="k")
                    plt.draw()
                    plt.show(block=False)
                    axm.set_title(ttl, fontsize=10)
                    if SAVEFIG_BINAVGMAPS:
                        figm.savefig(headout_figs_binavgmap + "binavg_frame-%s.png"%str(count_fig).zfill(3), bbox_inches="tight")
                        count_fig += 1
                    # _ = input("Press any key")
                    # plt.pause(0.5)
                    pp[0].set_visible(False)
                    pc[0].set_visible(False)
                    pc2[0].set_visible(False)
                    ppin[0].set_visible(False)

            # Extract profiles within strip.
            if nnear>0:
                dsin = dsikk.isel(n=np.where(fnear)[0])
            else:
                dsin = dsikk.mean("n").expand_dims("n").transpose() * np.nan # No profiles within strip (or less profiles than ), add a NaN.

        else: # Get profiles within circle with search radius maxdy.
            xys_en4 = zip(lonk0, latk0)
            drs = np.array([distance([xik, xx], [yik, yy]) for xx, yy in xys_en4]).squeeze()*1e-3 # [km]
            fnear = drs<maxdy
            nnear = fnear.sum()
            if nnear==0:
                continue
                print("***No profiles within %d km radius of lon, lat (%1.1f, %1.1f). Skipping..............."%(maxdy, xik, yik))
            else:
                print("Profiles within %d km radius of lon, lat (%1.1f, %1.1f): %d"%(maxdy, xik, yik, nnear))
                dsin = dsikk.isel(n=np.where(fnear)[0])

        Nnear_all += nnear

        # Remove depths with less than a minimum number of valid points before averaging.
        SAin, CTin = dsin["SA"], dsin["CT"]
        SAg, CTg = np.isfinite(SAin).sum(dim="n"), np.isfinite(CTin).sum(dim="n")
        SAin[SAg<ming, :] = np.nan
        CTin[CTg<ming, :] = np.nan

        # Remove outliers using a local range test.
        dsin["SA"] = SAin
        dsin["CT"] = CTin
        dsinmed, dsinstd_bounds = dsin.median(dim="n"), dsin.std(dim="n")*nstd_thresh
        dsinlo, dsinhi = dsinmed - dsinstd_bounds, dsinmed + dsinstd_bounds
        SAlo, SAhi  = dsinlo["SA"], dsinhi["SA"]
        CTlo, CThi  = dsinlo["CT"], dsinhi["CT"]
        fbadSA = np.logical_or(SAin<SAlo, SAin>SAhi)
        fbadCT = np.logical_or(CTin<CTlo, CTin>CThi)

        if fbadSA.any():
            SAin[np.where(fbadSA)] = np.nan
            print("************** %d outlier data points found in SA."%fbadSA.sum())
            dsin["SA"] = SAin

        if fbadCT.any():
            CTin[np.where(fbadCT)] = np.nan
            print("************** %d outlier data points found in CT."%fbadCT.sum())
            dsin["CT"] = CTin

        # Remove outliers using a local density inversion test.
        # 0.03 kg/m3 per 20 db = 0.015 kg/m3 per 10 db.
        # Density inversion QC step following Jones et al. (2023):
        # 1) Interpolate to 20 db and calculate density.
        # 2) Check if there are density inversions greater than 0.03 kg/m3 in each profile
        # 3) Remove entire profiles where one (or more) density inversions is found.
        if GLOBAL_DENSINV_TEST:
            dsin_aux = dsin.interp(z=zin_density_inversion_test, method="linear")
            rhoaux = sigma0(dsin_aux["SA"], dsin_aux["CT"]).drop_attrs()
            ddzrhoaux = rhoaux.diff("z")
            badprofs_densinv = 0
            iprof_bad_densinv = []
            for iprof in range(nnear):
                ddzrhoi = ddzrhoaux.isel(n=iprof)<density_inversion_threshold
                ninversions = ddzrhoi.values.sum()
                if ninversions>0:
                    badprofs_densinv += 1
                    iprof_bad_densinv.append(iprof)

            if badprofs_densinv>0:
                dsin = dsin.drop_isel(n=iprof_bad_densinv) # Drop profiles with density inversions.

                badprof_frac = 100*badprofs_densinv/nnear
                nnear -= badprofs_densinv
                print("+++++ %d profiles (%1.1f%% of total) with density inversions > %1.2f kg/m3 rejected in this subset."%(badprofs_densinv, badprof_frac, density_inversion_threshold))
                total_badprofs_densinv += badprofs_densinv

        if SORT_DENS_PROFILES: # Use isntead of interpolating over density inversions. This assumes the overturns are due to local turbulence only.
            sig0_aux = sigma0(dsin["SA"], dsin["CT"])
            for naux in range(dsin.sizes["n"]):
                sig0_auxn = sig0_aux.isel(n=naux)
                if np.isnan(sig0_auxn.values).all():
                    continue
                fg = np.isfinite(sig0_auxn.values)
                fbad_orig = ~fg
                zauxn = sig0_auxn["z"].values
                if fg.sum()>0:
                    sig0_auxn_filled = np.interp(zauxn, zauxn[fg], sig0_auxn.values[fg]) # Just for getting sorting indices with the NaNs in place.
                else:
                    sig0_auxn_filled = sig0_auxn.values

                idx = np.argsort(sig0_auxn_filled)
                saaux = dsin["SA"].values[:, naux]
                ctaux = dsin["CT"].values[:, naux]
                fg = np.isfinite(saaux)
                saaux = np.interp(zauxn, zauxn[fg], saaux[fg])
                fg = np.isfinite(ctaux)
                ctaux = np.interp(zauxn, zauxn[fg], ctaux[fg])
                saaux = saaux[idx]
                ctaux = ctaux[idx]
                saaux[fbad_orig] = np.nan # Add original NaNs back.
                ctaux[fbad_orig] = np.nan # Add original NaNs back.
                dsin["SA"].values[:, naux] = saaux
                dsin["CT"].values[:, naux] = ctaux

                if PLOT_RANDOM_PROFILES_SORT_DENS:
                    if np.random.randint(0, nrand_plot_densinv)==1:
                        fig, axrand = plt.subplots()
                        sig0_auxn.attrs["units"] = "kg/m3"
                        sig0_auxn.attrs["long_name"] = "Potential density"
                        sig0_auxn_sorted = sig0_auxn.copy()
                        sig0_auxn_sorted.values = sigma0(dsin["SA"].values[:, naux], dsin["CT"].values[:, naux])
                        sig0_auxn.plot(ax=axrand, y="z", yincrease=False)
                        sig0_auxn_sorted.plot(ax=axrand, y="z", yincrease=False)
                        fsorted = sig0_auxn.values != sig0_auxn_sorted.values
                        sig0_pltaux = sig0_auxn_sorted[fsorted]
                        if np.isfinite(sig0_pltaux.values).any():
                            sig0_pltaux.plot(ax=axrand, y="z", yincrease=False, linestyle="none", marker="o", ms=4, mfc="m", mec="m")
                        axrand.axis("tight")
                        rand_count += 1
                        if rand_count>max_rand_plots:
                            plt.show()
                            _ = input("Press any key.")


        # Extrapolate shallowest value towards surface.
        if EXTRAPNN_SURFACE:
            for naux in range(dsin.sizes["n"]):
                dsin_aux = dsin.isel(n=naux)

                SAaux, CTaux, zaux = dsin_aux["SA"].values.copy(), dsin_aux["CT"].values.copy(), dsin_aux["z"].values.copy()
                fsurf = zaux<=znearsurf
                SAtop, CTtop = SAaux[fsurf], CTaux[fsurf]

                # Find shallowest valid data point and extrapolate.
                if np.isfinite(SAtop).any():
                    ftop = np.where(np.isfinite(SAtop))[0][0]
                    if ftop>0:
                        SAaux[:ftop] = SAtop[ftop]

                if np.isfinite(CTtop).any():
                    ftop = np.where(np.isfinite(CTtop))[0][0]
                    if ftop>0:
                        CTaux[:ftop] = CTtop[ftop]

                dsin["SA"].values[:, naux] = SAaux
                dsin["CT"].values[:, naux] = CTaux

        if DISTANCE_WEIGHTING and nnear>0:
            lonsin, latsin = dsin["lon"].values, dsin["lat"].values

            # Ensure sum of weights adds up to zero at every depth (NaNs may be in different positions)
            dr = np.array([distance([xik, xx], [yik, yy]) for xx, yy in zip(lonsin, latsin)]).squeeze()*1e-3 # [km]
            weights0 = gaussr(dr, maxdy, a=gaussr_frac)

            nz = dsin["z"].size
            w = np.ones((nz, nnear))*weights0.copy()
            wCT, wSA = w.copy(), w.copy()
            fnanCT, fnanSA = np.isnan(dsin["CT"].values), np.isnan(dsin["SA"].values)
            wCT[fnanCT] = np.nan
            wSA[fnanSA] = np.nan
            for klev in range(nz):
                wCT[klev, :] /= np.nansum(wCT[klev, :])
                wSA[klev, :] /= np.nansum(wSA[klev, :])

            dsinavg = dsin.mean(dim="n")
            fnanCT = np.isnan(dsinavg["CT"].values)
            fnanSA = np.isnan(dsinavg["SA"].values)
            dsinavg["CT"].values = np.nansum(dsin["CT"].values*wCT, axis=1)
            dsinavg["SA"].values = np.nansum(dsin["SA"].values*wSA, axis=1)
            dsinavg["CT"][fnanCT] = np.nan
            dsinavg["SA"][fnanSA] = np.nan
        else:
            dsinavg = dsin.mean(dim="n")

        # Standard deviation and number of observations per depth level.
        dsinstd = dsin.std(dim="n")
        dsinn = np.isfinite(dsin).sum(dim="n")

        if nnear>0:
            lonin_in_strip = np.append(lonin_in_strip, dsin["lon"].values)
            latin_in_strip = np.append(latin_in_strip, dsin["lat"].values)

        xin = Variable("x", [di[k]], attrs=dict(long_name="Along-isobath distance", units="km"))
        lonin, latin = Variable("x", [xi[k]], attrs=dict(long_name="Longitude", units="Degrees east")), Variable("x", [yi[k]], attrs=dict(long_name="Latitude", units="Degrees north"))

        coordsin = dict(x=xin, lon=lonin, lat=latin)
        dsinmk = dsinavg.assign_coords(coordsin)
        dsinmk_std = dsinstd.assign_coords(coordsin)
        dsinmk_n = dsinn.assign_coords(coordsin)

        if k==kstart:
             dsiso = dsinmk.copy()
             dsiso_std = dsinmk_std.copy()
             dsiso_n = dsinmk_n.copy()
        else:
             dsiso = concat((dsiso, dsinmk), dim="x")
             dsiso_std = concat((dsiso_std, dsinmk_std), dim="x")
             dsiso_n = concat((dsiso_n, dsinmk_n), dim="x")

    dsiso_std, dsiso_n = dsiso_std.transpose(), dsiso_n.transpose()
    dsiso_std["z"], dsiso_n["z"] = -dsiso_std["z"], -dsiso_n["z"]

    # Derive potential density.
    dsiso = dsiso.transpose()
    rhoiso = sigma0(dsiso["SA"], dsiso["CT"])
    rhoiso.name = "rho"
    rhoiso.attrs["units"] = "kg/m3"
    rhoiso.attrs["long_name"] = "Potential density referenced to 0 dbar"
    dsiso["rho"] = rhoiso
    dsiso["z"] = -dsiso["z"]

    ssn = Variable("t", [tspan], attrs=dict(long_name="Season"))
    coordsint = dict(t=ssn)
    dsiso = dsiso.assign_coords(coordsint)
    dsiso_std = dsiso_std.assign_coords(coordsint)
    dsiso_n = dsiso_n.assign_coords(coordsint)

    if i==0:
        dsisos = dsiso
        dsisos_std = dsiso_std
        dsisos_n = dsiso_n
    else:
        dsisos = concat((dsisos, dsiso), dim="t")
        dsisos_std = concat((dsisos_std, dsiso_std), dim="t")
        dsisos_n = concat((dsisos_n, dsiso_n), dim="t")

    # Save (lon, lat) of points in strip only.
    np.savez("en4_xyNWES-" + tspan + ".npz", lon=lonin_in_strip, lat=latin_in_strip)

frac = 100*total_badprofs_densinv/Nnear_all
print("")
print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
print("++++++++++++++ Total %d profiles (%1.1f%%) with density inversions rejected."%(total_badprofs_densinv, frac))
print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
print("")

dsisos["CT"].attrs = dict(long_name="Conservative Temperature", units="Degrees Celsius")
dsisos["SA"].attrs = dict(long_name="Absolute Salinity", units="g/kg")
dsisos["rho"].attrs = dict(long_name="Potential density referenced to 0 dbar", units="kg/m3")

# Merge std and number of valid observations with mean profiles.
dsisos_std = dsisos_std.rename_vars(dict(CT="CT_std", SA="SA_std"))
dsisos_std["CT_std"].attrs = dict(long_name="Conservative Temperature standard deviation", units="Degrees Celsius")
dsisos_std["SA_std"].attrs = dict(long_name="Absolute Salinity standard deviation", units="g/kg")

dsisos_n = dsisos_n.rename_vars(dict(CT="CT_n", SA="SA_n"))
dsisos_n["CT_n"].attrs = dict(long_name="Number of valid Conservative Temperature observations", units="")
dsisos_n["SA_n"].attrs = dict(long_name="Number of valid Absolute Salinity observations", units="")
dsisos = merge((dsisos, dsisos_std, dsisos_n), compat="identical")

dims = ("t", "z", "x")
coords = dsisos.coords

# Get ADT interpolated to isobath and add seasonal ug(y, z).
ug_long_name = "Absolute geostrophic cross-isobath velocity"
fname_adt_ssnavg = "seasonal_altimetry_adt_200misob.nc"
adt_ssnavg = open_dataset(fname_adt_ssnavg)["adt"]

# Calculate ug from geopotential anomaly.
assert dsisos.sizes["x"] == adt_ssnavg.sizes["x"], "Different along-isobath grid points"
xx, zz, lat = dsisos.x.values, dsisos.z.values, dsisos.lat.values
g = grav(lat, 0)  # [m/s^2]
f = fcor(lat)  # [1/s]

PLOT_Phi_ug = False
ssns = dsisos["t"].values
uxis = []
for ssn in ssns:
    print("Season: ", ssn)
    if ssn == "annual":
        adt = adt_ssnavg.mean("season").values
    else:
        adt = adt_ssnavg.sel(season=ssn).values
    SA, CT = dsisos["SA"].sel(t=ssn).values, dsisos["CT"].sel(t=ssn).values

    # Replace first row (z = 0 m, all NaN), with the second row (z = -5 m) just for this ug calculation.
    SA[0, :] = SA[1, :]
    CT[0, :] = CT[1, :]

    p = p_from_z(zz[:, np.newaxis], lat)
    Phi = geo_strf_dyn_height(SA, CT, p, p_ref=0)
    ugbc = np.gradient(Phi, xx*1e3, axis=1)/f  # [m/s]
    ug0 = - g*np.gradient(adt, xx*1e3)/f # Surface geostrophic velocity [m/s].
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

uxis = DataArray(np.array(uxis)*1e2, dims=dims, coords=coords, attrs=dict(long_name=ug_long_name, units="cm/s"))
dsisos["ug"] = uxis

dsisos.to_netcdf(fout)
