import matplotlib as mpl
import numpy as np
from scipy.interpolate import griddata
# uncomment for command line use, comment for notebook
mpl.use("Agg")
import matplotlib.pyplot as plt
import visu_ramses
from matplotlib.colors import LogNorm

# NOTE: this test is designed to run on 2 proc
# on a different number of proc, it will fail unless tolerance are significantly increased
# this is because of the random number generator use in the star formation routine


X_H, M_H, K_B = 0.76, 1.6605390e-24, 1.3806490e-16   # amr/constants.f90

# Load RAMSES output
data = visu_ramses.load_snapshot(3,read_hydro=True,read_cr=True)
x      = data["data"]["x"]
y      = data["data"]["y"]
z      = data["data"]["z"]
rho    = data["data"]["density"]
p      = data["data"]["pressure"]
metals = data["data"]["metallicity"]
ecr    = data["data"]["CRegy_01"]
unit_d = data["data"]["unit_d"]
unit_v = data["data"]["unit_l"]/data["data"]["unit_t"]

# 10 kpc wide slices through the box centre, sampled on npix^2 pixels
centre = np.array([25., 25., 25.])   # kpc
half   = 5.
npix   = 2**7
s      = np.linspace(-half, half, npix, endpoint=False) + half/npix   # pixel centres
cut    = centre + s[npix//2]         # slice just off the centre, avoiding cell-face ties
h, v   = np.meshgrid(s, s)

fields = np.transpose([rho*unit_d*X_H/M_H,           # n_H [cm^-3]
                       p/rho*unit_v**2*M_H/K_B,      # T/mu [K]
                       metals,                       # Z
                       ecr*unit_d*unit_v**2])        # E_CR [erg cm^-3]
points = np.transpose([x, y, z])
planes = [(0, 1, 2), (1, 0, 2), (2, 0, 1)]   # (normal, horizontal, vertical) axis
slices = []
for normal, hor, ver in planes:
    q = np.empty((npix*npix, 3))
    q[:, normal] = cut[normal]
    q[:, hor]    = centre[hor] + h.ravel()
    q[:, ver]    = centre[ver] + v.ravel()
    slices.append(griddata(points, fields, q, method='nearest').reshape(npix, npix, -1))

# (row title, colorbar label, colormap, max decades shown or None for the full range)
rows = [("Gas density",               r"$n_{\rm H}$ [cm$^{-3}$]",      "viridis", None),
        ("Gas temperature",           r"$T/\mu$ [K]",                  "inferno", None),
        ("Metallicity",               r"$Z$ (metal mass fraction)",    "cividis", None),
        ("Cosmic-ray energy density", r"$E_{\rm CR}$ [erg cm$^{-3}$]", "magma",   6)]
axname = "xyz"

fig = plt.figure(figsize=(11, 13), layout="constrained")
for k, (sf, (title, label, cmap, decades)) in enumerate(zip(fig.subfigures(len(rows), 1), rows)):
    imgs = [sl[:, :, k] for sl in slices]
    pos  = np.concatenate([im[im > 0] for im in imgs])
    if pos.size == 0:
        pos = np.ones(1)
    vmax = pos.max()
    vmin = pos.min() if decades is None else max(pos.min(), vmax*10.**(-decades))
    cm = mpl.colormaps[cmap].copy()
    cm.set_bad(cm(0.)); cm.set_under(cm(0.))           # zero and clipped values at the low end
    norm = LogNorm(vmin=vmin, vmax=vmax)               # one scale shared by the row
    axs = sf.subplots(1, 3)
    for ax, img, (normal, hor, ver) in zip(axs, imgs, planes):
        im = ax.imshow(np.ma.masked_less_equal(img, 0.), origin="lower", cmap=cm, norm=norm,
                       interpolation='none',
                       extent=[centre[hor]-half, centre[hor]+half, centre[ver]-half, centre[ver]+half])
        ax.set_title("%s = %.2f kpc" % (axname[normal], cut[normal]), fontsize=9)
        ax.set_xlabel(axname[hor]+" [kpc]")
        ax.set_ylabel(axname[ver]+" [kpc]")
    sf.colorbar(im, ax=axs, label=label, extend="min" if vmin > pos.min() else "neither")
    sf.suptitle(title, fontweight="bold")

fig.savefig('g8-lowres.pdf',bbox_inches='tight')

#to_check={}
#to_check = data["particle"]
#to_check["pressure"]=z5.flatten()  # data["data"]["pressure"]

to_check = data["data"]
for key in data["particle"].keys():
    to_check['particle '+key] = data["particle"][key]

# default test nproc is 2
nproc=2

# tweak tolerance to allow for 2p vs 4p/8p/12p deviations:
tolerance={}
if(nproc!=2):
    tolerance={"density":1.0e-4,\
               "pressure":1.e-4,\
               "temperature":1.e-4,\
               "metallicity":1.e-4,\
               "CRegy_01":1.e-4,\
               "CRflx_01_x":1.e-4,\
               "CRflx_01_y":1.e-4,\
               "CRflx_01_z":1.e-4,\
}

# and this is the "default tolerance"
if(nproc==2):
     tolerance={"all":1.0e-13,\
               }
visu_ramses.check_solution(to_check, 'g8-lowres',tolerance=tolerance,overwrite=True)
