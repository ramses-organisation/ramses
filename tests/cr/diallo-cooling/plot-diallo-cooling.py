import matplotlib as mpl
mpl.use("Agg")
import os, glob, sys
import numpy as np
import matplotlib.pyplot as plt
import visu_ramses

title = "diallo-cooling"

# Must match amr/constants.f90, hydro/cooling_module.f90, cr_condinit 'diallo_314' and the namelist
X_H = 0.76
M_H = 1.6605390e-24                      # g
PC = 3.0856776e18                        # cm
MYR = 3.156e13                           # s
E_MAX_CGS = 2.5e-10                      # erg cm^-3
SIGMA_PC = 10.0
ZETA_CR, CR_NE, CR_FNEUT = 7.51e-16, 1e-3, 0.875
LAMBDA_CR = ZETA_CR * (1.0 + 0.22 * CR_NE + 0.125 * CR_FNEUT)   # cm^3 s^-1
TOL_PROFILE = 1e-3
TOL_TOTAL = 1e-9


def e0_profile(x, unit_d, unit_l, unit_t):
    """Initial CR energy density [code units] of cr_condinit_kind='diallo_314'."""
    e_max = E_MAX_CGS / (unit_d * (unit_l / unit_t) ** 2)
    sigma = SIGMA_PC * PC / unit_l
    x0 = 0.5 * (x.min() + x.max())       # box centre on a uniform periodic grid
    return e_max * np.exp(-0.5 * ((x - x0) / sigma) ** 2)


def analytic(d):
    """E(x,t) = E0(x) exp(-lambda_cr n_H t)."""
    nH = d["density"] * d["unit_d"] * X_H / M_H
    decay = np.exp(-LAMBDA_CR * nH * d["time"] * d["unit_t"])
    e0 = e0_profile(d["x"], d["unit_d"], d["unit_l"], d["unit_t"])
    return e0, e0 * decay, decay


def analytic_errors(d):
    e0, e_an, decay = analytic(d)
    e = d["CRegy_01"]
    err_profile = np.max(np.abs(e - e_an)) / np.max(e_an)
    err_total = abs(np.sum(e) / np.sum(e0) / np.mean(decay) - 1.0)
    return err_profile, err_total


if __name__ == "__main__":
    snaps = sorted(int(os.path.basename(p).split("_")[-1])
                   for p in glob.glob("output_?????") if os.path.isdir(p))
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(5.2, 5.0), sharex=True)
    cols = plt.cm.viridis(np.linspace(0.0, 0.82, len(snaps)))
    failed = []
    for i, si in enumerate(snaps):
        data = visu_ramses.load_snapshot(si, read_cr=True)
        d = data["data"]
        e0, e_an, decay = analytic(d)
        err_profile, err_total = analytic_errors(d)
        t_myr = float(d["time"]) * d["unit_t"] / MYR
        print("t=%g Myr: profile error %.3e, total-energy error %.3e" % (t_myr, err_profile, err_total))
        if err_profile > TOL_PROFILE or err_total > TOL_TOTAL:
            failed.append(t_myr)
        x = d["x"]; srt = np.argsort(x)
        ax1.plot(x[srt], d["CRegy_01"][srt], "o", ms=3, color=cols[i], label="t=%g Myr" % t_myr)
        ax1.plot(x[srt], e_an[srt], "-", lw=1, color=cols[i])
        ax2.plot(x[srt], (d["CRegy_01"][srt] - e_an[srt]) / np.max(e_an), color=cols[i])
    x0 = 0.5 * (x.min() + x.max())
    ax2.set_xlim(x0 - 60.0, x0 + 60.0)
    ax1.set_ylabel(r"$E_{\rm CR}$")
    ax2.set_ylabel(r"$(E-E_{\rm an})/\max E_{\rm an}$")
    ax2.set_xlabel("x [pc]")
    ax1.legend(fontsize=8, framealpha=0.5)
    ax1.set_title(title)
    fig.tight_layout(); plt.subplots_adjust(hspace=0.08)
    fig.savefig(title + ".pdf", bbox_inches="tight")

    if failed:
        print("FAILED: analytic CR cooling check at t =", failed, "Myr")
        sys.exit(1)

    visu_ramses.check_solution(data["data"], title, tolerance={"all": 1e-12}, overwrite=False)
