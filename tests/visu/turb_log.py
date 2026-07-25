"""
Helpers to recover the turbulence history of a single test from the RAMSES
test suite log (tests/test_suite.log).

The log holds every test that was run, so the block belonging to one test has
to be isolated before anything is parsed out of it. A block opens with e.g.
'Test 1/1: turb/driving' and closes with 'Test turb/driving passed/failed'.
Note that run_test_suite.sh calls the plotting scripts *before* it writes that
closing line, so in practice a block normally runs to the end of the file; the
header of the next test and EOF are used as fallbacks.

Within a block, update_time writes

    Main step=  6 mcons=... ekin= 1.23E+00 eint= 1.02E+00
    Fine step=  6 t= 1.75446E-02 dt= 2.924E-03 ...
     Current turbulent rms:    8.5578662989441927

in that order, and only advances t *after* the step lines. So ekin belongs to
time t, while the rms that follows belongs to time t+dt - both times are read
off the same 'Fine step=' line.
"""

import os
import re
import numpy as np

_STEP_RE = re.compile(r'Fine step=\s*\d+\s+t=\s*(\S+)\s+dt=\s*(\S+)')
_RMS_RE = re.compile(r'Current turbulent rms:\s*(\S+)')
_EKIN_RE = re.compile(r'Main step=.*\sekin=\s*(\S+)')


def _to_float(text):
    """Fortran may write exponents with D rather than E."""
    return float(text.replace('D', 'E').replace('d', 'e'))


def read_turb_history(logfile, testname):
    """
    Parameters
    ----------
    logfile: str
        Path to the test suite log.
    testname: str
        Test to isolate, as it appears in the log, e.g. 'turb/driving'.

    Returns
    -------
    dict with keys
        't'         time of each step
        'ekin'      total kinetic energy at 't' (NaN where not reported)
        't_rms'     time the turbulent rms applies to, i.e. t+dt
        'rms'       turbulent forcing rms at 't_rms'
        't_restart' time the run was restarted at (suite run with -r), else None
    All arrays are empty if the log or the test block cannot be found, so that
    a missing log degrades to an empty plot rather than breaking the test.
    """
    empty = {'t': np.array([]), 'ekin': np.array([]), 't_rms': np.array([]),
             'rms': np.array([]), 't_restart': None}

    if not os.path.exists(logfile):
        return empty

    start_re = re.compile(r'^\s*Test\s+\d+/\d+:\s*' + re.escape(testname) + r'\s*$')
    other_re = re.compile(r'^\s*Test\s+\d+/\d+:')
    end_re = re.compile(r'^\s*Test\s+' + re.escape(testname) + r'\s+(passed|failed)')

    t, ekin, t_rms, rms = [], [], [], []
    t_restart = None
    inside = False
    at_restart = False
    pending_ekin = None
    pending_trms = None

    with open(logfile) as f:
        for line in f:
            if not inside:
                inside = bool(start_re.match(line))
                continue
            if other_re.match(line) or end_re.match(line):
                break

            # Second half of a restart test starts here
            if 'Restart: step 2' in line:
                at_restart = True
                continue

            match = _EKIN_RE.search(line)
            if match:
                pending_ekin = _to_float(match.group(1))
                continue

            match = _STEP_RE.search(line)
            if match:
                step_t = _to_float(match.group(1))
                step_dt = _to_float(match.group(2))
                if at_restart:
                    t_restart = step_t
                    at_restart = False
                t.append(step_t)
                ekin.append(pending_ekin if pending_ekin is not None else np.nan)
                pending_ekin = None
                pending_trms = step_t + step_dt
                continue

            match = _RMS_RE.search(line)
            if match and pending_trms is not None:
                t_rms.append(pending_trms)
                rms.append(_to_float(match.group(1)))
                pending_trms = None

    return {'t': np.array(t), 'ekin': np.array(ekin),
            't_rms': np.array(t_rms), 'rms': np.array(rms),
            't_restart': t_restart}


def _pad_ylim(axis, lo, hi):
    """Set y limits leaving a band clear at the bottom for a legend."""
    span = hi - lo
    # Keep at least a few percent of the scale in view. A series that is
    # constant, or nearly so, would otherwise land on a wildly magnified axis
    # with an offset label - which is exactly the expected case for the static
    # field of turb_type=2.
    floor = 0.02*max(abs(lo), abs(hi))
    if floor <= 0.0:
        floor = 0.1
    span = max(span, floor)
    axis.set_ylim(lo - 0.6*span, hi + 0.15*span)


def plot_history(ax, hist, target_rms=None, ndim=3, show_rms=True,
                 show_ekin=True, rms_label='Current turbulent rms',
                 missing='No history found in ../../test_suite.log'):
    """
    Draw the turbulent rms and/or kinetic energy history returned by
    read_turb_history onto `ax`.

    The rms goes on `ax` itself and the kinetic energy on a twinned right axis,
    so both can be read against a common time axis. Kinetic energy is the more
    sensitive tracer of the two: swapping the interpolation weights between the
    two bracketing fields changes the direction of the forcing but barely its
    magnitude, so such a bug moves the energy while leaving the rms flat.

    Returns the axis the kinetic energy was drawn on, or None.
    """
    t_rms, rms = hist['t_rms'], hist['rms']
    t, ekin = hist['t'], hist['ekin']
    good = np.isfinite(ekin)

    have_rms = show_rms and len(t_rms) > 0
    have_ekin = show_ekin and np.any(good)

    ax.set_xlabel('Time')

    if not (have_rms or have_ekin):
        ax.text(0.5, 0.5, missing, ha='center', va='center', transform=ax.transAxes)
        return None

    handles = []
    ax_ekin = None

    if have_rms:
        line, = ax.plot(t_rms, rms, marker='o', ms=3, lw=1, color='C0',
                        label=rms_label)
        handles.append(line)
        lo, hi = rms.min(), rms.max()
        if target_rms is not None:
            target = np.sqrt(ndim)*target_rms
            handles.append(ax.axhline(target, color='k', ls='--', lw=1,
                           label=r'$\sqrt{N_{dim}}\times$turb_rms = %.3f' % target))
            lo, hi = min(lo, target), max(hi, target)
        _pad_ylim(ax, lo, hi)
        ax.set_ylabel('Turbulent rms')

    if have_ekin:
        ax_ekin = ax.twinx() if have_rms else ax
        line, = ax_ekin.plot(t[good], ekin[good], marker='s', ms=3, lw=1,
                             color='C1', label='Kinetic energy')
        handles.append(line)
        _pad_ylim(ax_ekin, ekin[good].min(), ekin[good].max())
        ax_ekin.set_ylabel('Kinetic energy')
        if have_rms:
            ax.tick_params(axis='y', labelcolor='C0')
            ax_ekin.tick_params(axis='y', labelcolor='C1')
            ax_ekin.yaxis.label.set_color('C1')
            ax.yaxis.label.set_color('C0')

    if hist['t_restart'] is not None:
        handles.append(ax.axvline(hist['t_restart'], color='C3', ls=':', lw=1.5,
                       label='restart at t = %.4f' % hist['t_restart']))

    ax.legend(handles=handles, loc='lower left', fontsize='small',
              ncol=len(handles))
    return ax_ekin


def read_target_rms(nmlfile):
    """
    Forcing rms requested in a namelist, or None if it cannot be read.

    This is the per-component amplitude, whereas current_turb_rms sums the
    squares of all ndim components, so the rms reported in the log should sit
    at sqrt(ndim)*turb_rms.
    """
    if not os.path.exists(nmlfile):
        return None
    with open(nmlfile) as f:
        for line in f:
            match = re.match(r'\s*turb_rms\s*=\s*([0-9eEdD.+-]+)', line)
            if match:
                return _to_float(match.group(1))
    return None
