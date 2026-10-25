"""Theoretical invasion velocities, Eqs (7), (8) and (11) of the main text.

All functions are vectorized over a_fs and a_sf and return NaN where a velocity does not apply,
so curves can be passed to matplotlib directly; regime boundaries then appear as gaps.
Sign convention: u > 0 when the fast type (or the mixture) advances into the slow type,
u < 0 when the slow type advances into the fast type.
"""
import numpy as np


def _arrays(afs, asf):
    return np.broadcast_arrays(np.asarray(afs, float), np.asarray(asf, float))


def u_exclusion(afs, asf, dr, D=1.0):
    """Eq (8): invasion velocity outside the coexistence region (NaN inside it)."""
    afs, asf = _arrays(afs, asf)
    coex = (afs > -dr) & (asf > dr)
    with np.errstate(invalid="ignore"):
        u = np.where(asf > dr,
                     -2 * np.sqrt(D * (asf - dr)),                     # dent
                     2 * np.sqrt(D * np.maximum(dr, dr + afs)))        # circular / composite bulge
    return np.where(coex, np.nan, u)


def u_special_case(a, dr, D=1.0):
    """Eq (7): the line a_fs = -a_sf = a (Fig 2A)."""
    return u_exclusion(a, -np.asarray(a, float), dr, D)


def velocities(afs, asf, dr, D=1.0):
    """All four velocities plotted in Fig 3C, as a dict of arrays:
    u_fs (fast into slow), u_sf (slow into fast, as a positive speed),
    u_ms (mixture into slow), u_mf (mixture into fast), Eq (11)."""
    afs, asf = _arrays(afs, asf)
    coex = (afs > -dr) & (asf > dr)
    dent = (asf > dr) & (afs <= -dr)
    fast = asf <= dr
    with np.errstate(invalid="ignore", divide="ignore"):
        geo_mix = 2 * np.sqrt(D * asf * (dr + afs) / (afs + asf))
        fisher_fast = 2 * np.sqrt(D * (dr + afs))
        fisher_slow = 2 * np.sqrt(D * (asf - dr))
        return dict(
            u_fs=np.where(fast, 2 * np.sqrt(D * np.maximum(dr, dr + afs)), np.nan),
            u_sf=np.where(dent, fisher_slow, np.nan),
            u_ms=np.where(coex, np.where(afs <= 0, geo_mix, fisher_fast), np.nan),
            u_mf=np.where(coex, fisher_slow, np.nan),
        )


def f_mix(afs, asf, dr):
    """Eq (9): composition of the mixed state (NaN outside the coexistence region)."""
    afs, asf = _arrays(afs, asf)
    coex = (afs > -dr) & (asf > dr)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(coex, (dr + afs) / (afs + asf), np.nan)
