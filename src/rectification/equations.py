from rectification.utils import unitcheck
import numpy as np


def operating_line(x, R, xp):
    return R * x / (R + 1) + xp / (R + 1)


@unitcheck(L="kg/s", rhol="kg/m**3", Ft="m**2", res_unit="m**3/(m**2*s)")
def calc_Ucoef(L, rhol, Ft):
    return L / (rhol * Ft)
