from rectification.core import ureg
from rectification.utils import isclose
import rectification.equations as eqs


def test_calc_Ucoef():
    L = 500 * 3600 * ureg.kg / ureg.hour
    rhol = 1000
    Ft = 1.57 * ureg.m**2
    U_real = 0.318 * ureg.m**3 / (ureg.m**2 * ureg.s)
    # U_real = ureg.Quantity(0.318, "m**3 / (m**2 * s)")
    U_calc = eqs.calc_Ucoef(L, rhol, Ft)
    assert isclose(U_real, U_calc, atol=1e-3)
