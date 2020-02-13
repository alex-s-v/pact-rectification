"""
Module with core fuctionality of this package.
"""
import warnings
import pint


ureg = pint.UnitRegistry()

# Disable new NumPy - pint interaction warning
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    ureg.Quantity([])


propSI = {
    "molar mass": ureg.g / ureg.mol,
    "viscosity": ureg.Pa * ureg.s,
    "surface tension": ureg.N / ureg.m,
    "molar volume": ureg.m**3 / ureg.mol,
    "density": ureg.kg / ureg.m**3,
    "mass enthalpy": ureg.J / ureg.kg,
    "thermal conductivity": ureg.W / ureg.m / ureg.K,
    "heat capacity": ureg.J / ureg.kg / ureg.K
}
