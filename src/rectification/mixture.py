import numpy as np
from thermo import Mixture as TMixture, Chemical
from rectification.core import ureg, propSI


class Mixture(TMixture):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @property
    def uMWs(self):
        return np.array(self.MWs) * propSI["molar mass"]

    @property
    def uMW(self):
        return self.MW * propSI["molar mass"]

    @property
    def uTb(self):
        return self.Tb * ureg.K

    @property
    def usigma(self):
        return super().sigma * propSI["surface tension"]

    @property
    def urhol(self):
        return super().rhol * propSI["density"]

    @property
    def urhog(self):
        return super().rhog * propSI["density"]

    @property
    def uVmls(self):
        return np.array(super().Vmls) * propSI["molar volume"]

    @property
    def umuls(self):
        return np.array(super().muls) * propSI["viscosity"]

    @property
    def umugs(self):
        return np.array(super().mugs) * propSI["viscosity"]

    @property
    def umul(self):
        return super().mul * propSI["viscosity"]

    @property
    def umug(self):
        return super().mug * propSI["viscosity"]

    def __repr__(self):
        return super().__repr__()

    def __str__(self):
        return super().__str__()


class Water(Chemical):

    def __init__(self, *args, **kwargs):
        super().__init__("water", *args, **kwargs)

    def calculate(self, T=None, P=None):
        return super().calculate(T=T, P=P)

    @property
    def uT(self):
        return self.T * ureg.K

    @property
    def urhol(self):
        return super().rhol * propSI["density"]

    @property
    def urhog(self):
        return super().rhog * propSI["density"]

    @property
    def usigma(self):
        return super().sigma * propSI["surface tension"]

    @property
    def uHvapm(self):
        return super().Hvapm * propSI["mass enthalpy"]

    @property
    def umul(self):
        return super().mul * propSI["viscosity"]

    @property
    def ukl(self):
        return super().kl * propSI["thermal conductivity"]

    @property
    def ukg(self):
        return super().kg * propSI["thermal conductivity"]

    @property
    def uCpl(self):
        return super().Cpl * propSI["heat capacity"]

    @property
    def uCpg(self):
        return super().Cpg * propSI["heat capacity"]
