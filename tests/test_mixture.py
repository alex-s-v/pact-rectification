from rectification.mixture import Mixture, Water
from rectification.utils import isclose


def test_mix_basic():
    mix = Mixture(["benzen", "toluene"], zs=[0.2, 0.8])
    assert [0.2, 0.8] == mix.zs


def test_water_basic():
    water = Water()
    assert isclose(298.15, water.uT)
