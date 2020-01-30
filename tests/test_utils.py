from rectification.utils import unitcheck
from rectification.core import ureg


def test_unitcheck_full_kwargs():

    @unitcheck(b="s", c="l", res_unit="km*hour*m^3")
    def func1(a, b, c):
        return a * b * c

    a = 5e3 * ureg.m
    b = 30 * ureg.min
    c = 4
    res = func1(a, b, c)
    assert str(res.u) == "hour * kilometer * meter ** 3"
    assert res.m == 3.6e7
