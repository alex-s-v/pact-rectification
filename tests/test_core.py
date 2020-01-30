from rectification.core import ureg


def test_ureg():
    a = 5 * ureg.km
    b = 7.2 * ureg.m
    c = (a + b).to("m")
    assert c.m == 5007.2
