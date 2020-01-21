"""
Module with useful functions and others.
"""
import inspect
import functools

from rectification.core import ureg


def isunit(value, ureg=ureg):
    """Check if the value is the `pint.Quantity`

    Parameters
    ----------
    value : any
        Any value to check if it is of the `pint.Quantity` type.

    Returns
    -------
    bool
        True if value is the `pint.Quantity`, False otherwise.

    Notes
    -----
        For correct operation values must be of the same Unit Registry.
    """
    return isinstance(value, ureg.Quantity)


def convunits(units, values):
    """Convert values to the desired units.
    Parameters
    ----------
    units : dict[str, str]
        Dictionaty of name - unit pairs, where name is the key from
        values dict and unit is the desired unit of its value.
    values : dict[str, Any]
        Dictionary of name - value paires, where name is the name of the
        value which will be converted or assigned to the unit from
        the units dict if possible.
    Returns
    -------
    dict[str, Any]
        New dictionary with all possible values converted to the
        specified units.
    """
    nv = {}
    for p, v in values.items():
        if isunit(v):
            if p in units:
                v = v.to(units[p])
        else:
            if p in units:
                v = ureg.Quantity(v, units[p])
        nv[p] = v
    return nv


def convparams(f, units, args, kwargs):
    """Convert all parameters of the function to the desired units.

    Parameters
    ----------
    f : callable
        Function parameters of which should be converted.
    units : dict[str, str]
        Dictionary of name - unit pairs, where name is the
        name of the parameter which should be set or converted to
        the unit.
    args : list of Any
        List of the positional paramerers values.
    kwargs : dict[str, Any]
        Dictionary of the parameter name - value pairs.

    Returns
    -------
    list of Any
        List of the positional paramerers values.
    dict[str, Any]
        Dictionary of the parameter name - value pairs.
    """
    func_args = inspect.getcallargs(f, *args, **kwargs)
    la = len(args)
    na, nk = [], {}
    fixed_args = convunits(units, func_args)
    for i, (p, v) in enumerate(fixed_args.items()):
        if i < la:
            na.append(v.m)
        else:
            nk[p] = v.m
    return na, nk


def unitcheck(units, res_unit=None):
    """
    Convinience function to check and convert parameters
    of the function to the desired units.

    Intended to use as a function or method decorator.

    Parameters
    ----------
    units : dict[str, str]
        Dictionary of name - unit pairs, where name is the
        name of the parameter which should be set or converted to
        the unit.
    res_unit : str, optional
        If specified will determine units of the output value.

    Returns
    -------
    callable
        Wrapper for the function.
    """
    def decorator(f):
        @functools.wraps(f)
        def wrapper(*args, **kwargs):
            args, kwargs = convparams(f, units, args, kwargs)
            res = f(*args, **kwargs)
            if res_unit is not None:
                res = ureg.Quantity(res, res_unit)
            return res
        return wrapper
    return decorator
