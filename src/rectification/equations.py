from rectification.utils import unitcheck
from scipy.constants import g
import numpy as np


def operating_line(x, R, xp):
    return R * x / (R + 1) + xp / (R + 1)


@unitcheck(L="kg/s", rhol="kg/m**3", Ft="m**2", res_unit="m**3/(m**2*s)")
def calc_Ucoef(L, rhol, Ft):
    return L / (rhol * Ft)


@unitcheck(Massl="g/mol", Massh="g/mol", mu_solv="Pa*s", nul="sm**3/mol", nuh="sm**3/mol", res_unit="m**2/s")
def calc_Diffcoef20(Massl, Massh, A , B, mu_solv, nul, nuh):
    """
    Calculates the diffusion coefficient at 20 degrees celcium.
    Parameters
    ----------
    Massl : float
        The molar mass of low-boilling component, [g/mol]
    Massh : float
        The molar mass of high-boilling component, [g/mol]
    A : float
        The correction coefficient depending on the properties solute, [dimensionless]
    B : float
        The correction coefficient depending on the properties solvent, [dimensionless]
    mu_solv : float
        The viscocity of solvent liquid, [Pa/s]
    nul : float
        The molar volume of solute, [sm**3/s]
    nuh : float
        The molar volume of solvent, [sm**3/s]
    Returns
    -------
    calc_Diffcoef20 : float
        The diffusion coefficient at 20 degrees celcium, [m**2/s]
    References
    ----------
    Романков, страница 289, формула 6.22
    """
    return 1e-6 * ((1/Massl) + (1/Massh))**0.5 / (A * B * mu_solv**0.5 * ((nul)**0.66 + (nuh)*0.66)**2)


def b(mul_20, rhol_20):
    """
    Calculates the temperature coefficient.
    Parameters
    ----------
    mul_20 : float
        The viscocity of low-boilling component of liquid at 20 degrees celcium, [Pa/s]
    rhol_20 : float
        The destinity of low-boilling component of liquid at 20 degrees celcium, [kg/m**3]
    Returns
    -------
    b : float
        The temperature coefficient, [dimensionless]
    References
    ----------
    Романков, страница 289, формула 6.24
    """
    return mul_20**0.5/rhol_20**0.66


def calc_Diffliq(calc_Diffcoef20, b, t_boil):
    """
    Calculates the diffusion coefficient of liquid phaze.
    Parameters
    ----------
    calc_Diffcoef20 : float
        The diffusion coefficient at 20 degrees celcium, [m**2/s]
    b : float
        The temperature coefficient, [dimensionless]
    t_boil : float
        The temperature of low-boilling component of liquid, [degrees celcium]
    Returns
    -------
    calc_Diffliq : float
        The diffusion coefficient of liquid phaze.
    References
    ----------
    Романков, страница 289, формула 6.23
    """
    return calc_Diffcoef20 * (1 + b * (t_boil - 20))


@unitcheck(t_boil="K", Massl="g/mol", Massh="g/mol", P_abs="Pa", nul="sm**3/mol", nuh="sm**3/mol", res_unit="m**2/s")
def calc_Diffvapor(t_boil, P_abs, Massl, Massh, nul, nuh):
    """
    Calculates the diffusion coefficient of vapor.
    Parameters
    ----------
    Massl : float
        The molar mass of the low-boilling component, [g/mol]
    Massh : float
        The molar mass of the high-boilling component, [g/mol]
    t_boil : float
        The boiling temperature of the low-boiling component, [K]
    P : float
        The absolute pressure of the column, [Pa]
    nul : float
        The molar volume of solute, [sm**3/s]
    nuh : float
        The molar volume of solvent, [sm**3/s]
    Returns
    -------
    calc_Diffvapor : float
        The diffusion coefficient of vapor, [m**2/s]
    References
    ----------
    Романков, страница 234, формула 6.25
    """
    return 4.22e-2 * t_boil**(3/2) * ((1/Massl) + (1/Massh))**0.5 / (P_abs * ((nul)**0.66 + (nuh)*0.66)**2)


@unitcheck(q_liq="m**2/s", h_septum="m", w_oper="m/s", mu_mix="Pa/s", sigma_mix="N/m", sigma_water="N/m", res_unit="m")
def heigth_layer(q_liq, h_septum, w_oper, m_coef, mu_mix, sigma_mix, sigma_water):
    """
    Calculates the heigth ligth layer of  the liquid.
    Parameters
    ----------
    q_liq : float
        The specific flow rate of the liquid for 1 m of drain septum, [m**2/s]
    h_septum : float
        The heigth of drain septum, [m]
    w_oper : float
        The operating speed in the column, [m/s]
    mu_mix : float
        The viscocity of mix [Pa/s]
    m_coef : float
    The specific coefficient for this equation [dimensionless]
    sigma_mix : float
        The surface tension of mix [N/m]
    sigma_water : float
        The surface tension of water [N/m]
    Returns
    -------
    heigth_layer : float
        The heigth ligth layer of  the liquid, [m]
    References
    ----------
    Дытнерский, страница 239, формула 6.39
    """
    return 0.787 * q_liq**0.2 * h_septum**0.56 * w_oper**m_coef * (1 - 0.31 * np.exp(-0.11 * mu_mix)) * (sigma_mix/sigma_water)**0.09


def m_coef(h_septum):
    """
    Calculates the specific coefficient for calculation the heigth ligth layer of liquid equation.
    Parameters
    ----------
    h_septum : float
        The heigth of drain septum, [m]
    Returns
    -------
    m_coef : float
        The specific coefficient for this equation [dimensionless]
    References
    ----------
    Дытнерский, страница 239, формула 6.39
    """
    return 0.05 - h_septum*4.6


@unitcheck(rho_mix="kg/m**3", flate_liq="kg/s", L_septum="m", res_unit="m**2/s")
def q_liq(rho_mix, L_septum, flate_liq):
    """
    Calculates the specific flow rate of the liquid for 1 m of drain septum
    Parameters
    ----------
    rho_mix : float
        The destiny of mix, [kg/m**3]
    flate_liq : float
        The flow rate of liquid [kg/s]
    L_septum : float
        The length of drain septum [m]
    Returns
    -------
    q_liq : float
        The specific flow rate of the liquid for 1 m of drain septum, [m**2/s]
    References
    ----------
    Дытнерский, страница 239, формула 6.39
    """
    return flate_liq / (rho_mix * L_septum)


def epsi_vapor(Fr):
    """
    Calculates the vapor content of bubble layer
    Parameters
    ----------
    Fr : float
        The Frudo criterion, [dimensionless]
    Returns
    -------
    epsi_vapor : float
        The vapor content of bubble layer, [dimensionless]
    References
    ----------
    Дытнерский, страница 207, формула 5.47
    """
    return Fr**0.5 / (1 + Fr**0.5) 


def Fr(w_oper, heigth_layer):
    """
    Calculates the Frudo criterion
    Parameters
    ----------
    w_oper : float
        The operating speed in the column, [m/s]
    heigth_layer : float
        The heigth ligth layer of  the liquid, [m]
    Returns
    -------
    Fr : float
        The Frudo criterion, [dimensionless]
    References
    ----------
    Дытнерский, страница 240
    """
    return w_oper**2 / (g * heigth_layer)


@unitcheck(H_bwplate="m", h_bubble="m", res_unit="m")
def H_separate(H_bwplate, h_bubble):
    """
    Calculates the heigth of separation space.

    Parameters
    ----------
    H_bwplate : float
        The heigth of between plates, [m]
    h_bubble : float
        The heigth of bubble layer, [m]
    Returns
    -------
    H_separate : float
        The heigth of separation space. [m]
    References
    ----------
    Дытнерский, страница 242, формула 6.42
    """    
    return H_bwplate - h_bubble


def h_bubble(heigth_layer, epsi_vapor):
    """
    Calculates the heigth of bubble layer.
    Parameters
    ----------
    epsi_vapor : float
        The vapor content of bubble layer, [dimensionless]
    heigth_layer : float
        The heigth ligth layer of  the liquid, [m]
    Returns
    -------
    h_bubble : float
        The heigth of of bubble layer. [m]
    References
    ----------
    Дытнерский, страница 242
    """
    return heigth_layer / (1 - epsi_vapor)
