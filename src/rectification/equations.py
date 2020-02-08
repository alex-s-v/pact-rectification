from rectification.utils import unitcheck
from scipy.constants import g
import numpy as np


def operating_line(x, R, xp):
    return R * x / (R + 1) + xp / (R + 1)


@unitcheck(L="kg/s", rhol="kg/m**3", Ft="m**2", res_unit="m**3/(m**2*s)")
def U_coef(L, rhol, Ft):
    return L / (rhol * Ft)


@unitcheck(Massl="g/mol", Massh="g/mol", mu_solv="Pa*s", nul="sm**3/mol", nuh="sm**3/mol", res_unit="m**2/s")
def Diff_20(Massl, Massh, A , B, mu_solv, nul, nuh):
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
    Diff_20 : float
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


def Diff_liq(Diff_20, b, t_boil):
    """
    Calculates the diffusion coefficient of liquid phaze.
    Parameters
    ----------
    calc_Diffcoef20 : float
        The diffusion coefficient at 20 degrees celcium, [m**2/s]
    b : float
        The temperature coefficient, [dimensionless]
    t_boil : float
        The boiling temperature of liquid, [degrees celcium]
    Returns
    -------
    Diff_liq : float
        The diffusion coefficient of liquid phaze.
    References
    ----------
    Романков, страница 289, формула 6.23
    """
    return Diff_20 * (1 + b * (t_boil - 20))


@unitcheck(t_boil="K", Massl="g/mol", Massh="g/mol", P_abs="Pa", nul="sm**3/mol", nuh="sm**3/mol", res_unit="m**2/s")
def Diff_vapor(t_boil, P_abs, Massl, Massh, nul, nuh):
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
    Diff_vapor : float
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


def w_oper(rho_mix, rho_vapor):
    """
    Calculates the operating speed in column.
    Parameters
    ----------
    rho_mix : float
        The destiny of mix, [kg/m**3]
    rho_vapor : float
        The destiny of vapor, [kg/m**3]
    Returns
    -------
    w_oper : float
        The operating speed in column. [m]
    References
    ----------
    Дытнерский, страница 205, формула 5.33
    """
    return 0.05 * rho_mix**0.5 / rho_vapor**0.5


@unitcheck(flate_vapor="kg/s", w_oper="m/s", rho_vapor="kg/m**3", res_unit="m")
def d_clm(flate_vapor, w_oper, rho_vapor):
    """
    Calculates the diameter of column.
    Parameters
    ----------
    flate_vapor : float
        The mass flow rate of vapor, [kg/s]
    w_oper : float
        The operating speed in column, [m]
    rho_vapor : 
        The destiny of vapor, [kg/m**3]
    Returns
    -------
    d_clm : float
        The diameter of column. [m]
    References
    ----------
    Дытнерский, страница 197, формула 5.10
    """    
    return (4 * flate_vapor / (np.pi * w_oper * rho_vapor))**0.5


def xfeed_molar(xfeed_mass, Massl, Massh):
    """
    Calculates the molar concetration feed of liquid.
    Parameters
    ----------
    xfeed_mass : float
        The mass concetration feed of liquid, [kg/kg]
    Massl : float
        The molar mass of low-boilling component, [g/mol]
    Massh : float
        The molar mass of high-boilling component, [g/mol]
    Returns
    -------
    xfeed_molar : float
        The molar concetration feed of liquid. [kmol/kmol]
    References
    ----------
    Романков, страница 282, таблица 6.1 / Дытнерский стр.228 формула 6.3
    """ 
    return xfeed_mass * Massh / ((xfeed_mass * Massh) + (Massl - Massl * xfeed_mass))


def xdist_molar(xdist_mass, Massl, Massh):
    """
    Calculates the molar concetration distilliat of liquid.
    Parameters
    ----------
    xdist_mass : float
        The mass concetration distilliat of liquid, [kg/kg]
    Massl : float
        The molar mass of low-boilling component, [g/mol]
    Massh : float
        The molar mass of high-boilling component, [g/mol]
    Returns
    -------
    xdist_molar : float
        The molar concetration distilliat of liquid. [kmol/kmol]
    References
    ----------
    Романков, страница 282, таблица 6.1 / Дытнерский стр.228 формула 6.3
    """ 
    return xdist_mass * Massh / ((xdist_mass * Massh) + (Massl - Massl * xdist_mass))


def xwaste_molar(xwaste_mass, Massl, Massh):
    """
    Calculates the molar concetration waste of liquid.
    Parameters
    ----------
    xwaste_mass : float
        The mass concetration waste of liquid, [kg/kg]
    Massl : float
        The molar mass of low-boilling component, [g/mol]
    Massh : float
        The molar mass of high-boilling component, [g/mol]
    Returns
    -------
    xwaste_molar : float
        The molar concetration waste of liquid. [kmol/kmol]
    References
    ----------
    Романков, страница 282, таблица 6.1 / Дытнерский стр.228 формула 6.3
    """ 
    return xwaste_mass * Massh / ((xwaste_mass * Massh) + (Massl - Massl * xwaste_mass))


def M_feed(xfeed_molar, Massh, Massl):
    """
    Calculates the molar mass feed of liquid.
    Parameters
    ----------
    xfeed_molar : float
        The molar concetration feed of liquid. [kmol/kmol]
    Massl : float
        The molar mass of low-boilling component, [g/mol]
    Massh : float
        The molar mass of high-boilling component, [g/mol]
    Returns
    -------
    M_feed : float
        The molar mass feed of liquid. [kg/kmol]
    References
    ----------
    Дытнерский, стр. 230, формула 6.6
    """ 
    return xfeed_molar * Massl + Massh - Massh * xfeed_molar


def M_dist(xdist_molar, Massh, Massl):
    """
    Calculates the molar mass distilliat of liquid.
    Parameters
    ----------
    xdist_molar : float
        The molar concetration distilliat of liquid. [kmol/kmol]
    Massl : float
        The molar mass of low-boilling component, [g/mol]
    Massh : float
        The molar mass of high-boilling component, [g/mol]
    Returns
    -------
    M_dist : float
        The molar mass distilliat of liquid. [kg/kmol]
    References
    ----------
    Дытнерский, стр. 230, формула 6.6
    """ 
    return xdist_molar * Massl + Massh - Massh * xdist_molar


def M_waste(xwaste_molar, Massh, Massl):
    """
    Calculates the molar mass waste of liquid.
    Parameters
    ----------
    xwaste_molar : float
        The molar concetration waste of liquid. [kmol/kmol]
    Massl : float
        The molar mass of low-boilling component, [g/mol]
    Massh : float
        The molar mass of high-boilling component, [g/mol]
    Returns
    -------
    M_waste : float
        The molar mass waste of liquid. [kg/kmol]
    References
    ----------
    Дытнерский, стр. 230, формула 6.6
    """ 
    return xwaste_molar * Massl + Massh - Massh * xwaste_molar


def W_mass(F_mass, xdist_mass, xfeed_mass, xwaste_mass):
    """
    Calculates the mass flow rate of waste.
    Parameters
    ----------
    F_mass : float
        The mass flow rate of feed [kg/s]
    xdist_mass : float
        The mass concetration distilliat of liquid, [kg/kg]
    xwaste_mass : float
        The mass concetration waste of liquid, [kg/kg]
    xfeed_mass : float
        The mass concetration feed of liquid, [kg/kg]
    Returns
    -------
    W_mass : float
        The mass flow rate of waste. [kg/s]
    References
    ----------
    Дытнерский, стр. 228, формула 6.1
    """ 
    return F_mass * (xdist_mass - xfeed_mass) / (xdist_mass - xwaste_mass)


def P_mass(F_mass, W_mass):
    """
    Calculates the mass flow rate of distilliat.
    Parameters
    ----------
    F_mass : float
        The mass flow rate of feed [kg/s]
    W_mass : float
        The mass flow rate of waste. [kg/s]
    Returns
    -------
    P_mass : float
        The mass flow rate of distilliat. [kg/s]
    References
    ----------
    Дытнерский, стр. 228, формула 6.1
    """    
    return F_mass - W_mass


def F_mol(F_mass, M_feed):
    """
    Calculates the molar flow rate of feed.
    Parameters
    ----------
    F_mass : float
        The mass flow rate of feed [kg/s]
    M_feed : float
        The molar mass feed of liquid. [kg/kmol]
    Returns
    -------
    F_mol : float
        The molar flow rate of feed, [kmol/s]
    References
    ----------
    ???
    """        
    return F_mass / M_feed


def P_mol(P_mass, M_dist):
    """
    Calculates the molar flow rate of dist.
    Parameters
    ----------
    P_mass : float
        The mass flow rate of distilliat, [kg/s]
    M_dist : float
        The molar mass  of distilliat, [kg/kmol]
    Returns
    -------
    P_mol : float
        The molar flow rate of distilliat, [kmol/s]
    References
    ----------
    ???
    """        
    return P_mass / M_dist


def W_mol(W_mass, M_waste):
    """
    Calculates the molar flow rate of waste.
    Parameters
    ----------
    P_mass : float
        The mass flow rate of waste, [kg/s]
    M_waste : float
        The molar mass  of waste, [kg/kmol]
    Returns
    -------
    P_mol : float
        The molar flow rate of waste, [kmol/s]
    References
    ----------
    ???
    """        
    return W_mass / M_waste


def beta_liq(Diff_liq, epsi_vapor, U_coef, heigth_layer, mu_vapor, mu_mix):
    """
    Calculates the coefficient masstransfer of liquid.
    Parameters
    ----------
    Diff_liq : float
        The diffusion coefficient of liquid phaze.
    epsi_vapor : float
        The vapor content of bubble layer, [dimensionless]
    heigth_layer : float
        The heigth ligth layer of  the liquid, [m]
    mu_mix : float
        The mix viscocity of liquid, [Pa/s]
    mu_vapor : float
        The mix viscocity of vapor, [Pa/s]
    U_coef : float
        The specific coefficient for beta_liq equation.
    Returns
    -------
    beta_liq : float
        The coefficient masstransfer of liquid, [m/s]
    References
    ----------
    Дытнерский, формула 6.37, стр.239
    """              
    return 6.24e+5 * (Diff_liq**0.5) * heigth_layer * ((U_coef/(1-epsi_vapor))**0.5) * (mu_vapor / (mu_vapor + mu_mix))**0.5


def beta_vapor(Diff_vapor, w_oper, epsi_vapor, heigth_layer, Fc, mu_vapor, mu_mix):
    """
    Calculates the coefficient masstransfer of vapor.
    Parameters
    ----------
    Diff_vapor : float
        The diffusion coefficient of vapor phaze, []
    epsi_vapor : float
        The vapor content of bubble layer, [dimensionless]
    heigth_layer : float
        The heigth ligth layer of  the liquid, [m]
    mu_mix : float
        The mix viscocity of liquid, [Pa/s]
    mu_vapor : float
        The mix viscocity of vapor, [Pa/s]
    Fc : float
        The free section of a plate, [dimensionless]
    Returns
    -------
    beta_vapor : float
        The coefficient masstransfer of vapor, [m/s]
    References
    ----------
    Дытнерский, формула 6.38, стр.239
    """
    return 6.24e+5 * Diff_vapor**0.5 * ((w_oper/epsi_vapor)**0.5) * heigth_layer * Fc * ((mu_vapor / (mu_vapor + mu_mix))**0.5)