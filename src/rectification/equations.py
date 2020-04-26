from rectification.utils import unitcheck
from scipy.constants import g
import numpy as np
Vm = 22.4


def operating_line(x, R, xp):
    return R * x / (R + 1) + xp / (R + 1)


@unitcheck(L="kg/s", rhol="kg/m**3", Ft="m**2", res_unit="m**3/(m**2*s)")
def U_coef(L, rhol, Ft):
    return L / (rhol * Ft)

#region Calculating the coefficients of diffusion
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


#endregion
#region Calculating of bubble layer
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


#endregion
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


#region Calculating flow rate, concentration, mass of distilliat, feed and waste

def xfeed_molar(xfeed_mass, Massl, Massh):
    """
    Calculates the molar concetration feed of liquid.
    Parameters
    ----------
    xfeed_mass : float
        The mass concetration feed of liquid, [kg/kg]
    Massl : float
        The molar mass of low-boilling component, [kg/kmol]
    Massh : float
        The molar mass of high-boilling component, [kg/kmol]
    Returns
    -------
    xfeed_molar : float
        The molar concetration feed of liquid. [kmol/kmol]
    References
    ----------
    Романков, страница 282, таблица 6.1 / Дытнерский стр.228 формула 6.3
    """ 
    return xfeed_mass * Massh / ((xfeed_mass * Massh) + (Massl - Massl * xfeed_mass))


@unitcheck(xdist_mass="kg/kg", Massl="kg/kmol", Massh="kg/kmol", res_unit="kmol/kmol")
def xdist_molar(xdist_mass, Massl, Massh):
    """
    Calculates the molar concetration distilliat of liquid.
    Parameters
    ----------
    xdist_mass : float
        The mass concetration distilliat of liquid, [kg/kg]
    Massl : float
        The molar mass of low-boilling component, [kg/kmol]
    Massh : float
        The molar mass of high-boilling component, [kg/kmol]
    Returns
    -------
    xdist_molar : float
        The molar concetration distilliat of liquid. [kmol/kmol]
    References
    ----------
    Романков, страница 282, таблица 6.1 / Дытнерский стр.228 формула 6.3
    """ 
    return xdist_mass * Massh / ((xdist_mass * Massh) + (Massl - Massl * xdist_mass))


@unitcheck(xwaste_mass="kg/kg", Massl="kg/kmol", Massh="kg/kmol", res_unit="kmol/kmol")
def xwaste_molar(xwaste_mass, Massl, Massh):
    """
    Calculates the molar concetration waste of liquid.
    Parameters
    ----------
    xwaste_mass : float
        The mass concetration waste of liquid, [kg/kg]
    Massl : float
        The molar mass of low-boilling component, [kg/kmol]
    Massh : float
        The molar mass of high-boilling component, [kg/kmol]
    Returns
    -------
    xwaste_molar : float
        The molar concetration waste of liquid. [kmol/kmol]
    References
    ----------
    Романков, страница 282, таблица 6.1 / Дытнерский стр.228 формула 6.3
    """ 
    return xwaste_mass * Massh / ((xwaste_mass * Massh) + (Massl - Massl * xwaste_mass))





@unitcheck(F_mass="kg/s", Ь_ауув="kg/kmol", res_unit="kg/s")
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


#endregion
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


#region The heat exchanger of feed
def deltaT_diff(deltaT_larger, deltaT_less):
    """
    Calculates the difference of temperatures.
    Parameters
    ----------
    deltaT_larger : float
        The difference temperatures between the temperature of vapor and the lowest temperature of liquid , [degrees celcium]
    deltaT_less : float
        The difference temperatures between the temperature of vapor and the biggest temperature of liquid , [degrees celcium]
    Returns
    -------
    deltaT_diff : float
        The coefficient difference of temperatures, [degrees celcium]
    References
    ----------
    Романков, формула 4.78, стр.169
    """   
    if deltaT_larger/deltaT_less < 2:
        return (deltaT_larger + deltaT_less) / 2
    if deltaT_larger/deltaT_less > 2:
        return (deltaT_larger - deltaT_less) / np.log(deltaT_larger / deltaT_less)


def deltaT_larger(t_vapor, tinit_mix):
    """
    Calculates the difference temperatures between the temperature of vapor and the lowest temperature of liquid.
    Parameters
    ----------
    t_vapor : float
        The temperature of vapor, [degrees celcium]
    tinit_mix : float
        The initial temperature of mix, [degrees celcium]
    Returns
    -------
    deltaT_larger : float
        The difference temperatures between the temperature of vapor and the lowest temperature of liquid , [degrees celcium]
    References
    ----------
    Романков, формула 4.78, стр.169
    """       
    return t_vapor - tinit_mix


def deltaT_less(t_vapor, tboil_mix):
    """
    Calculates the difference temperatures between the temperature of vapor and the lowest temperature of liquid.
    Parameters
    ----------
    t_vapor : float
        The temperature of vapor, [degrees celcium]
    tboil_mix : float
        The boilling temperature of mix, [degrees celcium]
    Returns
    -------
    deltaT_less : float
        The difference temperatures between the temperature of vapor and the biggest temperature of liquid, [degrees celcium]
    References
    ----------
    Романков, формула 4.78, стр.169
    """       
    return t_vapor - tboil_mix


@unitcheck(deltaT="degrees celcium", F_mass="kg/s", Cp="J/(kg * degrees celcium", res_unit="W")
def Qload_feed(deltaT, F_mass, Cp, phi_vapor, Feed_vaporazation):
    """
    Calculates the heat load of heat exchanger.
    Parameters
    ----------
    deltaT : float
        The difference temperatures between the initial and the ultimate temperatures of mix , [degrees celcium]
    F_mass : float
        The mass flow rate of feed [kg/s]
    Cp : float
        The heat capacity of mix [J/(kg * degrees C)]
    phi_vapor: float
        The part of vapor in feed, [dimensionless]
    Feed_vaporazation : float
        The heat vaporazation of mix, [J/kg]
    Returns
    -------
    Qload_feed : float
        The heat load of heat exchanger, [W] , [J/s]
    References
    ----------
    Дытнерский, формула 2.2, стр.45
    """   
    return F_mass * (Cp * deltaT + phi_vapor * Feed_vaporazation)


def deltaT(tinit_mix, tboil_mix):
    """
    Calculates the difference temperatures between the initial and the ultimate temperatures of mix , [degrees celcium]
    Parameters
    ----------
    tinit_mix : float
        The initial temperature of mix, [degrees celcium]
    tboil_mix : float
        The boilling temperature of mix, [degrees celcium]
    Returns
    -------
    deltaT : float
        The difference temperatures between the initial and the ultimate temperatures of mix , [degrees celcium]
    References
    ----------
    Дытнерский, формула 2.2, стр.45
    """       
    return tinit_mix - tboil_mix


@unitcheck(Qload_feed="W", deltaT_diff="degrees celcium",  Kt_approx="W/(m**2 * degrees celcium", res_unit="m**2")
def A_approx(Qload_feed, deltaT_diff, Kt_approx):
    """
    Calculates the approximate heatransfer area.
    Parameters
    ----------
    Qload_feed : float
        The heat load of heat exchanger, [W] , [J/s]
    deltaT_diff : float
        The coefficient difference of temperatures, [degrees celcium]
    Kt_approx : float
        The heatransfer coefficient [W/(m**2 * degrees celcium)]
    Returns
    -------
    A_approx : float
        The approximate heatransfer area, [m**2]
    References
    ----------
    Романков, формула 4.72, стр.168
    """           
    return Qload_feed / (deltaT_diff * Kt_approx)


def flatesteam_feed(Qload_feed, vaporazation_steam):
    """
    Calculates the flow rate steam of boiler.
    Parameters
    ----------
    Qload_feed : float
        The heat load feed of heat exchanger, [W] [J/s]
    vaporazation_steam : float
        The heat vaporazation of dist [J/kg]
    Returns
    -------
    flatesteam_feed : float
        The flow rate steam of feed heat exchanger, [W] [J/s]
    References
    ----------
    Дытнерский, формула 2.3, стр.45
    """               
    return Qload_feed / vaporazation_steam
    

def Re(F_mass, z_way, d_inner, n_pipe, mu_mix):
    """
    Calculates the Reynold criterion.
    Parameters
    ----------
    F_mass : float
        The mass flow rate of feed [kg/s]
    z_way : float
        The number of ways in heat exchanger [dimensionless]
    d_inner : float
        The diametr of inner pipe, [m]
    n_pipe : float
        The number of pipes in heat exchanger, [dimensionless]
    mu_mix : float
        The mix viscocity of liquid, [Pa/s]
    Returns
    -------
    Re : float
        The Reynold criterion, [dimensionless]
    References
    ----------
    &&&&&&
    """         
    return 0.785 * F_mass * z_way / (d_inner * n_pipe * mu_mix)


def Pr(C_capacity, mu_mix, lyambda_feed):
    """
    Calculates the Prandtl criterion.
    Parameters
    ----------
    C_capacity : float
        The heat capacity of mix [J/(kg * degrees celcium)]
    mu_mix : float
        The mix viscocity of liquid, [Pa/s]
    lyambda_feed : float
        The thermal conductivity of feed, [W / (m * degreec celcium)]
    Returns
    -------
    Pr : float
        The Prandtl criterion, [dimensionless]
    References
    ----------
    Романков, формула 4.12, стр.151
    """        
    return C_capacity * mu_mix / lyambda_feed


def Nu(Re, Pr):
    """
    Calculates the Nusselt criterion.
    Parameters
    ----------
    Re : float
        The Reynold criterion, [dimensionless]
    Pr : float
        The Prandtl criterion, [dimensionless]
    Returns
    -------
    Nu : float
        The Nusselt criterion, [dimensionless]
    References
    ----------
    Романков, формула 4.17, стр.152
    """      
    return 0.021 * (Re**0.8) * (Pr**0.4)


def alpha_liq(Nu, lyambda_feed, d_inner):
    """
    Calculates the coefficent of heat transfer(alpha) from liquid to wall of pipe.
    Parameters
    ----------
    Nu : float
        The Nusselt criterion, [dimensionless]
    lyambda_feed : float
        The thermal conductivity of feed, [W / (m * degreec celcium)]
    d_inner : float
        The diametr of inner pipe, [m]
    Returns
    -------
    alpha_liq : float
        The coefficent of heat transfer(alpha), [W / (m**2 * degrees celcium)]
    References
    ----------
    Романков, формула 4.11, стр.150
    """          
    return Nu * lyambda_feed / d_inner


def alpha_vap(lyambda_cond, rho_cond, mu_cond, flatesteam_feed, n_pipe, d_outside):
    """
    Calculates the coefficent of heat transfer(alpha) from steam to wall of pipe.
    Parameters
    ----------
    lyambda_cond : float
        The thermal conducivity of condensate, [W / (m * degrees celcium)]
    rho_cond : float
        The destiny of condensate, [kg / m**3]
    mu_cond : float
        The viscosity of condensate, [Pa / s]
    flatesteam_feed : float
        The flow rate steam of feed heat exchanger, [W] [J/s]
    n_pipe : float
        The number of pipes in heat exchanger, [dimensionless]
    d_outside : float
        The outside diameter of pipe, [m]    
    Returns
    -------
    alpha_vap : float
        The coefficent of heat transfer(alpha) from steam to wall of pipe, [W / (m**2 * degrees celcium)]
    References
    ----------
    Дытнерский, формула 2.24, стр.53
    """              
    return lyambda_cond * 3.78 * ((rho_cond**2)* n_pipe * d_outside / (mu_cond * flatesteam_feed))**(1/3)


def sigma_thermpollution(pollution_1, pollution_2, sigma, lyambda_wall):
    """
    Calculates the sum of thermal pollutions.
    Parameters
    ----------
    pollution_1 : float
        The thermal pollution of the first coolant, [m**2 * degrees celcium / W]
    pollution_2 : float
        The thermal pollution of the second coolant, [m**2 * degrees celcium / W]
    sigma : float
        The thickness of pipe wall, [m]
    lyambda_wall : float
        The thermal conducivity of wall, [W / (m * degrees celcium)]  
    Returns
    -------
    sigma_thermpollution : float
        The the sum of thermal pollutions, [m**2 * degrees celcium / W]
    References
    ----------
    &&&&&
    """                  
    return (sigma / lyambda_wall) + (1 / pollution_1) + (1 / pollution_2)


def Kt_real(alpha_liq, alpha_vap, sigma_thermpollution):
    """
    Calculates the coefficient of heat transfer (Kt).
    Parameters
    ----------
    alpha_liq : float
        The coefficent of heat transfer(alpha), [W / (m**2 * degrees celcium)]
    alpha_vap : float
        The coefficent of heat transfer(alpha) from steam to wall of pipe, [W / (m**2 * degrees celcium)]
    sigma_thermpollution : float
        The the sum of thermal pollutions, [m**2 * degrees celcium / W]  
    Returns
    -------
    Kt_real : float
        The coefficient of heat transfer (Kt), [W / (m**2 * degrees celcium)]
    References
    ----------
    Романков, формула 4.74, стр. 168
    """      
    return ((1 / alpha_liq) + (1 / alpha_vap) + (sigma_thermpollution))**-1

@unitcheck(Qload_feed="W", deltaT_diff="degrees celcium",  Kt_approx="W/(m**2 * degrees celcium", res_unit="m**2")
def A_real(Qload_feed, Kt_real, deltaT_diff):
    """
    Calculates the real heatransfer area.
    Parameters
    ----------
    Qload_feed : float
        The heat load of heat exchanger, [W] , [J/s]
    deltaT_diff : float
        The coefficient difference of temperatures, [degrees celcium]
    Kt_real : float
        The heat ransfer coefficient [W/(m**2 * degrees celcium)]
    Returns
    -------
    A_real : float
        The real heat ransfer area, [m**2]
    References
    ----------
    Романков, формула 4.72, стр.168
    """      
    return Qload_feed / (Kt_real * deltaT_diff)


def surface_margin (A_approx, A_real):
    """
    Calculates the surface margin.
    Parameters
    ----------
    A_approximate : float
        The approximate heat ransfer area, [m**2]
    A_real : float
        The real heat transfer area, [m**2]
    Returns
    -------
    surface_margin : float
        The surface margin [%]
    References
    ----------
    &&&&
    """          
    return (A_approx - A_real) * 100 / A_approx
#endregion


#region The boiler
def Qload_boiler(W_mass, Cw, tw, P_mass, R, dist_vaporazation, F_mass, Cf, tf, Cp, tp, Q_loss):
    """
    Calculates the heat load of boiler.
    Parameters
    ----------
    W_mass : float
        The mass flow rate of waste, [kg/s]
    Cw : float
        The heat capacity of waste [J/(kg * degrees C)]
    tw : float
        The boiling temperature of waste [degrees celcium]
    P_mass : float
        The mass flow rate of waste, [kg/s]
    R : float
        The number of reflux [dimensionless]
    dist_vaporazition : float
        The heat vaporazation of dist [J/kg]
    F_mass : float
        The mass flow rate of feed, [kg/s]
    Cf : float
        The heat capacity of feed, [J/(kg * degrees C)]
    tf : float
        The boiling temperature of feed, [degrees celcium]    
    Cp : float
        The heat capacity of mix [J/(kg * degrees C)]
    tp : float
        The boiling temperature of dist, [degrees celcium]     
    Q_loss : float
        The heat loss, [W] [J/s]
    Returns
    -------
    Qload_boiler : float
        The heat load of boiler, [W] [J/s]
    References
    ----------
    14.	Дытнерский Ю.И. Процессы и аппараты химической технологии. Учебник для вузов. Часть 2. стр 123-124
    """     
    return W_mass * Cw * tw + P_mass * (R + 1) * dist_vaporazation - F_mass * Cf * tf + P_mass * Cp * tp + Q_loss


def flatesteam_boil(Qload_boiler, vaporazation_steam):
    """
    Calculates the flow rate steam of boiler.
    Parameters
    ----------
    Qload_boiler : float
        The heat load of boiler, [W] [J/s]
    vaporazation_steam : float
        The heat vaporazation of dist [J/kg]
    Returns
    -------
    flatesteam_boil : float
        The flow rate steam of boiler, [W] [J/s]
    References
    ----------
    Дытнерский, формула 2.3, стр.45
    """               
    return Qload_boiler / vaporazation_steam


def deltaT_boil(t_boil, t_cond):
    """
    Calculates the temperature difference of boiler.
    Parameters
    ----------
    t_boil : float
        The boiling temperature of liquid, [С]
    t_cond : float
        The condensation temperature of steam, [C]
    Returns
    -------
    deltaT_boil : float
        The temperature difference of boiler, [C]
    References
    ----------
    &&&&&
    """     
    return t_cond - t_boil


@unitcheck(Qload_boil="W", deltaT_boil="degrees celcium",  Kt_approx="W/(m**2 * degrees celcium", res_unit="m**2")
def Aboil_approx(Qload_boiler, deltaT_boil, Kt_approx):
    """
    Calculates the approximate heatransfer area.
    Parameters
    ----------
    Qload_boiler : float
        The heat load of boiler, [W] , [J/s]
    deltaT_boil : float
        The temperature difference of boiler, [C]
    Kt_approx : float
        The heatransfer coefficient [W/(m**2 * degrees celcium)]
    Returns
    -------
    Aboil_approx : float
        The approximate heatransfer area of boiler, [m**2]
    References
    ----------
    Романков, формула 4.72, стр.168
    """           
    return Qload_boiler / (deltaT_boil * Kt_approx)


#endregion


#region Calculating the dephlegmator
def deltaT_dephleg(deltaT_larger_deph, deltaT_less_deph):
    """
    Calculates the difference of temperatures.
    Parameters
    ----------
    deltaT_larger_deph : float
        The difference temperatures between the temperature of vapor and the lowest temperature of liquid , [degrees celcium]
    deltaT_less_deph : float
        The difference temperatures between the temperature of vapor and the biggest temperature of liquid , [degrees celcium]
    Returns
    -------
    deltaT_depleg : float
        The coefficient difference of temperatures, [degrees celcium]
    References
    ----------
    Романков, формула 4.78, стр.169
    """   
    if deltaT_larger_deph / deltaT_less_deph < 2:
        return (deltaT_larger_deph + deltaT_less_deph) / 2
    if deltaT_larger_deph / deltaT_less_deph > 2:
        return (deltaT_larger_deph - deltaT_less_deph) / np.log(deltaT_larger_deph / deltaT_less_deph)


def deltaT_larger_deph (t_cond, t_dist):
    """
    Calculates the difference temperatures between the temperature of vapor and the lowest temperature of liquid.
    Parameters
    ----------
    t_cond : float
        The temperature of condensation, [C]
    tinit_water : float
        The initial temperature of cool water, [C]
    Returns
    -------
    deltaT_larger_deph : float
        The difference temperatures between the temperature of condensation and the lowest temperature of water , [C]
    References
    ----------
    Романков, формула 4.78, стр.169
    """       
    return t_cond - t_water


def deltaT_less_deph(t_cond, tboil_mix):
    """
    Calculates the difference temperatures between the temperature of condensation and the highest temperature of water.
    Parameters
    ----------
    t_cond : float
        The temperature condensation of dist, [C]
    tulti_water : float
        The ultimate temperature of cool water, [C]
    Returns
    -------
    deltaT_less_deph : float
        The difference temperatures between the temperature condensation of dist and the highest temperature of water, [C]
    References
    ----------
    Романков, формула 4.78, стр.169
    """       
    return t_cond - tulti_water


@unitcheck(P_mass="kg/s", rdist_vaporazation="J/kg", res_unit="W")
def Qload_deph(P_mass, rdist_vaporazation, R):
    """
    Calculates the heat load of dephlegmator.
    Parameters
    ----------
    P_mass : float
        The mass flow rate of dist , [kg/s]
    R : float
        The reflux number [dimensionless]
    rdist_vaporazation : float
        The heat vaporazation of dist, [J/kg]
    Returns
    -------
    Qload_deph : float
        The heat load of dephlegmator, [W] , [J/s]
    References
    ----------
    Дытнерский, формула 2.2, стр.45
    """   
    return P_mass * (R + 1) * rdist_vaporazation


def deltaT_deph(tinit_water, tulti_water):
    """
    Calculates the difference temperatures between the initial and the ultimate temperatures of cool water, [C]
    Parameters
    ----------
    tinit_water : float
        The initial temperature of cool water, [C]
    tulti_water : float
        The ultimate temperature of cool water, [C]
    Returns
    -------
    deltaT_deph : float
        Calculates the difference temperatures between the initial and the ultimate temperatures of cool water, [C]
    References
    ----------
    Дытнерский, формула 2.2, стр.45
    """       
    return tinit_water - tulti_water


@unitcheck(Qload_feed="W", deltaT_diff="degrees celcium",  Kt_approx="W/(m**2 * degrees celcium", res_unit="m**2")
def A_approx_deph(Qload_deph, deltaT_diff_deph, Kt_approx):
    """
    Calculates the approximate heatransfer area.
    Parameters
    ----------
    Qload_deph : float
        The heat load of dephlegmator, [W] , [J/s]
    deltaT_diff_deph : float
        The coefficient difference of temperatures, [C]
    Kt_approx : float
        The heatransfer coefficient [W/(m**2 * degrees celcium)]
    Returns
    -------
    A_approx_deph : float
        The approximate heatransfer area of dephlegmator, [m**2]
    References
    ----------
    Романков, формула 4.72, стр.168
    """           
    return Qload_deph / (deltaT_diff_deph * Kt_approx)


def flatewater_deph(Qload_deph, Cp, deltaT_deph):
    """
    Calculates the flow rate steam of boiler.
    Parameters
    ----------
    Qload_depg : float
        The heat load feed of dephlegmator, [W] [J/s]
    Cp : float
        The heat capacity of mix [J/(kg * degrees C)]
    deltaT_deph : float
        Calculates the difference temperatures between the initial and the ultimate temperatures of cool water, [C]
    Returns
    -------
    flatewater_deph : float
        The flow rate cool water of dephlegmator, [W] [J/s]
    References
    ----------
    Дытнерский, формула 2.3, стр.45
    """               
    return Qload_deph / Cp * deltaT_deph
    

def Re_deph(flatewater_deph, z_way, d_inner, n_pipe, mu_cool_water):
    """
    Calculates the Reynold criterion.
    Parameters
    ----------
    flatewater_deph : float
        The flow rate cool water of dephlegmator, [W] [J/s]
    z_way : float
        The number of ways in dephlegmator [dimensionless]
    d_inner : float
        The diametr of inner pipe, [m]
    n_pipe : float
        The number of pipes in heat dephlegmator, [dimensionless]
    mu_cool_water : float
        The viscocity of cool water, [Pa/s]
    Returns
    -------
    Re : float
        The Reynold criterion, [dimensionless]
    References
    ----------
    &&&&&&
    """         
    return 0.785 * flatewater_deph * z_way / (d_inner * n_pipe * mu_cool_water)


def Pr_deph(C_capacity_water, mu_water, lyambda_water):
    """
    Calculates the Prandtl criterion.
    Parameters
    ----------
    C_capacity_water : float
        The heat capacity of cool water [J/(kg * degrees celcium)]
    mu_water : float
        The viscocity of cool water, [Pa/s]
    lyambda_water : float
        The thermal conductivity of cool water, [W / (m * degreec celcium)]
    Returns
    -------
    Pr_deph : float
        The Prandtl criterion, [dimensionless]
    References
    ----------
    Романков, формула 4.12, стр.151
    """        
    return C_capacity_water * mu_cool_water / lyambda_water


def Nu_deph(Re_deph, Pr_deph):
    """
    Calculates the Nusselt criterion.
    Parameters
    ----------
    Re_deph : float
        The Reynold criterion, [dimensionless]
    Pr_deph : float
        The Prandtl criterion, [dimensionless]
    Returns
    -------
    Nu_deph : float
        The Nusselt criterion, [dimensionless]
    References
    ----------
    Романков, формула 4.17, стр.152
    """      
    return 0.021 * (Re_deph**0.8) * (Pr_deph**0.4)


def alpha_liq_deph(Nu_deph, lyambda_water, d_inner):
    """
    Calculates the coefficent of heat transfer(alpha) from liquid to wall of pipe.
    Parameters
    ----------
    Nu_deph : float
        The Nusselt criterion, [dimensionless]
    lyambda_water : float
        The thermal conductivity of cool water, [W / (m * degreec celcium)]
    d_inner : float
        The diametr of inner pipe, [m]
    Returns
    -------
    alpha_liq_deph : float
        The coefficent of heat transfer(alpha), [W / (m**2 * degrees celcium)]
    References
    ----------
    Романков, формула 4.11, стр.150
    """          
    return Nu_deph * lyambda_water / d_inner


def alpha_cond_deph(lyambda_cond_dist, rho_cond_dist, mu_cond_dist, P_mass, n_pipe_deph, L_pipe_deph):
    """
    Calculates the coefficent of heat transfer(alpha) from steam to wall of pipe.
    Parameters
    ----------
    lyambda_cond_dist : float
        The thermal conducivity condensate of dist , [W / (m * degrees celcium)]
    rho_cond_dist : float
        The destiny condensate of dist, [kg / m**3]
    mu_cond_dist : float
        The viscosity condensate of dist, [Pa / s]
    P_mass : float
        The mass flow rate of distilliat, [kg/s]
    n_pipe_deph : float
        The number of pipes in dephlegmator, [dimensionless]
    L_pipe_deph : float
        The length of pipes, [m]    
    Returns
    -------
    alpha_cond_deph : float
        The coefficent of heat transfer(alpha) from steam to wall of pipe, [W / (m**2 * degrees celcium)]
    References
    ----------
    Дытнерский, формула 2.24, стр.53
    """        
        if n < 100:      
            return lyambda_cond_water * 3.78 * ((rho_cond_water**2)* n_pipe_deph * L_pipe_deph / (mu_cond_dist * P_mass))**(1/3)
        if n > 100:
            return 0.6 * lyambda_cond_water * 3.78 * ((rho_cond_water**2)* n_pipe_deph * L_pipe_deph / (mu_cond_dist * P_mass))**(1/3)


def sigma_thermpollution_deph(pollution_1_deph, pollution_2_deph, sigma_deph, lyambda_wall_deph):
    """
    Calculates the sum of thermal pollutions.
    Parameters
    ----------
    pollution_1_deph : float
        The thermal pollution of the first coolant, [m**2 * degrees celcium / W]
    pollution_2_deph : float
        The thermal pollution of the second coolant, [m**2 * degrees celcium / W]
    sigma_deph : float
        The thickness of pipe wall, [m]
    lyambda_wall_deph : float
        The thermal conducivity of wall, [W / (m * degrees celcium)]  
    Returns
    -------
    sigma_thermpollution_deph : float
        The the sum of thermal pollutions, [m**2 * degrees celcium / W]
    References
    ----------
    &&&&&
    """                  
    return (sigma_deph / lyambda_wall_deph) + (1 / pollution_1_deph) + (1 / pollution_2_deph)


def Kt_real_deph(alpha_liq_deph, alpha_cond_deph, sigma_thermpollution_deph):
    """
    Calculates the coefficient of heat transfer (Kt).
    Parameters
    ----------
    alpha_liq_deph : float
        The coefficent of heat transfer(alpha), [W / (m**2 * degrees celcium)]
    alpha_cond_deph : float
        The coefficent of heat transfer(alpha) from steam to wall of pipe, [W / (m**2 * degrees celcium)]
    sigma_thermpollution_deph : float
        The the sum of thermal pollutions, [m**2 * degrees celcium / W]  
    Returns
    -------
    Kt_real : float
        The coefficient of heat transfer (Kt), [W / (m**2 * degrees celcium)]
    References
    ----------
    Романков, формула 4.74, стр. 168
    """      
    return ((1 / alpha_liq_deph) + (1 / alpha_cond_deph) + (sigma_thermpollution_deph))**-1

@unitcheck(Qload_deph="W", deltaT_diff_deph="degrees celcium",  Kt_real_deph="W/(m**2 * degrees celcium", res_unit="m**2")
def A_real_deph(Qload_deph, Kt_real_deph, deltaT_diff_deph):
    """
    Calculates the real heatransfer area.
    Parameters
    ----------
    Qload_deph : float
        The heat load of dephlegmator, [W] , [J/s]
    deltaT_diff_deph : float
        The coefficient difference of temperatures, [degrees celcium]
    Kt_real_deph : float
        The heat ransfer coefficient [W/(m**2 * degrees celcium)]
    Returns
    -------
    A_real_deph : float
        The real heat ransfer area, [m**2]
    References
    ----------
    Романков, формула 4.72, стр.168
    """      
    return Qload_deph / (Kt_real_deph * deltaT_diff_deph)


def surface_margin_deph (A_approx_deph, A_real_deph):
    """
    Calculates the surface margin.
    Parameters
    ----------
    A_approximate_deph : float
        The approximate heat ransfer area, [m**2]
    A_real_deph : float
        The real heat transfer area, [m**2]
    Returns
    -------
    surface_margin_deph : float
        The surface margin [%]
    References
    ----------
    &&&&
    """          
    return (A_approx_deph - A_real_deph) * 100 / A_approx_deph
#endregion

#region Calculating of pump for Feed
def flate_pump_feed(F_mass, rho_F):
    """
    Calculates the flow rate pump for Feed.
    Parameters
    ----------
    F_mass : float
        The flow rate of Feed, [kg / s]
    rho_F : float
        The density of feed, [kg / m**3]
    Returns
    -------
    flate_pump_feed : float
        The flow rate pump for Feed [m**3 / s]
    References
    ----------
    &&&&
    """       
    return (F_mass / rho_F)


def high_geometric(plate_low, lenght_bottom, lenght_support, H_bwplate):
    """
    Calculates the  geometric heigth.
    Parameters
    ----------
    H_bwplate : float
        The heigth of between plates, [m]
    plate_low : float
        The quantity of plates of low column, [dismensionless]
    lenght_bottom : float
        The length between bottom and lowest plate of column, [m]
    lenght_support : float
        The length of support, [m]    
    Returns
    -------
    high_geometric : float
        The geometric heigth, [m]
    References
    ----------
    &&&&
    """    
    return (lenght_plate_low * H_bwplate + lenght_bottom + lenght_support)

def head_pump(high_geometric, hydraulic_losses):
    """
    Calculates the  hydraulic head of pump.
    Parameters
    ----------
    high_geometric : float
        The geometric heigth, [m]
    hydraulic_losses : float
        The hydraulic losses  of pump, [m]
    Returns
    -------
    head_pump : float
        The  hydraulic head of pump, [m]
    References
    ----------
    &&&&
    """       
    return (high_geometric + high_losses)


def power_pump(flate_pump_feed, rho_F, g, head_pump, ECE_motor, ECE_trans):
    """
    Calculates the  power of pump.
    Parameters
    ----------
    flate_pump_feed : float
        The flow rate pump for Feed [m**3 / h]
    rho_F : float
        The density of feed, [kg / m**3]
    head_pump : float
        The  hydraulic head of pump, [m]
    ECE_motor : float
        The energy conversion efficiency of motor, [dismensionless]
    ECE_trans : float
        The energy conversion efficiency of transfer, [dismensionless]       
    Returns
    -------
    power_pump : float
        The  power of pump, [kW]
    References
    ----------
    &&&&
    """   
    return (flate_pump_feed * rho_F * g * head_pump / (ECE_motor * ECE_trans))

#endregion


#region The maximum theoretical suction height
def heigth_max_suction(Pa, rho_F, g, P_satur_vapor, speed_suction, hydraulic_losses_suct, heigth_cavitation):
    """
    Calculates the maximum theoretical suction height.
    Parameters
    ----------
    Pa : float
        The atmosphere pressure [Pa]
    rho_F : float
        The density of feed, [kg / m**3]
    P_satur_vapor : float
        The  pressure of saturated vapor, [Pa]
    speed_suction : float
        The speed of suction line , [m/s]
    hydraulic_losses_suct : float
        The hydraulic losses of suction line, [m]      
    heigth_cavitation : float
        The losses due to cavitation, [m]  
    Returns
    -------
    heigth_max_suction : float
        The maximum theoretical suction height, [m]
    References
    ----------
    &&&&
    """   
    return ((Pa/(rho_F * g) - ((P_satur_vapor)/(rho_F * g) + ((speed_suction / (2 * g))) + hydraulic_losses_satur + heigth_cavitation))


def heigth_cavitation(flate_pump_feed, n_turnover):
    """
    Calculates the losses due to cavitation.
    Parameters
    ----------
    flate_pump_feed : float
        The flow rate pump for Feed [m**3 / s]
    n_turnover : float
        The quantity turnover of pump, [ turnover / s] 
    Returns
    -------
    heigth_cavitation : float
        The losses due to cavitation, [m]
    References
    ----------
    &&&&
    """      
    return (0,3* (flate_pump_feed * n_turnover**2)**(2/3))


def hydraulic_losses_suct(dzeta_enter, dzeta_turn90, n_turn90, dzeta_ventil, n_ventil, g, speed_suction):
    """
    Calculates the hydraulic losses of suction line.
    Parameters
    ----------
    dzeta_enter : float
        The local resistance of tube enter, [m]
    dzeta_turn90 : float
        The local resistance of turn to 90 degrees, [m]
    n_turn90 : float
        The quantity of turn to 90 degrees, [dismensionless]
    dzeta_ventil : float
        The local resistance of ventil on sunction line, [m]
    n_ventil : float
        The quantity of ventil on suction line, [dismensionless]
    speed_suction : float
        The speed of suction line , [m/s] 
    Returns
    -------
    hydraulic_losses_suct : float
        The hydraulic losses of suction line, [m]
    References
    ----------
    &&&&
    """         
    return ((dzeta_enter + dzeta_turn90 + dzeta_ventil)*speed_suction/(2 * g))


#endregion


#region The material balance of masstransfer 
@unitcheck(xf_mass="kg/kg", M_lc="kg/kmol", M_hc="kg/kmol", res_unit="kmol/kmol")
def xf_mol(xf_mass, M_lc, M_hc):
    """
    Calculates the mol composition of feed.
    Parameters
    ----------
    xf_mass : float
        The mass concentration low component of feed, [kg/kg]
    M_lc : float
        The molar mass of low-boilling component, [kg/kmol]
    M_hc : float
        The molar mass of high-boilling component, [kg/kmol]
    Returns
    -------
    xf_mol : float
        The mol concentration low component of feed, [kmol/kmol]
    References
    ----------
    Романков, стр.283, табл. 6.2. 
    """      
    return (xf_mass * M_hc) / ((xf_mass * M_hc) + (M_lc - M_lc * xf_mass))


@unitcheck(xp_mass="kg/kg", M_lc="kg/kmol", M_hc="kg/kmol", res_unit="kmol/kmol")
def xp_mol(xp_mass, M_lc, M_hc):
    """
    Calculates the mol composition of distilliat.
    Parameters
    ----------
    xp_mass : float
        The mass concentration low component of distilliat, [kg/kg]
    M_lc : float
        The molar mass of low-boilling component, [kg/kmol]
    M_hc : float
        The molar mass of high-boilling component, [kg/kmol]
    Returns
    -------
    xp_mol : float
        The mol concentration low component of distilliat, [kmol/kmol]
    References
    ----------
    Романков, стр.283, табл. 6.2. 
    """      
    return (xp_mass * M_hc) / ((xp_mass * M_hc) + (M_lc - M_lc * xp_mass))


@unitcheck(xw_mass="kg/kg", M_lc="kg/kmol", M_hc="kg/kmol", res_unit="kmol/kmol")
def xw_mol(xw_mass, M_lc, M_hc):
    """
    Calculates the mol composition of waste.
    Parameters
    ----------
    xw_mass : float
        The mass concentration low component of waste, [kg/kg]
    M_lc : float
        The molar mass of low-boilling component, [kg/kmol]
    M_hc : float
        The molar mass of high-boilling component, [kg/kmol]
    Returns
    -------
    xw_mol : float
        The mol concentration low component of waste, [kmol/kmol]
    References
    ----------
    Романков, стр.283, табл. 6.2. 
    """      
    return (xw_mass * M_hc) / ((xw_mass * M_hc) + (M_lc - M_lc * xw_mass))


@unitcheck(xf_mol="kmol/kmol", M_lc="kg/kmol", M_hc="kg/kmol", res_unit="kg/kmol")
def M_feed(xf_mol, M_lc, M_hc):
    """
    Calculates the molar mass of feed.
    Parameters
    ----------
    xf_mol : float
        The mol concentration low component of feed, [kmol/kmol]
    M_lc : float
        The molar mass of low-boilling component, [kg/kmol]
    M_hc : float
        The molar mass of high-boilling component, [kg/kmol]
    Returns
    -------
    M_feed : float
        The molar mass of feed, [kg/kmol]
    References
    ----------
    Дытнерский, стр. 230, формула 6.6
    """ 
    return (M_lc * xf_mol + M_hc * (1 - xf_mol))


@unitcheck(xp_mol="kmol/kmol", M_lc="kg/kmol", M_hc="kg/kmol", res_unit="kg/kmol")
def M_dist(xp_mol, M_lc, M_hc):
    """
    Calculates the molar mass of distilliat.
    Parameters
    ----------
    xp_mol : float
        The mol concentration low component of distilliat, [kmol/kmol]
    M_lc : float
        The molar mass of low-boilling component, [kg/kmol]
    M_hc : float
        The molar mass of high-boilling component, [kg/kmol]
    Returns
    -------
    M_dist : float
        The molar mass of distilliat, [kg/kmol]
    References
    ----------
    Дытнерский, стр. 230, формула 6.6
    """ 
    return (M_lc * xp_mol + M_hc * (1 - xp_mol))


@unitcheck(xw_mol="kmol/kmol", M_lc="kg/kmol", M_hc="kg/kmol", res_unit="kg/kmol")
def M_waste(xw_mol, M_lc, M_hc):
    """
    Calculates the molar mass of waste.
    Parameters
    ----------
    xw_mol : float
        The mol concentration low component of waste, [kmol/kmol]
    M_lc : float
        The molar mass of low-boilling component, [kg/kmol]
    M_hc : float
        The molar mass of high-boilling component, [kg/kmol]
    Returns
    -------
    M_waste : float
        The molar mass of waste, [kg/kmol]
    References
    ----------
    Дытнерский, стр. 230, формула 6.6
    """ 
    return (M_lc * xw_mol + M_hc * (1 - xw_mol))


@unitcheck(F_mass="kg/s", xp_mass="kg/kg", xf_mass="kg/kg", x_w="kg/kg", res_unit="kg/s")
def W_mass(F_mass, xp_mass, xfeed_mass, xw_mass):
    """
    Calculates the mass flow rate of waste.
    Parameters
    ----------
    F_mass : float
        The mass flow rate of feed [kg/s]
    xp_mass : float
        The mass concetration distilliat, [kg/kg]
    xw_mass : float
        The mass concetration waste, [kg/kg]
    xf_mass : float
        The mass concetration feed, [kg/kg]
    Returns
    -------
    W_mass : float
        The mass flow rate of waste. [kg/s]
    References
    ----------
    Дытнерский, стр. 228, формула 6.1
    """ 
    return F_mass * (xp_mass - xf_mass) / (xp_mass - xw_mass)


@unitcheck(F_mass="kg/s", W_mass="kg/s", res_unit="kg/s")
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


@unitcheck(F_mass="kg/s", M_feed="kg/kmol", res_unit="kmol/s")
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


@unitcheck(P_mass="kg/s", M_dist="kg/kmol", res_unit="kmol/s")
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


@unitcheck(W_mass="kg/s", M_waste="kg/kmol", res_unit="kmol/s")
def W_mol(W_mass, M_waste):
    """
    Calculates the molar flow rate of waste.
    Parameters
    ----------
    W_mass : float
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


#output for mat.balance : P_mol, F_mol, W_mol, xf_mol, xp_mol, xw_mol, W_mass, P_mass
#endregion


#region Reflux and mass ways of liq and vapor
def R_min (xp_mol, y*_mol_xf, xf_mol):
    """
    Calculates the minimal reflux number.
    Parameters
    ----------
    xp_mol : float
        The mol concentration low component of distilliat, [kmol/kmol]
    xf_mol : float
        The mol concentration low component of feed, [kmol/kmol]
    y*_mol_xf : float
        The equilibrium concentration low component of vapor at concentration of feed (liquid), [kmol/kmol]    
    Returns
    -------
    R_min : float
        The minimal reglux number, [dismensionless]
    References
    ----------
    Дытнерский, стр.228, формула 6.2
    """     
    return (xp_mol - y*_mol_xf) / (y*_mol_xf  - xf_mol)


def R (Rmin, beta):
    """
    Calculates the actual reflux number.
    Parameters
    ----------
    R_min : float
        The minimal reglux number, [dismensionless]
    beta : float
        The overflow of reflux, [dismensionless]    
    Returns
    -------
    R : float
        The actual reflux number, [dismensionless]
    References
    ----------
    Романков, стр.321, формула 7.11
    """         
    return Rmin * beta





def x_aver_top(xp_mol, xpf_mol):
    """
    Calculates the average mol concentration at the top of column.
    Parameters
    ----------
    xp_mol : float
        The mol concentration  of distilliat, [kmol/kmol]
    xpf_mol : float
        The mol concentration of point of feed, [kmol/kmol]  
    Returns
    -------
    x_aver_top : float
        The average mol concentration at top of column, [kmol/kmol]
    References
    ----------
    Дытнерский, стр.230, формула 6.6
    """            
    return (xp_mol + xpf_mol) / 2


def x_aver_bot(xw_mol, xpf_mol):
    """
    Calculates the average mol concentration  at bottom of column.
    Parameters
    ----------
    xw_mol : float
        The mol concentration  of waste, [kmol/kmol]
    xpf_mol : float
        The mol concentration of point of feed, [kmol/kmol]  
    Returns
    -------
    x_aver_top : float
        The average mol concentration at top of column, [kmol/kmol]
    References
    ----------
    Дытнерский, стр.230, формула 6.6
    """            
    return (xw_mol + xpf_mol) / 2


def M_top(M_lc, x_aver_top, M_hc):
    """
    Calculates the molar mass at top of column.
    Parameters
    ----------
    x_aver_top : float
        The average mol concentration at top of column, [kmol/kmol]
    M_lc : float
        The molar mass of low-boilling component, [kg/kmol]
    M_hc : float
        The molar mass of high-boilling component, [kg/kmol]
    Returns
    -------
    M_top : float
        The molar mass at top of column, [kg/kmol]
    References
    ----------
    Дытнерский, стр. 230, формула 6.6
    """ 
    return (M_lc * x_aver_top + M_hc * (1 - x_aver_top))


def M_bot(M_lc, x_aver_bot, M_hc):
    """
    Calculates the molar mass at bottom of column.
    Parameters
    ----------
    x_aver_bot : float
        The average mol concentration at bottom of column, [kmol/kmol]
    M_lc : float
        The molar mass of low-boilling component, [kg/kmol]
    M_hc : float
        The molar mass of high-boilling component, [kg/kmol]
    Returns
    -------
    M_bot : float
        The molar mass at bottom of column, [kg/kmol]
    References
    ----------
    Дытнерский, стр. 230, формула 6.6
    """ 
    return (M_lc * x_aver_bot + M_hc * (1 - x_aver_bot))


def L_top(P_mass, R, M_top, M_dist):
    """
    Calculates the flow rate liquid at the top of column.
    Parameters
    ----------
    P_mass : float
        The mass flow rate of distilliat, [kg/s]
    M_dist : float
        The molar mass  of distilliat, [kg/kmol]
    R : float
        The actual reflux number, [dismensionless]
    M_top : float
        The molar mass at top of column, [kg/kmol]    
    Returns
    -------
    L_top : float
        The  flow rate liquid at the top of column, [kg/s]
    References
    ----------
    Дытнерский, стр. 229, формула 6.4
    """ 
    return P_mass * R * M_top / M_dist


def L_bot(P_mass, R, M_bot, M_dist, F_mass, M_feed, phi):
    """
    Calculates the flow rate liquid at the bottom of column.
    Parameters
    ----------
    P_mass : float
        The mass flow rate of distilliat, [kg/s]
    F_mass : float
        The mass flow rate of feed, [kg/s]
    M_dist : float
        The molar mass  of distilliat, [kg/kmol]
    M_feed : float
        The molar mass  of feed, [kg/kmol]
    phi : float
        The fraction of vapor at the feed point
    R : float
        The actual reflux number, [dismensionless]
    M_bot : float
        The molar mass at bottom of column, [kg/kmol]    
    Returns
    -------
    L_bot : float
        The  flow rate liquid at the bottom of column, [kg/s]
    References
    ----------
    Дытнерский, стр. 229, формула 6.5
    """ 
    return (P_mass * R * M_bot / M_dist) + (F_mass * M_bot * (1 - phi) / M_feed)


def y_aver_top(yp_mol, ypf_mol):
    """
    Calculates the average mol concentration at the top of column.
    Parameters
    ----------
    yp_mol : float
        The mol concentration  of distilliat, [kmol/kmol]
    ypf_mol : float
        The mol concentration of point of feed, [kmol/kmol]  
    Returns
    -------
    y_aver_top : float
        The average mol concentration at top of column, [kmol/kmol]
    References
    ----------
    Дытнерский, стр.230, формула 6.8
    """            
    return (yp_mol + ypf_mol) / 2


def y_aver_bot(yw_mol, ypf_mol):
    """
    Calculates the average mol concentration  at bottom of column.
    Parameters
    ----------
    yw_mol : float
        The mol concentration  of waste, [kmol/kmol]
    ypf_mol : float
        The mol concentration of point of feed, [kmol/kmol]  
    Returns
    -------
    y_aver_top : float
        The average mol concentration at top of column, [kmol/kmol]
    References
    ----------
    Дытнерский, стр.230, формула 6.8
    """            
    return (yw_mol + ypf_mol) / 2


def M_top_vap(M_lc, y_aver_top, M_hc):
    """
    Calculates the molar mass at top of column.
    Parameters
    ----------
    y_aver_top : float
        The average mol concentration at top of column, [kmol/kmol]
    M_lc : float
        The molar mass of low-boilling component, [kg/kmol]
    M_hc : float
        The molar mass of high-boilling component, [kg/kmol]
    Returns
    -------
    M_top : float
        The molar mass at top of column, [kg/kmol]
    References
    ----------
    Дытнерский, стр. 230, формула 6.8
    """ 
    return (M_lc * y_aver_top + M_hc * (1 - y_aver_top))


def M_bot_vap(M_lc, y_aver_bot, M_hc):
    """
    Calculates the molar mass at bottom of column.
    Parameters
    ----------
    y_aver_bot : float
        The average mol concentration at bottom of column, [kmol/kmol]
    M_lc : float
        The molar mass of low-boilling component, [kg/kmol]
    M_hc : float
        The molar mass of high-boilling component, [kg/kmol]
    Returns
    -------
    M_bot : float
        The molar mass at bottom of column, [kg/kmol]
    References
    ----------
    Дытнерский, стр. 230, формула 6.8
    """ 
    return (M_lc * y_aver_bot + M_hc * (1 - y_aver_bot))


def G_top(P_mass, R, M_top_vap, M_dist):
    """
    Calculates the flow rate liquid at the top of column.
    Parameters
    ----------
    P_mass : float
        The mass flow rate of distilliat, [kg/s]
    M_dist : float
        The molar mass  of distilliat, [kg/kmol]
    R : float
        The actual reflux number, [dismensionless]
    M_top_vap : float
        The molar mass at top of column, [kg/kmol]    
    Returns
    -------
    G_top : float
        The  flow rate of the vapor at the top of column, [kg/s]
    References
    ----------
    Дытнерский, стр. 229, формула 6.7
    """ 
    return P_mass * (R + 1) * M_top_vap / M_dist


def G_bot(P_mass, R, M_bot_vap, M_dist, F_mass, M_feed, phi):
    """
    Calculates the flow rate liquid at the bottom of column.
    Parameters
    ----------
    P_mass : float
        The mass flow rate of distilliat, [kg/s]
    F_mass : float
        The mass flow rate of feed, [kg/s]
    M_dist : float
        The molar mass  of distilliat, [kg/kmol]
    M_feed : float
        The molar mass  of feed, [kg/kmol]
    phi : float
        The fraction of vapor at the feed point
    R : float
        The actual reflux number, [dismensionless]
    M_bot_vap : float
        The molar mass of the vapor at bottom of column, [kg/kmol]    
    Returns
    -------
    G_bot : float
        The  flow rate of the vapor at the bottom of column, [kg/s]
    References
    ----------
    Дытнерский, стр. 229, формула 6.7
    """ 
    return (P_mass * R * M_bot / M_dist) - (F_mass * M_bot_vap * phi/ M_feed) ## 2018 и 1880 проверка знака второго слагаемого


def x_aver_top_mass(xp_mass, xpf_mass):
    """
    Calculates the average mass concentration at the top of column.
    Parameters
    ----------
    xp_mass : float
        The mass concentration  of distilliat, [kg/kg]
    xpf_mass : float
        The mass concentration of point of feed, [kg/kg]  
    Returns
    -------
    x_aver_top_mass : float
        The average mass concentration at top of column, [kg/kg]
    References
    ----------
    Дытнерский, стр.230, формула 6.8
    """            
    return (xp_mass + xpf_mass) / 2


def x_aver_bot_mass(xw_mol, xpf_mass):
    """
    Calculates the average mass concentration  at bottom of column.
    Parameters
    ----------
    xw_mass : float
        The mass concentration  of waste, [kg/kg]
    xpf_mass : float
        The mass concentration of point of feed, [kg/kg]  
    Returns
    -------
    x_aver_bot_mass : float
        The average mass concentration at bot of column, [kg/kg]
    References
    ----------
    Дытнерский, стр.230, формула 6.8
    """            
    return (xw_mass + xpf_mass) / 2

def rho_top_liq(x_aver_top_mass, rho_lc_x_aver_top, rho_hc_x_aver_top):
    """
    Calculates the destiny of liquid at top of column.
    Parameters
    ----------
    x_aver_top_mass : float
        The average mass concentration at top of column, [kg/kg]
    rho_lc_x_aver_top : float
        The destiny of low-boilling component for mass average concentration at the top of column, [kg/m**3]
    rho_hc_x_aver_top : float
        The destiny of high-boilling component for mass average concentration at the top of column, [kg/m**3] 
    Returns
    -------
    rho_top_liq : float
        The destiny of liquid at top of column, [kg/m**3]
    References
    ----------
    Романков, стр.12, формула 1.3
    """            
    return ((x_aver_top_mass / rho_lc_x_aver_top) + ((1 - x_aver_top_mass) / rho_hc_x_aver_top))


def rho_bot_liq(x_aver_bot_mass, rho_lc_x_aver_bot, rho_hc_x_aver_bot):
    """
    Calculates the destiny of liquid at the bottom of column.
    Parameters
    ----------
    x_aver_bot_mass : float
        The average mass concentration at the bottom of column, [kg/kg]
    rho_lc_x_aver_bot : float
        The destiny of low-boilling component for mass average concentration at the bottom of column, [kg/m**3]
    rho_hc_x_aver_bot : float
        The destiny of high-boilling component for mass average concentration at the bottom of column, [kg/m**3]   
    Returns
    -------
    rho_bot_liq : float
        The destiny of liquid at the bottom of column, [kg/m**3]
    References
    ----------
    Романков, стр.12, формула 1.3
    """            
    return ((x_aver_bot_mass / rho_lc_x_aver_bot) + ((1 - x_aver_bot_mass) / rho_hc_x_aver_bot))


def rho_top_vap(t_boil_y_top, T0, Vm, M_top_vap):
    """
    Calculates the destiny of vapor at the top of column.
    Parameters
    ----------
    t_boil_y_top : float
        The boilling temperature of low-bolling component of vapor at the top of column, [C]
    T0 : float
        The initial temperature, [K]
    M_top_vap : float
        The molar mass at top of column, [kg/kmol]    
    Returns
    -------
    rho_top_vap : float
        The destiny of vapor at top of column, [kg/m**3]
    References
    ----------
    Романков, стр.13, формула 1.5
    """            
    return (M_top_vap / Vm) * (T0 / (T0 + t_boil_y_top))


def rho_bot_vap(t_boil_y_bot, T0, Vm, M_bot_vap):
    """
    Calculates the destiny of vapor at the bottom of column.
    Parameters
    ----------
    t_boil_y_bot : float
        The boilling temperature of low-bolling component of vapor at the bottom of column, [C]
    T0 : float
        The initial temperature, [K]
    M_bot_vap : float
        The molar mass at bottom of column, [kg/kmol]    
    Returns
    -------
    rho_bot_vap : float
        The destiny of vapor at bottom of column, [kg/m**3]
    References
    ----------
    Романков, стр.13, формула 1.5
    """            
    return (M_bot_vap / Vm) * (T0 / (T0 + t_boil_y_bot))  

#operating line #N(R+1) N
#output: G_top, G_bot, L_bot, L_top, Rmin
#endregion


#region The speed and diameter of column
def speed_limit_top(rho_top_liq, rho_top_vap):
    """
    Calculates the limit speed of the vapor in the column.
    Parameters
    ----------
    rho_top_vap : float
        The destiny of vapor at the top of column, [kg/m**3]
    rho_top_liq : float
        The destiny of liquid at top of column, [kg/m**3] 
    Returns
    -------
    speed_limit_top : float
        The limit speed of the vapor at the top of  column, [m/s]
    References
    ----------
    Дытнерский, стр.205, формула 5.33
    """   
    return (rho_top_liq / rho_bot_vap)^0.5


def speed_limit_bot(rho_bot_liq, rho_bot_vap):
    """
    Calculates the limit speed of the vapor in the column.
    Parameters
    ----------
    rho_bot_vap : float
        The destiny of vapor at the bottom of column, [kg/m**3]
    rho_bot_liq : float
        The destiny of liquid at the bottom of column, [kg/m**3] 
    Returns
    -------
    speed_limit_bot : float
        The limit speed of the vapor at the bottom in  column, [m/s]
    References
    ----------
    Дытнерский, стр.205, формула 5.33
    """   
    return (rho_top_liq / rho_bot_vap)^0.5


def D_top(G_top, pi, speed_limit_top, rho_top_vap):
    """
    Calculates the top diameter of column.
    Parameters
    ----------
    rho_top_vap : float
        The destiny of vapor at the top of column, [kg/m**3]
    speed_limit_top : float
        The limit speed of the vapor at the top of  column, [m/s]
    G_top : float
        The  flow rate of the vapor at the top of column, [kg/s]    
    Returns
    -------
    D_top : float
        The top diameter of column, [m]
    References
    ----------
    &&&&&
    """       
    return (4 * G_top / (pi * speed_limit_top * rho_top_vap))^0.5


def D_bot(G_bot, pi, speed_limit_bot, rho_bot_vap):
    """
    Calculates the bottom diameter of column.
    Parameters
    ----------
    rho_bot_vap : float
        The destiny of vapor at the bottom of column, [kg/m**3]
    speed_limit_bot : float
        The limit speed of the vapor at the bottom in  column, [m/s]
    G_bot : float
        The  flow rate of the vapor at the bottom of column, [kg/s]    
    Returns
    -------
    D_bot : float
        The bottom diameter of column, [m]
    References
    ----------
    &&&&&
    """       
    return (4 * G_bot / (pi * speed_limit_bot * rho_bot_vap))^0.5


def speed_section_top(speed_limit_top, D_top, D):
    """
    Calculates the section speed of vapor
    Parameters
    ----------
    speed_limit_top : float
        The limit speed of the vapor at the top of  column, [m/s]
    D_top : float
        The calculating top diameter of column, [m]
    D : float
        The choosing diameter of column, [m]   
    Returns
    -------
    speed_section_top : float
        The section speed of vapor, [m/s]
    References
    ----------
    &&&&&
    """         
    return speed_limit_top * (D/top / D)**2


def speed_section_bot(speed_limit_bot, D_bot, D):
    """
    Calculates the section speed of vapor at the bottom
    Parameters
    ----------
    speed_limit_bot: float
        The limit speed of the vapor at the bottom of  column, [m/s]
    D_bot : float
        The calculating bottom diameter of column, [m]
    D : float
        The choosing diameter of column, [m]   
    Returns
    -------
    speed_section_bot : float
        The section speed of vapor at the bottom, [m/s]
    References
    ----------
    &&&&&
    """         
    return speed_limit_bot * (D/bot / D)**2


def speed_operating_section_top(speed_section_top, D, ft):
    """
    Calculates the operating section speed of vapor at the top
    Parameters
    ----------
    speed_section_top : float
        The section speed of vapor at the top, [m]
    D : float
        The choosing diameter of column, [m]  
    ft : float
        The operating section of plate, [m**2] 
    Returns
    -------
    speed_operating_section_top : float
        The operating section speed of vapor at the top, [m/s]
    References
    ----------
    &&&&&
    """         
    return speed_section_top * (0.785 * D**2 / ft)


def speed_operating_section_bot(speed_section_bot, D, ft):
    """
    Calculates the operating section speed of vapor at the bottom
    Parameters
    ----------
    speed_section_bot : float
        The section speed of vapor at the bottom, [m]
    D : float
        The choosing diameter of column, [m]  
    ft : float
        The operating section of plate, [m**2] 
    Returns
    -------
    speed_operating_section_bot : float
        The operating section speed of vapor at the bottom, [m/s]
    References
    ----------
    &&&&&
    """         
    return speed_section_bot * (0.785 * D**2 / ft)

    #output D_top, D_bot, D, w_operating_section_top, w_operating_section_bot
    #endregion 