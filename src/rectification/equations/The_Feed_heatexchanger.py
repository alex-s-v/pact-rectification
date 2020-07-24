from rectification.utils import unitcheck
import numpy as np



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


@unitcheck(deltaT="degrees celcium", F_mass="kg/s", Cp="J/(kg * degrees celcium", r_feed="J/kg", res_unit="W")
def Q_feed(deltaT, F_mass, Cp, phi_vapor, r_feed):
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
    r_feed : float
        The heat vaporazation of mix, [J/kg]
    Returns
    -------
    Q_feed : float
        The heat load of heat exchanger, [W] , [J/s]
    References
    ----------
    Дытнерский, формула 2.2, стр.45
    """   
    return F_mass * (Cp * deltaT + phi_vapor * r_feed)


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



@unitcheck(Q_feed="W", deltaT_diff="degrees celcium",  Kt_approx="W/(m**2 * degrees celcium", res_unit="m**2")
def A_approx(Q_feed, deltaT_diff, Kt_approx):
    """
    Calculates the approximate heatransfer area.
    Parameters
    ----------
    Q_feed : float
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
    return Q_feed / (deltaT_diff * Kt_approx)


@unitcheck(Q_feed="J/s", r_steam="J/kg", res_unit="kg/s")
def m_steam_feed(Q_feed, r_steam):
    """
    Calculates the flow rate steam of boiler.
    Parameters
    ----------
    Q_feed : float
        The heat load feed of heat exchanger, [W] [J/s]
    r_steam : float
        The heat vaporazation of dist [J/kg]
    Returns
    -------
    m_steam_feed : float
        The flow rate steam of feed heat exchanger, [kg/s]
    References
    ----------
    Дытнерский, формула 2.3, стр.45
    """               
    return Q_feed / r_steam
    

@unitcheck(F_mass="kg/s", d_inner="m", mu_feed="Pa/s")
def Re_feed(F_mass, z_way, d_inner, n_pipe, mu_feed):
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
    mu_feed : float
        The mix viscocity of liquid, [Pa/s]
    Returns
    -------
    Re_feed : float
        The Reynold criterion, [dimensionless]
    References
    ----------
    &&&&&&
    """         
    return 0.785 * F_mass * z_way / (d_inner * n_pipe * mu_feed)


@unitcheck(C_feed="J/(kg * degrees celcium)", mu_feed="Pa/s", lyambda_feed="W / (m * degreec celcium)")
def Pr(C_feed, mu_feed, lyambda_feed):
    """
    Calculates the Prandtl criterion.
    Parameters
    ----------
    C_feed : float
        The heat capacity of mix [J/(kg * degrees celcium)]
    mu_feed : float
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
    return C_feed * mu_feed / lyambda_feed


def Nu(Re_feed, Pr):
    """
    Calculates the Nusselt criterion.
    Parameters
    ----------
    Re_feed : float
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
    return 0.021 * (Re_feed**0.8) * (Pr**0.4)


@unitcheck(lyambda_feed="W / (m * degreec celcium)", d_inner="m", res_unit="W / (m**2 * degrees celcium)")
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


@unitcheck(lyambda_cond="W / (m * degrees celcium)", rho_cond="kg / m**3",  mu_cond="Pa / s", m_steam_feed="kg/s", d_outside="m", res_unit="W / (m**2 * degrees celcium)")
def alpha_vap(lyambda_cond, rho_cond, mu_cond, m_steam_feed, n_pipe, d_outside):
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
    m_steam_feed : float
        The flow rate steam of feed heat exchanger, [kg/s]
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
    return lyambda_cond * 3.78 * ((rho_cond**2)* n_pipe * d_outside / (mu_cond * m_steam_feed))**(1/3)


@unitcheck(pollution_1="m**2 * degrees celcium / W", pollution_2="m**2 * degrees celcium / W", sigma="m",  lyambda_wall="W / (m * degrees celcium)", res_unit="m**2 * degrees celcium / W")
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


@unitcheck(alpha_liq="W / (m**2 * degrees celcium)", alpha_vap="W / (m**2 * degrees celcium)",  sigma_thermpollution="W/(m**2 * degrees celcium", res_unit="W/(m**2 * degrees celcium")
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


@unitcheck(Q_feed="W", deltaT_diff="degrees celcium",  Kt_approx="W/(m**2 * degrees celcium", res_unit="m**2")
def A_real(Q_feed, Kt_real, deltaT_diff):
    """
    Calculates the real heatransfer area.
    Parameters
    ----------
    Ql_feed : float
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
    return Q_feed / (Kt_real * deltaT_diff)


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
        The surface margin, [%]
    References
    ----------
    &&&&
    """          
    return (A_approx - A_real) * 100 / A_approx