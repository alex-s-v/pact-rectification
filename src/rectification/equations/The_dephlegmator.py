from rectification.utils import unitcheck
import numpy as np


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


def deltaT_larger_deph(t_cond, t_dist):
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
    return t_cond - t_dist


def deltaT_less_deph(t_cond, tulti_water):
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


@unitcheck(P_mass="kg/s", r_dist="J/kg", res_unit="W")
def Q_deph(P_mass, r_dist, R):
    """
    Calculates the heat load of dephlegmator.
    Parameters
    ----------
    P_mass : float
        The mass flow rate of dist , [kg/s]
    R : float
        The reflux number [dimensionless]
    r_dist : float
        The heat vaporazation of dist, [J/kg]
    Returns
    -------
    Q_deph : float
        The heat load of dephlegmator, [W] , [J/s]
    References
    ----------
    Дытнерский, формула 2.2, стр.45
    """   
    return P_mass * (R + 1) * r_dist


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


@unitcheck(Q_deph="W", deltaT_diff="degrees celcium",  Kt_approx="W/(m**2 * degrees celcium", res_unit="m**2")
def A_approx_deph(Q_deph, deltaT_diff_deph, Kt_approx):
    """
    Calculates the approximate heatransfer area.
    Parameters
    ----------
    Q_deph : float
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
    return Q_deph / (deltaT_diff_deph * Kt_approx)


@unitcheck(Q_deph="J/s", Cp="J/(kg * degrees C)",  deltaT_deph="degrees C", res_unit="J/s")
def m_coolwater_deph(Q_deph, Cp, deltaT_deph):
    """
    Calculates the flow rate steam of boiler.
    Parameters
    ----------
    Q_deph : float
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
    return Q_deph / Cp * deltaT_deph
    

@unitcheck(m_coolwater_deph="kg/s", d_inner_deph="m",  mu_coolwater="Pa / s")
def Re_deph(m_coolwater_deph, z_way, d_inner_deph, n_pipe, mu_coolwater):
    """
    Calculates the Reynold criterion.
    Parameters
    ----------
    m_coolwater_deph : float
        The flow rate cool water of dephlegmator, [W] [J/s]
    z_way : float
        The number of ways in dephlegmator [dimensionless]
    d_inner_deph : float
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
    return 0.785 * m_coolwater_deph * z_way / (d_inner_deph * n_pipe * mu_coolwater)


@unitcheck(C_coolwater="J/(kg * degrees celcium)", mu_water="Pa/s",  lyambda_coolwater="W / (m * degreec celcium)")
def Pr_deph(C_coolwater, mu_coolwater, lyambda_coolwater):
    """
    Calculates the Prandtl criterion.
    Parameters
    ----------
    C_capacity_water : float
        The heat capacity of cool water [J/(kg * degrees celcium)]
    mu_water : float
        The viscocity of cool water, [Pa/s]
    lyambda_coolwater : float
        The thermal conductivity of cool water, [W / (m * degreec celcium)]
    Returns
    -------
    Pr_deph : float
        The Prandtl criterion, [dimensionless]
    References
    ----------
    Романков, формула 4.12, стр.151
    """        
    return C_coolwater * mu_coolwater / lyambda_coolwater


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


@unitcheck(lyambda_coolwater="W / (m * degreec celcium)", d_inner_deph="m", res_unit="W / (m**2 * degrees celcium)")
def alpha_liq_deph(Nu_deph, lyambda_coolwater, d_inner_deph):
    """
    Calculates the coefficent of heat transfer(alpha) from liquid to wall of pipe.
    Parameters
    ----------
    Nu_deph : float
        The Nusselt criterion, [dimensionless]
    lyambda_coolwater : float
        The thermal conductivity of cool water, [W / (m * degreec celcium)]
    d_inner_deph : float
        The diametr of inner pipe, [m]
    Returns
    -------
    alpha_liq_deph : float
        The coefficent of heat transfer(alpha), [W / (m**2 * degrees celcium)]
    References
    ----------
    Романков, формула 4.11, стр.150
    """          
    return Nu_deph * lyambda_coolwater / d_inner_deph


@unitcheck(lyambda_cond_dist="W / (m * degrees celcium)", rho_cond_dist="kg / m**3",  mu_cond_dist="Pa / s", P_mass="kg/s", res_unit="W / (m**2 * degrees celcium)")
def alpha_cond_deph(lyambda_cond_dist, rho_cond_dist, mu_cond_dist, P_mass, n_pipe_deph, L_pipe_deph):
    """
    Calculates the coefficent of heat transfer(alpha) from steam to wall of pipe.
    Parameters
    ----------
    lyambda_cond_dist : float
        The thermal conducivity condensate of distilliat , [W / (m * degrees celcium)]
    rho_cond_dist : float
        The destiny condensate of distilliat, [kg / m**3]
    mu_cond_dist : float
        The viscosity condensate of distilliat, [Pa / s]
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
    if n_pipe_deph < 100:      
        return lyambda_cond_dist * 3.78 * ((rho_cond_dist**2)* n_pipe_deph * L_pipe_deph / (mu_cond_dist * P_mass))**(1/3)
    if n_pipe_deph > 100:
        return 0.6 * lyambda_cond_dist * 3.78 * ((rho_cond_dist**2)* n_pipe_deph * L_pipe_deph / (mu_cond_dist * P_mass))**(1/3)


@unitcheck(pollution_1_deph="m**2 * degrees celcium / W", pollution_2_deph="m**2 * degrees celcium / W", sigma_deph="m",  lyambda_wall_deph="W / (m * degrees celcium)", res_unit="m**2 * degrees celcium / W")
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


@unitcheck(alpha_liq_deph="W / (m**2 * degrees celcium)", alpha_cond_deph="W / (m**2 * degrees celcium)",  sigma_thermpollution_deph="W/(m**2 * degrees celcium", res_unit="W/(m**2 * degrees celcium")
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


@unitcheck(Q_deph="W", deltaT_diff_deph="degrees celcium",  Kt_real_deph="W/(m**2 * degrees celcium", res_unit="m**2")
def A_real_deph(Q_deph, Kt_real_deph, deltaT_diff_deph):
    """
    Calculates the real heatransfer area.
    Parameters
    ----------
    Q_deph : float
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
    return Q_deph / (Kt_real_deph * deltaT_diff_deph)


def surface_margin_deph(A_approx_deph, A_real_deph):
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