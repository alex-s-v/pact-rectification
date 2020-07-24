from rectification.utils import unitcheck
import numpy as np




def deltaT_diff_dist(deltaT_larger_dist, deltaT_less_dist):
    """
    Calculates the difference of temperatures.
    Parameters
    ----------
    deltaT_larger_dist : float
        The difference temperatures between the temperature of vapor and the lowest temperature of liquid , [degrees celcium]
    deltaT_less_dist : float
        The difference temperatures between the temperature of vapor and the biggest temperature of liquid , [degrees celcium]
    Returns
    -------
    deltaT_diff_dist : float
        The coefficient difference of temperatures, [degrees celcium]
    References
    ----------
    Романков, формула 4.78, стр.169
    """   
    if deltaT_larger_dist/deltaT_less_dist < 2:
        return (deltaT_larger_dist + deltaT_less_dist) / 2
    if deltaT_larger_dist/deltaT_less_dist > 2:
        return (deltaT_larger_dist - deltaT_less_dist) / np.log(deltaT_larger_dist / deltaT_less_dist)


def deltaT_larger_dist(tp, t_coolwater_exit):
    """
    Calculates the difference temperatures between the temperature of boiling distilliat and the exit temperature of coolwater.
    Parameters
    ----------
    tp : float
        The temperature of boiling dist, [degrees celcium]
    t_coolwater_exit : float
        The exit temperature of coolwater, [degrees celcium]
    Returns
    -------
    deltaT_larger_dist : float
        The difference temperatures between the temperature of boiling dist and the exit temperature of coolwater, [degrees celcium]
    References
    ----------
    Романков, формула 4.78, стр.169
    """       
    return tp - t_coolwater_exit


def deltaT_less_dist(tp, t_coolwater_enter):
    """
    Calculates the difference temperatures between the temperature of boiling distilliat and the enter temperature of coolwater.
    Parameters
    ----------
    tw : float
        The temperature of boiling distilliat, [degrees celcium]
    t_coolwater_enter : float
        The enter temperature of coolwater, [degrees celcium]
    Returns
    -------
    deltaT_less_dist : float
        the difference temperatures between the temperature of boiling distilliat and the enter temperature of coolwater, [degrees celcium]
    References
    ----------
    Романков, формула 4.78, стр.169
    """       
    return tp - t_coolwater_enter


@unitcheck(deltaT_dist="degrees celcium", P_mass="kg/s", Cp="J/(kg * degrees celcium", r_feed="J/kg", res_unit="W")
def Q_distcooler(P_mass, Cp, t_coolwater_exit, tp):
    """
    Calculates the heat load of distilliat cooler.
    Parameters
    ----------
    W_mass : float
        The mass flow rate of distilliat, [kg/s]
    Cp : float
        The heat capacity of distilliat, [J/(kg * degrees C)]
    t_coolwater_exit: float
        The end temperature of cool distilliat, [degrees C]
    tp : float
        The temperature of boiling distilliat, [degrees celcium]
    Returns
    -------
    Q_distcooler : float
        The heat load of distilliat cooler, [W] , [J/s]
    References
    ----------
    Дытнерский, формула 2.2, стр.45
    """   
    return P_mass * Cp * (tp - t_coolwater_exit)


@unitcheck(Q_distcooler="W", deltaT_diff_dist="degrees celcium",  Kt_approx_dist="W/(m**2 * degrees celcium", res_unit="m**2")
def A_approx_dist(Q_distcooler, deltaT_diff_dist, Kt_approx_dist):
    """
    Calculates the approximate heatransfer area.
    Parameters
    ----------
    Q_distcooler : float
        The heat load of distilliat cooler, [W] , [J/s]
    deltaT_diff_dist : float
        The coefficient difference of temperatures, [degrees celcium]
    Kt_approx_dist : float
        The heatransfer coefficient [W/(m**2 * degrees celcium)]
    Returns
    -------
    A_approx_dist : float
        The approximate heatransfer area, [m**2]
    References
    ----------
    Романков, формула 4.72, стр.168
    """           
    return Q_distcooler / (deltaT_diff_dist * Kt_approx_dist)


@unitcheck(Q_distcooler="J/s", r_steam="J/kg", res_unit="kg/s")
def m_coolwater_dist(Q_distcooler, Cp_distcooler, t_coolwater_exit, t_coolwater_enter):
    """
    Calculates the flow rate steam of distilliat cooler.
    Parameters
    ----------
    Q_distcooler : float
        The heat load feed of distilliat cooler, [W] [J/s]
    Cp_distcooler : float
        The heat capacity of water, [J/(kg * degrees C)]
    t_coolwater_exit: float
        The end temperature of cool distilliat, [degrees C]
    t_coolwater_enter: float
        The start temperature of cool distilliat, [degrees C]    
    Returns
    -------
    m_coolwater_dist : float
        The flow rate steam of dist cooler, [kg/s]
    References
    ----------
    Дытнерский, формула 2.3, стр.45
    """               
    return Q_distcooler / (Cp_distcooler * (t_coolwater_exit - t_coolwater_enter))
    

@unitcheck(W_mass="kg/s", d_inner_dist="m", mu_dist_avrg="Pa/s")
def Re_intube_dist(W_mass, z_way_dist, d_inner_dist, n_pipe_dist, mu_dist_avrg):
    """
    Calculates the Reynold criterion.
    Parameters
    ----------
    F_mass : float
        The mass flow rate of feed [kg/s]
    z_way_dist : float
        The number of ways in heat exchanger [dimensionless]
    d_inner_dist : float
        The diametr of inner pipe, [m]
    n_pipe_dist : float
        The number of pipes in heat exchanger, [dimensionless]
    mu_dist_avrg : float
        The mix viscocity of liquid, [Pa/s]
    Returns
    -------
    Re_intube_dist : float
        The Reynold criterion in tube(distilliat), [dimensionless]
    References
    ----------
    &&&&&&
    """         
    return 0.785 * W_mass * z_way_dist / (d_inner_dist * n_pipe_dist * mu_dist_avrg)


@unitcheck(m_coolwater_dist="kg/s", d_outer_dist="m", mu_feed="Pa/s")
def Re_outtube_dist(m_coolwater_dist, d_outer_dist, mu_coolwater_dist):
    """
    Calculates the Reynold criterion.
    Parameters
    ----------
    m_coolwater_dist : float
        The mass flow rate of feed [kg/s]
    d_outer_dist : float
        The diametr of outer pipe, [m]
    mu_coolwater_dist : float
        The mix viscocity of coolwater, [Pa/s]
    Returns
    -------
    Re_outtube_dist : float
        The Reynold criterion in outer tube(distilliat), [dimensionless]
    References
    ----------
    &&&&&&
    """         
    return 0.785 * m_coolwater_dist * d_outer_dist / ( * mu_coolwater_dist)


@unitcheck(Cp_distcooler="J/(kg * degrees celcium)", mu_coolwater_dist="Pa/s", lyambda_coolwater_dist="W / (m * degreec celcium)")
def Pr_coolwater_dist(Cp_distcooler, mu_coolwater_dist, lyambda_coolwater_dist):
    """
    Calculates the Prandtl criterion of coolwater.
    Parameters
    ----------
    Cp_distcooler : float
        The heat capacity of water, [J/(kg * degrees C)]
    mu_coolwater_dist : float
        The mix viscocity of coolwater, [Pa/s]
    lyambda_coolwater_dist : float
        The thermal conductivity of coolwater, [W / (m * degreec celcium)]
    Returns
    -------
    Pr_coolwater_dist : float
        The Prandtl criterion of coolwater, [dimensionless]
    References
    ----------
    Романков, формула 4.12, стр.151
    """        
    return Cp_distcooler * mu_coolwater_dist / lyambda_coolwater_dist


def Nu_dist(Re_outtube_dist, Pr_coolwater_dist):
    """
    Calculates the Nusselt criterion of coolwater.
    Parameters
    ----------
    Re_outtube_dist : float
        The Reynold criterion, [dimensionless]
    Pr_coolwater_dist : float
        The Prandtl criterion, [dimensionless]
    Returns
    -------
    Nu_dist : float
        The Nusselt criterion, [dimensionless]
    References
    ----------
    Романков, формула 4.17, стр.152
    """      
    return 0.24 * (Re_outtube_dist**0.6) * (Pr_coolwater_dist**0.36)*0,93


@unitcheck(lyambda_coolwater_dist="W / (m * degreec celcium)", d_outer_dist="m", res_unit="W / (m**2 * degrees celcium)")
def alpha_coolwater_dist(Nu_dist, lyambda_coolwater_dist, d_outer_dist):
    """
    Calculates the coefficent of heat transfer(alpha) from liquid to wall of pipe.
    Parameters
    ----------
    Nu_dist : float
        The Nusselt criterion, [dimensionless]
    lyambda_coolwater_dist : float
        The thermal conductivity of coolwater, [W / (m * degreec celcium)]
    d_outer : float
        The diametr of outer pipe, [m]
    Returns
    -------
    alpha_coolwater_dist : float
        The coefficent of heat transfer(alpha), [W / (m**2 * degrees celcium)]
    References
    ----------
    Романков, формула 4.11, стр.150
    """          
    return Nu_dist * lyambda_coolwater_dist / d_outer_dist


@unitcheck(lyambda_cond_dist_avrg="W / (m * degrees celcium)", rho_cond_dist_avrg="kg / m**3",  mu_cond_dist_avrg="Pa / s", W_mass="kg/s", L_dist="m", res_unit="W / (m**2 * degrees celcium)")
def alpha_dist(lyambda_cond_dist_avrg, rho_cond_dist_avrg, mu_cond_dist_avrg, W_mass, n_pipe_dist, L_dist):
    """
    Calculates the coefficent of heat transfer(alpha) from steam to wall of pipe.
    Parameters
    ----------
    lyambda_cond_dist_avrg : float
        The thermal conducivity of condensate, [W / (m * degrees celcium)]
    rho_cond_dist_avrg : float
        The destiny of condensate, [kg / m**3]
    mu_cond_dist_avrg : float
        The viscosity of condensate, [Pa / s]
    W_mass : float
        The flow rate steam of feed heat exchanger, [kg/s]
    n_pipe_dist : float
        The number of pipes in heat exchanger, [dimensionless]
    L_dist : float
        The length of tubes, [m]    
    Returns
    -------
    alpha_dist : float
        The coefficent of heat transfer(alpha) from steam to wall of pipe, [W / (m**2 * degrees celcium)]
    References
    ----------
    Дытнерский, формула 2.24, стр.53
    """              
    return lyambda_cond_dist_avrg * 2.02 * ((rho_cond_dist_avrg**2)* L_dist * n_pipe_dist / (mu_cond_dist_avrg * W_mass))**(1/3)


@unitcheck(pollution_1_dist="m**2 * degrees celcium / W", pollution_2_dist="m**2 * degrees celcium / W", sigma="m",  lyambda_wall="W / (m * degrees celcium)", res_unit="m**2 * degrees celcium / W")
def sigma_thermpollution_dist(pollution_1_dist, pollution_2_dist, sigma, lyambda_wall):
    """
    Calculates the sum of thermal pollutions.
    Parameters
    ----------
    pollution_1_dist : float
        The thermal pollution of the first coolant, [m**2 * degrees celcium / W]
    pollution_2_dist : float
        The thermal pollution of the second coolant, [m**2 * degrees celcium / W]
    sigma : float
        The thickness of pipe wall, [m]
    lyambda_wall : float
        The thermal conducivity of wall, [W / (m * degrees celcium)]  
    Returns
    -------
    sigma_thermpollution_dist : float
        The the sum of thermal pollutions, [m**2 * degrees celcium / W]
    References
    ----------
    &&&&&
    """                  
    return (sigma / lyambda_wall) + (1 / pollution_1_dist) + (1 / pollution_2_dist)


@unitcheck(alpha_dist="W / (m**2 * degrees celcium)", alpha_coolwater_dist="W / (m**2 * degrees celcium)",  sigma_thermpollution_dist="W/(m**2 * degrees celcium", res_unit="W/(m**2 * degrees celcium")
def Kt_real_dist(alpha_dist, alpha_coolwater_dist, sigma_thermpollution_dist):
    """
    Calculates the coefficient of heat transfer (Kt).
    Parameters
    ----------
    alpha_dist : float
        The coefficent of heat transfer(alpha), [W / (m**2 * degrees celcium)]
    alpha_coolwater_dist : float
        The coefficent of heat transfer(alpha) from steam to wall of pipe, [W / (m**2 * degrees celcium)]
    sigma_thermpollution_dist : float
        The the sum of thermal pollutions, [m**2 * degrees celcium / W]  
    Returns
    -------
    Kt_real_dist : float
        The coefficient of heat transfer (Kt), [W / (m**2 * degrees celcium)]
    References
    ----------
    Романков, формула 4.74, стр. 168
    """      
    return ((1 / alpha_dist) + (1 / alpha_coolwater_dist) + (sigma_thermpollution_dist))**-1


@unitcheck(Q_distcooler="W", deltaT_diff_dist="degrees celcium",  Kt_real_dist="W/(m**2 * degrees celcium", res_unit="m**2")
def A_real_dist(Q_distcooler, Kt_real_dist, deltaT_diff_dist):
    """
    Calculates the real heatransfer area.
    Parameters
    ----------
    Q_distcooler : float
        The heat load of distilliat cooler, [W] , [J/s]
    deltaT_diff_dist : float
        The coefficient difference of temperatures, [degrees celcium]
    Kt_real_dist : float
        The heat transfer coefficient [W/(m**2 * degrees celcium)]
    Returns
    -------
    A_real_dist : float
        The real heat ransfer area, [m**2]
    References
    ----------
    Романков, формула 4.72, стр.168
    """      
    return Q_distcooler / (Kt_real_dist * deltaT_diff_dist)


def surface_margin_dist (A_approx_dist, A_real_dist):
    """
    Calculates the surface margin.
    Parameters
    ----------
    A_approximate_dist : float
        The approximate heat ransfer area, [m**2]
    A_real_dist : float
        The real heat transfer area, [m**2]
    Returns
    -------
    surface_margin : float
        The surface margin, [%]
    References
    ----------
    ???
    """          
    return (A_approx_dist - A_real_dist) * 100 / A_approx_dist