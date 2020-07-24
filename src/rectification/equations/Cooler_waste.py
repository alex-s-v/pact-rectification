from rectification.utils import unitcheck
import numpy as np




def deltaT_diff_waste(deltaT_larger_waste, deltaT_less_waste):
    """
    Calculates the difference of temperatures.
    Parameters
    ----------
    deltaT_larger_waste : float
        The difference temperatures between the temperature of vapor and the lowest temperature of liquid , [degrees celcium]
    deltaT_less_waste : float
        The difference temperatures between the temperature of vapor and the biggest temperature of liquid , [degrees celcium]
    Returns
    -------
    deltaT_diff_waste : float
        The coefficient difference of temperatures, [degrees celcium]
    References
    ----------
    Романков, формула 4.78, стр.169
    """   
    if deltaT_larger_waste/deltaT_less_waste < 2:
        return (deltaT_larger_waste + deltaT_less_waste) / 2
    if deltaT_larger_waste/deltaT_less_waste > 2:
        return (deltaT_larger_waste - deltaT_less_waste) / np.log(deltaT_larger_waste / deltaT_less_waste)


def deltaT_larger_waste(tw, t_coolwater_exit):
    """
    Calculates the difference temperatures between the temperature of boiling waste and the exit temperature of coolwater.
    Parameters
    ----------
    tw : float
        The temperature of boiling waste, [degrees celcium]
    t_coolwater_exit : float
        The exit temperature of coolwater, [degrees celcium]
    Returns
    -------
    deltaT_larger_waste : float
        The difference temperatures between the temperature of boiling waste and the exit temperature of coolwater, [degrees celcium]
    References
    ----------
    Романков, формула 4.78, стр.169
    """       
    return tw - t_coolwater_exit


def deltaT_less_waste(tw, t_coolwater_enter):
    """
    Calculates the difference temperatures between the temperature of boiling waste and the enter temperature of coolwater.
    Parameters
    ----------
    tw : float
        The temperature of boiling waste, [degrees celcium]
    t_coolwater_enter : float
        The enter temperature of coolwater, [degrees celcium]
    Returns
    -------
    deltaT_less_waste : float
        the difference temperatures between the temperature of boiling waste and the enter temperature of coolwater, [degrees celcium]
    References
    ----------
    Романков, формула 4.78, стр.169
    """       
    return tw - t_coolwater_enter


@unitcheck(deltaT_waste="degrees celcium", W_mass="kg/s", Cp="J/(kg * degrees celcium", r_feed="J/kg", res_unit="W")
def Q_wastecooler(W_mass, Cw, t_coolwater_exit, tw):
    """
    Calculates the heat load of waste cooler.
    Parameters
    ----------
    W_mass : float
        The mass flow rate of waste, [kg/s]
    Cw : float
        The heat capacity of waste, [J/(kg * degrees C)]
    t_coolwater_exit: float
        The end temperature of cool waste, [degrees C]
    tw : float
        The temperature of boiling waste, [degrees celcium]
    Returns
    -------
    Q_wastecooler : float
        The heat load of waste cooler, [W] , [J/s]
    References
    ----------
    Дытнерский, формула 2.2, стр.45
    """   
    return W_mass * Cw * (tw - t_coolwater_exit)


@unitcheck(Q_wastecooler="W", deltaT_diff_waste="degrees celcium",  Kt_approx_waste="W/(m**2 * degrees celcium", res_unit="m**2")
def A_approx_waste(Q_wastecooler, deltaT_diff_waste, Kt_approx_waste):
    """
    Calculates the approximate heatransfer area.
    Parameters
    ----------
    Q_wastecooler : float
        The heat load of waste cooler, [W] , [J/s]
    deltaT_diff_waste : float
        The coefficient difference of temperatures, [degrees celcium]
    Kt_approx_waste : float
        The heatransfer coefficient [W/(m**2 * degrees celcium)]
    Returns
    -------
    A_approx_waste : float
        The approximate heatransfer area, [m**2]
    References
    ----------
    Романков, формула 4.72, стр.168
    """           
    return Q_wastecooler / (deltaT_diff_waste * Kt_approx_waste)


@unitcheck(Q_wastecooler="J/s", r_steam="J/kg", res_unit="kg/s")
def m_coolwater_waste(Q_wastecooler, Cp_wastecooler, t_coolwater_exit, t_coolwater_enter):
    """
    Calculates the flow rate steam of waste cooler.
    Parameters
    ----------
    Q_wastecooler : float
        The heat load feed of waste cooler, [W] [J/s]
    Cp_wastecooler : float
        The heat capacity of water, [J/(kg * degrees C)]
    t_coolwater_exit: float
        The end temperature of cool waste, [degrees C]
    t_coolwater_enter: float
        The start temperature of cool waste, [degrees C]    
    Returns
    -------
    m_coolwater_waste : float
        The flow rate steam of waste cooler, [kg/s]
    References
    ----------
    Дытнерский, формула 2.3, стр.45
    """               
    return Q_wastecooler / (Cp_wastecooler * (t_coolwater_exit - t_coolwater_enter))
    

@unitcheck(W_mass="kg/s", d_inner_waste="m", mu_waste_avrg="Pa/s")
def Re_intube_waste(W_mass, z_way_waste, d_inner_waste, n_pipe_waste, mu_waste_avrg):
    """
    Calculates the Reynold criterion.
    Parameters
    ----------
    F_mass : float
        The mass flow rate of feed [kg/s]
    z_way_waste : float
        The number of ways in heat exchanger [dimensionless]
    d_inner_waste : float
        The diametr of inner pipe, [m]
    n_pipe_waste : float
        The number of pipes in heat exchanger, [dimensionless]
    mu_waste_avrg : float
        The mix viscocity of liquid, [Pa/s]
    Returns
    -------
    Re_intube_waste : float
        The Reynold criterion in tube(waste), [dimensionless]
    References
    ----------
    &&&&&&
    """         
    return 0.785 * W_mass * z_way_waste / (d_inner_waste * n_pipe_waste * mu_waste_avrg)


@unitcheck(m_coolwater_waste="kg/s", d_outer_waste="m", mu_feed="Pa/s")
def Re_outtube_waste(m_coolwater_waste, d_outer_waste, mu_coolwater_waste):
    """
    Calculates the Reynold criterion.
    Parameters
    ----------
    m_coolwater_waste : float
        The mass flow rate of feed [kg/s]
    d_outer_waste : float
        The diametr of outer pipe, [m]
    mu_coolwater_waste : float
        The mix viscocity of coolwater, [Pa/s]
    Returns
    -------
    Re_outtube_waste : float
        The Reynold criterion in outer tube(waste), [dimensionless]
    References
    ----------
    &&&&&&
    """         
    return 0.785 * m_coolwater_waste * d_outer_waste / ( * mu_coolwater_waste)


@unitcheck(Cp_wastecooler="J/(kg * degrees celcium)", mu_coolwater_waste="Pa/s", lyambda_coolwater_waste="W / (m * degreec celcium)")
def Pr_coolwater_waste(Cp_wastecooler, mu_coolwater_waste, lyambda_coolwater_waste):
    """
    Calculates the Prandtl criterion of coolwater.
    Parameters
    ----------
    Cp_wastecooler : float
        The heat capacity of water, [J/(kg * degrees C)]
    mu_coolwater_waste : float
        The mix viscocity of coolwater, [Pa/s]
    lyambda_coolwater_waste : float
        The thermal conductivity of coolwater, [W / (m * degreec celcium)]
    Returns
    -------
    Pr_coolwater_waste : float
        The Prandtl criterion of coolwater, [dimensionless]
    References
    ----------
    Романков, формула 4.12, стр.151
    """        
    return Cp_wastecooler * mu_coolwater_waste / lyambda_coolwater_waste


def Nu_waste(Re_outtube_waste, Pr_coolwater_waste):
    """
    Calculates the Nusselt criterion of coolwater.
    Parameters
    ----------
    Re_outtube_waste : float
        The Reynold criterion, [dimensionless]
    Pr_coolwater_waste : float
        The Prandtl criterion, [dimensionless]
    Returns
    -------
    Nu_waste : float
        The Nusselt criterion, [dimensionless]
    References
    ----------
    Романков, формула 4.17, стр.152
    """      
    return 0.24 * (Re_outtube_waste**0.6) * (Pr_coolwater_waste**0.36)*0,93


@unitcheck(lyambda_coolwater_waste="W / (m * degreec celcium)", d_outer_waste="m", res_unit="W / (m**2 * degrees celcium)")
def alpha_coolwater_waste(Nu_waste, lyambda_coolwater_waste, d_outer_waste):
    """
    Calculates the coefficent of heat transfer(alpha) from liquid to wall of pipe.
    Parameters
    ----------
    Nu_waste : float
        The Nusselt criterion, [dimensionless]
    lyambda_coolwater_waste : float
        The thermal conductivity of coolwater, [W / (m * degreec celcium)]
    d_outer : float
        The diametr of outer pipe, [m]
    Returns
    -------
    alpha_coolwater_waste : float
        The coefficent of heat transfer(alpha), [W / (m**2 * degrees celcium)]
    References
    ----------
    Романков, формула 4.11, стр.150
    """          
    return Nu_waste * lyambda_coolwater_waste / d_outer_waste


@unitcheck(lyambda_cond_waste_avrg="W / (m * degrees celcium)", rho_cond_waste_avrg="kg / m**3",  mu_cond_waste_avrg="Pa / s", W_mass="kg/s", L_waste="m", res_unit="W / (m**2 * degrees celcium)")
def alpha_waste(lyambda_cond_waste_avrg, rho_cond_waste_avrg, mu_cond_waste_avrg, W_mass, n_pipe_waste, L_waste):
    """
    Calculates the coefficent of heat transfer(alpha) from steam to wall of pipe.
    Parameters
    ----------
    lyambda_cond_waste_avrg : float
        The thermal conducivity of condensate, [W / (m * degrees celcium)]
    rho_cond_waste_avrg : float
        The destiny of condensate, [kg / m**3]
    mu_cond_waste_avrg : float
        The viscosity of condensate, [Pa / s]
    W_mass : float
        The flow rate steam of feed heat exchanger, [kg/s]
    n_pipe_waste : float
        The number of pipes in heat exchanger, [dimensionless]
    L_waste : float
        The length of tubes, [m]    
    Returns
    -------
    alpha_waste : float
        The coefficent of heat transfer(alpha) from steam to wall of pipe, [W / (m**2 * degrees celcium)]
    References
    ----------
    Дытнерский, формула 2.24, стр.53
    """              
    return lyambda_cond_waste_avrg * 2.02 * ((rho_cond_waste_avrg**2)* L_waste * n_pipe_waste / (mu_cond_waste_avrg * W_mass))**(1/3)


@unitcheck(pollution_1_waste="m**2 * degrees celcium / W", pollution_2_waste="m**2 * degrees celcium / W", sigma="m",  lyambda_wall="W / (m * degrees celcium)", res_unit="m**2 * degrees celcium / W")
def sigma_thermpollution_waste(pollution_1_waste, pollution_2_waste, sigma, lyambda_wall):
    """
    Calculates the sum of thermal pollutions.
    Parameters
    ----------
    pollution_1_waste : float
        The thermal pollution of the first coolant, [m**2 * degrees celcium / W]
    pollution_2_waste : float
        The thermal pollution of the second coolant, [m**2 * degrees celcium / W]
    sigma : float
        The thickness of pipe wall, [m]
    lyambda_wall : float
        The thermal conducivity of wall, [W / (m * degrees celcium)]  
    Returns
    -------
    sigma_thermpollution_waste : float
        The the sum of thermal pollutions, [m**2 * degrees celcium / W]
    References
    ----------
    &&&&&
    """                  
    return (sigma / lyambda_wall) + (1 / pollution_1_waste) + (1 / pollution_2_waste)


@unitcheck(alpha_waste="W / (m**2 * degrees celcium)", alpha_coolwater_waste="W / (m**2 * degrees celcium)",  sigma_thermpollution_waste="W/(m**2 * degrees celcium", res_unit="W/(m**2 * degrees celcium")
def Kt_real_waste(alpha_waste, alpha_coolwater_waste, sigma_thermpollution_waste):
    """
    Calculates the coefficient of heat transfer (Kt).
    Parameters
    ----------
    alpha_waste : float
        The coefficent of heat transfer(alpha), [W / (m**2 * degrees celcium)]
    alpha_coolwater_waste : float
        The coefficent of heat transfer(alpha) from steam to wall of pipe, [W / (m**2 * degrees celcium)]
    sigma_thermpollution_waste : float
        The the sum of thermal pollutions, [m**2 * degrees celcium / W]  
    Returns
    -------
    Kt_real_waste : float
        The coefficient of heat transfer (Kt), [W / (m**2 * degrees celcium)]
    References
    ----------
    Романков, формула 4.74, стр. 168
    """      
    return ((1 / alpha_waste) + (1 / alpha_coolwater_waste) + (sigma_thermpollution_waste))**-1


@unitcheck(Q_wastecooler="W", deltaT_diff_waste="degrees celcium",  Kt_real_waste="W/(m**2 * degrees celcium", res_unit="m**2")
def A_real_waste(Q_wastecooler, Kt_real_waste, deltaT_diff_waste):
    """
    Calculates the real heatransfer area.
    Parameters
    ----------
    Q_wastecooler : float
        The heat load of wastecooler, [W] , [J/s]
    deltaT_diff_waste : float
        The coefficient difference of temperatures, [degrees celcium]
    Kt_real_waste : float
        The heat transfer coefficient [W/(m**2 * degrees celcium)]
    Returns
    -------
    A_real_waste : float
        The real heat ransfer area, [m**2]
    References
    ----------
    Романков, формула 4.72, стр.168
    """      
    return Q_wastecooler / (Kt_real_waste * deltaT_diff_waste)


def surface_margin_waste (A_approx_waste, A_real_waste):
    """
    Calculates the surface margin.
    Parameters
    ----------
    A_approximate_waste : float
        The approximate heat ransfer area, [m**2]
    A_real_waste : float
        The real heat transfer area, [m**2]
    Returns
    -------
    surface_margin : float
        The surface margin, [%]
    References
    ----------
    &&&&
    """          
    return (A_approx_waste - A_real_waste) * 100 / A_approx_waste