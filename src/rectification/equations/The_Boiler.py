from rectification.utils import unitcheck
import numpy as np




def Q_boiler(W_mass, Cw, tw, P_mass, Q_deph,  F_mass, Cf, tf, Cp, tp, Q_loss):
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
    Q_boiler : float
        The heat load of boiler, [W] [J/s]
    References
    ----------
    14.	Дытнерский Ю.И. Процессы и аппараты химической технологии. Учебник для вузов. Часть 2. стр 123-124
    """     
    return W_mass * Cw * tw + Q_deph - F_mass * Cf * tf + P_mass * Cp * tp + Q_loss


def m_steam_boil(Q_boiler, r_steam):
    """
    Calculates the flow rate steam of boiler.
    Parameters
    ----------
    Q_boiler : float
        The heat load of boiler, [W] [J/s]
    r_steam : float
        The heat vaporazation of dist [J/kg]
    Returns
    -------
    m_steam_boil : float
        The flow rate steam of boiler, [kg/s]
    References
    ----------
    Дытнерский, формула 2.3, стр.45
    """               
    return Q_boiler / r_steam


def deltaT_boil(tw, t_vapor):
    """
    Calculates the temperature difference of boiler.
    Parameters
    ----------
    tw : float
        The boiling temperature of liquid, [С]
    t_vapor : float
        The condensation temperature of steam, [C]
    Returns
    -------
    deltaT_boil : float
        The temperature difference of boiler, [C]
    References
    ----------
    &&&&&
    """     
    return t_vapor - tw


@unitcheck(Q_boil="W", deltaT_boil="degrees celcium",  Kt_approx_boiler="W/(m**2 * degrees celcium", res_unit="m**2")
def A_approx_boiler(Q_boiler, deltaT_boil, Kt_approx_boiler):
    """
    Calculates the approximate heatransfer area.
    Parameters
    ----------
    Qload_boiler : float
        The heat load of boiler, [W] , [J/s]
    deltaT_boil : float
        The temperature difference of boiler, [C]
    Kt_approx_boiler : float
        The heatransfer coefficient [W/(m**2 * degrees celcium)]
    Returns
    -------
    A_approx_boiler : float
        The approximate heatransfer area of boiler, [m**2]
    References
    ----------
    Романков, формула 4.72, стр.168
    """           
    return Q_boiler / (deltaT_boil * Kt_approx_boiler)


def A_boiler(rho_cond, r_steam, g, mu_cond, L_tube_boiler, lyambda_cond):
    """
    Calculates the coefficent A of boiler.
    Parameters
    ----------
    lyambda_cond : float
        The thermal conducivity of condensate, [W / (m * degrees celcium)]
    rho_cond : float
        The destiny of condensate, [kg / m**3]
    mu_cond : float
        The viscosity of condensate, [Pa / s]
    r_steam : float
        The vaporazation of vapor, [kg/s]
    L_tube_boiler : float
        The tube length of boiler, [m]    
    Returns
    -------
    A_boiler : float
        The coefficent A of boiler, [dismensionless]
    References
    ----------
    Дытнерский, формула 2.23, стр.53
    """                  
    return 1.21 * lyambda_cond * ((rho_cond**2) * r_steam * g / (mu_cond * L_tube_boiler))**(1/3)


def B_boiler(rho_W_boil, rho_W_vapor, sigma_W_boil, r_W_boil, Cw_boil, mu_W_boil, lyambda_W_boil):
    """
    Calculates the coefficent B of boiler.
    Parameters
    ----------
    lyambda_W_boil : float
        The thermal conducivity of waste at the boilling temperature, [W / (m * degrees celcium)]
    rho_W_boil : float
        The destiny of waste at boilling the temperature, [kg / m**3]
    rho_W_vapor : float
        The destiny of vapor at boilling the temperature, [kg / m**3]
    sigma_W_boil : float
        The surface tension of waste at the boilling temperature, [N/m]
    r_W_boil : float
        The heat vaporazation of waste at the boilling temperature, [N/m]    
    Cw_boiler : float
        The heat capacity of waste at the boiling temperature, [J/(kg * degrees C)]
    mu_W_boiler : float
        The viscocity of waste at the boiling temperature, [Pa/s]      
    Returns
    -------
    B_boiler : float
        The coefficent B of boiler, [dismensionless]
    References
    ----------
    Дытнерский, формула 2.28, стр.54
    """                  
    return 780 * lyambda_W_boil^(1.3) * rho_W_boil^(0.5) * rho_W_vapor^(0.06) / (sigma_W_boil^(0.5) * r_W_boil^(0.6) * rho_W_vapor^(0.66) * Cw_boil^(0.3) * mu_W_boil^(0.3))

@unitcheck(pollution_1_boiler="m**2 * degrees celcium / W", pollution_2_boiler="m**2 * degrees celcium / W", sigma_boiler="m",  lyambda_wall_boiler="W / (m * degrees celcium)", res_unit="m**2 * degrees celcium / W")
def sigma_thermpollution_boiler(pollution_1_boiler, pollution_2_boiler, sigma_boiler, lyambda_wall_boiler):
    """
    Calculates the sum of thermal pollutions for boiler.
    Parameters
    ----------
    pollution_1_boiler : float
        The thermal pollution of the first coolant, [m**2 * degrees celcium / W]
    pollution_2_boiler : float
        The thermal pollution of the second coolant, [m**2 * degrees celcium / W]
    sigma_boiler : float
        The thickness of pipe wall, [m]
    lyambda_wall_boiler : float
        The thermal conducivity of wall, [W / (m * degrees celcium)]  
    Returns
    -------
    sigma_thermpollution_boiler : float
        The the sum of thermal pollutions, [m**2 * degrees celcium / W]
    References
    ----------
    &&&&&
    """                  
    return (sigma_boiler / lyambda_wall_boiler) + (1 / pollution_1_boiler) + (1 / pollution_2_boiler)


def A_real_boiler(Q_boiler, q_boiler):
    """
    Calculates the boiler's real heatransfer area.
    Parameters
    ----------
    Q_boiler : float
        The heat load of boiler, [W]
    q_boiler : float
        The unit heat load of boiler, [W/m**2]
    Returns
    -------
    A_real_boiler : float
        The boiler's real heatransfer area, [m**2]
    References
    ----------
    &&&&&
    """     
    return Q_boiler/q_boiler


def surface_margin_boiler(A_approx_boiler, A_real_boiler):
    """
    Calculates the boiler's surface margin .
    Parameters
    ----------
    A_approximate_boiler : float
        The approximate heat ransfer area, [m**2]
    A_real_boiler : float
        The real heat transfer area, [m**2]
    Returns
    -------
    surface_margin_boiler: float
        The surface margin [%]
    References
    ----------
    &&&&
    """          
    return (A_approx_boiler - A_real_boiler) * 100 / A_approx_boiler