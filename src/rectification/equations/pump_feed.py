from rectification.utils import unitcheck
from scipy.constants import g



@unitcheck(F_mass="kg/s", rho_F="kg/m**3", res_unit="m**3/h")
def Q_volume_feed(F_mass, rho_F):
    """
    Calculates the volume flow rate of feed.
    Parameters
    ----------
    F_mass : float
        The mass flow rate of feed, [kg/s]
    rho_F : float
        The destiny of feed, [kg/m**3]
    Returns
    -------
    Q_volume : float
        The  volume flow rate of feed, [m**3/h]
    References
    ----------
    ????? 
    """          
    return F_mass/rho_F


@unitcheck(lyambda_rubbery="m", L_tube_feed="m", d_enter_feed="m", w_liq_real_enter_feed="m/s", epsi_local_resistance_sum_feed="m", res_unit="m")
def H_losses_feed(lyambda_rubbery, L_tube_feed, d_enter_feed, w_liq_real_enter_feed, epsi_local_resistance_sum_feed):
    """
    Calculates the hydraulic losses.
    Parameters
    ----------
    lyambda_rubbery : float
        The coefficient of rubbery, [m]
    L_tube : float
        The length of tube, [m]
    d_enter_feed : float
        The tube's diameter of enter to the heat exchanger of Feed, [m]
    w_liq_real_enter_feed : float
        The real speed of liquid at the tube, [m/s]
    epsi_local_resistance_sum: float
        The sum of local resistance on tube, [m]       
    Returns
    -------
    H_losses : float
        The hydraulic losses, [m]
    References
    ----------
    ????? 
    """        
    return (lyambda_rubbery_feed*L_tube_feed/d_enter_feed)*(w_liq_real_enter_feed / 2*g) + epsi_local_resistance_sum_feed*w_liq_real_enter_feed / 2*g


def lyambda_rubbery_feed(lyambda_injection_feed, lyambda_heatexchanger_feed):
    """
    Calculates the coefficient of rubbery.
    Parameters
    ----------
    lyambda_injection : float
        The coefficient of rubbery at the injection tube, [m]
    lyambda_heatexchanger : float
        The coefficient of rubbery at the line of heat exchanger, [m]       
    Returns
    -------
    lyambda_rubbery : float
        The coefficient of rubbery, [m]
    References
    ----------
    ????? 
    """       
    return lyambda_injection_feed + lyambda_heatexchanger_feed



def lyambda_heatexchanger_feed(Re_feed, e_roughness_heatexchanger_feed, d_inner):
    """
    Calculates the coefficient of rubbery at the injection tube.
    Parameters
    ----------
    Re_feed : float
        The Re coefficient of rubbery at the injection tube, [dimensionless]
    e_roughness_heatexchanger_feed : float
        The roughness of rubbery at the injection tube, [m]
    d_inner : float
        The tube's diameter of heat exchanger of Feed, [m]       
    Returns
    -------
    lyambda_heatexchanger_feed : float
        The coefficient of rubbery at the line of heat exchanger, [m]
    References
    ----------
    ????? 
    """       
    return ((- 2 * log((e_roughness/d_inner)/3,7 + (6,81/Re_feed)^0,9))^(-1))^2


def lyambda_injection_feed(Re_injection_feed, e_roughness_injection_feed, d_enter_feed_real):
    """
    Calculates the coefficient of rubbery at the injection tube.
    Parameters
    ----------
    Re_injection_feed : float
        The Re coefficient of rubbery at the injection tube, [dimensionless]
    e_roughness : float
        The roughness of rubbery at the injection tube, [m]
    d_enter_feed_real : float
        The tube's diameter of enter to the heat exchanger of Feed, [m]       
    Returns
    -------
    lyambda_injection_feed : float
        The coefficient of rubbery at the injection tube, [m]
    References
    ----------
    ????? 
    """       
    return - 2 * log((e_roughness/d_enter_feed_real)/3,7 + (6,81/Re_injection_feed)^0,9)


def Re_injection_feed(w_liq_real_enter_feed, rho_F_20, d_enter_feed_real, mu_F_20):
    """
    Calculates the Re criterion .
    Parameters
    ----------
    w_liq_real_enter_feed : float
        The real speed of liquid at the tube, [m/s]
    rho_F_20 : float
        The density of feed at 20 degrees, [kg/m**3]
    d_enter_feed_real : float
        The real tube's diameter, [m]
    mu_F_20 : float
        The viscosity of feed at 20 degrees, [Pa * s]        
    Returns
    -------
    Re_injection_feed : float
        The Re criterion, [dimensionless]
    References
    ----------
    ????? 
    """       
    return w_liq_real_enter_feed * rho_F_20 * d_enter_feed_real / mu_F_20


def H_hydrohead_feed(H_losses_feed, H_geometric_high_feed):
    """
    Calculates the head of pump.
    Parameters
    ----------
    H_losses_feed : float
        The hydraulic losses, [m]
    H_geometric_high_feed : float
        The geometric high supply of pump, [m]     
    Returns
    -------
    H_hydrohead_feed : float
        The head of pump, [m]
    References
    ----------
    ????? 
    """        
    return H_geometric_high_feed + H_losses_feed


def N_power_feed(Q_volume_feed, rho_F_avrg, g, H_hydrohead_feed_real, nu_motor_efficiency, nu_supply_efficiency):
    """
    Calculates the power of pump.
    Parameters
    ----------
    H_losses_feed_real : float
        The hydraulic losses, [m]
    Q_volume_feed : float
        The  volume flow rate of feed, [m**3/s]
    rho_F_avrg : float
        The averange density of feed at the line, [kg/m**3]
    nu_motor_efficiency: float
        The motor efficiency of pump, [dismensionless]
    nu_supply_efficiency: float
        The supply efficiency of pump, [dismensionless]     
    Returns
    -------
    N_power_feed : float
        The power of pump, [W]
    References
    ----------
    ????? 
    """    
    return Q_volume_feed * rho_F_avrg * g * H_hydrohead_feed_real / (nu_motor_efficiency * nu_supply_efficiency)


def hydraulic_losses_suct_feed(dzeta_enter_feed, dzeta_turn90_feed, n_turn90_feed, dzeta_ventil_feed, n_ventil_feed, g, w_liq_real_enter_feed):
    """
    Calculates the hydraulic losses of suction line.
    Parameters
    ----------
    dzeta_enter_feed : float
        The local resistance of tube enter, [m]
    dzeta_turn90_feed : float
        The local resistance of turn to 90 degrees, [m]
    n_turn90_feed : float
        The quantity of turn to 90 degrees, [dismensionless]
    dzeta_ventil_feed : float
        The local resistance of ventil on sunction line, [m]
    n_ventil_feed : float
        The quantity of ventil on suction line, [dismensionless]
    speed_suction : float
        The speed of suction line , [m/s] 
    Returns
    -------
    hydraulic_losses_suct_feed : float
        The hydraulic losses of suction line, [m]
    References
    ----------
    &&&&
    """         
    return ((dzeta_enter_feed + dzeta_turn90_feed + dzeta_ventil_feed) * w_liq_real_enter_feed/(2 * g))


def heigth_cavitation_feed(Q_volume_feed, n_turnover):
    """
    Calculates the losses due to cavitation.
    Parameters
    ----------
    Q_volume_feed : float
        The flow rate pump for Feed [m**3/s]
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
    return 0.3 * (Q_volume_feed * n_turnover**2)**(2/3)

def heigth_max_suction(Pa, rho_F_20, g, P_satur_vapor_feed, w_liq_real_enter_feed, hydraulic_losses_suct_feed, heigth_cavitation_feed):
    """
    Calculates the maximum theoretical suction height.
    Parameters
    ----------
    Pa : float
        The atmosphere pressure [Pa]
    rho_F_20 : float
        The density of feed, [kg/m**3]
    P_satur_vapor : float
        The  pressure of saturated vapor, [Pa]
    w_liq_real_enter_feed : float
        The speed of suction line , [m/s]
    hydraulic_losses_suct : float
        The hydraulic losses of suction line, [m]      
    heigth_cavitation_feed : float
        The losses due to cavitation, [m]  
    Returns
    -------
    heigth_max_suction : float
        The maximum theoretical suction height, [m]
    References
    ----------
    &&&&
    """   
    return ((Pa/(rho_F_20 * g) - ((P_satur_vapor_feed)/(rho_F_20 * g) + ((w_liq_real_enter_feed / (2 * g))) + hydraulic_losses_satur_feed + heigth_cavitation_feed)))