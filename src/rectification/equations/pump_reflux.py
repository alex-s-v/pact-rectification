from rectification.utils import unitcheck
from scipy.constants import g



@unitcheck(reflux_mass="kg/s", rho_F="kg/m**3", res_unit="m**3/h")
def Q_volume_reflux(reflux_mass, rho_reflux):
    """
    Calculates the volume flow rate of reflux.
    Parameters
    ----------
    reflux_mass : float
        The mass flow rate of reflux, [kg/s]
    rho_reflux : float
        The destiny of reflux, [kg/m**3]
    Returns
    -------
    Q_volume : float
        The  volume flow rate of reflux, [m**3/h]
    References
    ----------
    ????? 
    """          
    return reflux_mass/rho_reflux

def reflux_mass(P_mass, R):
    """
    Calculates the flow rate of reflux.
    Parameters
    ----------
    P_mass : float
        The mass flow rate of distilliat, [kg/s]
    R : float
        The reflux number, [dismensionless]
    Returns
    -------
    reflux_mass : float
        The flow rate of reflux, [kg/s]
    References
    ----------
    ????? 
    """
    return P_mass * R


@unitcheck(lyambda_rubbery_reflux="m", L_tube_reflux="m", d_enter_reflux="m", w_liq_real_enter_reflux="m/s", epsi_local_resistance_sum_reflux="m", res_unit="m")
def H_losses_reflux(lyambda_rubbery_reflux, L_tube_reflux, d_enter_reflux, w_liq_real_enter_reflux, epsi_local_resistance_sum_reflux):
    """
    Calculates the hydraulic losses.
    Parameters
    ----------
    lyambda_rubbery : float
        The coefficient of rubbery, [m]
    L_tube : float
        The length of tube, [m]
    d_enter_reflux : float
        The tube's diameter of enter to the heat exchanger of Feed, [m]
    w_liq_real_enter_reflux : float
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
    return (lyambda_rubbery_reflux*L_tube_reflux/d_enter_reflux)*(w_liq_real_enter_reflux / 2*g) + epsi_local_resistance_sum_reflux*w_liq_real_enter_reflux / 2*g


def lyambda_rubbery_reflux(Re_reflux, e_roughness_reflux, d_enter_reflux_real):
    """
    Calculates the coefficient of rubbery.
    Parameters
    ----------
    Re_reflux : float
        The Re coefficient of rubbery at the injection tube, [dimensionless]
    e_roughness_reflux : float
        The roughness of rubbery at the tube, [m]
    d_enter_reflux_real : float
        The tube's diameter, [m]        
    Returns
    -------
    lyambda_rubbery_reflux : float
        The coefficient of rubbery, [m]
    References
    ----------
    ????? 
    """       
    return - 2 * log((e_roughness_reflux/d_enter_reflux_real)/3,7 + (6,81/Re_reflux)^0,9)


def Re_reflux(w_liq_real_enter_reflux, rho_reflux, d_enter_reflux_real, mu_reflux):
    """
    Calculates the Re criterion .
    Parameters
    ----------
    w_liq_real_enter_reflux : float
        The real speed of liquid at the tube, [m/s]
    rho_reflux : float
        The density of reflux, [kg/m**3]
    d_enter_reflux_real : float
        The real tube's diameter, [m]
    mu_reflux : float
        The viscosity of Reflux, [Pa * s]        
    Returns
    -------
    Re_reflux : float
        The Re criterion, [dimensionless]
    References
    ----------
    ????? 
    """       
    return w_liq_real_enter_reflux * rho_reflux * d_enter_reflux_real / mu_reflux


def H_hydrohead_reflux(H_losses_reflux, H_geometric_high_reflux):
    """
    Calculates the head of pump.
    Parameters
    ----------
    H_losses_reflux : float
        The hydraulic losses, [m]
    H_geometric_high_reflux : float
        The geometric high supply of pump, [m]     
    Returns
    -------
    H_hydrohead_reflux : float
        The head of pump, [m]
    References
    ----------
    ????? 
    """        
    return H_geometric_high_reflux + H_losses_reflux


def N_power_reflux(Q_volume_reflux, rho_reflux_avrg, g, H_hydrohead_reflux_real, nu_motor_efficiency, nu_supply_efficiency):
    """
    Calculates the power of pump.
    Parameters
    ----------
    H_losses_reflux_real : float
        The hydraulic losses, [m]
    Q_volume_reflux : float
        The  volume flow rate of reflux, [m**3/s]
    rho_reflux : float
        The density of reflux at the line, [kg/m**3]
    nu_motor_efficiency: float
        The motor efficiency of pump, [dismensionless]
    nu_supply_efficiency: float
        The supply efficiency of pump, [dismensionless]     
    Returns
    -------
    N_power_reflux : float
        The power of pump, [W]
    References
    ----------
    ????? 
    """    
    return Q_volume_reflux * rho_reflux * g * H_hydrohead_reflux_real / (nu_motor_efficiency * nu_supply_efficiency)


def hydraulic_losses_suct_reflux(dzeta_enter_reflux, dzeta_turn90_reflux, n_turn90_reflux, dzeta_ventil_reflux, n_ventil_reflux, g, w_liq_real_enter_reflux):
    """
    Calculates the hydraulic losses of suction line.
    Parameters
    ----------
    dzeta_enter_reflux : float
        The local resistance of tube enter, [m]
    dzeta_turn90_reflux : float
        The local resistance of turn to 90 degrees, [m]
    n_turn90_reflux : float
        The quantity of turn to 90 degrees, [dismensionless]
    dzeta_ventil_reflux : float
        The local resistance of ventil on sunction line, [m]
    n_ventil_reflux : float
        The quantity of ventil on suction line, [dismensionless]
    speed_suction : float
        The speed of suction line , [m/s] 
    Returns
    -------
    hydraulic_losses_suct_reflux : float
        The hydraulic losses of suction line, [m]
    References
    ----------
    &&&&
    """         
    return ((dzeta_enter_reflux + dzeta_turn90_reflux + dzeta_ventil_reflux) * w_liq_real_enter_reflux/(2 * g))


def heigth_cavitation_reflux(Q_volume_reflux, n_turnover):
    """
    Calculates the losses due to cavitation.
    Parameters
    ----------
    Q_volume_reflux : float
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
    return 0.3 * (Q_volume_reflux * n_turnover**2)**(2/3)

def heigth_max_suction(Pa, rho_reflux_20, g, P_satur_vapor_reflux, w_liq_real_enter_reflux, hydraulic_losses_suct_reflux, heigth_cavitation_reflux):
    """
    Calculates the maximum theoretical suction height.
    Parameters
    ----------
    Pa : float
        The atmosphere pressure [Pa]
    rho_reflux_20 : float
        The density of reflux, [kg/m**3]
    P_satur_vapor : float
        The  pressure of saturated vapor, [Pa]
    w_liq_real_enter_reflux : float
        The speed of suction line , [m/s]
    hydraulic_losses_suct : float
        The hydraulic losses of suction line, [m]      
    heigth_cavitation_reflux : float
        The losses due to cavitation, [m]  
    Returns
    -------
    heigth_max_suction : float
        The maximum theoretical suction height, [m]
    References
    ----------
    &&&&
    """
    return ((Pa/(rho_reflux_20 * g) - ((P_satur_vapor_reflux)/(rho_reflux_20 * g) + ((w_liq_real_enter_reflux / (2 * g))) + hydraulic_losses_satur_reflux + heigth_cavitation_reflux)))