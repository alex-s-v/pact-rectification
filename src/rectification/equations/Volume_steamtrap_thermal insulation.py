from rectification.utils import unitcheck
import pandas as pd





@unitcheck(F_mass="kg/s", rho_F_20="kg/m**3", tau="s", res_unit="m**3")
def V_tank_F(F_mass, tau, rho_F_20, dzeta_reserve):
    """
    Calculates the tank for feed.
    Parameters
    ----------
    F_mass : float
        The mass flowrate of feed, [kg/s]
    tau : float
        The time, [s]
    rho_f_20 : float
        The destiny of feed for 20 degrees celcium, [kg/m**3]
    dzeta_reserve : float
        The coeefificent of reserve, [dismensionless]
    Returns
    -------
    V_tank_F : float
        The tank for feed, [m**3]
    References
    ----------
    &
    """     
    return F_mass * tau * dzeta_reserve / rho_F_20


@unitcheck(P_mass="kg/s", rho_P_20="kg/m**3", tau="s", res_unit="m**3")
def V_tank_P(P_mass, tau, rho_P_20, dzeta_reserve):
    """
    Calculates the tank for distilliat.
    Parameters
    ----------
    P_mass : float
        The mass flowrate of distilliat, [kg/s]
    tau : float
        The time, [s]
    rho_P_20 : float
        The destiny of distilliat for 20 degrees celcium, [kg/m**3]
    dzeta_reserve : float
        The coeefificent of reserve, [dismensionless]
    Returns
    -------
    V_tank_P : float
        The tank for distilliat, [m**3]
    References
    ----------
    &
    """     
    return P_mass * tau * dzeta_reserve / rho_P_20


@unitcheck(W_mass="kg/s", rho_W_20="kg/m**3", tau="s", res_unit="m**3")
def V_tank_W(W_mass, tau, rho_W_20, dzeta_reserve):
    """
    Calculates the tank for waste.
    Parameters
    ----------
    W_mass : float
        The mass flowrate of waste, [kg/s]
    tau : float
        The time, [s]
    rho_W_20 : float
        The destiny of waste for 20 degrees celcium, [kg/m**3]
    dzeta_reserve : float
        The coeefificent of reserve, [dismensionless]
    Returns
    -------
    V_tank_W : float
        The tank for waste, [m**3]
    References
    ----------
    &
    """     
    return W_mass * tau * dzeta_reserve / rho_W_20


@unitcheck(Reflux_mass="kg/s", rho_Reflux_20="kg/m**3", tau="s", res_unit="m**3")
def V_tank_Reflux(Reflux_mass, tau, rho_Reflux_20, dzeta_reserve):
    """
    Calculates the tank for waste.
    Parameters
    ----------
    Reflux_mass : float
        The mass flowrate of Reflux, [kg/s]
    tau : float
        The time, [s]
    rho_Reflux_20 : float
        The destiny of waste for 20 degrees celcium, [kg/m**3]
    dzeta_reserve : float
        The coefificent of reserve, [dismensionless]
    Returns
    -------
    V_tank_Reflux : float
        The tank for Reflux, [m**3]
    References
    ----------
    &&&&&&&&&&&&
    """     
    return Reflux_mass * tau * dzeta_reserve / rho_Reflux_20


@unitcheck(m_steam_boil="t/h", delta_P_boiler="MPa", res_unit="t/h")
def k_steamtrap_boiler(m_steam_boil, delta_P_boiler):
    """
    Calculates the steam trap for boiler.
    Parameters
    ----------
    m_steam_boil : float
        The flow rate steam of boiler, [t/h]
    delta_P_boiler : float
        The differential pressure between steam pressure and atmospheric pressure, [MPa]
    Returns
    -------
    k_steamtrap_boiler : float
        The steam trap for boiler, [t/h]
    References
    ----------
    &&&&&&&&&&&&
    """        
    return 0,575 * m_steam_boil / (delta_P_boiler)^(0,5)


@unitcheck(m_steam_feed="t/h", delta_P_feed="MPa", res_unit="t/h")
def k_steamtrap_feed(m_steam_feed, delta_P_feed):
    """
    Calculates the steam trap for feed.
    Parameters
    ----------
    m_steam_feed : float
        The flow rate steam of feed, [t/h]
    delta_P_feed : float
        The differential pressure between steam pressure and atmospheric pressure, [MPa]
    Returns
    -------
    k_steamtrap_boiler : float
        The steam trap for boiler, [t/h]
    References
    ----------
    &&&&&&&&&&&&
    """        
    return 0,575 * m_steam_feed / (delta_P_feed)^(0,5)


