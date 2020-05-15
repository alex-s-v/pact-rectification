from rectification.utils import unitcheck





@unitcheck(F_mass="kg/s", rho_F_20="kg/m**3", w_liq="m/s", res_unit="m")
def d_enter_feed(F_mass, rho_F_20, w_liq):
    """
    Calculates the tube's diameter of enter to the heat exchanger of Feed.
    Parameters
    ----------
    F_mass : float
        The mass flow rate of feed [kg/s]
    rho_F_20 : float
        The density of feed at 20 degrees, [kg/m**3]
    w_liq :float
        The speed of liquid at the tube, [m/s]
    Returns
    -------
    d_enter_feed : float
        The tube's diameter of enter to the heat exchanger of Feed, [m]
    References
    ----------
    &&&
    """  
    return F_mass/(0,785*rho_F_20*w_liq)


@unitcheck(F_mass="kg/s", rho_F_20="kg/m**3", d_enter_feed_real="m", res_unit="m/s")
def w_liq_real_enter_feed(F_mass, d_enter_feed_real, rho_F_20):
    """
    Calculates the real speed of liquid at the tube.
    Parameters
    ----------
    F_mass : float
        The mass flow rate of feed [kg/s]
    rho_F_20 : float
        The density of feed at 20 degrees, [kg/m**3]
    d_enter_feed_real : float
        The real tube's diameter, [m]
    Returns
    -------
    w_liq_real_enter_feed : float
        The real speed of liquid at the tube, [m/s]
    References
    ----------
    &&&
    """  
    return F_mass/((d_enter_feed_real^2)*rho_F_20*0,785)


@unitcheck(F_mass="kg/s", rho_F_tboil="kg/m**3", w_liq="m/s", res_unit="m")
def d_enter_feed_column(F_mass, rho_F_tboil_feed, w_liq):
    """
    Calculates the tube's diameter of enter to the column at the feed plate.
    Parameters
    ----------
    F_mass : float
        The mass flow rate of feed [kg/s]
    rho_F_tboil : float
        The density of feed at the temperature of feed boiling , [kg/m**3]
    w_liq :float
        The speed of liquid at the tube, [m/s]
    Returns
    -------
    d_enter_feed_column : float
        The tube's diameter of enter to the column at the feed plate, [m]
    References
    ----------
    &&&&
    """  
    return F_mass/(0,785*rho_F_tboil*w_liq)


@unitcheck(F_mass="kg/s", rho_F_tboil="kg/m**3", d_enter_feed_column_real="m", res_unit="m/s")
def w_liq_real_feed_column(F_mass, d_enter_feed_column, rho_F_tboil):
    """
    Calculates the real speed of liquid at the tube.
    Parameters
    ----------
    F_mass : float
        The mass flow rate of feed [kg/s]
    rho_F_tboil : float
        The density of feed at the temperature of feed boiling , [kg/m**3]
    d_enter_feed_column_real : float
        The real tube's diameter, [m]
    Returns
    -------
    w_liq_real_enter_feed_column : float
        The real speed of liquid at the tube, [m/s]
    References
    ----------
    &&&
    """  
    return F_mass/((d_enter_feed_column_real^2)*rho_F_tboil*0,785)


@unitcheck(G_mass="kg/s", rho_P_vapor="kg/m**3", w_vapor="m/s", res_unit="m")
def d_out_dist(G_mass, rho_P_vapor, w_vapor):
    """
    Calculates the tube's diameter of out vapor distilliat from column to dephlegmator.
    Parameters
    ----------
    G_mass : float
        The mass flow rate of vapor, [kg/s]
    rho_P_vapor : float
        The density of distilliat vapor at boilling temperature, [kg/m**3]
    w_vapor :float
        The speed of vapor at the tube, [m/s]
    Returns
    -------
    d_out_dist : float
        The tube's diameter of enter to the heat exchanger of Feed, [m]
    References
    ----------
    &&&
    """  
    return G_mass/(0,785*rho_P_vapor*w_vapor)


@unitcheck(G_mass="kg/s", rho_P_vapor="kg/m**3", d_out_dist_real="m", res_unit="m/s")
def w_vapor_real_out_dist(G_mass, d_out_dist_real, rho_P_vapor):
    """
    Calculates the real speed of liquid at the tube.
    Parameters
    ----------
    G_mass : float
        The mass flow rate of vapor, [kg/s]
    rho_P_vapor : float
        The density of feed at 20 degrees, [kg/m**3]
    d_out_dist_real : float
        The real tube's diameter, [m]
    Returns
    -------
    w_vapor_real_out_dist : float
        The real speed of vapor at the tube, [m/s]
    References
    ----------
    &&&
    """  
    return G_mass/((d_out_dist_real^2)*rho_P_vapor*0,785)


@unitcheck(Reflux_mass="kg/s", rho_P_liq="kg/m**3", w_liq="m/s", res_unit="m")
def d_enter_reflux(Reflux_mass, rho_P_liq, w_liq):
    """
    Calculates the tube's diameter of out vapor distilliat from column to dephlegmator.
    Parameters
    ----------
    Reflux_mass : float
        The mass flow rate of reflux, [kg/s]
    rho_P_liq : float
        The density of distilliat liquid at boilling temperature, [kg/m**3]
    w_liq :float
        The speed of liquid at the tube, [m/s]
    Returns
    -------
    d_enter_reflux : float
        The tube's diameter of out vapor distilliat from column to dephlegmator, [m]
    References
    ----------
    &&&
    """  
    return Reflux_mass/(0,785*rho_P_liq*w_liq)


@unitcheck(Reflux_mass="kg/s", rrho_P_liq="kg/m**3", d_enter_reflux_real="m", res_unit="m/s")
def w_enter_reflux_real(Reflux_mass, rho_P_liq, d_enter_reflux_real):
    """
    Calculates the real speed of liquid at the tube.
    Parameters
    ----------
    Reflux_mass : float
        The mass flow rate of reflux, [kg/s]
    rho_P_liq : float
        The density of distilliat liquid at boilling temperature, [kg/m**3]
    d_enter_reflux_real : float
        The real tube's diameter, [m]
    Returns
    -------
    w_enter_reflux_real : float
        The real speed of liquid at the tube, [m/s]
    References
    ----------
    &&&
    """  
    return Reflux_mass/((d_enter_reflux_real^2)*rho_P_liq*0,785)


@unitcheck(W_mass="kg/s", rho_W_liq="kg/m**3", w_liq_drift="m/s", res_unit="m")
def d_enter_waste_boiler(W_mass, rho_W_liq, w_liq_drift):
    """
    Calculates the tube's diameter of enter waste to boiler from column.
    Parameters
    ----------
    W_mass : float
        The mass flow rate of waste, [kg/s]
    rho_W_liq : float
        The density of waste liquid at boilling temperature, [kg/m**3]
    w_liq_drift :float
        The drift speed of liquid at the tube, [m/s]
    Returns
    -------
    d_enter_reflux : float
        The tube's diameter of enter waste to boiler from column, [m]
    References
    ----------
    &&&
    """  
    return W_mass/(0,785*rho_W_liq*w_liq_drift)


@unitcheck(W_mass="kg/s", rho_W_liq="kg/m**3", d_enter_waste_boiler_real="m", res_unit="m/s")
def w_enter_waste_boiler_real(W_mass, rho_W_liq, d_enter_waste_boiler_real):
    """
    Calculates the real speed of liquid at the tube.
    Parameters
    ----------
    W_mass : float
        The mass flow rate of waste, [kg/s]
    rho_W_liq : float
        The density of waste liquid at boilling temperature, [kg/m**3]
    d_enter_waste_boiler_real : float
        The real tube's diameter, [m]
    Returns
    -------
    w_enter_waste_boiler_real : float
        The real speed of liquid at the tube, [m/s]
    References
    ----------
    &&&
    """  
    return W_mass/((d_enter_waste_boiler_real^2)*rho_W_liq*0,785)


@unitcheck(W_mass="kg/s", rho_W_vapor="kg/m**3", w_vapor="m/s", res_unit="m")
def d_out_waste_boiler(W_mass, rho_W_vapor, w_vapor):
    """
    Calculates the tube's diameter of out waste to column from boiler.
    Parameters
    ----------
    W_mass : float
        The mass flow rate of waste, [kg/s]
    rho_W_vapor : float
        The density of waste liquid at boilling temperature, [kg/m**3]
    w_vapor :float
        The speed of vapor at the tube, [m/s]
    Returns
    -------
    d_out_waste_boiler : float
        The tube's diameter of out waste to column from boiler, [m]
    References
    ----------
    &&&
    """  
    return W_mass/(0,785*rho_W_vapor*w_vapor)


@unitcheck(W_mass="kg/s", rho_W_vapor="kg/m**3", d_out_waste_boiler_real="m", res_unit="m/s")
def w_out_waste_boiler_real(W_mass, rho_W_vapor, d_out_waste_boiler_real):
    """
    Calculates the real speed of liquid at the tube.
    Parameters
    ----------
    W_mass : float
        The mass flow rate of waste, [kg/s]
    rho_W_vapor : float
        The density of waste vapor at boilling temperature, [kg/m**3]
    d_out_waste_boiler_real : float
        The real tube's diameter, [m]
    Returns
    -------
    w_out_waste_boiler_real : float
        The real speed of vapor at the tube, [m/s]
    References
    ----------
    &&&
    """  
    return W_mass/((d_out_waste_boiler_real^2)*rho_W_vapor*0,785)


@unitcheck(m_steam_boil="kg/s", rho_steam="kg/m**3", w_vapor="m/s", res_unit="m")
def d_enter_steam_boiler(m_steam_boil, rho_steam, w_vapor):
    """
    Calculates the tube's diameter of enter steam to boiler.
    Parameters
    ----------
    m_steam_boil : float
        The mass flow rate of steam, [kg/s]
    rho_steam : float
        The density of steam at boilling temperature, [kg/m**3]
    w_vapor :float
        The speed of steam at the tube, [m/s]
    Returns
    -------
    d_enter_steam_boiler : float
        The tube's diameter of out waste to column from boiler, [m]
    References
    ----------
    &&&
    """  
    return m_steam_boil/(0,785*rho_steam*w_vapor)


@unitcheck(m_steam_boil="kg/s", rho_steam_vapor="kg/m**3", d_enter_steam_boiler_real="m", res_unit="m/s")
def w_enter_steam_boiler(m_steam_boil, rho_steam_vapor, d_enter_steam_boiler_real):
    """
    Calculates the real speed of vapor at the tube.
    Parameters
    ----------
    m_steam_boil : float
        The mass flow rate of steam, [kg/s]
    rho_steam_vapor : float
        The density of waste vapor at boilling temperature, [kg/m**3]
    d_enter_steam_boiler_real : float
        The real tube's diameter, [m]
    Returns
    -------
    w_enter_steam_boiler_real : float
        The real speed of vapor at the tube, [m/s]
    References
    ----------
    &&&
    """  
    return m_steam_boil/((d_enter_steam_boiler_real^2)*rho_steam_vapor*0,785)


#@unitcheck(m_steam_boil="kg/s", rho_water_liq="kg/m**3", w_drift="m/s", res_unit="m")
def d_out_cond_boiler(m_steam_boil, rho_water_liq, w_drift):
    """
    Calculates the tube's diameter of out condensat from boiler.
    Parameters
    ----------
    m_steam_boil : float
        The mass flow rate of steam, [kg/s]
    rho_water_liq : float
        The density of liquid at boilling temperature, [kg/m**3]
    w_drift :float
        The speed of steam at the tube, [m/s]
    Returns
    -------
    d_out_cond_boiler : float
        The tube's diameter of out waste to column from boiler, [m]
    References
    ----------
    &&&
    """  
    return m_steam_boil/(0,785*rho_water_liq*w_drift)


#@unitcheck(m_steam_boil="kg/s", rho_water_liq="kg/m**3", d_out_cond_boiler_real="m", res_unit="m/s")
def w_enter_steam_boiler(m_steam_boil, rho_water_liq, d_out_cond_boiler_real):
    """
    Calculates the real speed of vapor at the tube.
    Parameters
    ----------
    m_steam_boil : float
        The mass flow rate of steam, [kg/s]
    rho_water_liq : float
        The density of liquid at boilling temperature, [kg/m**3]
    d_out_cond_boiler_real : float
        The real tube's diameter, [m]
    Returns
    -------
    w_enter_steam_boiler_real : float
        The real speed of vapor at the tube, [m/s]
    References
    ----------
    &&&
    """  
    return m_steam_boil/((d_enter_steam_boiler_real^2)*rho_steam_vapor*0,785)