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
def d_enter_feed_column(F_mass, rho_F_tboil, w_liq):
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
def w_liq_real_feed_column(F_mass, d_enter_feed_column_real, rho_F_tboil):
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
def w_enter_steam_boiler_real(m_steam_boil, rho_steam_vapor, d_enter_steam_boiler_real):
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


@unitcheck(m_steam_boil="kg/s", rho_water_liq="kg/m**3", w_drift="m/s", res_unit="m")
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


@unitcheck(m_steam_boil="kg/s", rho_water_liq="kg/m**3", d_out_cond_boiler_real="m", res_unit="m/s")
def w_out_cond_boiler_real(m_steam_boil, rho_water_liq, d_out_cond_boiler_real):
    """
    Calculates the real speed of condensat at the tube.
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
    w_out_cond_boiler_real : float
        The real speed of condensat at the tube, [m/s]
    References
    ----------
    &&&
    """  
    return m_steam_boil/((d_out_cond_boiler_real^2)*rho_water_liq*0,785)


@unitcheck(m_steam_feed="kg/s", rho_steam="kg/m**3", w_vapor="m/s", res_unit="m")
def d_enter_steam_feed(m_steam_feed, rho_steam, w_vapor):
    """
    Calculates the tube's diameter of enter steam to feed heat exchanger.
    Parameters
    ----------
    m_steam_feed : float
        The mass flow rate of steam, [kg/s]
    rho_steam : float
        The density of steam at boilling temperature, [kg/m**3]
    w_vapor :float
        The speed of steam at the tube, [m/s]
    Returns
    -------
    d_enter_steam_feed : float
        The tube's diameter of enter feed to column from heat exchanger, [m]
    References
    ----------
    &&&
    """  
    return m_steam_feed/(0,785*rho_steam*w_vapor)


@unitcheck(m_steam_feed="kg/s", rho_steam_vapor="kg/m**3", d_enter_steam_feed_real="m", res_unit="m/s")
def w_enter_steam_feed_real(m_steam_feed, rho_steam_vapor, d_enter_steam_feed_real):
    """
    Calculates the real speed of vapor at the tube.
    Parameters
    ----------
    m_steam_feed : float
        The mass flow rate of steam, [kg/s]
    rho_steam_vapor : float
        The density of waste vapor at boilling temperature, [kg/m**3]
    d_enter_steam_feed_real : float
        The real tube's diameter, [m]
    Returns
    -------
    w_enter_steam_feed_real : float
        The real speed of vapor at the tube, [m/s]
    References
    ----------
    &&&
    """  
    return m_steam_feed/((d_enter_steam_feed_real^2)*rho_steam_vapor*0,785)


@unitcheck(m_steam_feed="kg/s", rho_water_liq="kg/m**3", w_drift="m/s", res_unit="m")
def d_out_cond_feed(m_steam_feed, rho_water_liq, w_drift):
    """
    Calculates the tube's diameter of out condensat from heat exchanger.
    Parameters
    ----------
    m_steam_feed : float
        The mass flow rate of steam, [kg/s]
    rho_water_liq : float
        The density of liquid at boilling temperature, [kg/m**3]
    w_drift :float
        The speed of steam at the tube, [m/s]
    Returns
    -------
    d_out_cond_feed : float
        The tube's diameter of out waste to column from heat exchanger, [m]
    References
    ----------
    &&&
    """  
    return m_steam_feed/(0,785*rho_water_liq*w_drift)


@unitcheck(m_steam_feed="kg/s", rho_water_liq="kg/m**3", d_out_cond_feed_real="m", res_unit="m/s")
def w_out_cond_feed_real(m_steam_feed, rho_water_liq, d_out_cond_feed_real):
    """
    Calculates the real speed of condensat at the tube.
    Parameters
    ----------
    m_steam_feed : float
        The mass flow rate of steam, [kg/s]
    rho_water_liq : float
        The density of liquid at boilling temperature, [kg/m**3]
    d_out_cond_feed_real : float
        The real tube's diameter, [m]
    Returns
    -------
    w_out_cond_feed_real : float
        The real speed of condensat at the tube, [m/s]
    References
    ----------
    &&&
    """  
    return m_steam_feed/((d_out_cond_feed_real^2)*rho_water_liq*0,785)


@unitcheck(G_mass="kg/s", rho_dist="kg/m**3", w_drift="m/s", res_unit="m")
def d_out_cond_deph(G_mass, rho_dist, w_drift):
    """
    Calculates the tube's diameter of out distilliat from dephlegmator.
    Parameters
    ----------
    G_mass : float
        The mass flow rate of distilliat, [kg/s]
    rho_dist : float
        The density of distilliat at boilling temperature, [kg/m**3]
    w_drift :float
        The speed of steam at the tube, [m/s]
    Returns
    -------
    d_out_cond_deph : float
        The tube's diameter of out distilliat from dephlegmator, [m]
    References
    ----------
    &&&
    """  
    return G_mass/(0,785*rho_dist*w_drift)


@unitcheck(G_mass="kg/s", rho_dist="kg/m**3", d_out_deph_real="m", res_unit="m/s")
def w_out_deph_real(G_mass,rho_dist, d_out_deph_real):
    """
    Calculates the real speed of condensat at the tube.
    Parameters
    ----------
    G_mass : float
        The mass flow rate of distilliat, [kg/s]
    rho_dist : float
        The density of liquid at boilling temperature, [kg/m**3]
    d_out_deph_real : float
        The real tube's diameter, [m]
    Returns
    -------
    w_out_deph_real : float
        The real speed of condensat at the tube, [m/s]
    References
    ----------
    &&&
    """  
    return G_mass/((d_out_deph_real^2)*rho_dist*0,785)


@unitcheck(P_mass="kg/s", rho_dist="kg/m**3", w_drift="m/s", res_unit="m")
def d_enter_dist_cooler(P_mass, rho_dist, w_drift):
    """
    Calculates the tube's diameter of enter distilliat to distilliat cooler.
    Parameters
    ----------
    P_mass : float
        The mass flow rate of distilliat, [kg/s]
    rho_dist : float
        The density of liquid at boilling temperature, [kg/m**3]
    w_drift :float
        The speed of steam at the tube, [m/s]
    Returns
    -------
    d_enter_dist_cooler : float
        The tube's diameter of enter distilliat to distilliat cooler, [m]
    References
    ----------
    &&&
    """  
    return P_mass/(0,785*rho_dist*w_drift)


@unitcheck(P_mass="kg/s", rho_dist ="kg/m**3", d_enter_dist_cooler_real="m", res_unit="m/s")
def w_enter_dist_cooler_real(P_mass, rho_dist , d_enter_dist_cooler_real):
    """
    Calculates the real speed of liquid at the tube.
    Parameters
    ----------
    P_mass : float
        The mass flow rate of distilliat, [kg/s]
    rho_dist  : float
        The density of distilliat at boilling temperature, [kg/m**3]
    d_enter_dist_cooler_real : float
        The real tube's diameter, [m]
    Returns
    -------
    w_enter_dist_cooler_real : float
        The real speed of liquid at the tube, [m/s]
    References
    ----------
    &&&
    """  
    return P_mass/((d_enter_dist_cooler_real^2)*rho_dist *0,785)


@unitcheck(P_mass="kg/s", rho_dist_cool="kg/m**3", w_drift="m/s", res_unit="m")
def d_out_dist_cooler(P_mass, rho_dist_cool, w_drift):
    """
    Calculates the tube's diameter of out distilliat from distilliat cooler to distilliat volume.
    Parameters
    ----------
    P_mass : float
        The mass flow rate of distilliat, [kg/s]
    rho_dist_cool : float
        The density of liquid at cooling temperature, [kg/m**3]
    w_drift :float
        The speed of steam at the tube, [m/s]
    Returns
    -------
    d_out_dist_cooler : float
        The tube's diameter of out distilliat from distilliat cooler to distilliat volume, [m]
    References
    ----------
    &&&
    """  
    return P_mass/(0,785*rho_dist_cool*w_drift)


@unitcheck(P_mass="kg/s", rho_dist_cool ="kg/m**3", d_out_dist_cooler_real="m", res_unit="m/s")
def w_out_dist_cooler_real(P_mass, rho_dist_cool , d_out_dist_cooler_real):
    """
    Calculates the real speed of liquid at the tube.
    Parameters
    ----------
    P_mass : float
        The mass flow rate of distilliat, [kg/s]
    rho_dist_cool : float
        The density of distilliat at boilling temperature, [kg/m**3]
    d_out_dist_cooler_real : float
        The real tube's diameter, [m]
    Returns
    -------
    w_out_dist_cooler_real : float
        The real speed of liquid at the tube, [m/s]
    References
    ----------
    &&&
    """  
    return P_mass/((d_out_dist_cooler_real^2)*rho_dist_cool *0,785)


@unitcheck(W_mass="kg/s", rho_waste="kg/m**3", w_drift="m/s", res_unit="m")
def d_enter_waste_cooler(W_mass, rho_waste, w_drift):
    """
    Calculates the tube's diameter of enter waste to waste cooler.
    Parameters
    ----------
    W_mass : float
        The mass flow rate of waste, [kg/s]
    rho_waste : float
        The density of liquid at boilling temperature, [kg/m**3]
    w_drift :float
        The speed of steam at the tube, [m/s]
    Returns
    -------
    d_enter_waste_cooler : float
        The tube's diameter of enter waste to waste cooler, [m]
    References
    ----------
    &&&
    """  
    return W_mass/(0,785*rho_waste*w_drift)


@unitcheck(W_mass="kg/s", rho_waste ="kg/m**3", d_enter_waste_cooler_real="m", res_unit="m/s")
def w_enter_waste_cooler_real(W_mass, rho_waste , d_enter_waste_cooler_real):
    """
    Calculates the real speed of liquid at the tube.
    Parameters
    ----------
    W_mass : float
        The mass flow rate of waste, [kg/s]
    rho_waste  : float
        The density of waste at boilling temperature, [kg/m**3]
    d_enter_waste_cooler_real : float
        The real tube's diameter, [m]
    Returns
    -------
    w_enter_waste_cooler_real : float
        The real speed of liquid at the tube, [m/s]
    References
    ----------
    &&&
    """  
    return W_mass/((d_enter_waste_cooler_real^2)*rho_waste *0,785)


@unitcheck(W_mass="kg/s", rho_waste_cool="kg/m**3", w_drift="m/s", res_unit="m")
def d_out_waste_cooler(W_mass, rho_waste_cool, w_drift):
    """
    Calculates the tube's diameter of out waste from waste cooler to waste volume.
    Parameters
    ----------
    W_mass : float
        The mass flow rate of waste, [kg/s]
    rho_waste_cool : float
        The density of liquid at cooling temperature, [kg/m**3]
    w_drift :float
        The speed of steam at the tube, [m/s]
    Returns
    -------
    d_out_waste_cooler : float
        The tube's diameter of out waste from waste cooler to waste volume, [m]
    References
    ----------
    &&&
    """  
    return W_mass/(0,785*rho_waste_cool*w_drift)


@unitcheck(W_mass="kg/s", rho_waste_cool ="kg/m**3", d_out_waste_cooler_real="m", res_unit="m/s")
def w_out_waste_cooler_real(W_mass, rho_waste_cool , d_out_waste_cooler_real):
    """
    Calculates the real speed of liquid at the tube.
    Parameters
    ----------
    W_mass : float
        The mass flow rate of waste, [kg/s]
    rho_waste_cool : float
        The density of waste at boilling temperature, [kg/m**3]
    d_out_waste_cooler_real : float
        The real tube's diameter, [m]
    Returns
    -------
    w_out_waste_cooler_real : float
        The real speed of liquid at the tube, [m/s]
    References
    ----------
    &&&
    """  
    return W_mass/((d_out_waste_cooler_real^2)*rho_waste_cool *0,785)


@unitcheck(m_coolwater_dist="kg/s", rho_dist_coolwater="kg/m**3", w_liq="m/s", res_unit="m")
def d_dist_coolwater(m_coolwater_dist, rho_dist_coolwater, w_liq):
    """
    Calculates the tube's diameter of enter and out cooling water to distilliat cooler.
    Parameters
    ----------
    m_coolwater_dist : float
        The mass flow rate of cooling water, [kg/s]
    rho_dist_coolwater : float
        The density of cool water, [kg/m**3]
    w_liq :float
        The speed of liquid at the tube, [m/s]
    Returns
    -------
    d_dist_cooler : float
        The tube's diameter of enter and out cooling water to distilliat cooler, [m]
    References
    ----------
    &&&
    """  
    return m_coolwater_dist/(0,785*rho_dist_coolwater*w_liq)


@unitcheck(m_coolwater_dist="kg/s", rho_dist_coolwater ="kg/m**3", d_dist_coolwater_real="m", res_unit="m/s")
def w_dist_coolwater_real(m_coolwater_dist, rho_dist_coolwater , d_dist_coolwater_real):
    """
    Calculates the real speed of liquid at the tube.
    Parameters
    ----------
    m_coolwater_dist : float
        The mass flow rate of cooling water, [kg/s]
    rho_dist_coolwater : float
        The density of cool water, [kg/m**3]
    d_dist_coolwater_real : float
        The real tube's diameter, [m]
    Returns
    -------
    w_dist_coolwater_real : float
        The real speed of liquid at the tube, [m/s]
    References
    ----------
    &&&
    """  
    return m_coolwater_dist/((d_dist_coolwater_real^2)*rho_dist_coolwater *0,785)


@unitcheck(m_coolwater_waste="kg/s", rho_waste_coolwater="kg/m**3", w_liq="m/s", res_unit="m")
def d_waste_coolwater(m_coolwater_waste, rho_waste_coolwater, w_liq):
    """
    Calculates the tube's diameter of enter and out cooling water to waste cooler.
    Parameters
    ----------
    m_coolwater_waste : float
        The mass flow rate of cooling water, [kg/s]
    rho_waste_coolwater : float
        The density of cool water, [kg/m**3]
    w_liq :float
        The speed of liquid at the tube, [m/s]
    Returns
    -------
    d_waste_cooler : float
        The tube's diameter of enter and out cooling water to waste cooler, [m]
    References
    ----------
    &&&
    """  
    return m_coolwater_waste/(0,785*rho_waste_coolwater*w_liq)


@unitcheck(m_coolwater_waste="kg/s", rho_waste_coolwater ="kg/m**3", d_waste_coolwater_real="m", res_unit="m/s")
def w_waste_coolwater_real(m_coolwater_waste, rho_waste_coolwater , d_waste_coolwater_real):
    """
    Calculates the real speed of liquid at the tube.
    Parameters
    ----------
    m_coolwater_waste : float
        The mass flow rate of cooling water, [kg/s]
    rho_waste_coolwater : float
        The density of cool water, [kg/m**3]
    d_waste_coolwater_real : float
        The real tube's diameter, [m]
    Returns
    -------
    w_waste_coolwater_real : float
        The real speed of liquid at the tube, [m/s]
    References
    ----------
    &&&
    """  
    return m_coolwater_waste/((d_waste_coolwater_real^2)*rho_waste_coolwater *0,785)


@unitcheck(m_coolwater_deph="kg/s", rho_deph_coolwater="kg/m**3", w_liq="m/s", res_unit="m")
def d_deph_coolwater(m_coolwater_deph, rho_deph_coolwater, w_liq):
    """
    Calculates the tube's diameter of enter and out cooling water to dephlegmator.
    Parameters
    ----------
    m_coolwater_deph : float
        The mass flow rate of cooling water, [kg/s]
    rho_deph_coolwater : float
        The density of cool water, [kg/m**3]
    w_liq :float
        The speed of liquid at the tube, [m/s]
    Returns
    -------
    d_deph_cooler : float
        The tube's diameter of enter and out cooling water to dephlegmator, [m]
    References
    ----------
    &&&
    """  
    return m_coolwater_deph/(0,785*rho_deph_coolwater*w_liq)


@unitcheck(m_coolwater_deph="kg/s", rho_deph_coolwater ="kg/m**3", d_deph_coolwater_real="m", res_unit="m/s")
def w_deph_coolwater_real(m_coolwater_deph, rho_deph_coolwater , d_deph_coolwater_real):
    """
    Calculates the real speed of liquid at the tube.
    Parameters
    ----------
    m_coolwater_deph : float
        The mass flow rate of cooling water, [kg/s]
    rho_deph_coolwater : float
        The density of cool water, [kg/m**3]
    d_deph_coolwater_real : float
        The real tube's diameter, [m]
    Returns
    -------
    w_deph_coolwater_real : float
        The real speed of liquid at the tube, [m/s]
    References
    ----------
    &&&
    """  
    return m_coolwater_deph/((d_deph_coolwater_real^2)*rho_deph_coolwater *0,785)

