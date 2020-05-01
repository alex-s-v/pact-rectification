from rectification.utils import unitcheck




#region Reflux and mass ways of liq and vapor
@unitcheck(xp_mol="kmol/kmol", y_mol_xf_equil="kmol/kmol", xf_mol="kmol/kmol")
def R_min (xp_mol, y_mol_xf_equil, xf_mol):
    """
    Calculates the minimal reflux number.
    Parameters
    ----------
    xp_mol : float
        The mol concentration low component of distilliat, [kmol/kmol]
    xf_mol : float
        The mol concentration low component of feed, [kmol/kmol]
    y_mol_xf_equil : float
        The equilibrium concentration low component of vapor at concentration of feed (liquid), [kmol/kmol]    
    Returns
    -------
    R_min : float
        The minimal reglux number, [dismensionless]
    References
    ----------
    Дытнерский, стр.228, формула 6.2
    """     
    return (xp_mol - y_mol_xf_equil) / (y_mol_xf_equil  - xf_mol)


def R (Rmin, beta):
    """
    Calculates the actual reflux number.
    Parameters
    ----------
    R_min : float
        The minimal reglux number, [dismensionless]
    beta : float
        The overflow of reflux, [dismensionless]    
    Returns
    -------
    R : float
        The actual reflux number, [dismensionless]
    References
    ----------
    Романков, стр.321, формула 7.11
    """         
    return Rmin * beta


@unitcheck(xp_mol="kmol/kmol", xpf_mol="kmol/kmol")
def x_aver_top(xp_mol, xpf_mol):
    """
    Calculates the average mol concentration at the top of column.
    Parameters
    ----------
    xp_mol : float
        The mol concentration  of distilliat, [kmol/kmol]
    xpf_mol : float
        The mol concentration of point of feed, [kmol/kmol]  
    Returns
    -------
    x_aver_top : float
        The average mol concentration at top of column, [kmol/kmol]
    References
    ----------
    Дытнерский, стр.230, формула 6.6
    """            
    return (xp_mol + xpf_mol) / 2


@unitcheck(xw_mol="kmol/kmol", xpf_mol="kmol/kmol")
def x_aver_bot(xw_mol, xpf_mol):
    """
    Calculates the average mol concentration  at bottom of column.
    Parameters
    ----------
    xw_mol : float
        The mol concentration  of waste, [kmol/kmol]
    xpf_mol : float
        The mol concentration of point of feed, [kmol/kmol]  
    Returns
    -------
    x_aver_top : float
        The average mol concentration at top of column, [kmol/kmol]
    References
    ----------
    Дытнерский, стр.230, формула 6.6
    """            
    return (xw_mol + xpf_mol) / 2


@unitcheck(M_lc="kg/kmol", x_aver_top="kmol/kmol", M_hc="kg/kmol", res_unit="kg/kmol")
def M_top(M_lc, x_aver_top, M_hc):
    """
    Calculates the molar mass at top of column.
    Parameters
    ----------
    x_aver_top : float
        The average mol concentration at top of column, [kmol/kmol]
    M_lc : float
        The molar mass of low-boilling component, [kg/kmol]
    M_hc : float
        The molar mass of high-boilling component, [kg/kmol]
    Returns
    -------
    M_top : float
        The molar mass at top of column, [kg/kmol]
    References
    ----------
    Дытнерский, стр. 230, формула 6.6
    """ 
    return (M_lc * x_aver_top + M_hc * (1 - x_aver_top))


@unitcheck(M_lc="kg/kmol", x_aver_bot="kmol/kmol", M_hc="kg/kmol", res_unit="kg/kmol")
def M_bot(M_lc, x_aver_bot, M_hc):
    """
    Calculates the molar mass at bottom of column.
    Parameters
    ----------
    x_aver_bot : float
        The average mol concentration at bottom of column, [kmol/kmol]
    M_lc : float
        The molar mass of low-boilling component, [kg/kmol]
    M_hc : float
        The molar mass of high-boilling component, [kg/kmol]
    Returns
    -------
    M_bot : float
        The molar mass at bottom of column, [kg/kmol]
    References
    ----------
    Дытнерский, стр. 230, формула 6.6
    """ 
    return (M_lc * x_aver_bot + M_hc * (1 - x_aver_bot))


@unitcheck(P_mass="kg/s", M_top="kg/kmol", M_dist="kg/kmol", res_unit="kg/s")
def L_top(P_mass, R, M_top, M_dist):
    """
    Calculates the flow rate liquid at the top of column.
    Parameters
    ----------
    P_mass : float
        The mass flow rate of distilliat, [kg/s]
    M_dist : float
        The molar mass  of distilliat, [kg/kmol]
    R : float
        The actual reflux number, [dismensionless]
    M_top : float
        The molar mass at top of column, [kg/kmol]    
    Returns
    -------
    L_top : float
        The  flow rate liquid at the top of column, [kg/s]
    References
    ----------
    Дытнерский, стр. 229, формула 6.4
    """ 
    return P_mass * R * M_top / M_dist


@unitcheck(P_mass="kg/s", M_bot="kg/kmol", M_dist="kg/kmol", F_mass="kg/s", M_feed="kg/kmol", res_unit="kg/s")
def L_bot(P_mass, R, M_bot, M_dist, F_mass, M_feed, phi):
    """
    Calculates the flow rate liquid at the bottom of column.
    Parameters
    ----------
    P_mass : float
        The mass flow rate of distilliat, [kg/s]
    F_mass : float
        The mass flow rate of feed, [kg/s]
    M_dist : float
        The molar mass  of distilliat, [kg/kmol]
    M_feed : float
        The molar mass  of feed, [kg/kmol]
    phi : float
        The fraction of vapor at the feed point
    R : float
        The actual reflux number, [dismensionless]
    M_bot : float
        The molar mass at bottom of column, [kg/kmol]    
    Returns
    -------
    L_bot : float
        The  flow rate liquid at the bottom of column, [kg/s]
    References
    ----------
    Дытнерский, стр. 229, формула 6.5
    """ 
    return (P_mass * R * M_bot / M_dist) + (F_mass * M_bot * (1 - phi) / M_feed)


@unitcheck(yp_mol="kmol/kmol", ypf_mol="kmol/kmol")
def y_aver_top(yp_mol, ypf_mol):
    """
    Calculates the average mol concentration at the top of column.
    Parameters
    ----------
    yp_mol : float
        The mol concentration  of distilliat, [kmol/kmol]
    ypf_mol : float
        The mol concentration of point of feed, [kmol/kmol]  
    Returns
    -------
    y_aver_top : float
        The average mol concentration at top of column, [kmol/kmol]
    References
    ----------
    Дытнерский, стр.230, формула 6.8
    """            
    return (yp_mol + ypf_mol) / 2


@unitcheck(yw_mol="kmol/kmol", ypf_mol="kmol/kmol")
def y_aver_bot(yw_mol, ypf_mol):
    """
    Calculates the average mol concentration  at bottom of column.
    Parameters
    ----------
    yw_mol : float
        The mol concentration  of waste, [kmol/kmol]
    ypf_mol : float
        The mol concentration of point of feed, [kmol/kmol]  
    Returns
    -------
    y_aver_top : float
        The average mol concentration at top of column, [kmol/kmol]
    References
    ----------
    Дытнерский, стр.230, формула 6.8
    """            
    return (yw_mol + ypf_mol) / 2


@unitcheck(M_lc="kg/kmol", y_aver_top="kmol/kmol", M_hc="kg/kmol", res_unit="kg/kmol")
def M_top_vap(M_lc, y_aver_top, M_hc):
    """
    Calculates the molar mass at top of column.
    Parameters
    ----------
    y_aver_top : float
        The average mol concentration at top of column, [kmol/kmol]
    M_lc : float
        The molar mass of low-boilling component, [kg/kmol]
    M_hc : float
        The molar mass of high-boilling component, [kg/kmol]
    Returns
    -------
    M_top : float
        The molar mass at top of column, [kg/kmol]
    References
    ----------
    Дытнерский, стр. 230, формула 6.8
    """ 
    return (M_lc * y_aver_top + M_hc * (1 - y_aver_top))


@unitcheck(M_lc="kg/kmol", y_aver_bot="kmol/kmol", M_hc="kg/kmol", res_unit="kg/kmol")
def M_bot_vap(M_lc, y_aver_bot, M_hc):
    """
    Calculates the molar mass at bottom of column.
    Parameters
    ----------
    y_aver_bot : float
        The average mol concentration at bottom of column, [kmol/kmol]
    M_lc : float
        The molar mass of low-boilling component, [kg/kmol]
    M_hc : float
        The molar mass of high-boilling component, [kg/kmol]
    Returns
    -------
    M_bot : float
        The molar mass at bottom of column, [kg/kmol]
    References
    ----------
    Дытнерский, стр. 230, формула 6.8
    """ 
    return (M_lc * y_aver_bot + M_hc * (1 - y_aver_bot))


@unitcheck(P_mass="kg/s", M_top_vap="kg/kmol", M_dist="kg/kmol", res_unit="kg/s")
def G_top(P_mass, R, M_top_vap, M_dist):
    """
    Calculates the flow rate liquid at the top of column.
    Parameters
    ----------
    P_mass : float
        The mass flow rate of distilliat, [kg/s]
    M_dist : float
        The molar mass  of distilliat, [kg/kmol]
    R : float
        The actual reflux number, [dismensionless]
    M_top_vap : float
        The molar mass at top of column, [kg/kmol]    
    Returns
    -------
    G_top : float
        The  flow rate of the vapor at the top of column, [kg/s]
    References
    ----------
    Дытнерский, стр. 229, формула 6.7
    """ 
    return P_mass * (R + 1) * M_top_vap / M_dist


@unitcheck(P_mass="kg/s", M_bot_vap="kg/kmol", M_dist="kg/kmol", F_mass="kg/s", M_feed="kg/kmol", res_unit="kg/s")
def G_bot(P_mass, R, M_bot_vap, M_dist, F_mass, M_feed, phi):
    """
    Calculates the flow rate liquid at the bottom of column.
    Parameters
    ----------
    P_mass : float
        The mass flow rate of distilliat, [kg/s]
    F_mass : float
        The mass flow rate of feed, [kg/s]
    M_dist : float
        The molar mass  of distilliat, [kg/kmol]
    M_feed : float
        The molar mass  of feed, [kg/kmol]
    phi : float
        The fraction of vapor at the feed point
    R : float
        The actual reflux number, [dismensionless]
    M_bot_vap : float
        The molar mass of the vapor at bottom of column, [kg/kmol]    
    Returns
    -------
    G_bot : float
        The  flow rate of the vapor at the bottom of column, [kg/s]
    References
    ----------
    Дытнерский, стр. 229, формула 6.7
    """ 
    return (P_mass * R * M_bot / M_dist) - (F_mass * M_bot_vap * phi/ M_feed) ## 2018 и 1880 проверка знака второго слагаемого


@unitcheck(xp_mass="kg/kg", xpf_mass="kg/kg")
def x_aver_top_mass(xp_mass, xpf_mass):
    """
    Calculates the average mass concentration at the top of column.
    Parameters
    ----------
    xp_mass : float
        The mass concentration  of distilliat, [kg/kg]
    xpf_mass : float
        The mass concentration of point of feed, [kg/kg]  
    Returns
    -------
    x_aver_top_mass : float
        The average mass concentration at top of column, [kg/kg]
    References
    ----------
    Дытнерский, стр.230, формула 6.8
    """            
    return (xp_mass + xpf_mass) / 2


@unitcheck(xw_mass="kg/kg", xpf_mass="kg/kg")
def x_aver_bot_mass(xw_mass, xpf_mass):
    """
    Calculates the average mass concentration  at bottom of column.
    Parameters
    ----------
    xw_mass : float
        The mass concentration  of waste, [kg/kg]
    xpf_mass : float
        The mass concentration of point of feed, [kg/kg]  
    Returns
    -------
    x_aver_bot_mass : float
        The average mass concentration at bot of column, [kg/kg]
    References
    ----------
    Дытнерский, стр.230, формула 6.8
    """            
    return (xw_mass + xpf_mass) / 2


@unitcheck(x_aver_top_mass="kg/kg", rho_lc_x_aver_top="kg/m**3", rho_hc_x_aver_top="kg/m**3", res_unit="kg/m**3")
def rho_top_liq(x_aver_top_mass, rho_lc_x_aver_top, rho_hc_x_aver_top):
    """
    Calculates the destiny of liquid at top of column.
    Parameters
    ----------
    x_aver_top_mass : float
        The average mass concentration at top of column, [kg/kg]
    rho_lc_x_aver_top : float
        The destiny of low-boilling component for mass average concentration at the top of column, [kg/m**3]
    rho_hc_x_aver_top : float
        The destiny of high-boilling component for mass average concentration at the top of column, [kg/m**3] 
    Returns
    -------
    rho_top_liq : float
        The destiny of liquid at top of column, [kg/m**3]
    References
    ----------
    Романков, стр.12, формула 1.3
    """            
    return ((x_aver_top_mass / rho_lc_x_aver_top) + ((1 - x_aver_top_mass) / rho_hc_x_aver_top))


@unitcheck(x_aver_bot_mass="kg/kg", rho_lc_x_aver_bot="kg/m**3", rho_hc_x_aver_bot="kg/m**3", res_unit="kg/m**3")
def rho_bot_liq(x_aver_bot_mass, rho_lc_x_aver_bot, rho_hc_x_aver_bot):
    """
    Calculates the destiny of liquid at the bottom of column.
    Parameters
    ----------
    x_aver_bot_mass : float
        The average mass concentration at the bottom of column, [kg/kg]
    rho_lc_x_aver_bot : float
        The destiny of low-boilling component for mass average concentration at the bottom of column, [kg/m**3]
    rho_hc_x_aver_bot : float
        The destiny of high-boilling component for mass average concentration at the bottom of column, [kg/m**3]   
    Returns
    -------
    rho_bot_liq : float
        The destiny of liquid at the bottom of column, [kg/m**3]
    References
    ----------
    Романков, стр.12, формула 1.3
    """            
    return ((x_aver_bot_mass / rho_lc_x_aver_bot) + ((1 - x_aver_bot_mass) / rho_hc_x_aver_bot))


#unitcheck(!!!!!)
def rho_top_vap(t_boil_y_top, T0, Vm, M_top_vap):
    """
    Calculates the destiny of vapor at the top of column.
    Parameters
    ----------
    t_boil_y_top : float
        The boilling temperature of low-bolling component of vapor at the top of column, [C]
    T0 : float
        The initial temperature, [K]
    M_top_vap : float
        The molar mass at top of column, [kg/kmol]    
    Returns
    -------
    rho_top_vap : float
        The destiny of vapor at top of column, [kg/m**3]
    References
    ----------
    Романков, стр.13, формула 1.5
    """            
    return (M_top_vap / Vm) * (T0 / (T0 + t_boil_y_top))


#@unitcheck(!!!!!)
def rho_bot_vap(t_boil_y_bot, T0, Vm, M_bot_vap):
    """
    Calculates the destiny of vapor at the bottom of column.
    Parameters
    ----------
    t_boil_y_bot : float
        The boilling temperature of low-bolling component of vapor at the bottom of column, [C]
    T0 : float
        The initial temperature, [K]
    M_bot_vap : float
        The molar mass at bottom of column, [kg/kmol]    
    Returns
    -------
    rho_bot_vap : float
        The destiny of vapor at bottom of column, [kg/m**3]
    References
    ----------
    Романков, стр.13, формула 1.5
    """            
    return (M_bot_vap / Vm) * (T0 / (T0 + t_boil_y_bot))  

#operating line #N(R+1) N
#output: G_top, G_bot, L_bot, L_top, Rmin
#endregion


#region The speed and diameter of column
#@unitcheck(rho_top_liq="kg/m**3", rho_top_vap="kg/m**3", res_unit="m/s" !!!!!!!!!!!)
def speed_limit_top(rho_top_liq, rho_top_vap):
    """
    Calculates the limit speed of the vapor in the column.
    Parameters
    ----------
    rho_top_vap : float
        The destiny of vapor at the top of column, [kg/m**3]
    rho_top_liq : float
        The destiny of liquid at top of column, [kg/m**3] 
    Returns
    -------
    speed_limit_top : float
        The limit speed of the vapor at the top of  column, [m/s]
    References
    ----------
    Дытнерский, стр.205, формула 5.33
    """   
    return 0.05 * ((rho_top_liq / rho_bot_vap)^0.5)


#@unitcheck(rho_bot_liq="kg/m**3", rho_bot_vap="kg/m**3", res_unit="m/s" !!!!!!!!!!!)
def speed_limit_bot(rho_bot_liq, rho_bot_vap):
    """
    Calculates the limit speed of the vapor in the column.
    Parameters
    ----------
    rho_bot_vap : float
        The destiny of vapor at the bottom of column, [kg/m**3]
    rho_bot_liq : float
        The destiny of liquid at the bottom of column, [kg/m**3] 
    Returns
    -------
    speed_limit_bot : float
        The limit speed of the vapor at the bottom in  column, [m/s]
    References
    ----------
    Дытнерский, стр.205, формула 5.33
    """   
    return 0.05 * ((rho_top_liq / rho_bot_vap)^0.5)


@unitcheck(G_top="kg/s", speed_limit_top="m/s", rho_top_vap="kg/m**3", res_unit="m")
def D_top(G_top, pi, speed_limit_top, rho_top_vap):
    """
    Calculates the top diameter of column.
    Parameters
    ----------
    rho_top_vap : float
        The destiny of vapor at the top of column, [kg/m**3]
    speed_limit_top : float
        The limit speed of the vapor at the top of  column, [m/s]
    G_top : float
        The  flow rate of the vapor at the top of column, [kg/s]    
    Returns
    -------
    D_top : float
        The top diameter of column, [m]
    References
    ----------
    &&&&&
    """       
    return (4 * G_top / (pi * speed_limit_top * rho_top_vap))^0.5


@unitcheck(G_bot="kg/s", speed_limit_bot="m/s", rho_bot_vap="kg/m**3", res_unit="m")
def D_bot(G_bot, pi, speed_limit_bot, rho_bot_vap):
    """
    Calculates the bottom diameter of column.
    Parameters
    ----------
    rho_bot_vap : float
        The destiny of vapor at the bottom of column, [kg/m**3]
    speed_limit_bot : float
        The limit speed of the vapor at the bottom in  column, [m/s]
    G_bot : float
        The  flow rate of the vapor at the bottom of column, [kg/s]    
    Returns
    -------
    D_bot : float
        The bottom diameter of column, [m]
    References
    ----------
    &&&&&
    """       
    return (4 * G_bot / (pi * speed_limit_bot * rho_bot_vap))^0.5


@unitcheck(speed_limit_top="m/s", D_top="m", D="m", res_unit="m/s")
def speed_section_top(speed_limit_top, D_top, D):
    """
    Calculates the section speed of vapor
    Parameters
    ----------
    speed_limit_top : float
        The limit speed of the vapor at the top of  column, [m/s]
    D_top : float
        The calculating top diameter of column, [m]
    D : float
        The choosing diameter of column, [m]   
    Returns
    -------
    speed_section_top : float
        The section speed of vapor, [m/s]
    References
    ----------
    &&&&&
    """         
    return speed_limit_top * (D_top / D)**2


@unitcheck(speed_limit_bot="m/s", D_bot="m", D="m", res_unit="m/s")
def speed_section_bot(speed_limit_bot, D_bot, D):
    """
    Calculates the section speed of vapor at the bottom
    Parameters
    ----------
    speed_limit_bot: float
        The limit speed of the vapor at the bottom of  column, [m/s]
    D_bot : float
        The calculating bottom diameter of column, [m]
    D : float
        The choosing diameter of column, [m]   
    Returns
    -------
    speed_section_bot : float
        The section speed of vapor at the bottom, [m/s]
    References
    ----------
    &&&&&
    """         
    return speed_limit_bot * (D_bot / D)**2


@unitcheck(speed_section_top="m/s", D="m", ft="m**2", res_unit="m/s")
def speed_operating_section_top(speed_section_top, D, ft):
    """
    Calculates the operating section speed of vapor at the top
    Parameters
    ----------
    speed_section_top : float
        The section speed of vapor at the top, [m]
    D : float
        The choosing diameter of column, [m]  
    ft : float
        The operating section of plate, [m**2] 
    Returns
    -------
    speed_operating_section_top : float
        The operating section speed of vapor at the top, [m/s]
    References
    ----------
    &&&&&
    """         
    return speed_section_top * (0.785 * D**2 / ft)


@unitcheck(speed_section_bot="m/s", D="m", ft="m**2", res_unit="m/s")
def speed_operating_section_bot(speed_section_bot, D, ft):
    """
    Calculates the operating section speed of vapor at the bottom
    Parameters
    ----------
    speed_section_bot : float
        The section speed of vapor at the bottom, [m]
    D : float
        The choosing diameter of column, [m]  
    ft : float
        The operating section of plate, [m**2] 
    Returns
    -------
    speed_operating_section_bot : float
        The operating section speed of vapor at the bottom, [m/s]
    References
    ----------
    &&&&&
    """         
    return speed_section_bot * (0.785 * D**2 / ft)

    #output D_top, D_bot, D, w_operating_section_top, w_operating_section_bot
    #endregion 
