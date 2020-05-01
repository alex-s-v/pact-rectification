from rectification.utils import unitcheck



@unitcheck(xf_mass="kg/kg", M_lc="kg/kmol", M_hc="kg/kmol", res_unit="kmol/kmol")
def xf_mol(xf_mass, M_lc, M_hc):
    """
    Calculates the mol composition of feed.
    Parameters
    ----------
    xf_mass : float
        The mass concentration low component of feed, [kg/kg]
    M_lc : float
        The molar mass of low-boilling component, [kg/kmol]
    M_hc : float
        The molar mass of high-boilling component, [kg/kmol]
    Returns
    -------
    xf_mol : float
        The mol concentration low component of feed, [kmol/kmol]
    References
    ----------
    Романков, стр.283, табл. 6.2. 
    """      
    return (xf_mass * M_hc) / ((xf_mass * M_hc) + (M_lc - M_lc * xf_mass))


@unitcheck(xp_mass="kg/kg", M_lc="kg/kmol", M_hc="kg/kmol", res_unit="kmol/kmol")
def xp_mol(xp_mass, M_lc, M_hc):
    """
    Calculates the mol composition of distilliat.
    Parameters
    ----------
    xp_mass : float
        The mass concentration low component of distilliat, [kg/kg]
    M_lc : float
        The molar mass of low-boilling component, [kg/kmol]
    M_hc : float
        The molar mass of high-boilling component, [kg/kmol]
    Returns
    -------
    xp_mol : float
        The mol concentration low component of distilliat, [kmol/kmol]
    References
    ----------
    Романков, стр.283, табл. 6.2. 
    """      
    return (xp_mass * M_hc) / ((xp_mass * M_hc) + (M_lc - M_lc * xp_mass))


@unitcheck(xw_mass="kg/kg", M_lc="kg/kmol", M_hc="kg/kmol", res_unit="kmol/kmol")
def xw_mol(xw_mass, M_lc, M_hc):
    """
    Calculates the mol composition of waste.
    Parameters
    ----------
    xw_mass : float
        The mass concentration low component of waste, [kg/kg]
    M_lc : float
        The molar mass of low-boilling component, [kg/kmol]
    M_hc : float
        The molar mass of high-boilling component, [kg/kmol]
    Returns
    -------
    xw_mol : float
        The mol concentration low component of waste, [kmol/kmol]
    References
    ----------
    Романков, стр.283, табл. 6.2. 
    """      
    return (xw_mass * M_hc) / ((xw_mass * M_hc) + (M_lc - M_lc * xw_mass))


@unitcheck(xf_mol="kmol/kmol", M_lc="kg/kmol", M_hc="kg/kmol", res_unit="kg/kmol")
def M_feed(xf_mol, M_lc, M_hc):
    """
    Calculates the molar mass of feed.
    Parameters
    ----------
    xf_mol : float
        The mol concentration low component of feed, [kmol/kmol]
    M_lc : float
        The molar mass of low-boilling component, [kg/kmol]
    M_hc : float
        The molar mass of high-boilling component, [kg/kmol]
    Returns
    -------
    M_feed : float
        The molar mass of feed, [kg/kmol]
    References
    ----------
    Дытнерский, стр. 230, формула 6.6
    """ 
    return (M_lc * xf_mol + M_hc * (1 - xf_mol))


@unitcheck(xp_mol="kmol/kmol", M_lc="kg/kmol", M_hc="kg/kmol", res_unit="kg/kmol")
def M_dist(xp_mol, M_lc, M_hc):
    """
    Calculates the molar mass of distilliat.
    Parameters
    ----------
    xp_mol : float
        The mol concentration low component of distilliat, [kmol/kmol]
    M_lc : float
        The molar mass of low-boilling component, [kg/kmol]
    M_hc : float
        The molar mass of high-boilling component, [kg/kmol]
    Returns
    -------
    M_dist : float
        The molar mass of distilliat, [kg/kmol]
    References
    ----------
    Дытнерский, стр. 230, формула 6.6
    """ 
    return (M_lc * xp_mol + M_hc * (1 - xp_mol))


@unitcheck(xw_mol="kmol/kmol", M_lc="kg/kmol", M_hc="kg/kmol", res_unit="kg/kmol")
def M_waste(xw_mol, M_lc, M_hc):
    """
    Calculates the molar mass of waste.
    Parameters
    ----------
    xw_mol : float
        The mol concentration low component of waste, [kmol/kmol]
    M_lc : float
        The molar mass of low-boilling component, [kg/kmol]
    M_hc : float
        The molar mass of high-boilling component, [kg/kmol]
    Returns
    -------
    M_waste : float
        The molar mass of waste, [kg/kmol]
    References
    ----------
    Дытнерский, стр. 230, формула 6.6
    """ 
    return (M_lc * xw_mol + M_hc * (1 - xw_mol))


@unitcheck(F_mass="kg/s", xp_mass="kg/kg", xf_mass="kg/kg", x_w="kg/kg", res_unit="kg/s")
def W_mass(F_mass, xp_mass, xf_mass, xw_mass):
    """
    Calculates the mass flow rate of waste.
    Parameters
    ----------
    F_mass : float
        The mass flow rate of feed [kg/s]
    xp_mass : float
        The mass concetration distilliat, [kg/kg]
    xw_mass : float
        The mass concetration waste, [kg/kg]
    xf_mass : float
        The mass concetration feed, [kg/kg]
    Returns
    -------
    W_mass : float
        The mass flow rate of waste. [kg/s]
    References
    ----------
    Дытнерский, стр. 228, формула 6.1
    """ 
    return F_mass * (xp_mass - xf_mass) / (xp_mass - xw_mass)


@unitcheck(F_mass="kg/s", W_mass="kg/s", res_unit="kg/s")
def P_mass(F_mass, W_mass):
    """
    Calculates the mass flow rate of distilliat.
    Parameters
    ----------
    F_mass : float
        The mass flow rate of feed [kg/s]
    W_mass : float
        The mass flow rate of waste. [kg/s]
    Returns
    -------
    P_mass : float
        The mass flow rate of distilliat. [kg/s]
    References
    ----------
    Дытнерский, стр. 228, формула 6.1
    """    
    return F_mass - W_mass


@unitcheck(F_mass="kg/s", M_feed="kg/kmol", res_unit="kmol/s")
def F_mol(F_mass, M_feed):
    """
    Calculates the molar flow rate of feed.
    Parameters
    ----------
    F_mass : float
        The mass flow rate of feed [kg/s]
    M_feed : float
        The molar mass feed of liquid. [kg/kmol]
    Returns
    -------
    F_mol : float
        The molar flow rate of feed, [kmol/s]
    References
    ----------
    ???
    """        
    return F_mass / M_feed


@unitcheck(P_mass="kg/s", M_dist="kg/kmol", res_unit="kmol/s")
def P_mol(P_mass, M_dist):
    """
    Calculates the molar flow rate of dist.
    Parameters
    ----------
    P_mass : float
        The mass flow rate of distilliat, [kg/s]
    M_dist : float
        The molar mass  of distilliat, [kg/kmol]
    Returns
    -------
    P_mol : float
        The molar flow rate of distilliat, [kmol/s]
    References
    ----------
    ???
    """        
    return P_mass / M_dist


@unitcheck(W_mass="kg/s", M_waste="kg/kmol", res_unit="kmol/s")
def W_mol(W_mass, M_waste):
    """
    Calculates the molar flow rate of waste.
    Parameters
    ----------
    W_mass : float
        The mass flow rate of waste, [kg/s]
    M_waste : float
        The molar mass  of waste, [kg/kmol]
    Returns
    -------
    P_mol : float
        The molar flow rate of waste, [kmol/s]
    References
    ----------
    ???
    """        
    return W_mass / M_waste