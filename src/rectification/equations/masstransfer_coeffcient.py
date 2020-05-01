from rectification.utils import unitcheck



@unitcheck(sigma_lc_top="N/m", sigma_hc_top="N/m", x_aver_top_mass="kmol/kmol", res_unit="N/m")
def sigma_top(sigma_lc_top, sigma_hc_top, x_aver_top_mass):
    """
    Calculates the surface tension at the top of column.
    Parameters
    ----------
    sigma_lc_top : float
        The surface tension of low-boilling component at the top of column, [N / m]
    sigma_hc_top : float
        The surface tension of high-boilling component at the top of column, [N / m]
    x_aver_top_mass : float
        The average mass concentration at top of column, [kg/kg]
    Returns
    -------
    sigma_top : float
        The surface tension at the top of column, [N / m]
    References
    ----------
    &&&&&
    """       
    return (sigma_lc_top * x_aver_top_mass  + (1 - x_aver_top_mass) * sigma_hc_top)


@unitcheck(sigma_lc_bot="N/m", sigma_hc_bot="N/m", x_aver_bot_mass="kmol/kmol", res_unit="N/m")
def sigma_bot(sigma_lc_bot, sigma_hc_bot, x_aver_bot_mass):
    """
    Calculates the surface tension at the bottom of column.
    Parameters
    ----------
    sigma_lc_bot : float
        The surface tension of low-boilling component at the bottom of column, [N / m]
    sigma_hc_bot : float
        The surface tension of high-boilling component at the bottom of column, [N / m]
    x_aver_bot_mass : float
        The average mass concentration at bot of column, [kg/kg]
    Returns
    -------
    sigma_bot : float
        The surface tension at the  bottom of column, [N / m]
    References
    ----------
    &&&&&
    """       
    return (sigma_lc_bot * x_aver_bot_mass  + (1 - x_aver_bot_mass) * sigma_hc_bot)


