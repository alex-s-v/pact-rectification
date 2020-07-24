from rectification.utils import unitcheck
import numpy as np




def Ky(beta_liq, beta_vapor, m_distrib):
    """
    Calculates the masstransfer coefficent Ky.
    Parameters
    ----------
    beta_liq : float
        The masstransfer coefficent beta of liquid, [kmol / (m**2 * s)]
    beta_vapor : float
        The masstransfer coefficent beta of vapor, [kmol / (m**2 * s)]
    m_distrib : float
        The distribution coefficient, [dismensionless]
    Returns
    -------
    Ky : float
        The masstransfer coefficent Ky, [kmol / (m**2 * s)]
    References
    ----------
    Дытнерский, стр.194 формула 5.8
    """        
    return ((1 / beta_vapor) + (m_distrib / beta_liq))**(-1)


def noy(Ky, M, rho, w):
    """
    Calculates the number of units at the plate.
    Parameters
    ----------
    Ky : float
        The masstransfer coefficent Ky, [kmol / (m**2 * s)]
    M : float
        The molar mass, [kg / kmol]
    rho : float
        The destiny , [kg / m**3]
    w : float
        The speed [m / s]
    Returns
    -------
    noy : float
        The number of units at the plate, [dismensionless]
    References
    ----------
    Дытнерский, стр.239 формула 6.35
    """    
    return Ky * M / (rho * w)


def lyambda_fact_top(m_distrib, R):
    """
    Calculates the the factor of masstransfer.
    Parameters
    ----------
    m_distrib : float
        The distribution coefficient, [dismensionless]
    R : float
        The reflux number, [dismensionless]
    Returns
    -------
    lyambda_fact_top : float
        The factor of masstransfer, [dismensionless]
    References
    ----------
    Дытнерский, стр.239 формула 6.35
    """      
    return m_distrib * (R + 1) / R


def lyambda_fact_bot(m_distrib, R, P_mol, F_mol, phi):
    """
    Calculates the the factor of masstransfer.
    Parameters
    ----------
    m_distrib : float
        The distribution coefficient, [dismensionless]
    R : float
        The reflux number, [dismensionless]
    P_mol : float
        The flow rate of distilliat, [kmol/s]
    F_mol : float
        The flow rate of feed, [kmol/s]
    phi : float
        The fraction of vapor at the feed point
    Returns
    -------
    lyambda_fact_bot : float
        The factor of masstransfer, [dismensionless]
    References
    ----------
    Дытнерский, стр.239 формула 6.35
    """      
    return m_distrib * (P_mol * (R + 1) - F_mol*phi) / (P_mol * R + F_mol * (1- phi))


def Lt(D, Lc):
    """
    Calculates the length way of liquid
    Parameters
    ----------
    D : float
        The diameter of the column, [m]
    Lc : float
        The drain perimeter, [m]
    Returns
    -------
    Lt : float
        The length way of liquid, [m]
    References
    ----------
    &&&&
    """     
    return (D**2 - Lc**2)**(0.5)


def S(Lt, l):
    """
    Calculates the S coefficient
    Parameters
    ----------
    Lt : float
        The length way of liquid, [m]
    l : float
        The length way of liquid of one slot mix, [m]
    Returns
    -------
    Lt : float
        The S coefficient, [dismensionless]
    References
    ----------
    &&&&
    """     
    return Lt / l


def H_column(H_bwplate, N_plates, Zt, Zl, H_cap, H_support):
    """
    Calculates the S coefficient
    Parameters
    ----------
    H_bwplate : float
        The heigth of between plates, [m]  
    N_plates : float
        The number of plates, [dismensionless]
    Zt : float
        The separation space of top over the highest plate, [m]
    Zl : float
        The length between the lowest plate and bottom, [m] 
    H_cap : float
        The high of cap, [m]
    H_support : float
        The high of support, [m] 
    Returns
    -------
    Lt : float
        The S coefficient, [dismensionless]
    References
    ----------
    &&&&
    """       
    return ((N_plates - 1) * H_bwplate + Zt + Zl + H_cap + H_support)