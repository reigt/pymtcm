# Import math modules
from math import sqrt as sqrt

def compressive_height(
    d: float,
    w: float,
    c: float,
    As: float,
    As_c: float,
    Es: float,
    Ecm: float,
) -> float:
    """Compressive height of section based on State II. 

    Args:
        d (float): Effective height [mm]
        w (float): Width [mm]
        c (float): Cover [mm]
        As (float): Rebar area in tension [mm2]
        As_c (float): Rebar area in compression [mm2]
        Es (float): Young's modulus steel [MPa]
        Ecm (float): Young's modulus concrete

    Returns:
        float: Compressive height [mm]
    """
    
    # Geometry    
    rho_s = As/(w*d)
    rho_sc = As_c/(w*d)
    
    # Modular ratio
    alpha_e = Es/Ecm
    
    # Compressive height
    alpha = sqrt((alpha_e*rho_sc+alpha_e*rho_s)**2 + 2*alpha_e*(rho_sc*c/d+rho_s)) - (alpha_e*rho_sc + alpha_e*rho_s)
    xc = alpha*d
    
    return xc

def curv_factor(
    h: float,
    d: float,
    xc: float
) -> float:
    """Curvature factor

    Args:
        h (float): Height [mm]
        d (float): Effective height [mm]
        xc (float): Compressive height [mm]

    Returns:
        float: Curvature factor
    """
    
    return (h-xc)/(d-xc)