def mapc2p(xc,yc):
    """
    Specifies the mapping to curvilinear coordinates    
    """
    import numpy as np

    # Polar coordinates (x coordinate = radius,  y coordinate = theta)
    xp = xc * np.cos(yc)
    yp = xc * np.sin(yc)
    return xp,yp
