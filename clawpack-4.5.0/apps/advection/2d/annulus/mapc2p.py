

def mapc2p(xc,yc):
    """
    Specifies the mapping to curvilinear coordinates -- should be consistent
    with mapc2p.f
    """
    from numpy import sin, cos

    # Polar coordinates (xc = r,  yc = theta)
    xp = xc * cos(yc)
    yp = xc * sin(yc)
    return xp,yp
