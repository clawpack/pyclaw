def mapc2p_sphere(xc,yc):
    """
    Maps to points on a sphere of radius Rsphere. 
    """

    # Radius of the earth
    a = 6.37122e6
    r1 = a

    if (xc >= 1.0):
        xc = xc - 4.0
    
    if (xc <= -3.0):
        xc = xc + 4.0

    if (yc >= 1.0):
        yc = 2.0 - yc
        xc = -2.0 - xc

    if (yc <= -1.0):
        yc = -2.0 - yc
        xc = -2.0 - xc

    if (xc <= -1.0):
       # Points in [-3,-1] map to lower hemisphere - reflect about x=-1
       # to compute x,y mapping and set sgnz appropriately:
       xc = -2.0 - xc
       sgnz = -1.0
    else:
       sgnz = 1.0

    import math
    sgnxc = math.copysign(1.0,xc)
    sgnyc = math.copysign(1.0,yc)

    xc1 = np.abs(xc)
    yc1 = np.abs(yc)
    d = np.maximum(np.maximum(xc1,yc1), 1.0e-10)     

    DD = r1*d*(2.0 - d) / np.sqrt(2.0)
    R = r1
    center = DD - np.sqrt(np.maximum(R**2 - DD**2, 0.0))
            
    xp = DD/d * xc1
    yp = DD/d * yc1

    if (yc1 >= xc1):
        yp = center + np.sqrt(np.maximum(R**2 - xp**2, 0.0))
    else:
        xp = center + np.sqrt(np.maximum(R**2 - yp**2, 0.0))

    # Compute physical coordinates
    zp = np.sqrt(np.maximum(r1**2 - (xp**2 + yp**2), 0.0))

    return xp,yp,zp

