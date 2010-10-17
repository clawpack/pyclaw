"""
Tools for interpolating data.

Includes:
    pwcubic:        piecewise cubic interpolation.
    pwlinear:       piecewise linear interpolation.
"""


def pwcubic(xi, zl, zr, slopel, sloper, x):
    """
    Construct a piecewise cubic function and evaluate at the points in x.
    The interpolation data consists of points xi, and values
        zl[i] = value of z in limit as x approaches xi[i] from the left,
        zr[i] = value of z in limit as x approaches xi[i] from the right,
        slopel[i] = value of z'(x) in limit as x approaches xi[i] from the left,
        sloper[i] = value of z'(x) in limit as x approaches xi[i] from the right.
    To the left of xi[0] the linear function based on 
       zl[0] and slopel[0] is used.
    To the right of xi[-1] the linear function based on 
       zr[-1] and sloper[-1] is used.

    In the interior intervals a cubic us used that interpolates 
    the 4 values at the ends of the interval.  

    Note that the function will be linear in the i'th interval if
        sloper[i]   = (zl[i+1] - zr[i]) / (xi[i+1] - xi[i])
        slopel[i+1] = (zl[i+1] - zr[i]) / (xi[i+1] - xi[i])

    A piecewise linear interpolant can be created by setting this everywhere:
        slopel[1:]  = (zl[1:]-zr[:-1]) / (xi[1:] - xi[:-1])
        sloper[:-1] = (zl[1:]-zr[:-1]) / (xi[1:] - xi[:-1])
    This can be done automatically with the pwlinear function defined below.

    Note that the pwcubic function is continuous if zl == zr and is C^1 if
    in addition slopel == sloper.
    """

    from numpy import where, zeros

    dx = xi[1:] - xi[:-1]

    s = (zl[1:] - zr[:-1]) / dx
    c2 = (s - sloper[:-1]) / dx
    c3 = (slopel[1:] - 2.*s + sloper[:-1]) / (dx**2)

    # set to linear function for x<xi[0] or x>= xi[-1]
    z = where(x < xi[0],  zl[0] + (x-xi[0]) * slopel[0], 
                          zr[-1] + (x-xi[-1]) * sloper[-1])

    # replace by appropriate cubic in intervals:
    for i in range(len(xi)-1):
        cubic = zr[i] + sloper[i]*(x - xi[i]) \
                    + (x-xi[i])**2 * (c2[i] + c3[i]*(x-xi[i+1]))
        z = where((x >= xi[i]) & (x < xi[i+1]), cubic, z)

    return z


def pwlinear(xi, zl, zr, x, extrap=0):
    """
    Construct a piecewise linear function and evaluate at the points in x.
    The interpolation data consists of points xi, and values
        zl[i] = value of z in limit as x approaches xi[i] from the left,
        zr[i] = value of z in limit as x approaches xi[i] from the right,
    Extrapolate outside the interval by:
        if extrap==0:  slope = 0 outside of intervals
        if extrap==1:  slope taken from neighboring interval.
    Note that the pwlinear function is continuous if zl == zr.
    """

    from numpy import zeros
    slopel = zeros(len(zl))
    sloper = zeros(len(zr))
    dx = xi[1:] - xi[:-1]
    slopel[1:]  = (zl[1:]-zr[:-1]) / dx
    sloper[:-1] = slopel[1:]
    if extrap==1:
        slopel[0] = sloper[0]
        sloper[-1] = slopel[-1]
    z = pwcubic(xi, zl, zr, slopel, sloper, x)
    return z

