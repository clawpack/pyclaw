
def test_acoustics():
    import sys
    sys.path.append('./apps/acoustics/1d/homogeneous')
    import acoustics

    error=acoustics.acoustics(makePlot=False)
    print error

    assert(abs(error-0.0044043)<1.e-5)
