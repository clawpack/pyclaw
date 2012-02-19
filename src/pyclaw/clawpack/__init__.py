__all__ = []
__all__.extend(['classic1','classic2','classic3'])

try:
    import classic1, classic2,classic3
except:
    print "Some Clawpack object files missing.  Please make in pyclaw/clawpack/ before using."
