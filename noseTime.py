
#This class compute the wall time (WT) of each the regression test



from time import time

import operator

import nose
from nose.plugins.base import Plugin as plg

class timer(plg):
 
    def WallTime(self):
        if hasattr(self, '_timer'):
            taken = time() - self._timer
        else:
            taken = 0.0
        return taken

if __name__ == '__main__':
    nose.main(addplugins=[timer()])


