import numpy as nmp

class SGrid:
    """ Equivalent of Delphi's SGrid <T> with some modifications"""
    def __init__(self, arg_nX, arg_nY, arg_nZ, arg_val, **kwargs):
        """ Constructor"""
        self.nX, self.nY, self.nZ = arg_nX, arg_nY, arg_nZ
        self.value = arg_val

        # set attributes from the received keywords
        for prop, val in kwargs.items():
            setattr(self, prop, val)
        
    def __str__(self):
        """ String manifestation """
        return 'Grid idx : {} {} {} with Phi = {:>10.4f} (eps = {:>10.4f} and deb = {:>10.4f})'.format(self.nX, self.nY, self.nZ, self.value, self.eps, self.deb)
    
    def __le__(self, arg_other):
        return self.value <= arg_other.value
    
    def __lt__(self, arg_other):
        return self.value < arg_other.value
    
    def __ge__(self, arg_other):
        return self.value >= arg_other.value
    
    def __gt__(self, arg_other):
        return self.value > arg_other.value
    
    def __eq__(self, arg_other):
        return self.value == arg_other.value
    
    def get_grid_indices(self):
        return [self.nX, self.nY, self.nZ]

    def distance2(self, arg_other, arg_gridRes):
        """ Distance squared of a grid point from another one computed as
        dist2 = sum(Delta(i)**2 times gridRes**2,
                    Delta(j)**2 times gridRes**2,
                    Delta(k)**2 times gridRes**2)"""
        
        del_I2 = nmp.power(self.nX - arg_other.nX, 2)
        del_J2 = nmp.power(self.nY - arg_other.nY, 2)
        del_K2 = nmp.power(self.nZ - arg_other.nZ, 2)
        
        dist2 = (del_I2 + del_J2 + del_K2)*nmp.power(arg_gridRes, 2)
        return dist2
        
    def distance(self, arg_other, arg_gridRes):
        return sqrt(self.distance2(arg_other, arg_gridRes))
        
