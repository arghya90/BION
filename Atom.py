""" A class named Atom that stores information about an atom in a structure
as is read from its corresponding PDB/GRO/PQR etc. file """

class Atom(object):
    
    def __init__(self, arg_serial, **kwargs):
        self.serial = arg_serial
        for attr,value in kwargs.items():
            setattr(self, attr, value)
#END


