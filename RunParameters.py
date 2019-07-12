class RunParameters(object):
    """ This is a class that has as its attributes all the parameters
    that will be used for running external or internal programs/functions."""

    def __init__(self):
        # PARAMTERS FOR DELPHIPKA RUN
        self.dPkaParams = {'Remove HETATM                      ': 'T',
                        'Remove water molecule              ': 'T',
                        'HETATM in PQR format               ': 'F',
                        'Calculate More Residues            ': 'F',
                        'Output PQR file (With Topology)    ': 'T',
                        'Gaussian surface                   ': '1',
                        'Surface cutoff                     ': '20',
                        'Variance of Gaussian distribution  ': '0.70',
                        'Internal Dielectric                ': '8.0',
                        'External Dielectric                ': '80.0',
                        'Cluster Delimitation Threshold (A) ': '12.0',
                        'Kmean++ cluster number             ': 'AUTO',
                        'Hydrogen of GLU Attached to Atom   ': 'OE1',
                        'Hydrogen of ASP Attached to Atom   ': 'OD1',
                        'pH Initial Value                   ': '0.0',
                        'pH End Value                       ': '14.0',
                        'pH Interval                        ': '0.5'}

        # PARAMTERS FOR DELPHIPKA RUN
        self.delphiParams = {'perfil':'70', 
                    'sigma':'0.93',
                    'srfcut' :'20', 
                    'bndcon':'2', 
                    'prbrad':'1.4', 
                    'linit':'800', 
                    'maxc':'0.001'}

        # OTHER PARAMETERS
        self.otherParams = {}

        return
        
    # A FUNCTION TO SET PARAMTERS PROVIDED BY INPUT OPTIONS TO BION
    def add_run_parameter(self, arg_paramName, arg_paramValue):
        setattr(self, arg_paramName, arg_paramValue)
        self.otherParams[arg_paramName] = arg_paramValue
        return

    # A FUNCTION THAT PRINTS OUT THE PARAMS ASSIGNED SO FAR
    def show_params(self):
        for par, val in self.dPkaParams.items():
            print('{:<15s}{:>35s} : {}'.format('DELPHIPKA>', par, val))
        for par, val in self.delphiParams.items():
            print('{:<15s}{:>35s} : {}'.format('DELPHI>', par, val))
        for par, val in self.otherParams.items():
            print('{:<15s}{:>35s} : {}'.format('OTHERS>', par, val))
        return





