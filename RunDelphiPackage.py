import os, sys, re

#############################
# run delphi
#############################
def delphi(arg_exec, arg_pqrFile, arg_paramFile, arg_energies, **arg_params):
    """ run delphi uisng delphi_exec on the input PQR file"""
    
    with open(arg_paramFile, 'w') as fout:
        for k,v in arg_params.items():
            fout.write('{}={}\n'.format(k,v))
            # print('{}={}\n'.format(k,v))
            
        
        # energy terms
        energies = ','.join(arg_energies)
        fout.write('energy({})\n'.format(energies))
        
        # input files
        fout.write('in(modpdb4, file="{}",format=pqr)\n'.format(arg_pqrFile))
        
        # phimap (cube format name derived from the PQR file)
        # will be written in the  directory of the PQR file
        phiName = re.sub("\.pqr$",".phicube", arg_pqrFile)
        phiName = os.path.basename(phiName)
        print('Phimap will be written to {}'.format(os.path.join(os.path.dirname(arg_pqrFile),phiName)))
        fout.write('out(phi,file="{}",format=cube)\n'.format(os.path.join(os.path.dirname(arg_pqrFile),phiName)))
        
        # epsmap (cube format name derived from the PQR file)
        # will be written in the  directory of the PQR file
        epsName = re.sub("\.pqr$",".epscube", arg_pqrFile)
        epsName = os.path.basename(epsName)
        print('Epsmap will be written to {}'.format(os.path.join(os.path.dirname(arg_pqrFile),epsName)))
        fout.write('out(eps,file="{}",format=cube)\n'.format(os.path.join(os.path.dirname(arg_pqrFile),epsName)))
        
        # debmap (cube format name derived from the PQR file)
        # will be written in the  directory of the PQR file
        debName = re.sub("\.pqr$",".debcube", arg_pqrFile)
        debName = os.path.basename(debName)
        print('Debmap will be written to {}'.format(os.path.join(os.path.dirname(arg_pqrFile),debName)))
        fout.write('out(deb,file="{}")\n'.format(os.path.join(os.path.dirname(arg_pqrFile),debName)))
        
		# log file (text format file name derived from the PQR file)
        # will be written in the  directory of the PQR file
        logName = re.sub("\.pqr$",".log", arg_pqrFile)
        logName = os.path.basename(logName)
        print('LOG will be written to {}'.format(os.path.join(os.path.dirname(arg_pqrFile),logName)))
        
        
    # run Delphi
    cmd = '{} {} > {}'.format(arg_exec, arg_paramFile, logName)
    ifSuccess = os.system(cmd)
    return ifSuccess
#END


#############################
# run delphipka
#############################
def delphiPka(arg_exec, arg_pdbFile, arg_crgFile, arg_sizFile, arg_topolFile, arg_paramFile, **arg_params):
	""" run DelphiPkA to protonate the input pdb file based on the crg/siz/topology files provided. 
	The protonation can be pH dependent in which case it will take longer. """

	# create a prm file
	with open(arg_paramFile, 'w') as fout:

		#pdb name
		fout.write('PDB file name                      : {}\n'.format(arg_pdbFile))

		# crg/siz/topology files
		fout.write('Charge parameter                   : {}\n'.format(arg_crgFile))
		fout.write('Radius parameter                   : {}\n'.format(arg_sizFile))
		fout.write('Topology parameter                 : {}\n'.format(arg_topolFile))

		# other parameters
		for k, v in arg_params.items():
			fout.write('{:<35} : {}\n'.format(k, v))

	# run delphiPka 
	logName = re.sub("\.pdb$","_dpka.log", os.path.basename(arg_pdbFile))
	cmd = '{} {} > {}'.format(arg_exec, arg_paramFile, logName)
	print('Running : {}'.format(cmd))

	ifSuccess = os.system(cmd)
	return ifSuccess
#END
