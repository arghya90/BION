import os
import re
import sys
import CustomFunctions as CF
import InputParser
import RunParameters
import Reader

sys.path.append("/common/compbio/ARGO/myPython_modules")
import RunDelphiPackage 

def BION():
    
    ################################################################################################
    # INITIALIZE AN INSTANCE OF RunParamters TO PROVIDE FOR A GLOBAL STORAGE OF PARAMTERS
    ################################################################################################
    params = RunParameters.RunParameters()

    ################################################################################################
    # INPUTS
    ################################################################################################
    
    # read input options
    options, args = InputParser.parse_input()

    # --- ion params
    params.add_run_parameter("ionQ", options.ionQ)
    params.add_run_parameter("ionR", options.ionR)
    
    params.add_run_parameter("numMaxPredictions", options.numMaxPredictions)
    params.add_run_parameter("cutOffDist",options.cutOffDist)

    # --- delphi run params
    doGaussian = options.doGaussian
    params.add_run_parameter("delphiExec", options.delphiExec)
    params.add_run_parameter("scale", options.scale)
    params.add_run_parameter("indi", options.indi)
    params.add_run_parameter("exdi", options.exdi)
    params.add_run_parameter("salt", options.salt)
    if doGaussian:
        params.add_run_parameter("doGaussian", 1)
    else:
        params.add_run_parameter("doGaussian", 0)
    

    # --- delphipka run params
    if options.noDelphiPka:
        runDelphiPka = False
    else:
        runDelphiPka = True

    pH = options.pH
    crgFile = options.crgFile
    sizFile = options.sizFile
    topolFile = options.topolFile

    # --- Files
    inStructFile = options.inStructFile


    # show params
    params.show_params()

    ################################################################################################
    # protonate using Delphipka
    ################################################################################################
    if runDelphiPka:
        if pH <= 0 or pH >= 14: 
            print('Energy/pKa calculations will not be performed because no or invalid pH was provided') 
            pH = 7.0
            doEnergy = 'F'
            doPkA = 'F'
            doPQRWithPkA = 'F'
            iFlag = 1
            print("DelphiPkA will be run to protonate at pH = {} without energy/pka calculations".format(pH))
        else:
            doEnergy = 'T'
            doPkA = 'T'
            doPQRWithPkA = 'T'
            iFlag = 2
            print("DelphiPkA will be run to protonate at pH = {} with energy/pka calculations".format(pH))


        otherParams = {'At given pH Value                  ': pH,
                       'Do Protonation                     ': 'T',
                       'Do Energy Calculation              ': doEnergy,
                       'Do pKa\'s  Calculation             ': doPkA,
                       'Output PQR file (With pKa result)  ': doPQRWithPkA,
                       'Salt Concentration                 ': params.salt}
        
        isSuccess = RunDelphiPackage.delphiPka(arg_exec = "/common/compbio/ARGO/DelphiPKA_github/DelphiPka/bin/delphiPKa",
                                              arg_pdbFile = inStructFile,
                                              arg_crgFile = crgFile, 
                                              arg_sizFile = sizFile,
                                              arg_topolFile = topolFile,
                                              arg_paramFile = "run.prm",
                                              **params.dPkaParams, **otherParams)

        if isSuccess != 0:
            sys.exit("DelphiPka finished with errors. Therefore exiting ...")

        ################################################################################################
        # additional check (see if Delphi pka produced the correct PQR file)
        ################################################################################################
        pqrBaseName = re.sub("\.pdb$","_{}.pqr".format(iFlag),os.path.basename(inStructFile))
        pqrFile = os.path.join(os.path.dirname(inStructFile), pqrBaseName)
        try:
            file.exists(pqrFile)

        except:
            FileNotFoundError(pqrFile)
    else:
        # input must be a PQR file
        print("DelphiPka was asked not to run. Hoping for the input file to be in PQR format")
        if CF.isPQRFile(inStructFile):
            pqrFile = inStructFile

        else:
            sys.exit("Input file is not in the PQR format. There is no pooint in continuing .... Exiting therefore")

    ################################################################################################
    # Use the PQR file atom records 
    # to generate a list of instances of Class Atom for each atom read
    ################################################################################################
    atomList = Reader.read_PQR_into_atomList(pqrFile)
    print('A total of {} atoms were created after propotonation using DelphiPka'.format(len(atomList)))

    ################################################################################################
    # run Delphi
    ################################################################################################
    otherParams = {'indi' : params.indi, 
                   'exdi' : params.exdi, 
                   'scale' :  params.scale, 
                   'gaussian' : params.doGaussian,
                   'ionrad' : params.ionR,
                   'salt' : params.salt}

    isSuccess = RunDelphiPackage.delphi(arg_exec=params.delphiExec, 
                  arg_pqrFile = pqrFile, 
                  arg_paramFile = "param.txt", 
                  arg_energies = ['g'], 
                  **params.delphiParams, **otherParams)
    if isSuccess != 0:
        sys.exit("Delphi finished with errors. Therefore exiting ...")

    
    ################################################################################################
    # read the cube format Phimap, EpsMap and DebMap
    ################################################################################################
    epsName, phiName, debName = [ re.sub("\.pqr", extn, os.path.basename(pqrFile)) for extn in [".epscube",".phicube", ".debcube"] ]
    
    epsFile = os.path.join(os.path.dirname(inStructFile), epsName)
    try:
        file.exists(epsFile)
    except:
        FileNotFoundError(epsFile)
    params.add_run_parameter("epsFile", epsFile)
    
    potFile = os.path.join(os.path.dirname(inStructFile), phiName)
    try:
        file.exists(potFile)
    except:
        FileNotFoundError(potFile)
    params.add_run_parameter("potFile", potFile)

    debFile = os.path.join(os.path.dirname(inStructFile), debName)
    try:
        file.exists(debFile)
    except:
        FileNotFoundError(debFile)
    params.add_run_parameter("debFile", debFile)

    # READ THE CONTENTS OF THE PHIMAP AND EPSMAP 
    epsMap = Reader.read_cube_file(params.epsFile)
    params.add_run_parameter("epsMap", epsMap)

    phiMap = Reader.read_cube_file(params.potFile)
    params.add_run_parameter("phiMap", phiMap)

    debMap = Reader.read_cube_file(params.debFile)
    params.add_run_parameter("debMap", debMap)
    
    ################################################################################################
    # FETCH COMPUTATIONAL BOX'S INFORMATION (USING PHIMAP CUBE FILE; CAN ALSO USE EPSMAP FILE AS LONG AS IT IS FROM THE SAME RUN)
    ################################################################################################
    gridOrigin, gridSizes = Reader.fetch_box_geometry(potFile)
    params.add_run_parameter("origin", gridOrigin)
    params.add_run_parameter("gsize", gridSizes)

    gridRes = [1.0/params.scale for i in range(3)]
    params.add_run_parameter("resolution", gridRes)
    
    # FETCH A PRIORITY QUEUE OF GRID POINTS ORDERED ACCORDING TO ITS POTENTIAL ENERGY.
    potential_PQ = CF.generate_potential_PQ(arg_phiMap = params.phiMap, 
                                            arg_epsMap = params.epsMap, 
                                            arg_debMap = params.debMap, 
                                            arg_ionQ = params.ionQ, 
                                            arg_ionR = params.ionR, 
                                            arg_gridSizes = params.gsize, 
                                            arg_exdi = params.exdi)

    # GLORIFY THE CANDIDATES (at most numMaxPredictions)
    gridSites = CF.glorify_candidates(arg_gridRes=params.resolution, 
                                      arg_origin = params.origin, 
                                      arg_minDist=6,
                                      arg_numCandidates=params.numMaxPredictions,
                                      arg_PQ=potential_PQ, 
                                      arg_structAtoms = atomList,
                                      arg_ionRadius = params.ionR,
                                      arg_cutOffDist = params.cutOffDist)
    
    # RETURN
    return

#END

if __name__=="__main__":
    BION()


#END OF PROGRAM
