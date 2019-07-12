import optparse 
import re
import sys, os

def parse_input():
    inputParser = optparse.OptionParser()

    # ------- FILE IO --------
    # PDB file (-i)
    inputParser.add_option("-i",
                           action="store",
                           dest  ="inStructFile",
                           type  ="string",
                           help  ="PDB file (.pdb) or PQR file (.pqr)")
    
    # ------- ION PARAMS --------
    # ION CHARGE (-q, --ionq)
    inputParser.add_option("-q","--ionq",
                           action="store",
                           dest  ="ionQ",
                           type  =float,
                           help  ="Charge of the ion type for which predictions are to be made")
    
    # ION RADIUS (-r, --ionr)
    inputParser.add_option("-r","--ionr",
                           action="store",
                           dest  ="ionR",
                           type  =float,
                           help  ="Radius (Å) of the ion type for which predictions are to be made")

    # MAX NUMBER OF PREDICTIONS (-n, --nump)
    inputParser.add_option("-n","--nump",
                           action="store",
                           dest  ="numMaxPredictions",
                           type  =int,
                           default = 10,
                           help  ="Maximum number of predictions for favorable sites to be made")

    # MAX ACCEPTABLE DISTANCE FROM A SOLUTE ATOM (-x)
    inputParser.add_option("-x",
                           action="store",
                           dest  ="cutOffDist",
                           type  =float,
                           default = 6,
                           help  ="Maximum acceptable distance of a predicted site from a solute atom")

    # ------- DELPHIPKA RUN PARAMS ------
    # TO RUN DELPHIPKA?
    inputParser.add_option("--nodpka",
                           action="store_true",
                           dest="noDelphiPka",
                           default=False,
                           help="Should DelphiPka be run? (Need not run if you have a PQR file)")

    # CRG file (-c, --crg)
    inputParser.add_option("-c","--crg",
                           action="store",
                           dest  ="crgFile",
                           type  ="string",
                           help  ="Charge parameter file (For protonation using Delphi PkA run)")
    
    # SIZ file (-s, --siz)
    inputParser.add_option("-s","--siz",
                           action="store",
                           dest  ="sizFile",
                           type  ="string",
                           help  ="Radius parameter file (For protonation using Delphi PkA run)")
    
    # TOPOLOGY file (-t, --topol)
    inputParser.add_option("-t","--topol",
                           action="store",
                           dest  ="topolFile",
                           type  ="string",
                           help  ="Topology file (For protonation using Delphi PkA run)")

    # pH for protonation (--ph)
    inputParser.add_option("--ph",
                           action="store",
                           dest  ="pH",
                           type  ="float",
                           default = -1,   # dummy that will decide if Delphipka needs to run energy/pka calcs or not
                           help  ="Value of pH to use for protonation")
    
    # ------- DELPHI RUN PARAMS --------
    # delphi executable
    inputParser.add_option("--delphiexec",
                           action="store",
                           dest  ="delphiExec",
                           type  ="string",
                           help  ="Path of the Delphi executable")

    # delphi scale
    inputParser.add_option("--scale",
                           action="store",
                           dest  ="scale",
                           type  = float,
                           default = 1,
                           help  ="Value of scale for Delphi run (grids/Å)")

    # delphi internal dielectric
    inputParser.add_option("--indi",
                           action="store",
                           dest  ="indi",
                           type  = float,
                           default = 2,
                           help  ="Value of the internal dielectric constant")

    # delphi external dielectric
    inputParser.add_option("--exdi",
                           action="store",
                           dest  ="exdi",
                           type  = float,
                           default = 80,
                           help  ="Value of the external or solvent dielectric constant")
    
    # delphi salt conc.
    inputParser.add_option("--salt",
                           action="store",
                           dest  ="salt",
                           type  = float,
                           default = 0.15,
                           help  ="Value of the salt concentration (M)")


    # delphi execute gaussian delphi 
    inputParser.add_option("--gaussian",
                           action="store_true",
                           dest  ="doGaussian",
                           default = False,
                           help  ="Run Delphi with Gaussian dielectric distribution (time consuming). By default, the standard 2-dielectric delphi is run.")



    options, args = inputParser.parse_args()
    return (options, args)

if __name__ == "__main__":
    parse_input()

