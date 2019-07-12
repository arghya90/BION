""" Module with functions required for reading operations on different kinds of
files that are being used throught the program"""

import numpy as nmp
import Atom
import PriorityQueue
import SGrid
import CustomFunctions as CF

def read_PQR_into_atomList(arg_pqrFile):
    """ Read lines from a PQR file and creates instances of Class Atom
    for each of them. Make a list of these instances and return it. """
    atomSerial = 0
    atomList = []
    with open(arg_pqrFile, 'r') as fPQR:
        for line in fPQR:
            line = line.strip()
            if line.startswith("ATOM") or line.startswith("HETATM"):
                atomSerial += 1
                atomArgs = {
                            'header'      : line[0:6].strip(),
                            'index'       : int(line[6:11]),
                            'atomName'    : line[12:16].strip(),
                            'resName'     : line[17:20].strip(),
                            'x'           : float(line[30:38]),
                            'y'           : float(line[38:46]),
                            'z'           : float(line[46:54]),
                            'charge'      : float(line[54:62]),
                            'radius'      : float(line[62:69])
                            }

                newAtom = Atom.Atom(arg_serial = atomSerial, **atomArgs)
                atomList.append(newAtom)
    return atomList
#END

def fetch_box_geometry(arg_cubeFile):
	""" Read the cube format file outout by Delphi and return the value of the grid size and origin of the computational box derived from it. Check to see if the origin is in Angsrtom or Bohr units (typically Bohr units when it is a Delphi output cube file). Report everything in Angstroms. The resolution of the box can be obtained from the scale value used for the delphi run."""

	Bohr2Ang = 0.529177
	with open(arg_cubeFile, "r") as fCube:
		cubeData = fCube.readlines()
	
	# read the number of atoms and origin of the voxel (line 3)
	numAtoms = int(cubeData[2].split()[0])	
	originX, originY, originZ = [Bohr2Ang * float(e) for e in cubeData[2].split()[1:4]]

	# read voxel information
	iLine = 3
	while (iLine < 6): 
		lineInfo = cubeData[iLine].split()
		if iLine == 3:
			gsizeX = int(lineInfo[0])
		elif iLine == 4:
			gsizeY = int(lineInfo[0])
		elif iLine == 5:
			gsizeZ = int(lineInfo[0])
		iLine += 1
	
	gridOrigin = [originX, originY, originZ]
	gridSizes = [gsizeX, gsizeY, gsizeZ]

	print("Origin: {:.2f} {:.2f} {:.2f}".format(*gridOrigin))
	print("Grid sizes: X = {}, Y = {}, Z = {}".format(*gridSizes))

	return [gridOrigin, gridSizes]
#END

def read_cube_file(arg_cubeFile):
	""" Read the cube-format file output by Delphi when asked to report the potential and dielectric
	value at each of the grid and returna linear array with values placed in the order of z-y-x."""

	with open(arg_cubeFile, "r") as fCube:
		cubeData = fCube.readlines()

	# read the number of atoms (line 3)
	numAtoms = int(cubeData[2].split()[0])	

	# gridSize information from the next 3 lines
	iLine = 3
	while (iLine < 6): 
		lineInfo = cubeData[iLine].split()
		if iLine == 3:
			gsizeX = int(lineInfo[0])
		elif iLine == 4:
			gsizeY = int(lineInfo[0])
		elif iLine == 5:
			gsizeZ = int(lineInfo[0])

		iLine += 1

	# initialize a return array of size gsizeX x gsizeY x gsizeZ
	retArray = nmp.zeros(gsizeX * gsizeY * gsizeZ)
	print("Initialized array of size = {} ({}, {}, {})".format(gsizeX * gsizeY * gsizeZ, gsizeX , gsizeY , gsizeZ))

	# populate the array with value of the quantities stored in the file
	idx = 0
	jdx = 0
	kdx = 0
	for line in cubeData[(6 + int(numAtoms)):]:
		lineInfo = line.strip().split()
		for val in lineInfo:
			vdx = kdx + (jdx * gsizeY) + (idx * gsizeX * gsizeX)
			retArray[vdx] = float(val)
			#print("{index} : {value}".format(index = vdx, value = val))
			#if vdx == 100:
			#	return	

			kdx += 1
			if kdx == gsizeZ:
				kdx = 0
				jdx += 1
				if jdx == gsizeY:
					jdx = 0
					idx += 1

	return retArray
#END
			

def read_binary_files(arg_runParamInstance):
    """ *** OUT OF USE ***
	
	In a custom defined way, read the custom-desgined epsmap, phimap and debmap files in the binary format
    and return 3-tuple of grid-resolution, 3-tuple of grid-origin, a dict indicating 
    which atoms span which grids and a priority queue of grid points ordered according to the potential.
    
    Use the instance of Class RunParameters to get access to other parameters defined since the program has run.
	
	*** OUT OF USE ***
	*** OUT OF USE ***	"""

    debFile = arg_runParamInstance.debFile
    potFile = arg_runParamInstance.potFile
    epsFile = arg_runParamInstance.epsFile
    endianness = arg_runParamInstance.endianness

    # read the first block of debmap
    # which contains the information about which atoms
    # span a certain grid point.
    fDebBin = open(debFile,'rb')

    # GRID SIZES (First of the two occureneces of grid-size in debmap)
    # the two redundant occurences of grid-size data
    # come from the two blocks of the deb map.
    # the first bloack is simply a data about which atoms span a grid point
    # the second is a binary map consistent with the styles used for 
    # potential and epsioon maps.
    fmt = '{}i i i'.format(endianness)
    gsizeX, gsizeY, gsizeZ = CF.decode(fDebBin, fmt)

    # for each block of data read (corresponding to each grid point) in the
    # format: i j k (number of atoms spanning grid i,j,k) [indices of atoms spanning grid i,j,k]
    # store the list of atom indices as values corresponding to a key defined by 
    # (i-1)*gsizeY*gsizeZ + (j-1)*gsizeY + (k - 1)
    atomsOnGrid = dict()

    for gridIdx in range(nmp.product((gsizeX, gsizeY, gsizeZ))):
        i, j, k = CF.decode(fDebBin, '{}i i i'.format(endianness))
        gridKey = (k-1)*gsizeY*gsizeZ + (j-1)*gsizeZ + (i-1)

        atomsOnGrid[gridKey] = []

        # number of atoms spanning this grid (unsigned int)
        numAtomsOnGrid = CF.decode(fDebBin, '{}I'.format(endianness))[0]
        
        # add the index of each atom spanning this grid to its value in the dict
        for iv in range(numAtomsOnGrid):
            atIdx = CF.decode(fDebBin, '{}i'.format(endianness))[0]
            atomsOnGrid[gridKey].append(atIdx)
            

    # Gather potential, epsilon and debmap-value in store them in HEAPS
    potential_PQ = PriorityQueue.PriorityQueue()

    # REAL STUFF WHERE INFORMATION IS READ THE POTENTIAL IS HEAPIFIED
    fPotBin = open(potFile,'rb')
    fEpsBin = open(epsFile,'rb')

    # GRID SIZES
    fmt = '{}i i i'.format(endianness)
    gsizePhiMap = CF.decode(fPotBin, fmt)
    gsizeEpsMap = CF.decode(fEpsBin, fmt)
    gsizeDebMap = CF.decode(fDebBin, fmt)   # the second occurence of grid size data in debmap
    try:
        assert(gsizePhiMap == gsizeEpsMap == gsizeDebMap)
        print( gsizePhiMap )
    except:
        print("Exiting because Grid sizes read from the phimap ({} x {} x {}), epsilon-map ({} x {} x {}) and debmap ({} x {} x {}) do not match".format(*gsizePhiMap, *gsizeEpsMap, *gsizeDebMap))


    # ORIGIN
    fmt = '{}d d d'.format(endianness)
    originPhiMap = CF.decode(fPotBin, fmt)
    originEpsMap = CF.decode(fEpsBin, fmt)
    originDebMap = CF.decode(fDebBin, fmt)
    try:
        assert(originPhiMap == originEpsMap == originDebMap)
    except:
        print("Exiting because the ORIGIN read from the phimap ({} x {} x {}), epsilon-map ({} x {} x {}) and debmap ({} x {} x {}) do not match".format(*originPhiMap, *originEpsMap, *originDebMap))

    # GRID RESOLUTION/SPACING
    fmt = '{}d d d'.format(endianness)
    gridResPhiMap = CF.decode(fPotBin, fmt)
    gridResEpsMap = CF.decode(fEpsBin, fmt)
    gridResDebMap = CF.decode(fDebBin, fmt)
    try:
        assert(gridResPhiMap == gridResEpsMap == gridResDebMap)
    except:
        print("Exiting because the GRID-RESOLUTION read from the phimap ({} x {} x {}), epsilon-map ({} x {} x {}) and debmap ({} x {} x {}) do not match".format(*gridResPhiMap, *gridResEpsMap, *gridResDebMap))


    fmtR = '{}iiid'.format(endianness)
    fmtB = '{}iiib'.format(endianness)

    for ig in range(nmp.product(gsizePhiMap)):
        idx, jdx, kdx, pot =  CF.decode(fPotBin, fmtR)
        eps =  CF.decode(fEpsBin, fmtR)[-1]
        deb =  CF.decode(fDebBin, fmtB)[-1]

        # desolvation penalty
        desolv = 561 * 0.5 * arg_runParamInstance.ionQ * arg_runParamInstance.ionQ * ((1/eps) - 1/arg_runParamInstance.exdi)/ arg_runParamInstance.ionR

        # push sites that have negative potential energy and 
        # lie outisde the vdw surface
        if pot * arg_runParamInstance.ionQ < 0 and deb == 1:
            gridKey = (idx-1)*gsizeY*gsizeZ + (jdx-1)*gsizeZ + (kdx-1)
            otherPropDict = {'eps':eps, 'deb': deb, 'spanningAtoms': atomsOnGrid[gridKey]}
            potential_PQ.push(SGrid.SGrid(idx, jdx, kdx, (arg_runParamInstance.ionQ * pot + desolv), **otherPropDict))

    return (gridResPhiMap, originPhiMap, atomsOnGrid, potential_PQ)
#END

   
if __name__ == "__main__":
	sampleCubeFile = "phimap.cube"
	fetch_box_geometry(sampleCubeFile)
	cubeData = read_cube_file(sampleCubeFile)
