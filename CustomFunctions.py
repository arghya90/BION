import struct
import numpy as nmp
import PriorityQueue as PQ
import SGrid

def isPQRFile(arg_fileName):
    """ A function that returns true if the ATOM/HETATM lines are in DELPHI'S MODPDB4 PQR format """

    with open(arg_fileName, 'r') as fin:
        for line in fin:
            if line.startswith("ATOM  ") or line.startswith("HETATM"):
                #check the line's XYZQR contents
                x = line[30:38]
                y = line[38:46]
                z = line[46:54]
                q = line[54:62]
                r = line[62:69]
                
                # check float-ness of X
                try:
                    float(x)
                except:
                    return False

                # check float-ness of Y
                try:
                    float(y)
                except:
                    return False

                # check float-ness of Z
                try:
                    float(z)
                except:
                    return False

                # check float-ness of Q
                try:
                    float(q)
                except:
                    return False

                # check float-ness of R
                try:
                    float(r) and float(r) >= 0.0
                except:
                    return False

    # if x y z q r fields are all floats and 'r' is a non-negative float
    return True

#END



def decode(arg_fHandler, arg_fmt):
    """ Decode/unpack a binary buffer based on the input format """
    buff = arg_fHandler.read(struct.calcsize(arg_fmt))
    return struct.unpack(arg_fmt, buff)
#END

def grid_to_coords(arg_origin, arg_gridRes, arg_gridIndex):
    """ For a grid-index descrobed a triplet of integer bound by the number of grid-points,
    return the cartesian coordinates of that grid point based on the origin and grid-resolution
    of the box """

    # checks
    try:
        assert(len(arg_origin) == 3)
    except:
        raise ValueError("Origin should be of length 3")

    try:
        assert(len(arg_gridRes) == 3)
    except:
        raise ValueError("Grid resolution should be of length 3")

    # assuming grid indices are 0-indexed
    coords = [arg_origin[ii] + (arg_gridIndex[ii] * arg_gridRes[ii]) for ii in range(3)]
    return coords
#END

def get_distance(arg_coord1, arg_coord2):
    """ Return the euclidean distance between two 3D euclidean points """
    vecsub = nmp.subtract( arg_coord1, arg_coord2 )
    dist2 = nmp.sum( [ nmp.power(x, 2) for x in vecsub ] )
    return nmp.sqrt(dist2)
#END

def generate_potential_PQ(arg_phiMap, arg_epsMap, arg_debMap, arg_ionQ, arg_ionR, arg_gridSizes, arg_exdi):
    """ Use the linear arrays containing the potential, dielectric and debmap data computed by Delphi to create a potetnial PQ and return it.
    Potential values will be used to determine the rank of each node. More negative, the liklier it is to be the best site for prediction.
    Due to the Gaussian-salt module, the potential will be influenced by a dielectric-dependent desolvation energy value. The selected points must also 
    have a deb value of 1 (outside the molecular surface)."""
    
    # initialize an empty PQ
    potential_PQ = PQ.PriorityQueue()
    
    # fetch grid sizes (should be the same on all the three directions)
    gX, gY, gZ = arg_gridSizes

    # loop through the indices of the phimap/epsmap simultaneously to compute
    # desolvation and total potential
    for idx in range(gX):
        for jdx in range(gY):
            for kdx in range(gZ):
                # linear index
                vdx = kdx + (jdx * gY) + (idx * gX * gX)
                deb = arg_debMap[vdx]

                # only bother with sites that have deb == 1 (outside mol. surface)
                if deb == 1:
                    pot = arg_phiMap[vdx]
                    eps = arg_epsMap[vdx]
                    
					# desolvation penalty
                    desolv = 561 * 0.5 * arg_ionQ * arg_ionQ * ((1/eps) - 1/arg_exdi)/ arg_ionR
                    
                    # potential modified by desolv
                    modPot = arg_ionQ * pot + desolv

                    # push sites that have negative potential energy and 
                    if modPot < 0:
                        otherPropDict = {'eps':eps}
                        potential_PQ.push(SGrid.SGrid(idx, jdx, kdx, modPot, **otherPropDict))

    return potential_PQ
#END

def glorify_candidates(arg_PQ, arg_gridRes, arg_origin, arg_numCandidates, arg_minDist, arg_structAtoms, arg_ionRadius, arg_cutOffDist):
    """ Return arg_numCandidates number of grid-points with:
    1. most favorable potential (more negative)
    2. no steric clashes with the solute atoms (i.e. distance < atomRadius + arg_ionRadius) but also no more than arg_cutOffDist away from any of the solute atoms.
    3. no less than arg_minDist of separation between predictions (if the space b/w them is't occupied by any solute atom)
    4. arg_minDist separation if the space b/w the sites is occupied by at least one(?) solute atom.
    
    """
    
    gloriousGridList = []
    minDist2 = nmp.power(arg_minDist, 2)
    gridRes = arg_gridRes
    
    # start popping the priority queue
    count = 0
    while count < arg_numCandidates:
        newGloriousGrid = arg_PQ.pop()
       
       # assume it is good to add
        isTooClose = False
        hasVdwClash = False

        # fetch its coordinates
        gIndex = newGloriousGrid.get_grid_indices()
        gCoord = grid_to_coords(arg_origin, arg_gridRes, gIndex)
        
        # (TEMPORARILY SUSPENDED) fetch its spanning atoms
        if False:
            gSpanningAtoms = newGloriousGrid.spanningAtoms
            for atIdx in gSpanningAtoms:
                # get their coordinates
                anAtom = arg_structAtoms[atIdx-1]
                atCoord = [ anAtom.x, anAtom.y, anAtom.z ]
                atRadius = anAtom.radius

                dist2Atom = get_distance(atCoord, gCoord)
                # check if its has vdw clash with solute atoms or is far away
                if dist2Atom < (atRadius + arg_ionRadius):
                    hasVdwClash = True
                    break
                elif dist2Atom > arg_cutOffDist:
                    # remove it form the site's list of spanning atoms
                    # because it may be too far away to coordinate with it
                    gSpanningAtoms.remove(atIdx)

        # calculate its distance from other glorious grids in the list
        # reject it if it is less than minDist from any of those grids
        for gGrid in gloriousGridList:
            dist2 = newGloriousGrid.distance2(gGrid, gridRes[0])
            if dist2 < minDist2:
                isTooClose = True
                break
                
        if not (isTooClose or hasVdwClash):
            gloriousGridList.append(newGloriousGrid)
            count += 1
            print("Prediction {:<5}".format(count))
            print("-"*16)
            print("{:<5}> Potential energy : {} kT".format(count, newGloriousGrid.value))
            print("{:<5}> Coordinates      : {:<8.3f}{:<8.3f}{:<8.3f}".format(count, *gCoord))
            print("{:<5}> Grid indices     : {:<8d}{:<8d}{:<8d}".format(count, *gIndex))
            print("{:<5}> Epsilon          : {:<6.2f}".format(count, newGloriousGrid.eps))
            #print("{:<5}> Nearby atoms     : ".format(count), end='')
            #print(('{} '*len(gSpanningAtoms)).format(*gSpanningAtoms))
            print("\n")
            
    return gloriousGridList
#END

