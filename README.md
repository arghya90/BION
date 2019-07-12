# BION
A python-based command line script that predicts the sites for non-specific binding of a certain type of ion around a biomolecule using electrostatics and geometrical factors

## Usage
In the command line, type:
```
python PATH/TO/BION2.py -h
```
to see the command line options.

To run BION2, one needs a PQR or a PDB file. If PQR file is being used, then DelphiPka need not be run (*--nodpka*). Otherwise,  a CRG, SIZ and topology file along with a PDB file must be provided to obtain a protonated PQR formatted file of the molecule. The protonation can be done at some pH value which can be provided using the *--pH* flag. 
