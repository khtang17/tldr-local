This software uses the Torsion Library (DOI: acs.jcim.5b00522) jointly developed by
the University of Hamburg, Center for Bioinformatics, Hamburg,
Germany and F. Hoffmann-La-Roche Ltd., Basel, Switzerland.

We make revised method for turning the frequency information 
for dihedral angles in a Torsion Library XML file into
energy information. The purpose is to then use the new XML file, which
preserves the original library's structure, to estimate the energy of a
dihedral angle in a given conformation for a molecule of interest.

To run the code, you need to install RDKit: https://www.rdkit.org/docs/Install.html
You can use this software like this in the commandline:
  $ Python3.7 Torsion_Strain.py {.mol2 or .db2 file name}

OUTPUT: This software outputs a .csv file. Each row corresponds to a different
molecular scructure in the .mol2 or .db2 file input. The first column is the
name of the structure, either the ZINC ID for each molecule in a .mol2 file
or the number of the pose for a .db2 file (which should contain multiple
poses for the same molecule). Following the name, we calculate the total energy
estimate and lower and upper bounds of the 95% confidence interval for the compound.
The remaining columns in the outputted .csv file contain information 
for each torsion pattern in that row's structure. For each torsion pattern, we present: 
the index of the torsion pattern, 
the energy estimate and 95% confidence interval (lower and upper bounds) for just that torsion pattern, 
the atom indices defining the torsion pattern, 
the dihedral angle in degrees, 
the SMARTS pattern for the torsion rule in the Torsion Library, 
the type of torsion rule (specific or general), 
the type of energy calculation method (exact or approximate), 
and a Boolean indicating whether or not the torsion pattern has an approximate method 
and a dihedral angle not within the second tolerance of any peak (whether it is "flagged").
Each torsion pattern is sorted by the energy estimate from largest to smallest. 
If torsion pattern is flagged, the flagged pattern is given an energy estimate of 100, 
and the total energy estimate is labeled as -1.