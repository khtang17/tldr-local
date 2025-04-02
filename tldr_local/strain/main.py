import sys
import csv
import argparse
import logging
import multiprocessing as mp
from STRAIN_FILTER.TL_Functions_new import db2MolSupplier, \
                                         Mol2MolSupplier, \
                                        TL_lookup
from STRAIN_FILTER.db2_to_mol2 import db2_file_like
from tqdm import tqdm



def calculate_strain_energy(mol_tupple):
    # mol_tupple is a tuple of (name, mol)
    name = mol_tupple[0]
    mol = mol_tupple[1]
    try:
        M = TL_lookup(mol) #Create a TP_list function
        mol_est = []
        mol_est += M.sum(0)
        mol_info = M.get_TPs() #The molecule's information
        bond_info = [item for sublist in sorted(mol_info, key = lambda l:l[1], reverse=True) for item in sublist]
        return [name] + mol_est + bond_info
    except:
        return [name] + ["NA"]


def worker(mol):
    return calculate_strain_energy(mol)

def main():
    parser = argparse.ArgumentParser(description="Calculate torsion strain of 3D molecules in mol2 or db2")
    parser.add_argument('input', help='Input file in mol2 or db2 format')
    parser.add_argument('output', help='Output file in csv format')
    parser.add_argument('--processes', type=int, default=1, help='Number of processes to use')
    args = parser.parse_args()

    input_file = args.input
    # read input files into a list of names and dictionary of molecules
    if input_file.endswith('.db2') or input_file.endswith('.db2.gz'):
        # Creates a string buffer object called "file" that looks like a .mol2 file
        input_file = db2_file_like(input_file)
        mol_list = db2MolSupplier(input_file)
    elif input_file.endswith('.mol2'):
        mol_list = Mol2MolSupplier(input_file)
    else:
        raise ValueError("Input file must be in mol2 or db2 format")
    
    logging.info(f"{len(mol_list)} molecules finished reading. Calculating strain energy...")

    with open(args.output, 'w', newline='') as out, mp.Pool(processes=args.processes) as pool:
        writer = csv.writer(out)
        for res in tqdm(pool.imap(worker, mol_list, chunksize=256), unit_scale=True, mininterval=1.0):
            writer.writerow(res)

if __name__ == '__main__':
    main()
