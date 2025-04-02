import sys
import argparse
import multiprocessing as mp
from rdkit import Chem as C
from tqdm import tqdm

#################################################
# Determine if the current compound contains, and is allowed to contain, multiple symmetric RC
def allowed_reactive_symm_check(compound_mol):
        symm = C.RemoveHs(compound_mol).GetSubstructMatches((C.RemoveHs(compound_mol)), useChirality =True, uniquify =False)
        return len(symm) == 2 or len(symm) == 4

def inclusion_rules_check(mol, inclusion_smarts):
    for pattern in inclusion_smarts:
        inSMARTS_pattern_mol = C.MolFromSmarts(pattern)
        if mol.HasSubstructMatch(inSMARTS_pattern_mol, useChirality=True):
            num_substructs = mol.GetSubstructMatches(inSMARTS_pattern_mol, uniquify=True)
            if len(num_substructs) == 1:
                return True
            else:
                return False

def exclusion_rules_check(mol, exclusion_smarts):
    for pattern in exclusion_smarts:
        exSMARTS_pattern_mol = C.MolFromSmarts(pattern)
        if mol.HasSubstructMatch(exSMARTS_pattern_mol, useChirality=True):
            return False
        else:
            return True

# Check
def filter_by_hac(bb_mol, hac_range):
    min_hac, max_hac = [int(num.replace('H','').lstrip()) for num in hac_range.split('-')]
    # failed_hac_smiles = open('failed_hac_filter.smi', 'w')
    bb_hac = bb_mol.GetNumHeavyAtoms()
    if bb_hac <= max_hac and bb_hac >= min_hac:
        return True
    else:
        return False



def bb_filter(bb_smiles, inSMARTS, exSMARTS, hac_range=None):
    if hac_range:
        hac_range = hac_range.strip()
    filtered = open('FilteredBB.smi', 'w')
    cmpd_mol = C.MolFromSmiles(bb_smiles)
    if not cmpd_mol:
        raise ValueError("Invalid SMILES")
    
    if filter_by_hac(cmpd_mol, hac_range) and inclusion_rules_check(cmpd_mol, inSMARTS) and exclusion_rules_check(cmpd_mol, exSMARTS):
        return True
    else:
        return False
    
def worker(line, inSMARTS, exSMARTS, hac_range):
    smiles, cid = line.strip().split()[:2]
    bb_filter(smiles, inSMARTS, exSMARTS, hac_range)
    return f"{smiles} {cid}\n"

def main():
    parser = argparse.ArgumentParser(description="Filter building blocks based on HAC, inclusion and exclusion rules")
    parser.add_argument('bb_file', help="Input file containing building blocks in SMILES format")
    parser.add_argument('output', help="Output file to write filtered building blocks")
    parser.add_argument('inclusion_smarts', help="Inclusion SMARTS file")
    parser.add_argument('exclusion_smarts', help="Exclusion SMARTS file")
    parser.add_argument('--hac-range', help="HAC range to filter by")
    args = parser.parse_args()

    infile = open(args.bb_file, 'r')
    input_list = [line for line in infile.readlines()]
    hac_range = args.hac_range
    insmart_f = open(args.inclusion_smarts, 'r')
    inSMARTS = [line.strip('\n') for line in  insmart_f.readlines()]
    exsmart_f = open(args.exclusion_smarts, 'r')
    exSMARTS = [line.strip('\n') for line in exsmart_f.readlines()]

    with open(args.bb_file) as f, open(args.output, 'w') as out, mp.Pool(processes=1) as pool:
        for res in tqdm(pool.imap(worker, input_list, inSMARTS, exSMARTS, hac_range, chunksize=256), unit_scale=True, mininterval=1.0):
            if res:
                out.write(res)
    print("Done")

if __name__ == "__main__":
    main()


