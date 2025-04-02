import sys
import argparse
import multiprocessing as mp
from rdkit.Chem import MolFromSmiles
from rdkit.Chem.Descriptors import MolLogP
from rdkit.Chem.SaltRemover import SaltRemover
from tqdm import tqdm

remover = SaltRemover()

def scale_logp_value(logp):
    if logp < -9.0:
        logp = -9.0
    elif logp > 9.0:
        logp = 9.0
    if logp < 0 or logp >= 5.0:
        logp = 100*int(logp)
    else:
        logp = 10*int(10*logp)
    return logp

def fine_tranche(smiles):
    smiles = smiles.strip()
    try:
        mol = MolFromSmiles(smiles)
        mol = remover.StripMol(mol)
        if mol is None:
            raise ValueError("Invalid SMILES")
        logp = MolLogP(mol)
        num_heavy_atoms = mol.GetNumHeavyAtoms()
        if num_heavy_atoms > 99:
            num_heavy_atoms = 99
        sign = 'M' if logp < 0 else 'P' # M for negative/minus logP and P for positive logP
        logp = scale_logp_value(logp)
        logp_str = f'{sign}{abs(logp):03d}'
        num_heavy_atoms_str = f'H{num_heavy_atoms:02d}'
        return logp_str, num_heavy_atoms_str
    except Exception as e:
        print(f"Error: {e}")
        return
    
def worker(line):
    smiles, cid = line.strip().split()[:2]
    logp, num_heavy_atoms = fine_tranche(smiles)
    return f"{smiles} {cid} {num_heavy_atoms}{logp}\n"

    
def main():
    parser = argparse.ArgumentParser(description="Generate fine-grained tranches for a list of SMILES")
    parser.add_argument('input', help="Input file containing SMILES")
    parser.add_argument('output', help="Output file to write fine-grained tranches")
    parser.add_argument('--processes', type=int, default=1, help="Number of processes to use")
    args = parser.parse_args()
    if args.processes:
        processes = args.processes
    else:
        processes = 1

    with open(args.input) as f, open(args.output, 'w') as out, mp.Pool(processes=processes) as pool:
        for res in tqdm(pool.imap(worker,f, chunksize=256), unit_scale=True, mininterval=1.0):
            if res:
                out.write(res)
    print("Done")


if __name__ == "__main__":
    main()