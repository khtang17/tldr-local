import os
import sys
import argparse
import multiprocessing as mp
import tqdm

from rdkit.Chem import (
    MolFromSmiles,
    MolToSmiles,
    AddHs
)

from rdkit.Chem.Descriptors import (
    MolLogP,
    MolWt
)
from rdkit.Chem.rdChemReactions import ReactionFromSmarts




def main():
    parser = argparse.ArgumentParser(description="Run bioisosteric expansion on a set of SMILES strings.")
    parser.add_argument('input', type=str, required=True, help="Input file with SMILES strings")
    parser.add_argument('--gen1', type=int, choices=[0,1,2,3,4], required=True, help="Conservative bioisosteric expansion iteration")
    parser.add_argument('--gen2', type=int, choices=[0,1,2,3,4], required=True, help="Medium bioisosteric expansion iteration")
    parser.add_argument('--gen3', type=int, choices=[0,1,2,3,4], required=True, help="Aggressive bioisosteric expansion iteration")


    parser.add_argument('--output', type=str, required=True, help="Output file for bioisosteric SMILES")

