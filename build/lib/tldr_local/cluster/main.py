import os
import sys
import argparse
import tempfile
from io import StringIO
from helpers.molgrep import run_similarity_clustering

def main():
    parser = argparse.ArgumentParser(description="Cluster molecules based on similarity")
    parser.add_argument('input', help="Input file containing SMILES")
    parser.add_argument('output', help="Output file to write clusters")
    