# TLDR-local
 Running TLDR modules locally
## Installation 
1. Setup Python virtual environment and install rdkit

Option 1: use python venv
 ```
    python3 -m venv <venv_name>
    pip install rdkit-pypi 
```
Option 2: conda env
```
    conda create -c conda-forge -n <venv_name> rdkit
    conda install -c conda-forge rdkit # for conda
```

2. Download `tldr-local` repo
```
 git clone https://github.com/khtang17/tldr-local.git
 cd tldr-local
python setup.py sdist 
 pip install dist/tldr_local-<version>.tar.gz
 ```
3. Check installation
 ```
 python -c 'import tldr_local'
 ```

## Usage
### Run existing scripts

It is currently support `bb_filter`, `strain` and `fine_tranche`

Example:
 ```
bb_filer --help
usage: bb_filter [-h] [--hac-range HAC_RANGE] bb_file output inclusion_smarts exclusion_smarts

Filter building blocks based on HAC, inclusion and exclusion rules

positional arguments:
  bb_file               Input file containing building blocks in SMILES format
  output                Output file to write filtered building blocks
  inclusion_smarts      Inclusion SMARTS file
  exclusion_smarts      Exclusion SMARTS file

options:
  -h, --help            show this help message and exit
  --hac-range HAC_RANGE
                        HAC range to filter by
```
 ### Import the module and functions in your own scripts

Example:
 ```
>>> from tldr_local.fine_tranche.main import fine_tranche
>>> fine_tranche('Cc1c(-c2ccccc2)nn(C)c1NS(=O)(=O)C1CCCCC1')
('P340', 'H23')
```



 


