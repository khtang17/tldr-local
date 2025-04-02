# TLDR-local
 Running TLDR modules locally
## Installation 
1. Setup Python virtual environment
 ```
    pip install rdkit-pypi # or pip
    or
    conda install conda-forge::rdkit # for conda
```

2. Download `tldr-local` repo
```
 git clone https://github.com/khtang17/tldr-local.git
 cd tldr-local
 pip install -e .
 ```
3. Check installation
 ```
 python -m 'import tldr-local
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
 >>> from tldr_local import fine_tranche
>>> from tldr_local.fine_tranche.main import scale_logp_value
>>> scale_logp_value(0.9)
90
```



 


