import os
import sys
import argparse


# test tldr_local.bbfilter
def test_bbfilter():
    from tldr_local.bbfilter import main as bbfilter_main
    bb_molecules = [
        'C#CC[C@@H](N)C(N)=O ZINC000078724393',
        'COc1ccc([C@@H]2CC(c3ccc4ccccc4c3)=NN2C(=O)CSc2nnc(-c3ccc(Cl)cc3Cl)n2-c2ccccc2C(C)C)cc1 ZINCNY0000000uYy',
        'Cc1cc(N=c2ccc3c([nH]2)Oc2[nH]c(=Nc4cc(C)nn4-c4ccccc4)ccc2C3c2ccc(C(C)C)cc2)n(-c2ccccc2)n1 ZINCNY0000000uah',
        'Cc1cc(-c2nc3cc(F)ccc3n2C(C)C)ccc1NC(=O)C/C(=C/c1cccc(OCc2ccccc2)c1)c1nc2ccccc2s1 ZINCNY0000000ubd',
        'Cc1nc(Oc2ccc(C(C)(C)c3ccc(Oc4nc(C)nc5sc(-c6ccccc6)cc45)cc3)cc2)c2cc(-c3ccccc3)sc2n1 ZINCNY0000000ukd',
        'COc1cc(/C=C2/CCCc3c2nc2ccccc2c3C(=O)O[C@@H](C)C(=O)N2c3ccccc3Sc3ccc(Cl)cc32)cc(OC)c1OC ZINCNY0000000uwF',
        'O=C(CCCCCCCC(=O)Nc1nc(-c2ccc(Oc3ccccc3)cc2)cs1)Nc1nc(-c2ccc(Oc3ccccc3)cc2)cs1 ZINCNY0000000uBY',
        'Cc1cc(N=c2ccc(C(c3ccc(Cl)c(Cl)c3)c3ccc(=Nc4cc(C)nn4-c4ccccc4)[nH]c3Cl)c(Cl)[nH]2)n(-c2ccccc2)n1 ZINCNY0000000uBL'
    ]
    inclusion_SMARTS = ['C#[CH]']
    exclusion_SMARTS = ['[Si,P,S,B][Cl,F,Br,I]', '[CH,CH2,SH,PH](!@=[O,N,S])', '[SiH,SiH2,SiH3,SiH4,PH,PH2,PH3,SH,SH2,BH,BH2,BH3,BH4]', '[Mg,Li,Zn][#6]', '[C,S,P](!@=[O,N,S])[F,Cl,Br]']
    result = []
    for bb_molecule in bb_molecules:
        smiles,cid = bb_molecule.split(' ')
        result.append(bbfilter_main.bb_filter(smiles, inclusion_SMARTS, exclusion_SMARTS))
    assert result[0] == True, f"Test failed for molecule: {bb_molecules[0]}"
    return "bbfilter test passed"
        

# test_parser = argparse.ArgumentParser(description="tldr_local unittest")
# test_parser.add_argument('--test', type=str, help="Test module to run")
# test_args = test_parser.parse_args()

# if test_args.test:
test_bbfilter()

# import and test each module


        
