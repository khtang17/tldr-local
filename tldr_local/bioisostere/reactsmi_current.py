#!/usr/bin/env python
from __future__ import print_function
import pdb
import itertools
import sys
import shlex
import subprocess

from rdkit.Chem import (
    MolFromSmiles,
    MolToSmiles,
    AddHs
)

from rdkit.Chem.Descriptors import MolLogP, MolWt
from rdkit.Chem.rdChemReactions import ReactionFromSmarts
try:
    from flask import (
        abort,
        Flask, 
        redirect,
        render_template,
        request,
        Response,
        url_for,
    )
except ImportError:
    Flask = None


DEBUG = True

class Unbuffered(object):
   def __init__(self, stream):
       self.stream = stream
   def write(self, data):
       self.stream.write(data)
       self.stream.flush()
   def writelines(self, datas):
       self.stream.writelines(datas)
       self.stream.flush()
   def __getattr__(self, attr):
       return getattr(self.stream, attr)

def convert_smarts_to_rxn_list(smarts_file):
    rxn_list = []
    for line in smarts_file:
        row_info = line.split(',')
        rxn_id = row_info[0]
        reactants = row_info[1]
        products = row_info[2]
        risk_level = row_info[3]
        rxn = '{0}>>{1}'.format(reactants,products)
        rxn_list.append((rxn,rxn_id))
    return rxn_list

def load_smiles_file(it):
    for line in it:
        #pdb.set_trace()
        #if statement skips the blank lines in the input
        if line.strip():
            smiles, cid = str(line).strip().split()[:2]
            mol = MolFromSmiles(smiles)
            if mol is not None:
                mol.SetProp('_Name', cid)
                yield mol


def reaction_matrix(rxn, *molsLists):
    for mols in itertools.product(*molsLists):
        products = rxn.RunReactants(mols)
        #print("Testing products print")
        #print(dir(products))
        for product in products:
            yield mols, product


def smiles_reaction_matrix(smarts, *sources, **kwargs):
    #pdb.set_trace()
    sep = kwargs.setdefault('sep', ' ')
    molValue = int(kwargs.get('molValue', 900))
    logValue = float(kwargs.get('logValue', 9.0))
    try:
        reaction = ReactionFromSmarts(smarts)
    except ValueError:
        #pdb.set_trace()
        with open('invalid_smarts.txt', 'a') as invalid_smarts:
            invalid_smarts.write('{0}\n'.format(smarts))
        return
    smilesLists = [load_smiles_file(source) for source in sources]
    #print('About to run reaction_matrix')
    print(smarts)
    debug_file.write('{0}\n'.format(smarts))
    products = reaction_matrix(reaction, *smilesLists)
    #print('outside of loop')
    for reactants, product in products:
        #print('do we have products?')
        #print("Reactants")
        #print(reactants)
        #print(products)
        cids = [r.GetProp("_Name") for r in reactants]
        product_id = '.'.join(cids)
        for mol in product:
            smiles = MolToSmiles(mol, isomericSmiles=True)
            mol.UpdatePropertyCache(strict=False)
            mh = AddHs(mol, addCoords=True)
            mwt = MolWt(mol)
            if mwt <= molValue:
                logp = MolLogP(mol)
                if logp < logValue:
                    yield sep.join((smiles, product_id, str(mwt), str(logp)))+"\n"

def generate_products_at_level(current_result_file, current_rxn_list, smiles_files):
    fs = []
    rxn_list = current_rxn_list
    result_file = current_result_file
    for element in rxn_list:
        #pdb.set_trace()
        for smiles_file in smiles_files:
            debug_file.write('Opening... {0}\n'.format(smiles_file))
            fs.append(open(smiles_file))
        #print('Reaction, ReactionID')
        #print('({0},{1})'.format(element[0],element[1]))
        #print(element[1])
        rxn = element[0]
        debug_file.write('Rxn ID: {0} \n'.format(element[1]))
        #print('Reaction: {0}'.format(rxn))
        for line in smiles_reaction_matrix(rxn, *fs):
            #print('Line: {0}'.format(line))
            #Line will be None if reaction fails
            if line is not None:
                result_line = line.split()
                #append rxnID for bioisosteric replacement just performed to the resulting molecule/product ID
                result_line[1] = '{0}_{1}'.format(result_line[1],element[1])
                #print(result_line)
                result_file.write('{0} {1} {2} {3} \n'.format(result_line[0],result_line[1],result_line[2], result_line[3]))
                print('{0} {1}'.format(result_line[0],result_line[1]))
                #result_file.write('{0} {1}\n'.format(result_line[0],result_line[1]))
            else:
                with open('invalid_smarts.txt', 'a') as invalid_smarts:
                    print('Invalid smarts: ')
                    invalid_smarts.write(' Rxn ID: {0} \n'.format(element[1]))
        for f in fs:
            if f and not f.closed:
                f.close()
        fs = []
    pass

def concatenate_smiles_files_to_result_file(smiles_file_list, result_file):
    #assumes outfile is already open but infile is not
    outfile = result_file
    for fname in smiles_file_list:
        with open(fname) as infile:
            for line in infile:
                outfile.write(line)
    pass

def unique_sort_on_column(result_file, sorted_result_file, column_number):
    print('sort -u -k{0},{0} {1}'.format(column_number,result_file))
    #cmd = shlex.split('sort -u -k{0},{0} {1} > {2}'.format(column_number,result_file,sorted_result_file))
    #print(cmd)
    #args = ['/bin/bash', '-c']
    #for comp in cmd:
    #    args.append(comp)
    #print(args)
    process = subprocess.Popen('sort -u -k{0},{0} {1} > {2}'.format(column_number,result_file,sorted_result_file), shell=True, close_fds=True)
    process.wait()
    #subprocess.Popen(cmd, shell=True, executable='/bin/bash')
    pass

if Flask is not None:
    app = Flask(__name__)
    app.config['DEBUG'] = DEBUG
    app.config['UPLOAD_FOLDER'] = '/tmp'
    
    
    @app.route('/')
    @app.route('/index', methods=['GET'])
    def index():
        smarts1 = request.args.get('smarts1', None)
        smarts2 = request.args.get('smarts2', None)
        smarts3 = request.args.get('smarts3', None)
        product = request.args.get('product', None)
        name = request.args.get('name', None)
        short_name = request.args.get('short_name', None)
        return render_template('index.html', smarts1=smarts1, smarts2=smarts2, smarts3=smarts3, product=product, name=name, short_name=short_name)
    
    
    @app.route('/react', methods=['GET', 'POST'])
    def react():
        product_smarts = request.form.get('product')
        smarts1 = request.form.get('smarts1', None)
        smarts2 = request.form.get('smarts2', None)
        smarts3 = request.form.get('smarts3', None)

        molValue = request.form.get('molW')
        logValue = request.form.get('logP')

        if smarts3:
            reactants = '.'.join((smarts1, smarts2, smarts3))
        else :
            reactants = '.'.join((smarts1, smarts2))
        reaction = '{}>>{}'.format(reactants, product_smarts)

        smiles1 = request.files.get('reactant1')
        smiles2 = request.files.get('reactant2')
        smiles3 = request.files.get('reactant3')

        if reaction and smiles1 and smiles2:
            reaction = str(reaction)
            smiles1 = list(smiles1)
            smiles2 = list(smiles2)
            if smiles3:
                smiles3 = list(smiles3)
                products = smiles_reaction_matrix(reaction, smiles1, smiles2, smiles3, molValue=molValue,
                                                  logValue=logValue)
            else:
                products = smiles_reaction_matrix(reaction, smiles1, smiles2, molValue=molValue, logValue=logValue)
            return Response(products, mimetype='chemical/x-daylight-smiles')
        else:
            return redirect(url_for('index'))
else:
    app = None


if __name__ == '__main__':
    if '--serve' in sys.argv:
        if app is None:
            print("The flask package is missing. Cannot run reaction server.", file=sys.stderr)
            print("Try `pip install flask`", file=sys.stderr)
            sys.exit(-1)
        else:
            sys.exit(app.run(host='0.0.0.0', port=8081, debug=DEBUG))
    else:
        sys.stdout = Unbuffered(sys.stdout)
        if len(sys.argv) < 8:
            print('python reactsmi_current.py conservative_table <# conservative> medium_table <# medium> aggressive_table <# aggressive> smiles')
            exit()
        #supports multiple smiles files
        #pdb.set_trace()
        smiles_files = sys.argv[7:]
        
        try:
            #combine user-provided SMILES into a single file for processing
            combined_original_smiles = open('gen0.ism', 'w+')
            concatenate_smiles_files_to_result_file(smiles_files, combined_original_smiles)
            combined_original_smiles.close()
        except:
            print('Could not open gen0.ism for writing')

        debug_file = open('debug.txt','w')
       
        #get info from command line args
        conservative_num_iterations = int(sys.argv[2])
        medium_num_iterations = int(sys.argv[4])
        aggressive_num_iterations = int(sys.argv[6])
        conservative_num_iterations_orig = conservative_num_iterations
        medium_num_iterations_orig = medium_num_iterations
        aggressive_num_iterations_orig = aggressive_num_iterations

        conservative_smarts_file = open(sys.argv[1], 'r')
        medium_smarts_file = open(sys.argv[3], 'r')
        aggressive_smarts_file = open(sys.argv[5], 'r')
      
        #generate lists of rxns
        conservative_rxn_list = convert_smarts_to_rxn_list(conservative_smarts_file)
        medium_rxn_list = convert_smarts_to_rxn_list(medium_smarts_file)
        aggressive_rxn_list = convert_smarts_to_rxn_list(aggressive_smarts_file)

        conservative_smarts_file.close()
        medium_smarts_file.close()
        aggressive_smarts_file.close()

        subprocess.Popen('touch gen0_conservative.ism', shell=True)
        subprocess.Popen('touch gen0_medium.ism', shell=True)
        subprocess.Popen('touch gen0_aggressive.ism', shell=True)

        print(conservative_rxn_list)
        print(medium_rxn_list)
        print(aggressive_rxn_list)
        fs = []
        try:
            #pdb.set_trace()
            #for smiles_file in smiles_files:
            #    fs.append(open(smiles_file))
            #for (rxn,rxn_id) in rxn_list:
            counter = 1
            while conservative_num_iterations > 0 or medium_num_iterations > 0 or aggressive_num_iterations > 0:
                output_file_string = 'gen{0}.ism'.format(counter)
                
                conservative_result_file = open('gen{0}_conservative.ism'.format(counter), 'w')
                if conservative_num_iterations != 0:
                    generate_products_at_level(conservative_result_file,conservative_rxn_list,smiles_files)
                    conservative_result_file.close()
                    #conservative_num_iterations = conservative_num_iterations - 1
                
                medium_result_file = open('gen{0}_medium.ism'.format(counter), 'w')
                if medium_num_iterations != 0:
                    generate_products_at_level(medium_result_file,medium_rxn_list,smiles_files)
                    medium_result_file.close()
                    #medium_num_iterations = medium_num_iterations - 1
                
                aggressive_result_file = open('gen{0}_aggressive.ism'.format(counter), 'w')
                if aggressive_num_iterations != 0:
                    generate_products_at_level(aggressive_result_file,aggressive_rxn_list,smiles_files)
                    aggressive_result_file.close()
                    #aggressive_num_iterations = aggressive_num_iterations - 1
                
                current_gen_list = ['gen{0}_conservative.ism'.format(counter), 'gen{0}_medium.ism'.format(counter), 'gen{0}_aggressive.ism'.format(counter)]
                if conservative_num_iterations == 0:
                    current_gen_list[0] = 'gen{0}_conservative.ism'.format(conservative_num_iterations_orig)
                if medium_num_iterations == 0:
                    current_gen_list[1] = 'gen{0}_medium.ism'.format(medium_num_iterations_orig)
                if aggressive_num_iterations == 0:
                    current_gen_list[2] = 'gen{0}_aggressive.ism'.format(aggressive_num_iterations_orig)
                print(current_gen_list)

                if conservative_num_iterations > 0:
                    conservative_num_iterations = conservative_num_iterations - 1
                if medium_num_iterations > 0:
                    medium_num_iterations = medium_num_iterations - 1
                if aggressive_num_iterations > 0:
                    aggressive_num_iterations = aggressive_num_iterations - 1

                #only one result file open at a time
                result_file = open(output_file_string, 'w')
                concatenate_smiles_files_to_result_file(current_gen_list, result_file)
                result_file.close()
                sorted_result_file = open('sorted_{0}'.format(output_file_string), 'w')
                #unique_sort_on_column(result_file, sorted_result_file, 1)
                unique_sort_on_column(output_file_string, 'sorted_{0}'.format(output_file_string), 1)
                #result_file.close()
                sorted_result_file.close()

                #reset smiles files input for the next iteration
                smiles_files = list()
                smiles_files.append('sorted_{0}'.format(output_file_string))
                #print('Output file string')
                #print(output_file_string)
                debug_file.write('New input file for next iteration: {0} \n'.format(smiles_files))
                counter = counter + 1
        finally:
            for f in fs:
                if f and not f.closed:
                    f.close()
            fs = []
        debug_file.close()
