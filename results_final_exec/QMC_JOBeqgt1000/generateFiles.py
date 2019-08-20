#!/usr/bin/python3

# Modeller imports
from modeller import *              # Load standard Modeller classes
from modeller.automodel import *    # Load the automodel class
from modeller.parallel import  *    # Load parallel Modeller class

# Python imports
import re
import os
import sys
import glob
import time
import shutil
import hashlib
import inspect
from pprint import pprint

# Biopython imports
from Bio import SeqIO
from Bio.PDB import PDBList
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB import PPBuilder
from Bio.PDB import CaPPBuilder

# Extract sequence from chain
def extract(chain):
    dict = {'ALA':'A',
            'ARG':'R',
            'ASN':'N',
            'ASP':'D',
            'CYS':'C',
            'GLU':'E',
            'GLN':'Q',
            'GLY':'G',
            'HIS':'H',
            'ILE':'I',
            'LEU':'L',
            'LYS':'K',
            'MET':'M',
            'PHE':'F',
            'PRO':'P',
            'SER':'S',
            'THR':'T',
            'TRP':'W',
            'TYR':'Y',
            'VAL':'V'
            }
    inseq = ''
    for residue in chain:
        resname = residue.get_resname()
        hetflag, resseq, icode = residue.get_id()
        if icode == ' ' and hetflag == ' ':
            inseq = inseq + dict[resname]
    return inseq

# Create seg file
def createSeg(elements, sequences, targetChain):
    f = open(basePath + elements[0] + ".seg", "w")
    f.write(">P1;" + elements[0] + "\n")
    f.write("sequence:" + elements[0] + ":FIRST:@:LAST:@: : : :\n")
    f.write(sequences[0] + "*")
    f.write("\n>P1;" + elements[1] + "\n")
    f.write("structureX:" + elements[1] + ":FIRST:" + targetChain + ":LAST:" + targetChain + ": : : :\n")
    f.write(sequences[1] + "*")
    f.close()

    
# Perform alignment with modeller salign
def doAlign(env, protein):
    env.io.atom_files_directory = ['.', './:../atom_files/']
    aln = alignment(env, file=basePath + protein + '.seg', align_codes='all')
    aln.salign(overhang=30, gap_penalties_1d=(-450, -50), alignment_type='pairwise', output='ALIGNMENT')
    aln.write(file=basePath + protein + '.ali', alignment_format='PIR')
    #aln.write(file=basePath + protein + '.pap', alignment_format='PAP')

# Create the models (same waht as MHOLline does)
def doModel(env, elements, start = 1, end = 50):
    # Start parallel job
    j = job(host='localhost')
    
    # Add 6 workers (one to each core)
    for processor in range(0, 10):
        j.append(local_slave())

    # Set some env values
    env.schedule_scale = physical.values(default=1.0, soft_sphere=0.7)
    env.io.hetatm = env.io.water = True
    env.io.atom_files_directory = ['.', '../atom_files']
    env.libs.topology.read(file='$(LIB)/top_heav.lib')
    env.libs.parameters.read(file='$(LIB)/par.lib')

    # Create automodel object
    a = automodel(env,
                  alnfile        = elements[0] + '.ali',     # alignment filename
                  knowns         = elements[1],              # codes of the templates
                  sequence       = elements[0],              # code of the target
                  assess_methods = (assess.DOPE,
                                    assess.DOPEHR,
                                    assess.normalized_dope,
                                    assess.GA341)            # Strcuture evaluation scores
                                    )
    a.starting_model= start                 # index of the first model
    a.ending_model  = end                   # index of the last model
                                            # (determines how many models to calculate)
                                            
    # Very trough variable target function method (VTFM) optimization
    a.library_schedule = autosched.slow
    a.max_var_iterations = 300
    
    a.md_level = refine.slow
    
    a.repeat_optimization = 2
    
    a.max_molpdf = 1e6
    
    a.use_parallel_job(j)               # Use the job for model building
    
    a.make()                            # do the actual comparative modeling
    
    # Create output file
    fo = open(elements[0] + '_analysis.out', 'w')
    
    ok_models = filter(lambda x: x['failure'] is None, a.outputs)
    
    lok_models = list(ok_models)
    
    fields = [x for x in lok_models[0].keys() if x.endswith(' score')]
    fields.sort()
    fields = ['molpdf'] + fields 
    
    header = '%-25s ' % 'Filename' + ' '.join(['%14s' % x for x in fields])
    
    fo.write('>> Summary of successfully produced model\n\n')
    fo.write(header)
    fo.write('-' * len(header))
    
    for mdl in lok_models:
        text = '\n%-25s' %mdl['name']
        for field in fields:
            if isinstance(mdl[field], (tuple, list)):
                text = text + ' %14.5f' % mdl[field][0]
            else:
                text = text + ' %14.5f' % mdl[field]
        fo.write(text)
    
    fo.write('\n\n>> Top model results:')
    
    keys = ['molpdf', 'DOPE score', 'GA341 score', 'Normalized DOPE score']

    for key in keys:
        m = sorted(lok_models, key = lambda i: i[key])[0]
        if key == 'GA341 score':
            fo.write("\n" + key + ': ' + str(m[key][0]) + ' (file: ' + m['name'] + ')')
        else:
            fo.write("\n" + key + ': ' + str(m[key]) + ' (file: ' + m['name'] + ')')
    
    fo.close()


def performProcesses (env, element, basePathLocal):
    createSeg([element["name"], element["targetName"]], [element["sequence"], element["targetSequence"]], element["targetChain"])
    time.sleep(1)
    doAlign(env, element["name"])
    time.sleep(1)
    doModel(env, [element["name"], element["targetName"]])
    time.sleep(1)
    if os.path.exists(basePathLocal + element["name"]):
        shutil.rmtree(basePathLocal + element["name"], ignore_errors=True)
        
    os.makedirs(basePathLocal + element["name"])
    
    for file in glob.glob(basePathLocal + element["name"] + "*"):
        if element["name"] in file or element["targetName"] in file:
            shutil.move(file, element["name"])
            
    shutil.move(basePathLocal + element["targetName"] + ".cif", element["name"])
    
    for file in glob.glob(basePathLocal + "generateFiles*"):
        if not '.py' in file:
            shutil.move(file, element["name"])
            
            
if len(sys.argv) > 1:
    fasta = sys.argv[1]
    basePath = sys.argv[2]
else:
    print("wrong inputs")
    exit()

log.verbose()    # request verbose output
env = environ()  # create a new MODELLER environment to build this model in

content = ""

# List of dicts
structures = []

# Target list (from model descending order)
targets = ["3iap_D", "2w54_H", "3b9j_K", "5j8v_D"]

errors = []

with open(basePath   + fasta, 'r') as content_file:
    content = content_file.read()

fastasSeq = content.split(">")

del fastasSeq[0]

for index, fastaSeq in enumerate(fastasSeq):
    inData = fastaSeq.split()
    
    if os.path.exists(basePath + inData[0]):
        errors.append("Protein " + inData[0] + " is already processed, I will not download its template file (" + targets[index] + ")")
        continue

    print(targets[index] + " " + str(index))
    targetSeq = targets[index].split("_")[0]
    ch = targets[index].split("_")[1]
    pdbl = PDBList()
    
    try:
        pdbl.retrieve_pdb_file(targetSeq, pdir=basePath)
    except:
        errors.append("Protein " + fastaSeq.split()[0] + " is experiencing problems on download.")
        continue
        
    parser = MMCIFParser()
    structure = parser.get_structure(targetSeq, basePath + targetSeq + ".cif")
    sequence = ''
    
    for model in structure:
       chain = model[ch]        
       sequence = sequence + extract(chain)
    
    structureList = {}
    structureList["name"] = inData[0]
    structureList["sequence"] = re.sub("(.{64})", "\\1\n", inData[1], 0, re.DOTALL)
    structureList["targetName"] = targetSeq
    structureList["targetChain"] = ch
    structureList["targetSequence"] = re.sub("(.{64})", "\\1\n", sequence, 0, re.DOTALL)
    structures.append(structureList)

#pprint(structures)

for index, data in enumerate(structures):

    element = structures[index]
    
    if os.path.exists(basePath + element["name"]):
        print("Protein " + element["name"] + " is already processed, skipping.")
        continue

    try:
        performProcesses(env, element, basePath)
    except:
        errors.append("Error while processing protein " + element["name"] + " against " + element["targetName"] + " template. Trying alternative another sequence getter. Method: cif-seqres")
        
        if os.path.exists(basePath + element["name"] + ".ali"):
            os.remove(basePath + element["name"] + ".ali")
        if os.path.exists(basePath + element["name"] + ".seg"):
            os.remove(basePath + element["name"] + ".seg")
        
        sequence = ''
        chs = SeqIO.parse(basePath + element["targetName"] + ".cif", "cif-seqres")
        
        for records in chs:
            if records.annotations["chain"] == element["targetChain"]:
                sequence = str(records.seq)

        element["targetSequence"] = re.sub("(.{64})", "\\1\n", sequence.replace('.', ''), 0, re.DOTALL)
        
        try:
            performProcesses(env, element, basePath)
        except:
            errors.append("Error while processing protein " + element["name"] + " against " + element["targetName"] + " template. Trying alternative another sequence getter. Method: PPBuilder (C-N)")
            
            if os.path.exists(basePath + element["name"] + ".ali"):
                os.remove(basePath + element["name"] + ".ali")
            if os.path.exists(basePath + element["name"] + ".seg"):
                os.remove(basePath + element["name"] + ".seg")
            
            try:
                structure = parser.get_structure(element["targetName"], basePath + element["targetName"] + ".cif")
            except:
                errors.append("Error while processing protein " + element["name"] + " against " + element["targetName"] + " template. Problems in processing structure from file. Skipping")
                continue
            
            sequence = ''
            
            try:
                ppb = PPBuilder()
            
                # Using C-N
                for pp in ppb.build_peptides(structure[0][element["targetChain"]]):
                    sequence = sequence + str(pp.get_sequence())
                    
                element["targetSequence"] = re.sub("(.{64})", "\\1\n", sequence.replace('.', ''), 0, re.DOTALL)
            
                performProcesses(env, element, basePath)
            except:
                errors.append("Error while processing protein " + element["name"] + " against " + element["targetName"] + " template. Trying alternative another sequence getter. Method: PPBuilder (CA-CA)")
            
                if os.path.exists(basePath + element["name"] + ".ali"):
                    os.remove(basePath + element["name"] + ".ali")
                if os.path.exists(basePath + element["name"] + ".seg"):
                    os.remove(basePath + element["name"] + ".seg")
                
                sequence = ''
                try:
                    ppb = CaPPBuilder()
                
                    # Using CA-CA
                    for pp in ppb.build_peptides(structure[0][element["targetChain"]]):
                        sequence = sequence + str(pp.get_sequence())
                        
                    element["targetSequence"] = re.sub("(.{64})", "\\1\n", sequence.replace('.', ''), 0, re.DOTALL)
                
                    performProcesses(env, element, basePath)
                except:
                    errors.append("Error while processing protein " + element["name"] + " against " + element["targetName"] + " template. Trying to download pdb file.")
                    
                    pdbl = PDBList()
                    
                    if os.path.exists(basePath + element["name"] + ".ali"):
                        os.remove(basePath + element["name"] + ".ali")
                    if os.path.exists(basePath + element["name"] + ".seg"):
                        os.remove(basePath + element["name"] + ".seg")
                    if os.path.exists(basePath + element["targetName"] + ".cif"):
                        os.remove(basePath + element["targetName"] + ".cif")
                    
                    try:
                        pdbl.retrieve_pdb_file(element["targetName"], file_format="pdb", pdir=basePath)
                    except:
                        print("Protein " + fastaSeq.split()[0] + " is experiencing problems on download. There is nothing more I can do... NEXT!")
                        continue
                    
                    parser = PDBParser()
                    try:
                        structure = parser.get_structure(element["targetName"], basePath + "pdb" + element["targetName"] + ".ent")
                    except:
                        try:
                            structure = parser.get_structure(element["targetName"], basePath + "pdb" + element["targetName"] + ".pdb")
                        except:
                            print("Protein " + fastaSeq.split()[0] + " is experiencing problems with its files. There is nothing more I can do... NEXT!")
                            continue
                        
                    sequence = ''
                    
                    for model in structure:
                       chain = model[element["targetChain"]]        
                       sequence = sequence + extract(chain)
                    
                    inData = fastaSeq.split()
                    element["targetSequence"] = re.sub("(.{64})", "\\1\n", sequence.replace('.', ''), 0, re.DOTALL)
                    
                    try:
                        performProcesses(env, element, basePath)
                    except:
                        errors.append("Error while processing protein " + element["name"] + " against " + element["targetName"] + " template. Trying alternative another sequence getter in pdb. Method: cif-seqres")
                        
                        if os.path.exists(basePath + element["name"] + ".ali"):
                            os.remove(basePath + element["name"] + ".ali")
                        if os.path.exists(basePath + element["name"] + ".seg"):
                            os.remove(basePath + element["name"] + ".seg")
                        
                        sequence = ''
                        try:
                            chs = SeqIO.parse(basePath + "pdb" + element["targetName"] + ".ent", "pdb-seqres")
                        except:
                            chs = SeqIO.parse(basePath + element["targetName"] + ".pdb", "pdb-seqres")
    
                        for records in chs:
                            if records.annotations["chain"] == element["targetChain"]:
                                sequence = str(records.seq)
                                
                        element["targetSequence"] = re.sub("(.{64})", "\\1\n", sequence.replace('.', ''), 0, re.DOTALL)
                        
                        try:
                            performProcesses(env, element, basePath)
                        except:
                            errors.append("Error while processing protein " + element["name"] + " against " + element["targetName"] + " template. Trying alternative another sequence getter in pdb. Method: PPBuilder (C-N)")
                            
                            if os.path.exists(basePath + element["name"] + ".ali"):
                                os.remove(basePath + element["name"] + ".ali")
                            if os.path.exists(basePath + element["name"] + ".seg"):
                                os.remove(basePath + element["name"] + ".seg")
                            
                            sequence = ''
                            try:
                                ppb = PPBuilder()
                            
                                # Using C-N
                                for pp in ppb.build_peptides(structure[0][element["targetChain"]]):
                                    sequence = sequence + str(pp.get_sequence())
                                    
                                element["targetSequence"] = re.sub("(.{64})", "\\1\n", sequence.replace('.', ''), 0, re.DOTALL)
                            
                                performProcesses(env, element, basePath)
                            except:
                                errors.append("Error while processing protein " + element["name"] + " against " + element["targetName"] + " template. Trying alternative another sequence getter in pdb. Method: PPBuilder (CA-CA)")
                                
                                if os.path.exists(basePath + element["name"] + ".ali"):
                                    os.remove(basePath + element["name"] + ".ali")
                                if os.path.exists(basePath + element["name"] + ".seg"):
                                    os.remove(basePath + element["name"] + ".seg")
                                
                                sequence = ''
                                try:
                                    ppb = CaPPBuilder()
                                
                                    # Using C-N
                                    for pp in ppb.build_peptides(structure[0][element["targetChain"]]):
                                        sequence = sequence + str(pp.get_sequence())
                                        
                                    element["targetSequence"] = re.sub("(.{64})", "\\1\n", sequence.replace('.', ''), 0, re.DOTALL)
                                
                                    performProcesses(env, element, basePath)
                                except:
                                    errors.append("Error while processing protein " + element["name"] + " against " + element["targetName"] + " template. I give up let's skip this thing...")
                                    
                                    if os.path.exists(basePath + element["name"] + ".ali"):
                                        os.remove(basePath + element["name"] + ".ali")
                                    if os.path.exists(basePath + element["name"] + ".seg"):
                                        os.remove(basePath + element["name"] + ".seg")
                                    try:
                                        os.remove(basePath + "pdb" + element["targetName"] + ".ent")
                                    except:
                                        os.remove(basePath + element["targetName"] + ".pdb")

with open('errorLog.txt', 'w') as f:
    for item in errors:
        f.write("%s\n" % item)
#pprint(errors)