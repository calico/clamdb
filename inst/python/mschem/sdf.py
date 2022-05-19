import os

import pandas as pd

from rdkit.Chem import *
from rdkit import Chem, Geometry
from rdkit.Chem import AllChem, Descriptors
from rdkit import RDConfig
from rdkit.Chem import rdDepictor

from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import MACCSkeys

def parse_sdf(sdf_data):
  
    """
    Parse SDF
    
    Extract structural information and identifiers from SDF data. Adapted from
      Gist (https://raw.githubusercontent.com/NLeSC/MAGMa/master/pubchem/process_hmdb.py)
      modified by Eugene Melamud
    
    Parameters:
    ----------
    sdf_data (_io.TextIOWrapper):
        An opened sdf file or streaming results
    
    Returns:
    ----------
    List of results or None if the query was invalid
    
    """
  
    sdf_entries = []
  
    line = '$$$$'
    while line != "":
        record = []
        amap = {}
        skip = False
        # read header:
        for x in range(4):
            line = sdf_data.readline()
            record.append(line)
        if line == "":
            continue
        natoms = int(record[-1][:3])
        nbonds = int(record[-1][3:6])
        bonds = 0
        y = 0

        for x in range(natoms):
            line = sdf_data.readline()
            #if line[31:33] == 'H ':
            #    # skip hydrogens
            #    continue
            y += 1
            amap[x + 1] = y
            if line[31:33] not in ['H ', 'C ', 'N ', 'O ', 'P ', 'S ', 'F ', 'Cl', 'Br', 'I ']:
                # filter non-organic compounds
                skip = True
            elif line[50:51] != '0':
                # this flag has something to do with polymeric structures
                # and resulted in deviation between calculated and given inchikeys, skip
                skip = True
            elif line[38:39] == '4':
                # radical, resulted in deviation between calculated and given inchikeys
                skip = True
            record.append(line[:42] + '\n')
        for x in range(nbonds):
            
            line = sdf_data.readline()
            a1 = int(line[:3])
            a2 = int(line[3:6])
            
            if a1 in amap and a2 in amap:
                bonds += 1
                record.append('%3i%3i%s%3i\n' %
                              (amap[a1], amap[a2], line[6:9], int(line[11:12])))
        while line.startswith('M  END') == False and line != '':
            line = sdf_data.readline()
            record.append(line)
            if line[:6] == 'M  ISO':
                skip = True
                print('Skipped isotopically labeled:', record[0][:-1])
        while line != "$$$$\n" and line != "":
            line = sdf_data.readline()
            if line.startswith("> <HMDB_ID>"):
                compound_id = str(sdf_data.readline()[:-1])
            if line.startswith(">  <MSMLS_ID>"):
                compound_id = str(sdf_data.readline()[:-1])
                molname = compound_id
            if line.startswith("> <LM_ID>"):
                compound_id = str(sdf_data.readline()[:-1])
            elif line.startswith("> <ChEBI ID>"):
                compound_id = str(sdf_data.readline()[:-1])
            elif line.startswith("> <DATABASE_ID>"):
                compound_id = str(sdf_data.readline()[:-1])
            elif line.startswith(">  <BARCODE>"):
                compound_id = str(sdf_data.readline()[:-1])
            #elif line.startswith(">  <MDLNUM>"):
            #    compound_id = str(sdf_data.readline()[:-1])
            elif line.startswith("> <PUBCHEM_COMPOUND_CID>"):
                compound_id = "PUBCHEM:" +  str(sdf_data.readline()[:-1])
            elif line.startswith(">  <KEGGID>"):
                compound_id = str(sdf_data.readline()[:-1])
                molname = compound_id
            elif line.startswith(">  <SUPPLIER_REGID>"):
                molname = str(sdf_data.readline()[:-1])
            elif line.startswith("> <COMMON_NAME>"):
                molname = str(sdf_data.readline()[:-1])
            elif line.startswith("> <GENERIC_NAME>"):
                molname = str(sdf_data.readline()[:-1])
            elif line.startswith("> <PUBCHEM_IUPAC_NAME>"):
                molname = str(sdf_data.readline()[:-1])
            elif line.startswith("> <PUBCHEM_TRADITIONAL_NAME>"):
                molname = str(sdf_data.readline()[:-1])
            elif line.startswith(">  <CHEM_NAME>"):
                molname = str(sdf_data.readline()[:-1])
            elif line.startswith("> <ChEBI Name>"):
                molname = str(sdf_data.readline()[:-1])
            elif line.startswith("> <PRODUCT>"):
                molname = str(sdf_data.readline()[:-1])
            elif line.startswith("> <SYSTEMATIC_NAME>"):
                molname = str(sdf_data.readline()[:-1])
            elif line.startswith("> <INCHI_KEY>"):
                inchi_key = sdf_data.readline()[9:-1]
            elif line.startswith("> <PUBCHEM_IUPAC_INCHIKEY>"):
                inchi_key = sdf_data.readline()[0:-1]
        if line != "" and skip == False:
            
            try:
                compound_id 
                molname
            except:
                print("Id or Name is Missing!.. skipping")
                continue
          
            # format a molecules raw attrbutes
            
            record[3] = repr(y).rjust(3) + repr(bonds).rjust(3) + record[3][6:]
            sdf_record = format_sdf_record(record)
            if sdf_record is not None:
              
                entry_dict = {
                  "id" : compound_id,
                  "name" : molname[:255]
                }
              
                entry_dict = {**entry_dict, **sdf_record}
              
                sdf_entries.append(entry_dict)

    return pd.DataFrame(sdf_entries)

def format_sdf_record(record):
  
    """
    Format SDF record
    
    Format an SDF record to extract and standardize structural features.
      Adapted from Gist (https://raw.githubusercontent.com/NLeSC/MAGMa/master/pubchem/process_hmdb.py)
      modified by Eugene Melamud
    
    Parameters:
    ----------
    record (list):
        A list containing an SDF's mol block
    
    Returns:
    ----------
    sdf_record (dict):
        Dictionary containing structral information or None if entry is invalid
    
    """
  
    # format a molecules raw attrbutes
    
    molblock = ''.join(record)
    mol = Chem.MolFromMolBlock(molblock)
  
    if mol == None or mol.GetNumAtoms() == 0: 
        print("Empty molecule..")
        return None
  
    fragments =  Chem.GetMolFrags(mol,asMols=True)
    if len(fragments) > 1:
        fragSize = 0
  
        for frag in fragments:
            if frag.GetNumAtoms() > fragSize:
                fragSize = frag.GetNumAtoms()
                mol = frag   #set to largest fragment
  
        if fragSize <= 1: 
            print("Skipping.. largest fragment is less then 1 atom")
            return None
            
    smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
    conf = mol.GetConformer(0)
    #molblock = zlib.compress(''.join(record))
    molform = Chem.rdMolDescriptors.CalcMolFormula(mol)
    mw = Chem.rdMolDescriptors.CalcExactMolWt(mol)
    natoms = mol.GetNumHeavyAtoms()
    logp = Chem.Crippen.MolLogP(mol)
  
    #finger prints
    fingerPrint2 = MACCSkeys.GenMACCSKeys(mol)
    #fingerPrint1 = FingerprintMols.FingerprintMol(mol)
    #fingerPrint3 = AllChem.GetMorganFingerprintAsBitVect(mol,2,nBits=1024)
    #print fingerPrint3.ToBase64()
  
    charge = 0
    if '-' in molform:
        if molform[-1] == '-':
            charge = int(-1)
        else:
            charge = -1 * int(molform[-1])
    elif '+' in molform:
        if molform[-1] == '+':
            charge = int(1)
        else:
            charge = int(molform[-1])
    
    inchikey=""
    inchi=""
    try:
        inchi = MolToInchi(mol)
        inchikey = InchiToInchiKey(inchi)[:14]
    except:
        print("error.. generating inchikey")
        return None
  
    # add a new record
    
    record_dict = {
      "mw" : mw,
      "charge" : charge,
      "natoms" : int(natoms),
      "molblock" : "",
      "inchikey" : inchikey,
      "inchi" : inchi,
      "smiles" : smiles,
      "molform" : molform,
      "logp" : float(logp),
      "fingerPrint" : fingerPrint2.ToBitString()
    }
    
    return record_dict

def sdf_to_tsv(sdf_path, save_path, overwrite = False):
  
    """
    SDF to TSV
    
    Read an .sdf file, parse the results and export a .tsv file
    
    Parameters:
    ----------
    record (list):
        A list containing an SDF's mol block
    
    Returns:
    ----------
    sdf_record (dict):
        Dictionary containing structral information or None if entry is invalid
    
    """
  
    assert type(sdf_path) == str
    assert type(save_path) == str
    assert type(overwrite) == bool
  
    if not os.path.exists(sdf_path):
        raise FileNotFoundError(sdf_path)
    
    if os.path.isfile(save_path) and not overwrite:
        raise FileExistsError(f"{save_path} exists and overwrite is False")
    
    if os.path.isfile(save_path):
        os.remove(save_path)

    with open(sdf_path) as sdf_data:
        sdf_results = parse_sdf(sdf_data)
    
    assert type(sdf_results) == pd.DataFrame
    
    sdf_results.to_csv(save_path, sep = "\t", index = False)
    
    return None
