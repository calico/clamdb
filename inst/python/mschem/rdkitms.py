from rdkit import Chem, rdBase
from rdkit.Chem.SaltRemover import SaltRemover
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

def canonicalize_input(
  input_list,
  input_type,
  output_types,
  remove_salts = True
  ):

  """
  Canonicalize Input
  
	Takes a list of SMILES or InCHI ids and returns a dictionary of input IDs -> canonicalized outputs: SMILES, InCHI, &/or InCHI-key.
	Standard InCHI, standard InCHI key and canonical SMILES will be the same for any rotation of a molecule which would otherwise generate a different ID.
	
	Parameters
	----------
	input_list: list
	  List of SMILES or InCHI IDs
	input_type: str
	  Type of input IDs: smiles or inchi
	output_types: str/list
	  Types of outputs to be created: smiles, inchi, inchikey, formula
	remove_salts: bool
	  Remove salts from a structure

  Returns
  -------
  output_dict: dict
    Dictionary of output types

	"""

  rdBase.DisableLog('rdApp.error') # prevent rdkit messages about parse errors since we don't want side effects
  rdBase.DisableLog('rdApp.warning')

  # test inputs
  if not isinstance(input_list, list):
    raise NameError('"input_list" must be a "list"')

  if not isinstance(input_type, str):
    raise NameError('"input_type" must be of type "str"')

  valid_inputs = ['smiles', 'inchi']
  if input_type not in valid_inputs:
    raise NameError('"input_type" is not a valid input; it needs to be one of:', valid_inputs)

  if not (isinstance(output_types, list) or isinstance(output_types, str)):
    raise NameError('"output_types" must be of type "str" or "list"')
  if isinstance(output_types, str):
    output_types = [output_types]

  valid_outputs = ['smiles', 'inchi', 'inchikey', 'formula']
  if not set(valid_outputs).issuperset(set(output_types)):
    raise NameError('"output_types" were not all valid; valid values are:', valid_outputs)

  if not isinstance(remove_salts, bool):
    raise NameError('"remove_salts" must be a bool')

  # define salts to use if remove_salts is True
  remover = SaltRemover(defnData="[Ba,Br,Ca,Cl,Co,Fe,Li,Na,O,K]")

  # define dictionaries for storing possible outputs

  smiles_dict = dict()
  inchi_dict = dict()
  inchikey_dict = dict()
  formula_dict = dict()

  for an_input in input_list:

    if input_type == 'smiles':
      mol = Chem.MolFromSmiles(an_input)
    elif input_type == 'inchi':
      mol = Chem.MolFromInchi(an_input)
    else:
      raise NameError('"input_types" import option not defined')
    if mol is None:
      continue

    if remove_salts:
      mol = remover.StripMol(mol)

    if 'smiles' in output_types:
      smiles_dict[an_input] = Chem.MolToSmiles(mol)
    if 'inchi' in output_types or 'inchikey' in output_types:
      inchi_value = Chem.MolToInchi(mol)

      if 'inchi' in output_types:
        inchi_dict[an_input] = inchi_value
      if 'inchikey' in output_types:
        inchikey_dict[an_input] = Chem.InchiToInchiKey(inchi_value)
    if 'formula' in output_types:
      formula_dict[an_input] = CalcMolFormula(mol)

  # define output
  output_dict = dict()
  if 'smiles' in output_types: output_dict['smiles'] = smiles_dict
  if 'inchi' in output_types: output_dict['inchi'] = inchi_dict
  if 'inchikey' in output_types: output_dict['inchikey'] = inchikey_dict
  if 'formula' in output_types: output_dict['formula'] = formula_dict
  
  return output_dict
