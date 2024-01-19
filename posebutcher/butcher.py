
from rdkit import Chem

class PoseButcher:

	"""Pose butcher segments a ligand into categories:

		- GOOD:

			* fragment space: in catalytic fragment space
			* pocket X: in desirable pocket X

		- BAD:
			
			* solvent space: Heading out of the protein/crystal
			* protein clash: Clashing with the protein

	Usage:

		1. Create the butcher:

			butcher = Butcher(protein, hits, pockets)

		2. Chop up a ligand pose:

			result = butcher.chop(pose)

	"""

	### DUNDERS

	def __init__(self, protein, hits, pockets):

		'''Create a Butcher with a protein structure and pocket information

		protein: path or MolParse.System object

		hits: list of MolParse.AtomGroups / rdkit.Mol / MolParse shape

		pocket: dictionary with pocket names as keys and values:

				- list of residue strings (MolParse naming)
					- 'n101' residue with number 101
					- 'TYR 101' tyrosine 101 
					- '101' residue with index 101
				
				- MolParse primitive or compound shape, convex hull

		'''

		self._pockets = {}
		self._protein = None

		self._parse_protein(protein)
		self._parse_pockets(pockets)

	def __call__(self, pose, **kwargs):
		# wrapper for Butcher.chop
		self.chop(pose, **kwargs)

	
	### PUBLIC METHODS

	def chop(self, pose, protein=None, base=None, draw='2d'):

		'''Butcher a pose:

		pose: path, MolParse.AtomGroup, or rdkit.Mol

		protein: optionally pass a different protein reference
		
		base: optionally pass a reference base compound

		draw: one of:

				- '2d' rdkit 2d drawing
				- '3d' rdkit 3d render
				- 'debug' plotly 3d graph showing the various hulls 

		Returns labelled atom indices: e.g.

		{
			0: 'GOOD fragment space',
			1: 'GOOD pocket ...',
			2: 'BAD solvent space'	,
			3: 'BAD protein/clash',
			4: 'BASE ignoring base',
			...
		}

		'''

		# parse arguments

		if draw:
			assert draw in ['2d', '3d', 'debug']

		if protein:
			protein = self._parse_protein(protein)
		else:
			protein = self.protein

		if isinstance(pose, Chem.Mol):
			atoms = ...
		else:
			assert hasattr(pose, 'atoms')
			atoms = pose.atoms

		# classify atoms in a loop

		output = {}

		for i,atom in enumerate(atoms):
			output[i] = self._classify_atom(atom)

		# render the result

		if draw:
			assert draw == 'debug'

			import plotly.graph_objects as go

			fig = go.Figure()

			group = 

		return output

		
	### PROPERTIES

	@property
	def protein(self):
		return self._protein

	@property
	def pockets(self):
		return self._pockets

	
	### INTERNAL METHODS

	def _parse_protein(self, protein):
		pass

	def _parse_pockets(self, pockets):
		pass

	def _classify_atom(self, atom):
		pass
