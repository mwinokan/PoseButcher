
from rdkit import Chem
from rdkit.Chem import PandasTools
import molparse as mp
from pathlib import Path
import mgo
import plotly.graph_objects as go

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
		self._hit_df = None
		self._hit_atomgroup = None

		self._parse_protein(protein)
		self._parse_hits(hits)
		self._parse_pockets(pockets)

		self._build_fragment_bolus()

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

		if isinstance(pose, Chem.rdchem.Mol):
			pose = mp.rdkit.mol_to_AtomGroup(pose)
			atoms = pose.atoms
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

			fig = go.Figure()

			fig = pose.plot3d(fig=fig,show=False)

			self._plot_fragment_bolus(fig=fig)

			fig.show()

			# group = 

		return output

		
	### PROPERTIES

	@property
	def protein(self):
		return self._protein

	@property
	def pockets(self):
		return self._pockets

	@property
	def hit_df(self):
		return self._hit_df

	@property
	def hit_atomgroup(self):
		return self._hit_atomgroup

	@property
	def hit_volume(self):
		return self._hit_volume

	
	### INTERNAL METHODS

	def _parse_protein(self, protein):
		
		if isinstance(protein,str) or isinstance(protein, Path):
			self._protein = mp.parse(protein)

		else:
			raise NotImplementedError

	def _parse_hits(self, hits):
		
		if isinstance(hits,str) or isinstance(hits, Path):
			self._hit_df = PandasTools.LoadSDF(hits)

		else:
			raise NotImplementedError

	def _parse_pockets(self, pockets):
		pass

	def _classify_atom(self, atom):
		
		if self.hit_volume.is_inside(atom.position):
			return 'GOOD: fragment space'
		else:
			return 'BAD: solvent space'

	def _build_fragment_bolus(self):

		atoms = []

		for name,mol in zip(self.hit_df['ID'],self.hit_df['ROMol']):

			for atom in mp.rdkit.mol_to_AtomGroup(mol).atoms:
				atom.residue = name
				atoms.append(atom)

		self._hit_atomgroup = mp.AtomGroup.from_any('Fragment Bolus', atoms)

		self._hit_volume = mp.monte.CompoundVolume()

		for atom in self.hit_atomgroup.atoms:

			vol = mp.monte.Sphere(atom.position, atom.vdw_radius)

			self.hit_volume.add_volume(vol)

		self.hit_volume.simplify()

		pass

	def _plot_fragment_bolus(self, fig=None):

		if not fig:
			fig = go.Figure()

		for vol in self.hit_volume.volumes:

			fig.add_trace(mgo.sphere_trace(vol.centre, vol.radius))

		return fig