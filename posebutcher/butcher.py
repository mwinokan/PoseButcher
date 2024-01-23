
from rdkit import Chem
from rdkit.Chem import PandasTools, AllChem, rdDepictor
import molparse as mp
from pathlib import Path
import mgo
import plotly.graph_objects as go
import copy
import mout

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
		self._hit_mesh = None
		
		self._protein_mesh = None
		self._protein_hull = None

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
				- 'render' open3d render

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
			assert draw in ['2d', '3d', 'debug', 'render']

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

		if draw == 'debug':

			fig = go.Figure()

			fig = pose.plot3d(fig=fig,show=False)

			self._plot_fragment_bolus(fig=fig)

			fig.show()

		elif draw == '2d':

			mout.header(pose.name)

			mol = mp.rdkit.mol_from_pdb_block(pose.pdb_block)
			rdDepictor.Compute2DCoords(mol)

			drawing = mp.rdkit.draw_highlighted_mol(mol,output_to_color_pairs(output))
			display(drawing)

		elif draw == 'render':

			from .o3d import render, mesh_from_AtomGroup
			render([self.protein_mesh, mesh_from_AtomGroup(pose, use_covalent=True)])

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
	def hit_mesh(self):
		if self._hit_mesh is None:
			from .o3d import mesh_from_AtomGroup, paint
			self._hit_mesh = mesh_from_AtomGroup(self.hit_atomgroup)
			# paint(self._hit_mesh, [1, 0.706, 0])
		return self._hit_mesh

	@property
	def protein_mesh(self):
		if self._protein_mesh is None:
			from .o3d import mesh_from_AtomGroup, paint
			self._protein_mesh = mesh_from_AtomGroup(self.protein)
			paint(self._protein_mesh, [0.098, 0.463, 0.824])
		return self._protein_mesh

	@property
	def protein_hull(self):
		if self._protein_hull is None:
			from .o3d import convex_hull
			mesh = copy.deepcopy(self.protein_mesh['geometry'])
			self._protein_hull = convex_hull(mesh)

		return self._protein_hull
	
	### INTERNAL METHODS

	def _parse_protein(self, protein):
		
		if isinstance(protein,str) or isinstance(protein, Path):
			self._protein = mp.parse(protein).protein_system

		else:
			raise NotImplementedError

	def _parse_hits(self, hits):
		
		if isinstance(hits,str) and hits.endswith('.sdf'):
			self._hit_df = PandasTools.LoadSDF(hits)
		
		if isinstance(hits, Path) and hits.is_dir():
			
			raise NotImplementedError

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
		
		from .o3d import is_point_in_mesh

		if is_point_in_mesh(self.hit_mesh, atom.position):
			return ('GOOD','fragment space')

		# elif is_point_in_mesh(self.protein_mesh, atom.position, within=atom.covalent_radius):
		elif is_point_in_mesh(self.protein_mesh, atom.position, within=atom.vdw_radius):
			return ('BAD','protein clash')

		# elif is_point_in_mesh(self.protein_hull, atom.position):
		# 	return ('GOOD','pocket?')

		else:
			return ('BAD','solvent space')

	def _build_fragment_bolus(self):

		atoms = []

		for name,mol in zip(self.hit_df['ID'],self.hit_df['ROMol']):

			for atom in mp.rdkit.mol_to_AtomGroup(mol).atoms:
				atom.residue = name
				atoms.append(atom)

		self._hit_atomgroup = mp.AtomGroup.from_any('Fragment Bolus', atoms)

	def _plot_fragment_bolus(self, fig=None):

		if not fig:
			fig = go.Figure()

		for vol in self.hit_volume.volumes:

			fig.add_trace(mgo.sphere_trace(vol.centre, vol.radius))

		return fig

	def _render_meshes(self, wireframe=False):
		from .o3d import render
		render([self.hit_mesh, self.protein_mesh], wireframe=wireframe)
		# render([self.hit_mesh], wireframe=wireframe)

def output_to_color_pairs(output):

	lookup = {
		'fragment space': None,
		'pocket?': (0,1,0),
		'protein clash': (1,0,0),
		'solvent space': (0,0,1),
	}

	pairs = []
	for k,v in output.items():

		c = lookup[v[1]]

		if c is None:
			continue

		pairs.append((k,c))

	return pairs