
from rdkit import Chem
from rdkit.Chem import PandasTools, AllChem, rdDepictor
import molparse as mp
from pathlib import Path
import mgo
import plotly.graph_objects as go
import copy
import mout
import mcol
import numpy as np

class PoseButcher:

	"""

	POSE BUTCHER
	============

	"A good butcher always trims the fat"
	
	Pose butcher segments a ligand into categories:

		- GOOD:

			* fragment space: within the fragment bolus
			* pocket X: in desirable pocket X

		- BAD:
			
			* solvent space: Heading out of the protein/crystal
			* protein clash: Clashing with the protein

	Usage:

		1. Create the butcher:

			butcher = Butcher(protein, hits, pockets)

		2. Chop up an rdkit.ROMol ligand (with a conformer):

			result = butcher.chop(mol, draw='2d')

	"""

	### DUNDERS

	def __init__(self, protein, hits, pockets):

		'''Create a Butcher with a protein structure and pocket information

		protein: pathlib.Path or str to a PDB/GRO file from which the protein will be extracted 

		hits: pathlib.Path or str to an SDF file containing the fragment hits

		pocket: dictionary with pocket names as keys and nested dictionary values:

			* type='sphere'
			
			* atoms: list of strings to define atoms. These must be formatted as "RES RNUM ATOM"

				- RES: residue name
				- RNUM: residue number
				- ATOM: atom name

			* radius: can be a fixed (float) number or 'mean'/'min'/'max' 

			e.g. 

			pockets = {
			    "P1":  dict(type='sphere', atoms=['GLY 127 O', 'PRO 107 CG', 'VAL 124 CG1'], radius='mean'),
			    "P1'": dict(type='sphere', atoms=['VAL 84 CG1', 'TYR 90 CD2', 'SER 87 CB'], radius='mean'),
			}

		'''

		mout.debug('PoseButcher')

		self._pockets = {}
		
		self._protein = None 		# molparse.System
		self._hit_df = None  		# pandas.DataFrame
		self._hit_atomgroup = None 	# molparse.AtomGroup
		self._hit_mesh = None 		# open3d.geometry.TriangleMesh
		
		self._protein_mesh = None 	# open3d.geometry.TriangleMesh
		self._protein_hull = None 	# open3d.geometry.TriangleMesh

		self._parse_protein(protein)
		self._parse_hits(hits)
		self._parse_pockets(pockets)

		self._build_fragment_bolus()

		self._protein_clash_function = lambda atom: atom.covalent_radius*1.5
		# self._protein_clash_function = lambda atom: atom.vdw_radius

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
				- '3d' open3d render

		Returns labelled atom indices: e.g.

		{
			0: ('GOOD', 'fragment space'),
			1: ('GOOD', 'pocket', 'X'),
			2: ('BAD', 'solvent space'),
			3: ('BAD', 'protein clash'),
			...
		}
		'''

		# parse arguments

		if draw:
			assert draw in ['2d', '3d']

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

		# classify atoms
		output = {}
		for i,atom in enumerate(atoms):
			output[i] = self._classify_atom(atom)

		# render the result
		if draw == '2d':
			mol = mp.rdkit.mol_from_pdb_block(pose.pdb_block)
			rdDepictor.Compute2DCoords(mol)

			for atom in mol.GetAtoms():
				atom.SetProp('atomNote',output_to_label(output, atom.GetIdx()))

			drawing = mp.rdkit.draw_highlighted_mol(
				mol,
				output_to_color_pairs(output),
				legend=pose.name,
			)

			display(drawing)

		elif draw == '3d':
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
			mout.out('Generating fragment bolus mesh...')
			from .o3d import mesh_from_AtomGroup, paint
			self._hit_mesh = mesh_from_AtomGroup(self.hit_atomgroup)
			paint(self._hit_mesh, [1, 0.706, 0])
		return self._hit_mesh

	@property
	def protein_mesh(self):
		if self._protein_mesh is None:
			mout.out('Generating protein mesh...')
			from .o3d import mesh_from_AtomGroup, paint
			self._protein_mesh = mesh_from_AtomGroup(self.protein)
			paint(self._protein_mesh, [0.098, 0.463, 0.824])
		return self._protein_mesh

	@property
	def pocket_meshes(self):
		return list(self.pockets.values())
	
	@property
	def pockets(self):
		return self._pockets

	@property
	def protein_hull(self):
		if self._protein_hull is None:
			mout.out('Generating protein convex hull..')
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
			mout.out(f'parsing {mcol.file}{hits}{mcol.clear} ...', end='')
			self._hit_df = PandasTools.LoadSDF(hits)
			mout.out('Done.')
			return
		
		elif isinstance(hits, Path) and hits.is_dir():
			raise NotImplementedError

		else:
			raise NotImplementedError

	def _parse_pockets(self, pockets):
		
		for name, d in pockets.items():

			assert d['type'] == 'sphere'
			assert 'atoms' in d

			radius = d['radius'] if 'radius' in d else None

			atoms = [self._get_protein_atom(s) for s in d['atoms']]

			self._spherical_pocket_from_atoms(name, atoms, radius)

	def _spherical_pocket_from_atoms(self, name, atoms, radius='mean'):

		from .o3d import sphere #, subtract_atoms

		# sphere centred between given atoms

		com = sum([a.np_pos for a in atoms])/len(atoms)

		if radius == 'mean':
			r = sum([np.linalg.norm(a.np_pos - com) for a in atoms])/len(atoms)
		elif radius == 'max':
			r = max([np.linalg.norm(a.np_pos - com) for a in atoms])
		elif radius == 'min':
			r = min([np.linalg.norm(a.np_pos - com) for a in atoms])
		else:
			r = float(radius)

		mout.header(f'Pocket "{name}":')
		mout.var('center', com, unit='Å')
		mout.var('radius', r, unit='Å')

		mesh = sphere(r, com)

		self._new_pocket(name, mesh)

	def _get_protein_atom(self, query: str):

		"""Query must be a string containing three white-space separated substrings:

			1. residue name
			2. residue number
			3. atom name

			e.g. "TYR 90 CD2"

		"""

		res_name, res_num, atom_name = query.split()

		res = self.protein.residues[f'n{res_num}']

		assert res.name == res_name, f"Residue number queried does not match the queried name: (protein='{res.name}', query='{res_name}')"

		return res.get_atom(atom_name)

	def _new_pocket(self, name, mesh):
		self._pockets[name] = {'name':name, 'geometry':mesh}

	def _classify_atom(self, atom):
		
		from .o3d import is_point_in_mesh

		if is_point_in_mesh(self.hit_mesh, atom.position):
			return ('GOOD','fragment space')

		# protein clash
		if is_point_in_mesh(self.protein_mesh, atom.position, within=self._protein_clash_function(atom)):
			return ('BAD','protein clash')

		# pockets
		for p_name, p_mesh in self.pockets.items():
			if is_point_in_mesh(p_mesh, atom.position):
				return ('GOOD', 'pocket', p_name)

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

	def _render_meshes(self, pockets=True, wireframe=False):
		from .o3d import render
		render([self.hit_mesh, self.protein_mesh] + self.pocket_meshes, wireframe=wireframe)
		# render([self.hit_mesh] + self.pocket_meshes, wireframe=wireframe)
		# render([self.hit_mesh], wireframe=wireframe)

def output_to_color_pairs(output):

	lookup = {
		'fragment space': None,
		'pocket': (0,1,0),
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

def output_to_label(output, index):

	output_tuple = output[index]

	if output_tuple[1] == 'pocket':
		return f'{output_tuple[2]}'

	return ''
