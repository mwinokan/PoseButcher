
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

		# self._protein_clash_function = lambda atom: atom.covalent_radius*1.5
		self._protein_clash_function = lambda atom: atom.vdw_radius*0.5

	def __call__(self, pose, **kwargs):
		# wrapper for Butcher.chop
		self.chop(pose, **kwargs)

	
	### PUBLIC METHODS

	def chop(self, pose, base=None, draw='2d', bolus=True):

		'''Butcher a pose:

		pose: path, MolParse.AtomGroup, or rdkit.Mol
		
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

		if isinstance(pose, Chem.rdchem.Mol):
			pose = mp.rdkit.mol_to_AtomGroup(pose)
			atoms = pose.atoms
		else:
			assert hasattr(pose, 'atoms')
			atoms = pose.atoms

		# classify atoms
		output = {}
		for i,atom in enumerate(atoms):
			output[i] = self._classify_atom(atom, bolus=bolus)

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
			from .o3d import mesh_from_AtomGroup
			# render([self.protein_mesh, mesh_from_AtomGroup(pose, use_covalent=True)])
			self._render_meshes(protein=True, pockets='hide', fragments='hide', hull='hide', extra=mesh_from_AtomGroup(pose, use_covalent=True))

		return output

	def tag(self, pose, draw='2d'):

		output = self.chop(pose, draw=draw, bolus=False)

		return list(set([d[2] for d in output.values() if d[1] == 'pocket']))
		
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

			# create fragment bolus PDB
			
			sys = mp.System('FragmentBolus')
			chain = mp.Chain('A')
			res = mp.Residue('LIG', 1, 1)

			for atom in self.hit_atomgroup.atoms:
				atom.heterogen = False
				res.add_atom(atom)

			chain.add_residue(res)
			sys.add_chain(chain)

			self._fragment_bolus_path = 'fragment_bolus.pdb'
			mp.writePDB(self._fragment_bolus_path, sys, shift_name=True)

			# create the mesh from the PDB
			from .o3d import mesh_from_pdb, paint
			self._hit_mesh = dict(
				name='fragments',
				geometry=mesh_from_pdb(self._fragment_bolus_path, gauss=False).to_legacy()
			)
			
			paint(self._hit_mesh, [1, 0.706, 0])

		return self._hit_mesh

	@property
	def protein_mesh(self):
		if self._protein_mesh is None:
			mout.out('Generating protein mesh...')
			# from .o3d import mesh_from_AtomGroup, paint
			# self._protein_mesh = mesh_from_AtomGroup(self.protein)

			from .o3d import mesh_from_pdb, paint
			self._protein_mesh = dict(
				name='protein',
				geometry=mesh_from_pdb(self._apo_protein_path).to_legacy()
			)

			# self._protein_mesh = self._protein_mesh.to_legacy()
			
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
			self._protein_hull = {'name':'protein hull', 'geometry':convex_hull(mesh)}

		return self._protein_hull
	
	### INTERNAL METHODS

	def _parse_protein(self, protein):
		
		if isinstance(protein,str) or isinstance(protein, Path):
			self._protein = mp.parse(protein).protein_system
			self._apo_protein_path = f'apo_template.pdb'
			mp.writePDB(self._apo_protein_path, self._protein, shift_name=True)

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
			assert 'atoms' in d or 'center' in d

			radius = d['radius'] if 'radius' in d else 'mean'

			if 'atoms' not in d:
				atoms = None
				center = d['center']
			else:
				atoms = [self._get_protein_atom(s) for s in d['atoms']]
				center = None

			if 'shift' in d:
				shift = d['shift']
			else:
				shift = None

			self._spherical_pocket_from_atoms(name, atoms, center=center, radius=radius, shift=shift)

		self._clip_pockets()

	def _spherical_pocket_from_atoms(self, name, atoms, center=None, radius='mean', shift=None, subtract_protein=False):

		import random
		from .o3d import sphere #, subtract_atoms

		# sphere centred between given atoms

		if center is None:
			com = sum([a.np_pos for a in atoms])/len(atoms)
			
			if radius == 'mean':
				r = sum([np.linalg.norm(a.np_pos - com) for a in atoms])/len(atoms)
			elif radius == 'max':
				r = max([np.linalg.norm(a.np_pos - com) for a in atoms])
			elif radius == 'min':
				r = min([np.linalg.norm(a.np_pos - com) for a in atoms])
			else:
				r = float(radius)
		
		else:
			com = center
			r = float(radius)

		if shift:
			com += np.array(shift)

		com_str = ', '.join([f'{v:.2f}' for v in com])
		mout.header(f'Pocket "{name}", radius={r:.2f}, center=[{com_str}]')
		# mout.var('center', com, unit='Å')
		# mout.var('radius', r, unit='Å')

		mesh = sphere(r, com)

		# if subtract_protein:
		# 	mout.out('subtracting protein...')
		# 	from open3d.t.geometry import TriangleMesh
		# 	protein = TriangleMesh.from_legacy(self.protein_mesh['geometry'])
		# 	mesh = mesh.boolean_difference(protein)

		from open3d.visualization.rendering import MaterialRecord
		mat = MaterialRecord()
		mat.base_color = [
			random.random(),
			random.random(),
			random.random(), 1.0
		]
		mat.shader = "defaultLit"

		self._new_pocket(name, mesh, mat, radius=r)

	def _clip_pockets(self, protein=True, pockets=True, hull=False, pocket_bisector=False):
		
		# from .o3d import convex_hull
		from open3d.t.geometry import TriangleMesh
		mout.out('Clipping pockets...')

		# clip the pockets to the protein
		if pockets:
			mout.out('pocket-pocket intersection...')

			for i,pocket1 in enumerate(self.pocket_meshes):

				for pocket2 in self.pocket_meshes[i+1:]:

					# mout.header(f"{pocket1['name']} {pocket2['name']}")

					center1 = pocket1['geometry'].get_center()
					center2 = pocket2['geometry'].get_center()

					distance = np.linalg.norm((center1 - center2).numpy())
					# mout.var('distance',distance)

					r_sum = pocket1['radius'] + pocket2['radius']
					# mout.var('r_sum',r_sum)

					if distance >= r_sum:
						continue

					if pocket_bisector:
						plane_center = (center1 + center2)/2
					else:
						x = (distance*distance - pocket2['radius']*pocket2['radius'] + pocket1['radius']*pocket1['radius'])/(2*distance)
						# mout.var('x',x)

						plane_center = center1 + x / distance * (center2 - center1).numpy()
					# mout.var('plane_center',plane_center)

					plane_normal = (center1 - center2)
					# mout.var('plane_normal',plane_normal)

					pocket1['geometry'] = pocket1['geometry'].clip_plane(plane_center, plane_normal)

					pocket2['geometry'] = pocket2['geometry'].clip_plane(plane_center, -plane_normal)

			mout.out('pocket convex hull...')
			for pocket in self.pocket_meshes:
				pocket['geometry'] = pocket['geometry'].compute_convex_hull()
				# pocket = convex_hull(pocket)

		# clip the pockets by their bisector planes
		if protein:
			mout.out('clipping pockets (protein)')
			
			protein = TriangleMesh.from_legacy(self.protein_mesh['geometry'])

			for mesh in self.pocket_meshes:
				mesh['geometry'] = mesh['geometry'].boolean_difference(protein)

		# clip the pockets to the convex hull of the protein
		if hull:
			mout.out('clipping pockets (protein hull)')

			protein = TriangleMesh.from_legacy(self.protein_hull)

			for mesh in self.pocket_meshes:
				mesh['geometry'] = mesh['geometry'].boolean_intersection(protein)

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

	def _new_pocket(self, name, mesh, material=None, **kwargs):
		self._pockets[name] = {'name':name, 'geometry':mesh}
		if material:
			self._pockets[name]['material'] = material

		if kwargs:
			for k,v in kwargs.items():
				self._pockets[name][k] = v

	def _classify_atom(self, atom, bolus=True):
		
		from .o3d import is_point_in_mesh

		if bolus and is_point_in_mesh(self.hit_mesh, atom.position):
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

	def _render_meshes(self, protein=True, pockets=True, fragments=True, hull=False, wireframe=False, extra=None):
		from .o3d import render

		meshes = []

		if protein:
			if protein == 'hide':
				self.protein_mesh['is_visible'] = False
			else:
				self.protein_mesh['is_visible'] = True
			
			meshes.append(self.protein_mesh)

		if fragments:
			if fragments == 'hide':
				self.hit_mesh['is_visible'] = False
			else:
				self.hit_mesh['is_visible'] = True

			meshes.append(self.hit_mesh)

		if hull:
			if hull == 'hide':
				self.protein_hull['is_visible'] = False
			else:
				self.protein_hull['is_visible'] = True

			meshes.append(self.protein_hull)

		if pockets:
			if pockets == 'hide':
				for pocket in self.pocket_meshes:
					pocket['is_visible'] = False
			else:
				for pocket in self.pocket_meshes:
					pocket['is_visible'] = True

			meshes += self.pocket_meshes

		if extra:
			if isinstance(extra, list):
				meshes += extra
			else:
				meshes.append(extra)

		render(meshes, wireframe=wireframe)
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
