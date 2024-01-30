
# MolParse
import molparse as mp
from molparse import AtomGroup

# for type hinting
from pathlib import Path
from rdkit.Chem.rdchem import Mol

# logging
import mcol
from mlog import setup_logger
logger = setup_logger('PoseButcher')

class PoseButcher:

	"""

	POSE BUTCHER
	============

	"A good butcher always trims the fat"

	PoseButcher is a tool for categorising and segmenting virtual hits with reference to experimental protein structures and (fragment) hits.
	
	Ligand atoms are tagged with categories:

		- GOOD:

			* fragment space: within the fragment bolus
			* pocket X: in a specified catalytic/allosteric pocket X

		- BAD:
			
			* protein clash: Clashing with the protein
			* solvent space: Heading out of the protein/crystal

	Usage
	-----

		1. Create the butcher (see PoseButcher.__init__):

			from posebutcher import PoseButcher
			butcher = PoseButcher(protein, hits, pockets)

		2. Chop up a posed virtual hit (rdkit.ROMol with a conformer):

			result = butcher.chop(mol)

		3. Tag a compound based on its pocket occupancy and clashes:

			tags = butcher.tag(mol)

		4. (Coming soon) Trim a parts of a compound that clash with a protein or leave the crystal

			mol = butcher.trim(mol)

		5. (Coming soon) Explore the expansion opportunities from a given atom in a virtual hit

			result = butcher.explore(mol, index, direction)

		6. (Coming soon) Score how well a virtual hit recapitulates shape and colour of the fragment bolus

			score: float = butcher.score(mol)

	"""

	### DUNDERS

	def __init__(self, 
		protein: str | Path, 
		fragments: str | Path, 
		pockets: dict [str, dict] | None = None,
	) -> None:

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

		logger.title('Creating PoseButcher')

		self._pockets = {}
		
		self._protein = None 				# molparse.System
		self._protein_mesh = None 			# open3d.geometry.TriangleMesh
		self._protein_hull = None 			# open3d.geometry.TriangleMesh

		self._fragment_df = None  			# pandas.DataFrame
		self._fragment_atomgroup = None 	# molparse.AtomGroup
		self._fragment_mesh = None 			# open3d.geometry.TriangleMesh

		self._parse_protein(protein)
		self._parse_fragments(fragments)
		
		if pockets:
			self._parse_pockets(pockets)

		self._build_fragment_bolus()

		# define atom clashes with the protein surface:
		self._protein_clash_function = lambda atom: atom.vdw_radius*0.5

		logger.success('PoseButcher initialised')

	
	### PUBLIC METHODS

	def chop(self, 
		pose: Mol | AtomGroup, 
		base: str | Mol | None = None, 
		draw: str | None | bool = '2d', 
		fragments: bool = False,
	) -> dict [str, tuple]:

		'''

		For each atom [1] in the provided molecule evaluates if it is:

		* in the fragment bolus [2]: 	('GOOD', 'fragment space')
		* in a defined pocket X: 		('GOOD', 'pocket', X)
		
		* clashing with the protein: 	('BAD', 'protein clash')
		* leaving the protein: 			('BAD', 'solvent space')

		[1] if a base is provided atoms in the MCS between [pose,base] are ignored
		[2] only if the fragments argument is truthy

		Draw must be one of:
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

		from rdkit import Chem

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
			output[i] = self._classify_atom(atom, fragments=fragments)

		if base:
			if isinstance(base,str):
				base = mp.rdkit.mol_from_smiles(base)

			mol = pose.rdkit_mol
			from rdkit.Chem import rdFMCS
			res = rdFMCS.FindMCS([mol, base])
			mcs_mol = Chem.MolFromSmarts(res.smartsString)
			matches = mol.GetSubstructMatch(mcs_mol)

			for i in matches:
				output[i] = ('BASE', 'atom in base')

		# render the result
		if draw == '2d':

			# get a flat depiction of the molecule
			from rdkit.Chem import rdDepictor
			mol = pose.rdkit_mol
			rdDepictor.Compute2DCoords(mol)

			# add labels to the atoms
			for atom in mol.GetAtoms():
				atom.SetProp('atomNote',output_to_label(output, atom.GetIdx()))

			# create the drawing
			drawing = mp.rdkit.draw_highlighted_mol(
				mol,
				self._output_to_color_pairs(output),
				legend=pose.name,
			)

			display(drawing)

		elif draw == '3d':

			from .o3d import mesh_from_AtomGroup
			self._render_meshes(
				protein=True, 
				pockets='hide', 
				fragments='hide', 
				hull='hide', 
				extra=mesh_from_AtomGroup(pose, use_covalent=True)
			)

		return output

	def tag(self, 
		pose: Mol | AtomGroup, 
		base: str | Mol | None = None, 
		draw: str | None | bool = False, 
		pockets_only: bool = False,
	) -> set:

		"""Return a set of string tags signifying which pockets and optionally 
		any protein/solvent clashes are relevant to this ligand.

		If pockets_only is truthy, protein/solvent clashes are not tagged.

		See the docstring for PoseButcher.chop() to see information about the other arguments."""

		output = self.chop(pose, base=base, draw=draw, fragments=False)

		if pockets_only:
			return set([d[2] for d in output.values() if d[1] == 'pocket'])
		else:
			return set([d[2] if d[1] == 'pocket' else d[1] for d in output.values()])
		
	def render(self, hull='hide', **kwargs):
		self._render_meshes(hull=hull, **kwargs)

	def trim(self):
		try: raise NotImplementedError; 
		except:
			logger.exception(f"{mcol.func}PoseButcher.trim{mcol.clear} method not implemented ")

	def explore(self):
		try: raise NotImplementedError; 
		except:
			logger.exception(f"{mcol.func}PoseButcher.explore{mcol.clear} method not implemented ")

	def score(self):
		try: raise NotImplementedError; 
		except:
			logger.exception(f"{mcol.func}PoseButcher.score{mcol.clear} method not implemented ")


	### PROPERTIES

	@property
	def protein(self):
		return self._protein

	@property
	def pockets(self):
		return self._pockets

	@property
	def fragment_df(self):
		return self._fragment_df

	@property
	def fragment_atomgroup(self):
		return self._fragment_atomgroup

	@property
	def fragment_mesh(self):
		if self._fragment_mesh is None:

			# create fragment bolus PDB
			
			sys = mp.System('FragmentBolus')
			chain = mp.Chain('A')
			res = mp.Residue('LIG', 1, 1)

			for atom in self.fragment_atomgroup.atoms:
				atom.heterogen = False
				res.add_atom(atom)

			chain.add_residue(res)
			sys.add_chain(chain)

			self._fragment_bolus_path = '_butcher_fragments.pdb'
			logger.writing(self._fragment_bolus_path)
			mp.writePDB(self._fragment_bolus_path, sys, shift_name=True, verbosity=0)

			# create the mesh from the PDB
			from .o3d import mesh_from_pdb, paint
			logger.warning('excuse the PyGAMer warnings... (they are safe to ignore)')
			self._fragment_mesh = dict(
				name='fragments',
				geometry=mesh_from_pdb(self._fragment_bolus_path, gauss=False).to_legacy()
			)
			
			paint(self._fragment_mesh, FRAGMENT_COLOR)

		return self._fragment_mesh

	@property
	def protein_mesh(self):
		if self._protein_mesh is None:
			logger.info('Generating protein mesh...')

			from .o3d import mesh_from_pdb, paint
			self._protein_mesh = dict(
				name='protein',
				geometry=mesh_from_pdb(self._apo_protein_path).to_legacy()
			)
			
			paint(self._protein_mesh, PROTEIN_COLOR)

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
			logger.debug('Generating protein convex hull...')
			from .o3d import convex_hull, paint
			from copy import deepcopy
			mesh = deepcopy(self.protein_mesh['geometry'])
			paint(mesh, PROTEIN_COLOR)
			self._protein_hull = {'name':'protein hull', 'geometry':convex_hull(mesh)}

		return self._protein_hull
	
	### INTERNAL METHODS

	def _parse_protein(self, protein):
		
		if isinstance(protein,str) or isinstance(protein, Path):
			logger.reading(protein)
			self._protein = mp.parse(protein, verbosity=0).protein_system
			self._apo_protein_path = f'_butcher_protein.pdb'
			logger.writing(self._apo_protein_path)
			mp.writePDB(self._apo_protein_path, self._protein, shift_name=True, verbosity=0)

		else:
			raise NotImplementedError

	def _parse_fragments(self, fragments):
		
		if isinstance(fragments,str) and fragments.endswith('.sdf'):
			logger.reading(fragments)
			from rdkit.Chem import PandasTools
			self._fragment_df = PandasTools.LoadSDF(fragments)
			return
		
		elif isinstance(fragments, Path) and fragments.is_dir():
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

			from numpy.linalg import norm
			
			if radius == 'mean':
				r = sum([norm(a.np_pos - com) for a in atoms])/len(atoms)
			elif radius == 'max':
				r = max([norm(a.np_pos - com) for a in atoms])
			elif radius == 'min':
				r = min([norm(a.np_pos - com) for a in atoms])
			else:
				r = float(radius)
		
		else:
			com = center
			r = float(radius)

		if shift:
			from numpy import array
			com += array(shift)

		com_str = ', '.join([f'{v:.2f}' for v in com])
		logger.header(f'Pocket "{name}", radius={r:.2f}, center=[{com_str}]')

		mesh = sphere(r, com)

		from open3d.visualization.rendering import MaterialRecord
		mat = MaterialRecord()
		color = POCKET_COLORS[len(self.pocket_meshes)]
		mat.base_color = [
			color[0],
			color[1],
			color[2],
			0.5,
		]
		mat.shader = "defaultLit"

		self._new_pocket(name, mesh, mat, radius=r, color=color)

	def _clip_pockets(self, protein=True, pockets=True, hull=False, pocket_bisector=False):
		
		from numpy.linalg import norm
		from open3d.t.geometry import TriangleMesh

		from open3d import utility
		utility.set_verbosity_level(utility.VerbosityLevel.Error)

		logger.info('Clipping pockets...')

		# clip the pockets to the protein
		if pockets:
			logger.debug('pocket-pocket intersection...')

			for i,pocket1 in enumerate(self.pocket_meshes):

				for pocket2 in self.pocket_meshes[i+1:]:

					center1 = pocket1['geometry'].get_center()
					center2 = pocket2['geometry'].get_center()

					distance = norm((center1 - center2).numpy())

					r_sum = pocket1['radius'] + pocket2['radius']

					if distance >= r_sum:
						continue

					if pocket_bisector:
						plane_center = (center1 + center2)/2
					else:
						x = (distance*distance - pocket2['radius']*pocket2['radius'] + pocket1['radius']*pocket1['radius'])/(2*distance)

						plane_center = center1 + x / distance * (center2 - center1).numpy()

					plane_normal = (center1 - center2)

					pocket1['geometry'] = pocket1['geometry'].clip_plane(plane_center, plane_normal)

					pocket2['geometry'] = pocket2['geometry'].clip_plane(plane_center, -plane_normal)

			logger.debug('pocket convex hull...')
			for pocket in self.pocket_meshes:
				pocket['geometry'] = pocket['geometry'].compute_convex_hull()

		# clip the pockets by their bisector planes
		if protein:
			logger.debug('clipping pockets (protein)')
			
			protein = TriangleMesh.from_legacy(self.protein_mesh['geometry'])

			for mesh in self.pocket_meshes:
				mesh['geometry'] = mesh['geometry'].boolean_difference(protein)

		# clip the pockets to the convex hull of the protein
		if hull:
			logger.debug('clipping pockets (protein hull)')

			protein = TriangleMesh.from_legacy(self.protein_hull)

			for mesh in self.pocket_meshes:
				mesh['geometry'] = mesh['geometry'].boolean_intersection(protein)

		logger.success('Pocket clipping complete')

		utility.set_verbosity_level(utility.VerbosityLevel.Warning)

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

	def _classify_atom(self, atom, fragments=True):
		
		from .o3d import is_point_in_mesh

		if fragments and is_point_in_mesh(self.fragment_mesh, atom.position):
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

		for name,mol in zip(self.fragment_df['ID'], self.fragment_df['ROMol']):

			for atom in mp.rdkit.mol_to_AtomGroup(mol).atoms:
				atom.residue = name
				atoms.append(atom)

		self._fragment_atomgroup = mp.AtomGroup.from_any('Fragment Bolus', atoms)

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
				self.fragment_mesh['is_visible'] = False
			else:
				self.fragment_mesh['is_visible'] = True

			meshes.append(self.fragment_mesh)

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

	def _output_to_color_pairs(self,output):

		pairs = []
		for k,v in output.items():

			if v[1] == 'pocket':
				c = self.pockets[v[2]]['color']

			elif v[1] == 'protein clash':
				c = PROTEIN_COLOR

			elif v[1] == 'solvent space':
				c = SOLVENT_COLOR

			elif v[0] == 'BASE':
				c = BASE_COLOR

			else:
				continue

			pairs.append((k,c))

		return pairs

def output_to_label(output, index):

	output_tuple = output[index]

	if output_tuple[1] == 'pocket':
		return f'{output_tuple[2]}'
	
	if output_tuple[1] == 'solvent space':
		return 'SOL.'
	
	if output_tuple[1] == 'protein clash':
		return 'PROT.'

	return ''
    
PROTEIN_COLOR = (0.8392156862745098, 0.15294117647058825, 0.1568627450980392)  # 'tab:red'
SOLVENT_COLOR = (0.12156862745098039, 0.4666666666666667, 0.7058823529411765)  # 'tab:blue'
FRAGMENT_COLOR = (1.0, 0.4980392156862745, 0.054901960784313725) 			   # 'tab:orange'
BASE_COLOR = (0.4980392156862745, 0.4980392156862745, 0.4980392156862745)   # 'tab:gray'

POCKET_COLORS = [
    (0.17254901960784313, 0.6274509803921569, 0.17254901960784313), # 'tab:green'
    (0.5803921568627451, 0.403921568627451, 0.7411764705882353),    # 'tab:purple'
    (0.5490196078431373, 0.33725490196078434, 0.29411764705882354), # 'tab:brown'
    (0.8901960784313725, 0.4666666666666667, 0.7607843137254902),   # 'tab:pink'
    (0.7372549019607844, 0.7411764705882353, 0.13333333333333333),  # 'tab:olive'
    (0.09019607843137255, 0.7450980392156863, 0.8117647058823529),  # 'tab:cyan'
    
    (0.4980392156862745, 0.4980392156862745, 0.4980392156862745),   # 'tab:gray'
    (1.0, 0.4980392156862745, 0.054901960784313725), 				# 'tab:orange'
    (0.8392156862745098, 0.15294117647058825, 0.1568627450980392),  # 'tab:red'
]
