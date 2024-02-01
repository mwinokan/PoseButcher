
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
		pocket_clip: bool = True,
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

		pocket_clip: disable this to improve startup performance for debugging. N.B. pockets will overlap with each other and the protein!

		'''

		logger.title('Creating PoseButcher')

		self._pockets = {}
		
		self._protein = None 				# molparse.System
		self._protein_mesh = None 			# open3d.geometry.TriangleMesh
		self._protein_hull = None 			# open3d.geometry.TriangleMesh

		self._fragment_df = None  			# pandas.DataFrame
		self._fragment_atomgroup = None 	# molparse.AtomGroup
		self._fragment_mesh = None 			# open3d.geometry.TriangleMesh

		self._pocket_clip = pocket_clip

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

	def explore(self,
		pose: Mol | AtomGroup,
		origin: int | list | None = None,
		direction: list | None = None,
		draw: str | bool | None = '3d',
		return_arrow_mesh: bool = False,
		color_by_distance: bool = False,
		samples: bool = False,
		sample_angle_max: float = 10,
		sample_angle_step: float = 2,
		sample_choose_best: bool = True,
		warnings: bool = True,
	):

		"""Explore an expansion vector(s).

		pose: rdkit.Chem.Mol or molparse.AtomGroup of the ligand

		origin: index of the atom from which to expand or None to try all atoms

		direction: direction of the expansion (optional)
		
		draw: Currently only '3d' and None/False are supported.
		
		return_arrow_mesh: Return open3d TriangleMesh objects for the arrow(s)
		
		color_by_distance: Color vectors by their length rather than their destination
		
		samples: Perform additional samples around the (guessed) direction

		sample_angle_max: Maximal angle of additional samples

		sample_angle_step: Angle distance between additional samples
		
		sample_choose_best: Return the best sample vector only
		
		warnings: Set False to silence warnings

		"""

		if isinstance(pose, Mol):
			mol = pose
			group = mp.rdkit.mol_to_AtomGroup(pose)
		else:
			mol = pose.rdkit_mol
			group = pose

		if origin is None:

			results = []
			arrows = []
			
			if not sample_choose_best:
				logger.warning('Setting sample_choose_best option to True when sampling every atom')

			for i,_ in enumerate(group.atoms):

				result = dict(atom_index=i)

				this_result, arrow = self.explore(group, 
					origin=i, 
					draw=None, 
					return_arrow_mesh=True, 
					samples=samples, 
					sample_angle_step=sample_angle_step, 
					sample_angle_max=sample_angle_max, 
					sample_choose_best=True,
				)

				# print(this_result)

				if this_result['success']:
					result |= this_result
					arrows.append(arrow)

				results.append(result)

			if samples:
				logger.success('Vector sampling complete')

			if draw:

				from .o3d import mesh_from_AtomGroup
				
				if draw:
				
					atoms_mesh = mesh_from_AtomGroup(group, use_covalent=True)
					
					if isinstance(atoms_mesh, list):
						extra = arrows + atoms_mesh
					else:
						extra = arrows + [atoms_mesh]

					self.render(extra=extra, pockets='hide', fragments='hide')

			return results

		elif samples:

			# perform multiple vector explorations and return the 'best' one
			logger.debug(f'Sampling {origin=} ...')

			results = []

			# guess the direction
			if not direction:
				direction = self._guess_vector_direction(origin, group, warnings=False)

			if direction is None:
				logger.warning(f'No successful vectors found for {origin=} (#bonds>2)')
				if return_arrow_mesh:
					return dict(success=False), None
				else:
					return dict(success=False)

			# create unit vectors perpendicular to the direction vector
			from .o3d import rotation_matrix_from_vectors
			from numpy import array, matmul, dot, tan
			rotation_matrix = rotation_matrix_from_vectors([0,0,1], direction)
			unit1 = matmul(rotation_matrix,array([1,0,0]))
			unit2 = matmul(rotation_matrix,array([0,1,0]))

			# first result along guessed direction
			result = self.explore(
				group, origin=origin, direction=direction, 
				draw=None, return_arrow_mesh=False, samples=False)
			result['sample_shift_x'] = 0
			result['sample_shift_y'] = 0
			results.append(result)

			# more results along shifted directions
			sample_coords = circular_samples(sample_angle_max, sample_angle_step)

			for i,(x,y) in enumerate(sample_coords):
				# if i%25 == 0:
				# 	logger.debug(f'Sampling {i=}/{len(sample_coords)}')
				new_direction = direction + x*unit1 + y*unit2
				result = self.explore(
					group, origin=origin, direction=new_direction, 
					draw=None, return_arrow_mesh=False, samples=False, warnings=False)
				result['sample_shift_x'] = x
				result['sample_shift_y'] = y
				results.append(result)

			# pick the best one
			sorted_results = sorted([r for r in results if r['success']], key=lambda x: x['last_intersection_distance'], reverse=True)
			if sorted_results:
				result = sorted_results[0]
			else:
				result = dict(success=False)

			if not result['success']:
				logger.warning(f'No successful vectors found for {origin=}')

			### visualise
			if draw or return_arrow_mesh:

				from .o3d import mesh_from_AtomGroup

				if isinstance(origin, int):
					orig_index = origin
					orig_atom = group.atoms[origin]
					origin = orig_atom.position
				else:
					orig_index = None
					orig_atom = None
				
				if result['success']:
					arrow_mesh = create_arrow_mesh(orig_atom, orig_index, origin, color_by_distance, result)
				else:
					arrow_mesh = None

				if draw:
				
					atoms_mesh = mesh_from_AtomGroup(group, use_covalent=True)
					
					if arrow_mesh is None:
						if isinstance(atoms_mesh, list):
							extra = atoms_mesh
						else:
							extra = [atoms_mesh]
					else:
						if isinstance(atoms_mesh, list):
							extra = [arrow_mesh] + atoms_mesh
						else:
							extra = [arrow_mesh, atoms_mesh]

					self.render(extra=extra, pockets='hide', fragments='hide')

			if return_arrow_mesh:
				if sample_choose_best:
					return result, arrow_mesh
				else:
					return results, arrow_mesh
			else:
				if sample_choose_best:
					return result
				else:
					return results

		else:

			if isinstance(origin, int):
				orig_index = origin
				orig_atom = group.atoms[origin]
				origin = orig_atom.position
			else:
				orig_index = None
				orig_atom = None

			import numpy as np
			origin = np.array(origin)

			# guess the direction from nearby atoms
			if direction is None:
				direction = self._guess_vector_direction(orig_index, group, warnings=False)

				if direction is None:
					if return_arrow_mesh:
						return dict(success=False), None
					else:
						return dict(success=False)

			assert len(direction) == 3

			direction = np.array(direction, dtype=np.float64)
			direction /= np.linalg.norm(direction)
			

			result = self._classify_vector(origin, direction, fragments=False, warnings=False)

			if not result['success']:
				if return_arrow_mesh:
					return result, None
				else:
					return result

			### visualise
			if draw or return_arrow_mesh:

				from .o3d import mesh_from_AtomGroup
				
				arrow_mesh = create_arrow_mesh(orig_atom, orig_index, origin, color_by_distance, result)
				
				if draw:
				
					atoms_mesh = mesh_from_AtomGroup(group, use_covalent=True)
					
					if isinstance(atoms_mesh, list):
						extra = [arrow_mesh] + atoms_mesh
					else:
						extra = [arrow_mesh, atoms_mesh]

					self.render(extra=extra, pockets='hide', fragments='hide')

			if return_arrow_mesh:
				return result, arrow_mesh
			else:
				return result

	def trim(self):
		try: raise NotImplementedError; 
		except:
			logger.exception(f"{mcol.func}PoseButcher.trim{mcol.clear} method not implemented ")

	def score(self):
		try: raise NotImplementedError; 
		except:
			logger.exception(f"{mcol.func}PoseButcher.score{mcol.clear} method not implemented ")

	def render(self, hull='hide', **kwargs):
		self._render_meshes(hull=hull, **kwargs)

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
			from .o3d import mesh_from_pdb, material
			logger.warning('excuse the PyGAMer warnings... (they are safe to ignore)')
			mesh = mesh_from_pdb(self._fragment_bolus_path, gauss=False).to_legacy()
			mesh.compute_vertex_normals()
			mat = material(FRAGMENT_COLOR, alpha=1.0)
			self._fragment_mesh = dict(
				name='fragments',
				geometry=mesh,
				material=mat,
			)			

		return self._fragment_mesh

	@property
	def protein_mesh(self):
		if self._protein_mesh is None:
			logger.info('Generating protein mesh...')

			from .o3d import mesh_from_pdb, material
			
			mesh = mesh_from_pdb(self._apo_protein_path).to_legacy()

			mesh.compute_vertex_normals()

			self._protein_mesh = dict(
				name='protein',
				geometry=mesh,
				material=material(PROTEIN_COLOR, alpha=1.0),
			)
			
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
			from .o3d import convex_hull, material
			from copy import deepcopy
			mesh = deepcopy(self.protein_mesh['geometry'])
			mesh = convex_hull(mesh)
			mesh.compute_vertex_normals()
			self._protein_hull = dict(
				name = 'protein hull', 
				geometry = mesh,
				material = material(PROTEIN_COLOR, alpha=1.0, shiny=False),
			)

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

		if self._pocket_clip:
			self._clip_pockets()

	def _spherical_pocket_from_atoms(self, name, atoms, center=None, radius='mean', shift=None):

		import random
		from .o3d import sphere, material, glass#, subtract_atoms

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

		color = POCKET_COLORS[len(self.pocket_meshes)]
		mat = material(color, alpha=0.7)
		# mat = glass()

		# mesh.compute_triangle_normals()
		mesh.compute_vertex_normals()

		self._new_pocket(name, mesh, mat, radius=r)

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
		clash_radius = self._protein_clash_function(atom)
		return self._classify_position(atom.position, fragments, clash_radius=clash_radius)

	def _classify_position(self, position, fragments=True, clash_radius=0.0):
		
		from .o3d import is_point_in_mesh

		if fragments and is_point_in_mesh(self.fragment_mesh, position):
			return ('GOOD','fragment space')

		# protein clash
		if is_point_in_mesh(self.protein_mesh, position, within=clash_radius):
			return ('BAD','protein clash')

		# pockets
		for p_name, p_mesh in self.pockets.items():
			if is_point_in_mesh(p_mesh, position):
				return ('GOOD', 'pocket', p_name)

		return ('BAD','solvent space')

	def _classify_vector(self, origin, direction, fragments=False, warnings=True):

		# logger.debug(f"butcher._classify_vector({origin=}, {direction=})")

		d = 0.0

		from numpy.linalg import norm
		
		direction /= norm(direction)

		origin_info = self._classify_position(origin, fragments=fragments)
		# logger.info(f'origin: {origin_info}')

		result = dict(origin=origin_info, direction=direction)

		skip_pockets = []
		if origin_info[1] == 'protein clash':
			if warnings:
				logger.warning('Skipping vector that begins within protein volume')
			result['success'] = False
			return result
		elif origin_info[1] == 'pocket':
			skip_pockets.append(origin_info[2])

		complete = False

		result['intersections'] = {}

		while not complete:

			id_lookup, ans = self._cast_ray(origin, direction, skip_pockets)

			hit_distance = ans['t_hit'][0]
			
			from numpy import inf
			if hit_distance == inf:
				logger.error('NO HIT')

			else:

				geo_id = ans['geometry_ids'].numpy()[0]
				object_name = id_lookup[geo_id]

				hit_distance = round(float(hit_distance.numpy()),3)
				
				match object_name:
					case 'hull':
						complete = True
						result['intersections'][hit_distance] = ('BAD', 'solvent space')

					case 'protein':
						complete = True
						result['intersections'][hit_distance] = ('BAD', 'protein clash')

					case _:
						pocket_name = object_name.split()[-1]
						skip_pockets.append(pocket_name)
						result['intersections'][hit_distance] = ('GOOD', 'pocket', pocket_name)

		result['first_intersection_distance'] = min(result['intersections'])
		result['new_pocket'] = 'pocket' in result['intersections'][result['first_intersection_distance']]
		result['last_intersection_distance'] = max(result['intersections'])
		result['destination'] = result['intersections'][result['last_intersection_distance']][1]
		result['max_atoms_added'] = self._num_heavy_atoms_from_distance(result['last_intersection_distance'])

		if result['last_intersection_distance'] < CC_DIST:
			if warnings:
				logger.warning('Vector ends within carbon-carbon bond distance')
			result['success'] = False
			return result

		result['success'] = True

		return result

	def _cast_ray(self, origin, direction, skip_pockets):

		# create a ray casting scene
		from open3d.t.geometry import RaycastingScene, TriangleMesh
		scene = RaycastingScene()

		# populate the ray casting scene
		id_lookup = {}

		geo_id = scene.add_triangles(TriangleMesh.from_legacy(self.protein_mesh['geometry']))
		id_lookup[geo_id] = 'protein'

		geo_id = scene.add_triangles(TriangleMesh.from_legacy(self.protein_hull['geometry']))
		id_lookup[geo_id] = 'hull'
		
		for pocket in self.pocket_meshes:

			# skip starting pocket
			if pocket['name'] in skip_pockets:
				continue

			geo_id = scene.add_triangles(pocket['geometry'])
			id_lookup[geo_id] = f"pocket {pocket['name']}"

		# cast the ray
		from open3d.core import Tensor, Dtype
		try:
			rays = Tensor([[*origin, *direction]], dtype=Dtype.Float32)
		except ValueError:
			logger.error(f'{origin=}')
			logger.error(f'{direction=}')
			logger.error(f'{[[*origin, *direction]]=}')
			raise

		return id_lookup, scene.cast_rays(rays)

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

	def _output_to_color_pairs(self, output):

		pairs = []
		for k,v in output.items():

			if v[1] == 'pocket':
				c = self.pockets[v[2]]['material'].base_color
				c = tuple([float(x) for x in c[:3]])

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

	def _num_heavy_atoms_from_distance(self, distance):

		from numpy import inf

		if distance <= CC_DIST + C_VDW_RADIUS/2:
			return 1

		elif distance <= 1.5*CC_DIST + C_VDW_RADIUS/2:
			return 3

		elif distance <= 2.5*CC_DIST + C_VDW_RADIUS/2:
			return 7

		elif distance <= 3.0*CC_DIST + C_VDW_RADIUS/2:
			return 14

		else:
			return inf

	def _guess_vector_direction(self, orig_index, group, warnings=True):

		assert orig_index is not None


		bonds = [pair for pair in group.guess_bonds() if orig_index in pair]
		bonded_atom_indices = {i for i in sum(bonds,[]) if i != orig_index}

		if len(bonded_atom_indices) > 2:
			if warnings:
				logger.warning(f'Skipping vector from atom with more than two bonds. i={orig_index}')
			return None

		nearby_atom_center = sum([group.atoms[i].np_pos for i in bonded_atom_indices])
		nearby_atom_center /= len(bonded_atom_indices)

		origin = group.atoms[orig_index].np_pos
		direction = origin - nearby_atom_center

		return direction

def output_to_label(output, index):

	output_tuple = output[index]

	if output_tuple[1] == 'pocket':
		return f'{output_tuple[2]}'
	
	if output_tuple[1] == 'solvent space':
		return 'SOL.'
	
	if output_tuple[1] == 'protein clash':
		return 'PROT.'

	return ''

def circular_samples(theta_max, theta_step):
	import numpy as np

	def radians(x):
		return np.pi*x/180
	
	coords = []
	for theta in np.arange(theta_step,theta_max+theta_step,theta_step):

		r = np.tan(radians(theta))
		c = 2 * np.pi * r
		N = round(c/np.tan(radians(theta_step)))
		phi_step = 360/N

		for phi in np.arange(0, 360, phi_step):
			x = r*np.cos(radians(phi))
			y = r*np.sin(radians(phi))
			coords.append([x,y])

	return coords
    
def create_arrow_mesh(orig_atom, orig_index, origin, color_by_distance, result):

	from .o3d import arrow

	direction = result['direction']
	dist_to_clash = result['last_intersection_distance']
	first_intersection = result['first_intersection_distance']

	detail_str = f'{orig_atom} [i={orig_index}]' if orig_index else origin

	if color_by_distance:

		if dist_to_clash <= CC_DIST + C_VDW_RADIUS/2: # one heavy atom
			color = [1,0,0]

		elif dist_to_clash <= 1.5*CC_DIST + C_VDW_RADIUS/2: # three heavy atoms
			color = [1,0.5,0]

		elif dist_to_clash <= 2.5*CC_DIST + C_VDW_RADIUS/2: # almost whole ring
			color = [1,1,0]

		elif dist_to_clash <= 3.0*CC_DIST + C_VDW_RADIUS/2: # whole ring
			color = [0.5,1,0]

		else: # more heavy atoms
			color = [0,1,0]

	else:

		# color by clash/destination

		if result['intersections'][first_intersection][1] == 'pocket':
			color = [0,1,0]

		elif result['intersections'][dist_to_clash][1] == 'solvent space':
			color = [0,0,1]

		elif result['intersections'][dist_to_clash][1] == 'protein clash':
			color = [1,0,0]

		else:
			color = [0,0,0]

	return dict(geometry=arrow(origin, direction, length=dist_to_clash, radius=0.2, color=color), name=f'Vector {detail_str}')

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

CC_DIST = 1.54
C_VDW_RADIUS = 1.7