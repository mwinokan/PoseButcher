
import open3d as o3d
import mout
import numpy as np

import mcol

import logging
logger = logging.getLogger("PoseButcher")

def sphere(radius=1.0, position=None, resolution=20, legacy=False):
	mesh = o3d.geometry.TriangleMesh.create_sphere(radius, resolution)
	if position is not None:
		mesh.translate(position)
	if not legacy:
		mesh = o3d.t.geometry.TriangleMesh.from_legacy(mesh)
	return mesh

def arrow(origin, direction, length, radius=0.5, color=[0,1,0]):

	mesh = o3d.geometry.TriangleMesh.create_arrow(
		cone_height=3*radius, 
		cylinder_height=length-3*radius, 
		cylinder_radius=radius, 
		cone_radius=1.5*radius,
	)

	# rotate
	rotation_matrix = rotation_matrix_from_vectors([0,0,1], direction)
	mesh.rotate(rotation_matrix, center=[0,0,0])

	# translate
	mesh.translate(origin)

	mesh.compute_triangle_normals()
	paint(mesh, color)

	return mesh

def cylinder(origin, direction, length=None, radius=0.25, legacy=False, resolution=20):

	origin = np.array(origin)
	direction = np.array(direction)
	# end = np.array(direction)

	# print(origin, direction)

	if length is not None:
		direction /= np.linalg.norm(direction)
		end = origin + direction*length
	else:
		end = np.array(direction)
		direction = end - origin
		length = np.linalg.norm(direction)
		direction /= length

	mesh = o3d.geometry.TriangleMesh.create_cylinder(
		radius=radius, 
		height=length,
		resolution=resolution,
	)

	# translate
	mesh.translate([0,0,length/2])

	# rotate
	rotation_matrix = rotation_matrix_from_vectors([0,0,1], direction)
	mesh.rotate(rotation_matrix, center=[0,0,0])

	# translate
	mesh.translate(origin)

	if not legacy:
		mesh = o3d.t.geometry.TriangleMesh.from_legacy(mesh)

	return mesh

def rotation_matrix_from_vectors(vec1, vec2):
    """ Find the rotation matrix that aligns vec1 to vec2
    :param vec1: A 3d "source" vector
    :param vec2: A 3d "destination" vector
    :return mat: A transform matrix (3x3) which when applied to vec1, aligns it with vec2.
    """
    a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)
    kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
    return np.array(rotation_matrix)

# render mesh or list of meshes, or list of dictionaries with keys 'name' and 'geometry'
def render(mesh, raw_mode=False, show_ui=True, wireframe=False):

	if wireframe:

		line_sets = []

		if not isinstance(mesh, list):
			mesh = [mesh]

		for m in mesh:
			if isinstance(m, dict):
				m = m['geometry']

			line_sets.append(o3d.geometry.LineSet.create_from_triangle_mesh(m))

		o3d.visualization.draw(line_sets)
	
	else:
		
		o3d.visualization.draw(mesh, raw_mode=raw_mode, show_ui=show_ui, show_skybox=False)

def paint(mesh, colour):
	if isinstance(mesh, dict):
		mesh['geometry'].paint_uniform_color(colour)
	else:
		mesh.paint_uniform_color(colour)

def material(color, alpha=1.0, shiny=False):

	mat = o3d.visualization.rendering.MaterialRecord()
	
	mat.base_color = [
		color[0],
		color[1],
		color[2],
		alpha,
	]
	
	
	if alpha < 1.0:
		mat.shader = 'defaultLitTransparency'
		mat.has_alpha = True

		# # non-rough materials are shiny
		# mat.base_roughness = 0.5

		# # low-reflectance == matte
		# mat.base_reflectance = 0.2
		
		# # # adds specular reflections
		# # mat.base_clearcoat = 0.1

		# # mat.thickness = 1.0
		# mat.transmission = 1.0 - alpha
		# mat.absorption_distance = 10
		# mat.absorption_color = color
	
	else:
		mat.shader = 'defaultLit'

	if shiny:
		mat.base_roughness = 0.0
		mat.base_reflectance = 0.0
		# mat.base_clearcoat = 1.0 - alpha
		mat.base_clearcoat = 0.1
		mat.thickness = 1.0
		mat.transmission = 1.0 - alpha
		mat.absorption_distance = 10
		mat.absorption_color = color

	return mat

def glass():
	mat = o3d.visualization.rendering.MaterialRecord()
	mat.shader = 'defaultLitSSR'
	mat.base_color = [0.467, 0.467, 0.467, 0.2]
	mat.base_roughness = 0.0
	mat.base_reflectance = 0.0
	mat.base_clearcoat = 1.0
	mat.thickness = 1.0
	mat.transmission = 1.0
	mat.absorption_distance = 10
	mat.absorption_color = [0.5, 0.5, 0.5]
	return mat

def compute_vertex_normals(mesh):
	if isinstance(mesh, dict):
		if hasattr(mesh['geometry'], 'compute_vertex_normals'):
			mesh['geometry'].compute_vertex_normals() 
	elif hasattr(mesh, 'compute_vertex_normals'):
			mesh.compute_vertex_normals()

def union(mesh1, mesh2):
	mesh = mesh1.boolean_union(mesh2)
	return mesh

def is_point_in_mesh(mesh, point, within=0.0, tolerance=0.1):
	return bool(signed_distance(mesh,point)[0] < tolerance + within)

def signed_distance(mesh, point):

	query_point = o3d.core.Tensor([point], dtype=o3d.core.Dtype.Float32)

	scene = o3d.t.geometry.RaycastingScene()

	if isinstance(mesh, dict):
		mesh = mesh['geometry']

	if not isinstance(mesh, o3d.t.geometry.TriangleMesh):
		mesh = o3d.t.geometry.TriangleMesh.from_legacy(mesh)
	_ = scene.add_triangles(mesh)

	return scene.compute_signed_distance(query_point)

def combine_many(meshes):
	combined = None
	for i,mesh in enumerate(meshes):
		
		if i == 0:
			combined = mesh
			continue
		
		combined += mesh

	return combined

def convex_hull(mesh):
	d = None
	if isinstance(mesh, dict):
		d = mesh
		mesh = d['geometry']

	if d:
		d['geometry'] = mesh.compute_convex_hull()[0]
		return d
	else:
		mesh = mesh.compute_convex_hull()[0]
		return mesh

def tensor_from_legacy(mesh):
	return o3d.t.geometry.TriangleMesh.from_legacy(mesh)

def subtract_atoms(mesh, group, r_scale=1.0, use_covalent=False):

	for i,atom in enumerate(group.atoms):

		mout.progress(i, group.num_atoms, prepend='subtracting')

		if use_covalent:
			r = atom.covalent_radius
		else:
			r = atom.vdw_radius

		mesh = mesh.boolean_difference(sphere(r, atom.position), tolerance=0.1)

	mout.finish()

	return mesh


### THIS IS NOT WORKING RELIABLY
def union_mesh_from_atoms(atoms, skip_hydrogen=True):

	raise NotImplementedError

	if skip_hydrogen:
		left = [a for a in atoms if a.species != 'H']
	else:
		left = list(atoms)
	used = []
	failed = []

	def p_r(atom):
		return atom.np_pos, atom.vdw_radius

	atom = left.pop()
	p,r = p_r(atom)
	combined = sphere(r,p, legacy=False)
	used.append((p,r))

	while left:

		mout.debug(f'#left={len(left)}')
		mout.debug(f'#used={len(used)}')

		found = False

		for atom in left:
			
			print(f'checking {atom}')

			p_test, r_test = p_r(atom)

			for (p,r) in used:
				
				print(f'checking used...')

				d = np.linalg.norm(p_test - p)

				diff = r_test + r - d

				if diff > 0.3:

					# intersection

					# left.pop(atom)
					mout.success(f'Found {atom}!')
					mout.var('p_test',p_test)
					mout.var('r_test',r_test)

					found = atom
					
					s = sphere(r_test, p_test, legacy=False)

					if len(used) == 5:
						render([combined, s], wireframe=False)

					try:
						new_combined = union(combined, s)

						if new_combined.to_legacy().has_triangles:
							combined = new_combined
							used.append((p_test, r_test))
						else:
							mout.error('Union has no triangles!')

					except IndexError as e:
						raise e
					
					break


			if found:
				mout.debug(f'Removing {atom}!')
				left = [a for a in left if a != found]
				break

	return combined

####

def mesh_from_AtomGroup(group, r_scale=1.0, use_covalent=False, licorice=True, licorice_radius=0.3):

	if licorice:
		
		from ase.data import vdw_radii, atomic_numbers
		from ase.data.colors import jmol_colors

		# params
		use_union=False
		resolution=40

		bond_pairs = group.guess_bonds()
		atom_positions = group.positions

		meshes = []

		symbols = group.present_symbols

		for symbol in symbols:
			
			symbol_meshes = []

			index_atom_dict = {i:a for i,a in enumerate(group.atoms) if a.symbol == symbol}

			atomic_number = atomic_numbers[symbol]
			color = jmol_colors[atomic_number]
			color = (color[0],color[1],color[2])

			for index, atom in index_atom_dict.items():
				
				atom_meshes = []

				# atom spheres
				mesh = sphere(licorice_radius, atom.position, legacy=not use_union, resolution=resolution)
				atom_meshes.append(mesh)

				# bonds
				bonds = [pair for pair in bond_pairs if index in pair]
				for i,j in bonds:
					a = atom_positions[i]
					b = atom_positions[j]
					center = (a+b)/2
					mesh = cylinder(center, atom.position, radius=licorice_radius, legacy=not use_union, resolution=resolution)
					atom_meshes.append(mesh)

				if use_union:
					# merge shapes with union
					combined = atom_meshes[0]
					for mesh in atom_meshes[1:]:
						combined = combined.boolean_union(mesh)
					combined = combined.to_legacy()
					combined.compute_triangle_normals()
					paint(combined,color)
					meshes.append(combined)
				else:
					symbol_meshes += atom_meshes

			if not use_union:
				combined = combine_many(symbol_meshes)
				combined.compute_triangle_normals()
				paint(combined,color)
				meshes.append(dict(geometry=combined,name=f'{group.name} {symbol}'))

		return meshes

	else:

		meshes = []
		for atom in group.atoms:
			if use_covalent:
				r = atom.covalent_radius
			else:
				r = atom.vdw_radius
			meshes.append(sphere(r*r_scale, atom.np_pos, legacy=True))

		combined = combine_many(meshes)

	return {'name':group.name, 'geometry':combined}

def mesh_from_pdb(path, gauss=False):

	import pygamer

	if gauss:
		logger.debug('pygamer.readPDB_gauss...')
		mesh = pygamer.readPDB_gauss(
				path, 
				blobbyness= -0.8,
				isovalue= 2.5,
			)

	else:
		logger.debug('pygamer.readPDB_molsurf...')
		mesh = pygamer.readPDB_molsurf(path)

	logger.debug('mesh.compute_orientation...')
	components, orientable, manifold = mesh.compute_orientation()
	mesh.correctNormals()

	logger.debug('creating tensors...')
	protverts, protedges, protfaces = mesh.to_ndarray()
	vertices = o3d.core.Tensor(protverts, dtype=o3d.core.Dtype.Float32)
	triangles = o3d.core.Tensor(protfaces, dtype=o3d.core.Dtype.Int32)

	logger.debug('creating mesh...')
	mesh = o3d.t.geometry.TriangleMesh()
	mesh.vertex.positions = vertices
	mesh.triangle.indices = triangles

	logger.success(f'Mesh {mcol.file}{path}{mcol.success} complete')
	return mesh

####

def test_sphere():
	mesh = sphere()
	render(mesh)

def test_2spheres():
	
	mesh1 = sphere()
	
	mesh2 = sphere(0.5, [0.75,0,0])

	mesh = union(mesh1, mesh2)

	render(mesh)

def test_3spheres():
	
	mesh1 = sphere()
	
	mesh2 = sphere(0.5, [0.75,0,0])

	mesh3 = sphere(0.5, [0,0.75,0])

	mesh = union_many([mesh1, mesh2, mesh3])

	render(mesh)

