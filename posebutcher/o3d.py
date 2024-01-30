
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
		
		if isinstance(mesh, list):
			for m in mesh:
				compute_vertex_normals(m)
		else:
			compute_vertex_normals(mesh)
		
		o3d.visualization.draw(mesh, raw_mode=raw_mode, show_ui=show_ui)

def paint(mesh, colour):
	if isinstance(mesh, dict):
		mesh['geometry'].paint_uniform_color(colour)
	else:
		mesh.paint_uniform_color(colour)

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
						# render([combined, sphere(r_test, p_test, legacy=False)], wireframe=True)
						# render([combined], wireframe=True)
						# render([s], wireframe=False)
						raise e
					
					break

			# if len(used) == 6:
			# 	left = 0
			# 	break

			if found:
				mout.debug(f'Removing {atom}!')
				left = [a for a in left if a != found]
				break

	return combined

####

def mesh_from_AtomGroup(group, r_scale=1.0, use_covalent=False):
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

