
import logging
logger = logging.getLogger("PoseButcher")

def read_ply(file):
	"""Read vtk polydata from a .ply file"""

	from vtk import vtkPolyDataMapper, vtkActor
	from vtkmodules.vtkIOPLY import vtkPLYReader

	reader = vtkPLYReader()
	reader.SetFileName(file)
	reader.Update()

	return reader.GetOutput()

def write_ply(file, polydata):
	"""Write vtk polydata to a .ply file"""
	from vtk.vtkIOPLY import vtkPLYWriter
	plyWriter = vtkPLYWriter()
	plyWriter.SetFileName(file)
	plyWriter.SetInputDataObject(polydata)
	plyWriter.Write()

def to_open3d(port):
	"""Convert some vtk poly data to an open3d triangle mesh"""

	import tempfile
	from .o3d import load_mesh

	with tempfile.NamedTemporaryFile(suffix='.ply', delete=False) as fp:
		write_ply(fp.name, port)
		mesh = load_mesh(fp.name, verbosity=False)

	return mesh

def from_open3d(o3d_mesh):
	"""Convert an open3d trianglemesh to some vtk poly data"""

	from .o3d import dump_mesh
	import tempfile

	with tempfile.NamedTemporaryFile(suffix='.ply', delete=False) as fp:
		dump_mesh(fp.name, o3d_mesh, write_ascii=True, compressed=True, verbosity=False)
		polydata = read_ply(fp.name)

	return polydata

def boolean_difference(m1, m2):
	"""Uses vtkBool to compute the boolean difference between two open3d triangle meshes"""

	import tempfile
	from .o3d import load_mesh
	from vtkbool.vtkBool import vtkPolyDataBooleanFilter
	from vtk.vtkIOPLY import vtkPLYWriter
	
	pd1 = from_open3d(m1)
	pd2 = from_open3d(m2)

	boolean = vtkPolyDataBooleanFilter()	
	boolean.SetInputDataObject(0, pd1)
	boolean.SetInputDataObject(1, pd2)
	boolean.SetOperModeToDifference()

	with tempfile.NamedTemporaryFile(suffix='.ply', delete=False) as fp:
		writer = vtkPLYWriter()
		writer.SetInputConnection(boolean.GetOutputPort())
		writer.SetFileName(fp.name)
		writer.Update()
		mesh = load_mesh(fp.name, verbosity=False)

	return mesh
