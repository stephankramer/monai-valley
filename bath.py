from thetis import *
import model_config

mesh = Mesh('mesh/MonaiValley_B.msh')
#mesh = RectangleMesh(50, 30, 5.488, 3.402)
V = FunctionSpace(mesh, "CG", 1)
bathymetry2d = Function(V, name='bathymetry2d')
model_config.interpolate_bathymetry(bathymetry2d)
File('bathymetry.pvd').write(bathymetry2d)
with CheckpointFile('monai_bathymetry_B.h5', 'w') as afile:
    afile.save_mesh(mesh)
    afile.save_function(bathymetry2d)

