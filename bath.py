from thetis import *
import model_config

#mesh = Mesh('mesh/MonaiValley_D.msh')
mesh = RectangleMesh(50, 30, 5.488, 3.402)
V = FunctionSpace(mesh, "CG", 1)
bathymetry2d = Function(V, name='bathymetry2d')
bathymetry_interpolator = model_config.BathymetryInterpolator(bathymetry2d)
bathymetry_interpolator.interpolate()
File('bathymetry.pvd').write(bathymetry2d)
with CheckpointFile('monai_bathymetry_F.h5', 'w') as afile:
    afile.save_mesh(mesh)
    afile.save_function(bathymetry2d)

