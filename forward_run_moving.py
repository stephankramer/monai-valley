from thetis import *
import time as time_mod
from model_config import *
import argparse
import os

from moving import MeshMovementCallback

# Parse user input
parser = argparse.ArgumentParser(
    description="Monai Valley tsunami propagation",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)
parser.add_argument("--suffix", type=str, default="")
parser.add_argument("--load", action="store_true")
args = parser.parse_args()
suffix = args.suffix
no_exports = os.getenv("THETIS_REGRESSION_TEST") is not None
pwd = os.path.abspath(os.path.dirname(__file__))
output_dir = f"{pwd}/outputs_move"
if suffix != "":
    output_dir = "_".join([output_dir, suffix])

mesh = RectangleMesh(50, 30, 5.488, 3.402)
V = FunctionSpace(mesh, "CG", 1)
bathymetry2d = Function(V, name="bathmetry")
bathymetry_interpolator = BathymetryInterpolator(bathymetry2d)

# Solve forward
solver_obj, update_forcings = construct_solver(
    bathymetry2d, 
    store_station_time_series=not no_exports,
    output_directory=output_dir,
    no_exports=no_exports,
    left_id=1,
)

mmc = MeshMovementCallback(solver_obj, interpolation_cb=bathymetry_interpolator.interpolate)
dexp = DepthExporter(solver_obj)
solver_obj.add_callback(mmc, "timestep")
solver_obj.add_callback(dexp, "export")
print_output(f"Exporting to {solver_obj.options.output_directory}")
tic = time_mod.perf_counter()
solver_obj.iterate(update_forcings=update_forcings)
toc = time_mod.perf_counter()
print_output(f"Total duration: {toc-tic:.2f} seconds")
