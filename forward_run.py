from thetis import *
import time as time_mod
from model_config import *
import argparse
import os


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
output_dir = f"{pwd}/outputsB_dt0.01_alpha0.01"
if suffix != "":
    output_dir = "_".join([output_dir, suffix])

# Read in bathymetry (created by running bath.py)
pwd = os.path.abspath(os.path.dirname(__file__))
with CheckpointFile(f"{pwd}/monai_bathymetry_B.h5", "r") as f:
    mesh2d = f.load_mesh("firedrake_default")
    bathymetry = f.load_function(mesh2d, "bathymetry2d")

# Solve forward
solver_obj, update_forcings = construct_solver(
    bathymetry, 
    store_station_time_series=not no_exports,
    output_directory=output_dir,
    no_exports=no_exports,
)
print_output(f"Exporting to {solver_obj.options.output_directory}")
tic = time_mod.perf_counter()
solver_obj.iterate(update_forcings=update_forcings)
toc = time_mod.perf_counter()
print_output(f"Total duration: {toc-tic:.2f} seconds")
