from thetis import *
import netCDF4
import os
import scipy.interpolate as si


def interpolate_bathymetry(bathymetry_2d, nc_file_name="raw_data/Bathymetry.grd"):
    """
    Interpolate a bathymetry field from some data set.

    :arg bathymetry_2d: :class:`Function` to store the data in
    :kwarg nc_file_name: netcdf file with bathymetry
    """
    mesh = bathymetry_2d.function_space().mesh()

    # Read data from file
    with netCDF4.Dataset(nc_file_name, "r") as nc:
        interp = si.RectBivariateSpline(
            nc.variables["y"][:],
            nc.variables["x"][:],
            nc.variables["z"][:, :],
        )

    # Interpolate at mesh vertices
    xy = mesh.coordinates.dat.data[:]
    bathymetry_2d.dat.data[:] = -interp(xy[:,1], xy[:,0], grid=False)


def construct_solver(bathymetry_2d, left_id=10, input_wave_csv='raw_data/InputWave.csv',
        store_station_time_series=False, **model_options):
    """
    Construct a *linear* shallow water equation solver for tsunami
    propagation modelling.
    """
    mesh2d = bathymetry_2d.function_space().mesh()
    #bathymetry_2d.assign(bathymetry_2d + Constant(0.2))

    t_end = 25
    u_mag = Constant(5.0)
    t_export = 0.1
    dt = 0.01
    if os.getenv("THETIS_REGRESSION_TEST") is not None:
        t_end = 5 * t_export

    # Create solver
    solver_obj = solver2d.FlowSolver2d(mesh2d, bathymetry_2d)
    options = solver_obj.options
    options.use_nonlinear_equations = True
    options.element_family = "dg-dg"
    options.simulation_export_time = t_export
    options.fields_to_export = ["elev_2d", "uv_2d"]
    options.fields_to_export_hdf5 = ["elev_2d", "uv_2d"]
    options.simulation_end_time = t_end
    options.horizontal_velocity_scale = u_mag
    options.swe_timestepper_type = "CrankNicolson"
    options.use_wetting_and_drying = True
    options.wetting_and_drying_alpha = Constant(0.01)
    options.manning_drag_coefficient = Constant(0.019)

    if not hasattr(options.swe_timestepper_options, "use_automatic_timestep"):
        options.timestep = dt
    #options.swe_timestepper_options.solver_parameters = {
    #    "snes_type": "newtonls",
    #    "ksp_type": "gmres",
    #    "pc_type": "bjacobi",
    #    "sub_pc_type": "ilu",
    #}
    options.update(model_options)


    # Set up gauges
    if store_station_time_series:
        for x, y, name in [
         [4.521, 1.196, 'gauge1'],
         [4.521, 1.696, 'gauge2'],
         [4.521, 2.196, 'gauge3']]:
            cb = TimeSeriesCallback2D(
                solver_obj,
                ["elev_2d"],
                x, y,
                name,
                append_to_log=False,
            )
            solver_obj.add_callback(cb)

    # Set boundary conditions
    elev_at_boundary = Constant(0.0)
    solver_obj.bnd_functions["shallow_water"] = {
        left_id: {"elev": elev_at_boundary}
    }
    input_data = np.loadtxt(input_wave_csv)
    intp = si.interp1d(input_data[:, 0], input_data[:, 1], fill_value='extrapolate')
    def update_forcings(t):
        elev_at_boundary.assign(intp(t))
    update_forcings(0)

    solver_obj.create_equations()
    solver_obj.assign_initial_conditions(uv=Constant((1e-20, 0.)))
    return solver_obj, update_forcings
