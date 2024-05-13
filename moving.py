from thetis import *
from movement import MongeAmpereMover
from animate import RiemannianMetric

class MeshMovementCallback:
    name = 'meshmovement'
    def __init__(self, solver, interpolation_cb=None):
        self.solver = solver
        self.interpolation_cb = interpolation_cb

        self.tpvd = File('tst.pvd')
        self.alpha = 4
        def monitor(mesh):
            V = FunctionSpace(mesh, "CG", 1)
            elev = Function(V, name='elev')
            elev.project(self.solver.fields.elev_2d)
            TV = TensorFunctionSpace(mesh, "CG", 1)
            hbc = DirichletBC(TV, 0, "on_boundary")
            Helev = RiemannianMetric(TV, name="Helev")
            Helev.compute_hessian(elev, method="L2")
            Hnorm = Function(V, name="Hnorm")
            Hnorm.interpolate(sqrt(inner(Helev, Helev)))
            Hnorm_max = Hnorm.dat.data.max()
            if Hnorm_max > 0.0:
                m = 1 + self.alpha * Hnorm / Hnorm_max
            else:
                m = 1
            mf = Function(V, name='monitor')
            mf.interpolate(m)
            self.tpvd.write(mf, elev, Helev)
            return m

        self.mover = MongeAmpereMover(self.solver.mesh2d, monitor, method="relaxation",
                rtol=1e-2)


    def evaluate(self, index):
        print("before move: xi[2]:", self.mover.xi.dat.data[2])
        self.mover.move()

        V = FunctionSpace(self.mover.mesh, self.solver.function_spaces.H_2d.ufl_element())
        elev = Function(V, name='elev')
        elev.project(self.solver.fields.elev_2d)
        self.solver.fields.elev_2d.dat.data[:] = elev.dat.data[:]

        W = FunctionSpace(self.mover.mesh, self.solver.function_spaces.U_2d.ufl_element())
        uv = Function(W, name='uv')
        uv.project(self.solver.fields.uv_2d)
        self.solver.fields.uv_2d.dat.data[:] = uv.dat.data[:]

        if self.interpolation_cb:
            self.interpolation_cb()

        self.solver.mesh2d.coordinates.dat.data[:] = self.mover.mesh.coordinates.dat.data[:]
        self.solver.mesh2d.clear_spatial_index()
