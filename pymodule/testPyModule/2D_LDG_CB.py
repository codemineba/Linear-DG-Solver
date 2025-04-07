import PY_DG_2D

mesh = PY_DG_2D.TriangleMesh()

path = '../../data/2D/CycleBoundary/'
file = "triangleMesh0.1.off"

mesh.read_off(path+file)
mesh.collect_edges()

DGSolver = PY_DG_2D.LinearDGSolver_2D_CycleBoundary(mesh)
DGSolver.computeTimeDiscretization(2)

print('iter over!')