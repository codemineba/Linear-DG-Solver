import PY_DG_3D

mesh = PY_DG_3D.TetrahedronMesh()

path = '../../data/3D/CycleBoundary'
file = "tetrahedronMesh0.1.off"

mesh.read_off(path+file)
mesh.collect_faces()

DGSolver = PY_DG_3D.LinearDGSolver_3D_CycleBoundary(mesh)
DGSolver.computeTimeDiscretization(2)

print('iter over!')