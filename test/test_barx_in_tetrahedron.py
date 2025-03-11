import gmsh
import os

def generate_tetrahedral_mesh(filename):
    # 获取当前脚本所在目录
    script_dir = os.path.dirname(os.path.fabspath(__file__))
    path = os.path.join(script_dir, filename)

    # 初始化 Gmsh
    gmsh.initialize()
    gmsh.model.add("TetrahedralMesh")

    # 创建立方体 [0, 2]^3
    box = gmsh.model.occ.addBox(0, 0, 0, 2, 2, 2)

    # 创建立方体的六个面
    gmsh.model.occ.synchronize()

    # 在六个面上设置一致分布
    divisions = 2  # 网格点划分数目，可调整

    # # 获取立方体的边并设置一致划分
    edges = gmsh.model.getEntities(dim=1)  # 获取所有 1D 实体（边）
    for edge in edges:
        gmsh.model.mesh.setTransfiniteCurve(edge[1], divisions)

    surfaces = gmsh.model.getEntities(dim=2)
    for surface in surfaces:
        gmsh.model.mesh.setTransfiniteSurface(surface[1], "Left")  # 设置平面为规则分布
    gmsh.model.mesh.setTransfiniteVolume(box)  # 设置整个立方体为结构化网格
    # # gmsh.model.mesh.setTransfiniteAutomatic()  # 自动优化一致性

    # 同步几何模型
    gmsh.model.occ.synchronize()


    # Define the points
    points = [
        # [0.33, 0.33, 1.33],
        # [0.33, 1.33, 0.33],
        # [1.33333, 0.333333, 0.333333],
        # [0.333333, 0.333333, 0],
        # [0.333333, 1.33333, 0],
        # [1.33333, 0.333333, 0],
        # [0.333333, 0, 1.33333],
        # [0.333333, 0, 0.333333],
        # [1.33333, 0, 0.333333],
        # [0, 0.333333, 1.33333],
        # [0, 0.333333, 0.333333],
        # [0, 1.33333, 0.333333],
        # [0.666667, 1.66667, 1.66667],
        # [1.66667, 0.666667, 1.66667],
        # [1.66667, 1.66667, 0.666667],
        # [0.666667, 1.66667, 2],
        # [1.66667, 0.666667, 2],
        # [1.66667, 1.66667, 2],
        # [2, 0.666667, 1.66667],
        # [2, 1.66667, 1.66667],
        # [2, 1.66667, 0.666667],
        [0.666667, 2, 1.66667],
        [1.66667, 2, 1.66667],
        [1.66667, 2, 0.666667],
        # [0.666667, 1.66667, 1.33333],
        # [1.66667, 0.666667, 0.333333],
        # [1.66667, 1.66667, 0.333333],
        # [0.333333, 1.66667, 1.33333],
        # [0.333333, 1.66667, 0.333333],
        # [1.33333, 0.666667, 0.333333],
        [0.333333, 2, 1.33333],
        [0.333333, 2, 0.333333],
        [1.33333, 2, 0.333333],
        # [0.666667, 1.66667, 0],
        # [1.66667, 0.666667, 0],
        # [1.66667, 1.66667, 0],
        # [0.666667, 1.33333, 1.66667],
        # [1.66667, 0.333333, 1.66667],
        # [1.66667, 0.333333, 0.666667],
        # [0.333333, 0.333333, 2],
        # [0.333333, 1.33333, 2],
        # [1.33333, 0.333333, 2],
        # [0.333333, 0.333333, 1.66667],
        # [0.333333, 1.33333, 1.66667],
        # [1.33333, 0.333333, 0.666667],
        # [0.666667, 0, 1.66667],
        # [1.66667, 0, 1.66667],
        # [1.66667, 0, 0.666667],
        # [0, 0.666667, 1.66667],
        # [0, 1.66667, 1.66667],
        # [0, 1.66667, 0.666667],
        # [2, 0.333333, 1.33333],
        # [2, 0.333333, 0.333333],
        # [2, 1.33333, 0.333333]
    ]

    # Add points to the gmsh model
    for i, (x, y, z) in enumerate(points, start=1):
        gmsh.model.geo.addPoint(x, y, z, 1.0, i)

    # gmsh.model.geo.addPoint(0.6, 0.6, 0.6, 1.0)
    # Synchronize and generate the model
    gmsh.model.geo.synchronize()

    # 生成三维四面体剖分网格
    gmsh.model.mesh.generate(3)
    
    # gmsh.write(path)

    # 如果需要直接在 GUI 中查看
    gmsh.fltk.run()

    # 关闭 Gmsh
    gmsh.finalize()

if __name__ == "__main__":
    generate_tetrahedral_mesh("tetrahedronMesh0.2.msh")
