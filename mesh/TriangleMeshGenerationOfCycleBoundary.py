import gmsh

gmsh.initialize()  # 初始化
gmsh.model.add("model1")  # 创建模型

lc = 0.1  # 设置网格尺寸值

# 创建点
gmsh.model.geo.addPoint(0, 0, 0, lc)
gmsh.model.geo.addPoint(2, 0, 0, lc)
gmsh.model.geo.addPoint(2, 2, 0, lc)
gmsh.model.geo.addPoint(0, 2, 0, lc)


# 创建线
gmsh.model.geo.addLine(1, 2, 1)
gmsh.model.geo.addLine(3, 2, 2)
gmsh.model.geo.addLine(3, 4, 3)
gmsh.model.geo.addLine(4, 1, 4)
# 创建边
gmsh.model.geo.addCurveLoop([4, 1, -2, 3], 1)
# 创建面
gmsh.model.geo.addPlaneSurface([1])

gmsh.model.geo.synchronize()  # 同步到模型

gmsh.model.mesh.generate(2)  # 生成网格

gmsh.write("triangleMesh0.1.off")  # 生成.msh文件


gmsh.fltk.run()  # 图形界面显示
gmsh.finalize()  # 关闭gmsh