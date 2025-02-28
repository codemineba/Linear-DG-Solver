import gmsh
import math

# 初始化 Gmsh
gmsh.initialize()
gmsh.model.add("Square with Center Circle")

# 定义全局网格尺寸因子
global_mesh_factor = 0.1  # 全局网格尺寸因子，调整这个数值可以缩放网格密度

# 参数设置
rect_size = 16                     # 矩形尺寸 (16x16)
circle_radius = 0.5               # 圆的半径
mesh_size_factor = global_mesh_factor  # 网格尺寸因子
num_points_on_circle = int(30 * (0.1 / global_mesh_factor))  # 根据 global_mesh_factor 调整圆上点的数量

# 定义矩形的四个顶点坐标（使用 OCC 内核）
half_size = rect_size / 2
points = [
    gmsh.model.occ.addPoint(-half_size, -half_size, 0,global_mesh_factor*4),
    gmsh.model.occ.addPoint(half_size, -half_size, 0,global_mesh_factor*4),
    gmsh.model.occ.addPoint(half_size, half_size, 0,global_mesh_factor*4),
    gmsh.model.occ.addPoint(-half_size, half_size, 0,global_mesh_factor*4)
]

# 创建矩形的四条边
lines = [
    gmsh.model.occ.addLine(points[0], points[1]),
    gmsh.model.occ.addLine(points[1], points[2]),
    gmsh.model.occ.addLine(points[2], points[3]),
    gmsh.model.occ.addLine(points[3], points[0])
]

# 创建矩形的闭合环面
rect_loop = gmsh.model.occ.addCurveLoop(lines)
rect_plane = gmsh.model.occ.addPlaneSurface([rect_loop])

# 创建近似圆的点
circle_points = []
for i in range(num_points_on_circle):
    angle = 2 * math.pi * i / num_points_on_circle  # 计算每个点的角度
    x = circle_radius * math.cos(angle)  # 计算点的 x 坐标
    y = circle_radius * math.sin(angle)  # 计算点的 y 坐标
    point = gmsh.model.occ.addPoint(x, y, 0)  # 在圆的边界上添加点
    circle_points.append(point)

# 通过点创建闭合圆形线段
circle_lines = []
for i in range(num_points_on_circle):
    start_point = circle_points[i]
    end_point = circle_points[(i + 1) % num_points_on_circle]
    line = gmsh.model.occ.addLine(start_point, end_point)
    circle_lines.append(line)

# 使用圆的线段生成封闭曲线环和面
circle_curve_loop = gmsh.model.occ.addCurveLoop(circle_lines)
circle_surface = gmsh.model.occ.addPlaneSurface([circle_curve_loop])

# 同步模型
gmsh.model.occ.synchronize()

# 定义圆与矩形的布尔减操作
result = gmsh.model.occ.cut([(2, rect_plane)], [(2, circle_surface)], removeTool=True)

# 同步模型
gmsh.model.occ.synchronize()

# 设置更密集的网格
gmsh.model.mesh.field.add("Distance", 1)
gmsh.model.mesh.field.setNumbers(1, "NodesList", circle_points)  # 使用圆周边界的点

gmsh.model.mesh.field.add("Threshold", 2)
gmsh.model.mesh.field.setNumber(2, "InField", 1)
gmsh.model.mesh.field.setNumber(2, "SizeMin", mesh_size_factor / 0.5)
gmsh.model.mesh.field.setNumber(2, "SizeMax", mesh_size_factor * 5)
gmsh.model.mesh.field.setNumber(2, "DistMin", circle_radius * 0.2)
gmsh.model.mesh.field.setNumber(2, "DistMax", circle_radius * 20)
gmsh.model.mesh.field.setAsBackgroundMesh(2)

# 生成网格
gmsh.model.mesh.generate(2)

# 可视化和保存
gmsh.write("rectangle_with_center_circle1.off")
gmsh.fltk.run()  # 运行图形界面（可选）

# 清理 Gmsh 资源
gmsh.finalize()
