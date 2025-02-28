import torch

# # 计算四面体体积的函数（与之前的C++代码类似）
# def directed_3D_tetrahedron_volume(v0, v1, v2, v3):
#     x0, y0, z0 = v0[0], v0[1], v0[2]
#     x1, y1, z1 = v1[0], v1[1], v1[2]
#     x2, y2, z2 = v2[0], v2[1], v2[2]
#     x3, y3, z3 = v3[0], v3[1], v3[2]

#     volume = (1.0 / 6.0) * (
#         x0 * (y1 * z2 + y2 * z3 + y3 * z1 - y1 * z3 - y2 * z1 - y3 * z2) -
#         y0 * (x1 * z2 + x2 * z3 + x3 * z1 - x1 * z3 - x2 * z1 - x3 * z2) +
#         z0 * (x1 * y2 + x2 * y3 + x3 * y1 - x1 * y3 - x2 * y1 - x3 * y2) -
#         (x1 * y2 * z3 + x2 * y3 * z1 + x3 * y1 * z2 - x1 * y3 * z2 - x2 * y1 * z3 - x3 * y2 * z1)
#     )
#     return -volume


# 计算带方向的二维三角形面积
def directed_2D_triangle_area(v0, v1, v2):
    x0, y0 = v0[0], v0[1]
    x1, y1 = v1[0], v1[1]
    x2, y2 = v2[0], v2[1]
    # 计算带方向的面积
    area = 0.5 * (x0 * (y1 - y2) + x1 * (y2 - y0) + x2 * (y0 - y1))
    return area

# 计算带方向的四面体体积（基于三角形面积的计算）
def directed_3D_tetrahedron_volume(v0, v1, v2, v3):
    v0_2d = v0[:2]
    v1_2d = v1[:2]
    v2_2d = v2[:2]
    v3_2d = v3[:2]
    
    z0, z1, z2, z3 = v0[2], v1[2], v2[2], v3[2]
    
    volume = (1.0 / 3.0) * (
        -z0 * directed_2D_triangle_area(v1_2d, v2_2d, v3_2d)
        +z1 * directed_2D_triangle_area(v0_2d, v2_2d, v3_2d)
        -z2 * directed_2D_triangle_area(v0_2d, v1_2d, v3_2d)
        +z3 * directed_2D_triangle_area(v0_2d, v1_2d, v2_2d)
    )
    return volume

# 创建顶点坐标的张量，并启用梯度
v0 = torch.tensor([0.0, 0.0, 0.0], dtype=torch.float32, requires_grad=True)
v1 = torch.tensor([1.0, 0.0, 0.0], dtype=torch.float32, requires_grad=True)
v2 = torch.tensor([0.0, 1.0, 0.0], dtype=torch.float32, requires_grad=True)
v3 = torch.tensor([0.0, 0.0, 1.0], dtype=torch.float32, requires_grad=True)

# v0 = torch.tensor([0.0, 0.0, 0.0], dtype=torch.float32, requires_grad=True)
# v1 = torch.tensor([100.0, 0.0, 0.0], dtype=torch.float32, requires_grad=True)
# v2 = torch.tensor([0.0, 100.0, 0.0], dtype=torch.float32, requires_grad=True)
# v3 = torch.tensor([0.0, 100.0, 1.0], dtype=torch.float32, requires_grad=True)

# v0 = torch.tensor([10.0, 10.0, 10.0], dtype=torch.float32, requires_grad=True)
# v1 = torch.tensor([20.0, 10.0, 10.0], dtype=torch.float32, requires_grad=True)
# v2 = torch.tensor([10.0, 20.0, 10.0], dtype=torch.float32, requires_grad=True)
# v3 = torch.tensor([10.0, 10.0, 30.0], dtype=torch.float32, requires_grad=True)

# 计算体积
volume = directed_3D_tetrahedron_volume(v0, v1, v2, v3)

# 计算梯度
volume.backward()

print(volume)
# 打印梯度（关于每个顶点的体积梯度）
print(f"Gradient w.r.t. v0: {v0.grad}")
print(f"Gradient w.r.t. v1: {v1.grad}")
print(f"Gradient w.r.t. v2: {v2.grad}")
print(f"Gradient w.r.t. v3: {v3.grad}")