import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# 2D 数据
grid_size_2d = np.array([0.1, 0.05, 0.025, 0.0125])
dt_2d = np.array([1214, 2604, 5412, 9909])
err_2d = np.array([2.91563e-05, 2.31078e-06, 2.88683e-07, 3.71515e-08])

# 3D 数据
grid_size_3d = np.array([0.4, 0.2, 0.1, 0.05])
dt_3d = np.array([2958, 3126, 7785, 19376])
err_3d = np.array([0.00360371, 0.000612675, 7.96446e-05, 9.55462e-06])

# 计算 x 轴 (grid size + dt)
x_values_2d = grid_size_2d + dt_2d
x_values_3d = grid_size_3d + dt_3d

# 绘制 2D 误差图
plt.figure(figsize=(8, 6))
plt.plot(x_values_2d, err_2d, marker='o', linestyle='-', color='b', label="2D Error")
plt.xlabel("Grid size + dt")
plt.xscale('log')
plt.ylabel("Error")
plt.yscale('log')  
plt.title("2D Error vs Grid size + dt (Second-Order Convergence)")
plt.legend()
plt.grid(True, which="both", linestyle='--', linewidth=0.5)
plt.show()

# 绘制 3D 误差图
plt.figure(figsize=(8, 6))
plt.plot(x_values_3d, err_3d, marker='o', linestyle='-', color='r', label="3D Error")
plt.xlabel("Grid size + dt")
plt.xscale('log')
plt.ylabel("Error")
plt.yscale('log')  
plt.title("3D Error vs Grid size + dt (Second-Order Convergence)")
plt.legend()
plt.grid(True, which="both", linestyle='--', linewidth=0.5)
plt.show()

# 计算 2D 误差阶数
data_2d = {
    "Grid Size (h)": grid_size_2d.tolist(),
    "Number of Cells": dt_2d.tolist(),
    "Error (e_n)": err_2d.tolist()
}
df_2d = pd.DataFrame(data_2d)
df_2d["Order of Accuracy (α)"] = [np.nan] + [
    np.log2(df_2d["Error (e_n)"].iloc[i-1] / df_2d["Error (e_n)"].iloc[i])
    for i in range(1, len(df_2d))
]
print("2D Error Table:")
print(df_2d.to_string(index=False))

# 计算 3D 误差阶数
data_3d = {
    "Grid Size (h)": grid_size_3d.tolist(),
    "Number of Cells": dt_3d.tolist(),
    "Error (e_n)": err_3d.tolist()
}
df_3d = pd.DataFrame(data_3d)
df_3d["Order of Accuracy (α)"] = [np.nan] + [
    np.log2(df_3d["Error (e_n)"].iloc[i-1] / df_3d["Error (e_n)"].iloc[i])
    for i in range(1, len(df_3d))
]
print("\n3D Error Table:")
print(df_3d.to_string(index=False))