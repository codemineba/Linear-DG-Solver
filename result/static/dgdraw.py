import matplotlib.pyplot as plt
import numpy as np

# 给定数据

grid_size = np.array([0.1, 0.05, 0.025, 0.0125])
dt = np.array([0.002, 0.001, 0.0005, 0.00025])
err = np.array([0.00166714, 0.00042522, 0.000108002, 0.0000270313])

# 计算 lc + dt 作为横坐标
x_values = grid_size + dt

# 创建图表
plt.figure(figsize=(8, 6))
plt.plot(x_values, err, marker='o', linestyle='-', color='b', label="Error")

# 设置坐标轴标签
plt.xlabel("Grid size + dt")
plt.xscale('log')
plt.ylabel("Error")
plt.yscale('log')  # 使用对数刻度展示二阶收敛精度
plt.title("Error vs Grid size  + dt (Second-Order Convergence)")

# 显示网格和图例
plt.grid(True, which="both", linestyle='--', linewidth=0.5)
plt.legend()
plt.show()


import pandas as pd
import numpy as np

# 给定数据
data = {
    "Grid Size (h)": [0.1, 0.05, 0.025, 0.0125],
    "Time Step (dt)": [0.002, 0.001, 0.0005, 0.00025],
    "Iterations (n)": [500, 1000, 2000, 4000],
    "Error (e_n)": [0.00166714, 0.00042522, 0.000108002, 0.0000270313]
}

# 创建 DataFrame
df = pd.DataFrame(data)

# 计算精度阶数 (α) 使用误差值
df["Order of Accuracy (α)"] = [np.nan] + [
    np.log2(df["Error (e_n)"].iloc[i-1] / df["Error (e_n)"].iloc[i])
    for i in range(1, len(df))
]

# 显示表格，保留 6 位小数
print(df.round(6))