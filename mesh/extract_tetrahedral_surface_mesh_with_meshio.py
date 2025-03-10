import meshio
import os

def convert_msh_to_off(file_name, output_off_filenamne):
    # 读取 .msh 文件
    script_dir = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(script_dir, file_name)
    mesh = meshio.read(path)

    # 提取点 (points) 和三角形面 (triangle)
    points = mesh.points
    triangles = mesh.cells_dict.get("triangle")

    if triangles is None:
        raise ValueError("输入网格文件中没有 'triangle' 表面数据！")

    output_off_dict = os.path.join(script_dir, output_off_filenamne)
    # 写入 .off 文件
    with open(output_off_dict, "w") as off_file:
        # 写入 OFF 文件头
        off_file.write("OFF\n")
        off_file.write(f"{len(points)} {len(triangles)} 0\n")

        # 写入顶点坐标
        for point in points:
            off_file.write(f"{point[0]} {point[1]} {point[2]}\n")

        # 写入三角形面
        for triangle in triangles:
            off_file.write(f"3 {triangle[0]} {triangle[1]} {triangle[2]}\n")

    print(f"表面网格文件已成功保存至：{output_off_file}")


if __name__ == "__main__":
    # 示例：输入 .msh 文件和输出 .off 文件路径
    input_msh_file = "tetrahedronMesh.msh"
    output_off_file = "surface_mesh.off"

    convert_msh_to_off(input_msh_file, output_off_file)