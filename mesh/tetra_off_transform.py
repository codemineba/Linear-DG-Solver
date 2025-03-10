import meshio
import os

def off_transform(file_name):

    script_dir = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(script_dir, file_name)
    mesh = meshio.read(path)

    vertexs = mesh.points
    tetra = mesh.cells_dict['tetra']    

    ver_num = len(vertexs)
    tetra_num = len(tetra)

    new_file_path = os.path.splitext(path)[0] + ".off"

    with open(new_file_path, 'w') as f:
        f.write('OFF\n')
        f.write('{} {} {}\n'.format(ver_num, tetra_num, 0))
        for i in range(ver_num):
            f.write('{} {} {}\n'.format(vertexs[i, 0], vertexs[i, 1], vertexs[i, 2]))
        for i in range(tetra_num):
            f.write('{} {} {} {} {}\n'.format(4, tetra[i, 0], tetra[i, 1], tetra[i, 2], tetra[i, 3]))

        f.close()


off_transform('tetrahedronMesh0.025.msh')