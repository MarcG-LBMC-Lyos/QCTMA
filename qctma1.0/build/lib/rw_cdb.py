import numpy as np


def read_cdbfile(path, type='Tet'):
    """
    Extract the elements number and their associated material from a cdb mesh file.
    :param path: Path to the cdb file.
    :return: Elements number array and associated materials, Nodes and associated coordinates x, y and z arrays.
    """

    # EXTRACTING THE TABLES OF CONNECTION AND COORDINATES
    materials = []  # Material property number
    elems = []  # Table of connection
    nodes = []  # Table of coordinates
    x = []  # x-coordinates
    y = []  # y-coordinates
    z = []  # z-coordinates
    with open(path, 'r', errors="ignore") as f:
        # Open the mesh_file_path file to extract the table of connection and the table of coordinates.
        extract_nodes = False
        extract_elems = False
        line = f.readline()
        while line:
            if extract_nodes:
                if line.find("-1,") != -1:
                    extract_nodes = False
                    line = f.readline()
                    continue
                nodes.append(int(line[:NODE_LEN]))
                x.append(float(line[NODE_PROPERTIES_NB * NODE_LEN:3 * NODE_LEN + COORD_LEN]))
                y.append(float(line[NODE_PROPERTIES_NB * NODE_LEN + COORD_LEN:3 * NODE_LEN + COORD_LEN * 2]))
                z.append(float(line[NODE_PROPERTIES_NB * NODE_LEN + COORD_LEN * 2:3 * NODE_LEN + COORD_LEN * 3]))
                line = f.readline()
                continue
            if extract_elems:
                if line.find("-1") != -1:
                    extract_elems = False
                    line = f.readline()
                    continue
                if IS2LINES:
                    line += f.readline()
                line = line.replace("\n", '')
                materials.append(int(line[:ELEM_LEN]))
                line_elem = [int(line[i * ELEM_LEN:(i + 1) * ELEM_LEN]) for i in
                             range(ELEM_START, len(line) // ELEM_LEN)]
                elems.append(line_elem)
                line = f.readline()
                continue

            # DETECTING TABLE OF COORDINATE (NODES)
            if line.find("NBLOCK") != -1:
                line = f.readline()
                NODE_LEN = int(line.split(',')[0].replace('(', '').split('i')[1])
                NODE_PROPERTIES_NB = 3  # line.split(',')[0].replace('(', '').split('i')[0]
                COORD_LEN = int(line.split(',')[1].split('.')[0].split('e')[1])
                extract_nodes = True
                line = f.readline()
                continue

            # DETECTING TABLE OF CONNECTION (ELEMENTS)
            if line.find("EBLOCK") != -1:
                line = f.readline()
                ELEM_LEN = int(line.split(',')[0].replace(')', '').split('i')[1])
                ELEM_START = 10  # Table of connection doesn't start before the 10th value
                extract_elems = True
                line = f.readline()
                IS2LINES = False
                curs_pos = f.tell()
                sec_line = f.readline()
                if len(line) != len(sec_line):
                    IS2LINES = True
                f.seek(curs_pos)
                continue

            line = f.readline()

    elems = np.array(elems)
    materials = np.array(materials)
    nodes = np.array(nodes)
    x = np.array(x)
    y = np.array(y)
    z = np.array(z)

    return elems, materials, nodes, x, y, z

def get_density(path):
    return

def get_E(path):
    return

def get_nu(path):
    return

def write_cdb_mat(self, source_mesh_path, save_mesh_path, matid, e_pool, density_pool):
    """
    Write a new CDB mesh file based on a source CDB file, with materials defined assigned to each element.
    :param save_mesh_path: Path to the CDB mesh file that will be saved.
    :param matid: Material number assigned to each element
    :param e_pool: Pool of Young's modulus
    :param density_pool: Pool of density
    :return: None.
    """
    with open(source_mesh_path, 'r') as f:
        with open(save_mesh_path, 'w') as f_out:
            # Open the mesh_file_path file to extract the table of connection and the table of coordinates.
            extract_elems = False
            extract_mat = False
            line = f.readline()
            while line:
                if extract_elems:
                    if line.find("-1") != -1:
                        extract_elems = False
                        extract_mat = True
                        f_out.write(line)
                        continue

                    elem = int(line[ELEM_START * ELEM_LEN:(ELEM_START + 1) * ELEM_LEN])
                    mat_nb = int(matid[elem - 1])
                    line = "{:{width}}".format(mat_nb, width=ELEM_LEN) + line[ELEM_LEN:]

                    f_out.write(line)
                    if IS2LINES:
                        f_out.write(f.readline())
                    line = f.readline()
                    continue

                if extract_mat:
                    for i in range(len(e_pool)):
                        f_out.write("MPDATA,R5.0, 1,EX,%6i, 1, %.8f    ,\n" % (i + 1, e_pool[i]))
                        f_out.write("MPDATA,R5.0, 1,NUXY,%6i, 1, %.8f    ,\n" % (i + 1, 0.3))
                        f_out.write("MPDATA,R5.0, 1,DENS,%6i, 1, %.8f    ,\n" % (i + 1, density_pool[i]))
                    f_out.write("\n")
                    f_out.write("/GO\n")
                    f_out.write("FINISH\n")
                    break

                # DETECTING TABLE OF CONNECTION (ELEMENTS)
                if line.find("EBLOCK") != -1:
                    f_out.write(line)
                    line = f.readline()
                    ELEM_LEN = int(line.split(',')[0].replace(')', '').split('i')[1])
                    ELEM_START = 10  # Table of connection doesn't start before the 10th value
                    extract_elems = True
                    f_out.write(line)
                    line = f.readline()
                    IS2LINES = False
                    curs_pos = f.tell()
                    sec_line = f.readline()
                    if len(line) != len(sec_line):
                        IS2LINES = True
                    f.seek(curs_pos)
                    continue

                f_out.write(line)
                line = f.readline()