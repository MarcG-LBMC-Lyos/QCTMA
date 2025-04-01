import pyvista as pv
import pydicom as pd
import numpy as np
import os
import rw_cdb

cdb_path = r""  # Path to the cdb mesh file
dicomDir_path = r""  # Path to the dicom directory

def cdb2pvmesh(cdb_path: str) -> pv.UnstructuredGrid:
    """
    Load a cdb file and return a pyvista mesh.
    :param cdb_path: Path to the cdb file.
    :return: Pyvista mesh.
    """
    # Read the cdb file
    elems, materials, nodes, x, y, z = rw_cdb.read_cdbfile(cdb_path)
    # Read the field file, except first line
    # Create the mesh
    densities = rw_cdb.get_density(cdb_path)
    Es = rw_cdb.get_E(cdb_path)
    nodes_nb = len(elems[0]) - 1
    if nodes_nb == 8:
        ELEM_TYPE = pv.CellType.HEXAHEDRON  # Hexahedron mesh
    elif nodes_nb == 20:
        ELEM_TYPE = pv.CellType.QUADRATIC_HEXAHEDRON  # Hexahedron mesh
    elif nodes_nb == 4:
        ELEM_TYPE = pv.CellType.TETRA  # Hexahedron mesh
    elif nodes_nb == 10:
        ELEM_TYPE = pv.CellType.QUADRATIC_TETRA  # Hexahedron mesh
    node_args = elems[:, 1:] - 1  # Nodes nb must start at 1 and be continuous.
    cells = np.hstack((nodes_nb * np.ones((len(elems), 1)), node_args)).astype(int)
    celltypes = ELEM_TYPE * np.ones(len(cells))
    node_coords = np.array((x, y, z)).T
    pv_mesh = pv.UnstructuredGrid(cells, celltypes, node_coords)
    return pv_mesh

def load_dicom(dicomDir_path):
    origins = []
    volume = []
    orientation_z = None
    last_z_pos = None
    z_coords = []
    for dcm_file in os.listdir(dicomDir_path):
        if not dcm_file.lower().endswith(".dcm"):
            continue
        dcm = pd.dcmread(os.path.join(dicomDir_path, dcm_file))
        z_coords.append(dcm.ImagePositionPatient[2])
    order_dcm = np.argsort(z_coords)
    thicknesses = np.diff(np.array(z_coords)[order_dcm])
    # Look for most frequent value of thickness
    unique, counts = np.unique(thicknesses, return_counts=True)
    thickness = unique[np.argmax(counts)]
    for dcm_file in np.array(os.listdir(dicomDir_path))[order_dcm]:
        dcm = pd.dcmread(os.path.join(dicomDir_path, dcm_file))
        volume.append(dcm.pixel_array)
        origins.append(dcm[0x20, 0x32].value)
        if last_z_pos is not None and not orientation_z:
            orientation_z = np.sign(float(dcm.ImagePositionPatient[2]) - last_z_pos)
        if not orientation_z and not last_z_pos:
            last_z_pos = float(dcm.ImagePositionPatient[2])
    orientation = np.array(
        [float(dcm.ImageOrientationPatient[0]), float(dcm.ImageOrientationPatient[4]), orientation_z])
    origins = np.array(origins)
    origin = np.array(np.min(origins.T, axis=1))
    # Check if slice thickness is given in dicom
    try:
        thickness = float(dcm[0x18, 0x50].value)
    except:
        print("/!\ Slice thickness not found in dicom.")
        pass
    # Check if spacing between slices is given in dicom
    try:
        thickness = float(dcm[0x18, 0x88].value)
    except:
        print("/!\ Spacing between slices not found in dicom.")
        pass
    spacing = np.abs(np.array((float(dcm.PixelSpacing[0]), float(dcm.PixelSpacing[1]), thickness)))
    volume = np.array(volume).T
    dimensions = np.array(volume.shape) + 1
    # Rectifying orientation and origin of the volume if dicom's orientation is negative
    if orientation[0] < 0:
        volume = np.flip(volume, axis=0)
        origin[0] -= dimensions[0] * spacing[0]
    if orientation[1] < 0:
        volume = np.flip(volume, axis=1)
        origin[1] -= dimensions[1] * spacing[1]
    if orientation[2] < 0:
        volume = np.flip(volume, axis=2)

    return volume, dimensions, origin, spacing

def ct_scan2ImageData(dicomDir_path: str):
    volume, dimensions, origin, spacing = load_dicom(dicomDir_path)
    volume[volume < 0] = 0  # For better contrast
    grid = pv.ImageData()
    grid.dimensions = dimensions
    grid.origin = origin
    grid.spacing = spacing
    grid.cell_data["values"] = volume.flatten(order="F")
    return grid

if __name__ == "__main__":
    pv_mesh = cdb2pvmesh(cdb_path)  # Load the cdb file
    elems, materials, nodes, x, y, z = rw_cdb.read_cdbfile(cdb_path)  # Just want to extract the materials associated to each element
    Es = rw_cdb.get_E(cdb_path)  # Get the Young's modulus associated to each material
    scalar_field = list(map(lambda x: Es[x - 1], materials))  # Get the densities associated to each element
    pv_mesh.cell_data["mat"] = scalar_field  # Assign the densities to the mesh

    grid = ct_scan2ImageData(dicomDir_path)

    barycenter = pv_mesh.center

    plotter = pv.Plotter(shape=(1, 2))  # Create a plotter with 1 row and 2 columns
    # Left subplot for the CT scan
    plotter.subplot(0, 0)
    plotter.add_mesh(grid.slice_orthogonal(x=barycenter[0], y=barycenter[1], z=barycenter[2]), name="CT scan", cmap="gray")
    plotter.add_mesh(pv_mesh, color="red", opacity=0.5)
    # Right subplot for the FEM mesh
    plotter.subplot(0, 1)
    plotter.add_mesh(pv_mesh.slice_orthogonal(x=barycenter[0], y=barycenter[1], z=barycenter[2]), name="FEM mesh", scalars="mat")
    plotter.link_views()

    plotter.show()
