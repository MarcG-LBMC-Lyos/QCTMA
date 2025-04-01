import numpy as np
import qctma


dicom_dir_path = r""  # Path to the folder containing the dicoms
source_mesh_path = r""  # Path to the source mesh
deltaE = 50  # Step of variation of Young's modulus
out_mesh_path = r""  # Path to the output mesh (must end with ".cdb")
smoothing = 0  # Gaussian smoothing parameter (0 for no smoothing)
window_smoothing = 0  # Kernel size for Gaussian smoothing
nb_process = 10  # Number of processes to use for processing

def gl2density(gl_array):
    """
    Function to convert gray levels to density.
    :param gl_array: List of gray levels
    :return: List of densities
    """
    rho = gl_array
    rho[rho <= 0] = 0.01
    return rho


def density2E(rho):
    """
    Function to convert density to Young's modulus.
    :param rho: List of densities
    :return: List of Young's moduli
    """
    E_mat = rho
    E_mat = np.array(E_mat)
    E_mat[E_mat < 0.01] = 0.01
    return E_mat


if __name__ == "__main__":  # !!! Do not forget to put the code in this condition to avoid errors due to multiprocessing
    qctma.qctma(dcm_path=dicom_dir_path, mesh_path=source_mesh_path, gl2density=gl2density,
                density2E=density2E, deltaE=deltaE, process=True, save_mesh_path=out_mesh_path,
                gaussian_smoothing=smoothing, window_gaussian_smoothing=window_smoothing,
                nb_process=nb_process)
    
    

