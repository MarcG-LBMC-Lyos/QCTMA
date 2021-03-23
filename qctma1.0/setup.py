from setuptools import setup


setup(name='qctma',
      version='1.0',
      description="Injects material (Young's modulus) to each element, based on a Dicom stack, and gray level to Young's"
                  "modulus relationships. Specifically designed to be used with Ansys .cdb meshes.",
      url='https://github.com/MarcG-LBMC-Lyos/QCTMA',
      author='Marc Gardegaront',
      author_email='m.gardegaront@gmail.com',
      license='GNU GPLv3',
      py_modules=['qctma', 'rw_cdb'],
      install_requires=["pydicom", "numpy", "scipy", "quadpy", "pyansys"],
      python_requires=">=3.6")