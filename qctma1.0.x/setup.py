from setuptools import setup


setup(name='qctma',
      version='1.0.11',
      description="Injects material (Young's modulus) to each element, based on a Dicom stack, and gray level to Young's"
                  "modulus relationships. Specifically designed to be used with Ansys .cdb meshes.",
      long_description="Injects material (Young's modulus) to each element, based on a Dicom stack, and gray level to Young's"
                       "modulus relationships. Specifically designed to be used with Ansys .cdb meshes.",
      url='https://github.com/MarcG-LBMC-Lyos/QCTMA',
      author='Marc Gardegaront',
      author_email='m.gardegaront@gmail.com',
      license='GNU GPLv3',
      py_modules=['qctma', 'rw_cdb'],
      install_requires=['matplotlib>=3.3.4', 'numpy>=1.20.1', 'pydicom>=2.1.2', 'quadpy>=0.16.6', 'scipy>=1.6.1',
                        'reportlab>=3.5.66'
                        ],
      python_requires=">=3.6")