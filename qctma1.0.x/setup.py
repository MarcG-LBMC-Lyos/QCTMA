from setuptools import setup


setup(name='qctma',
      version='1.0.12',
      description="Injects material (Young's modulus) to each element, based on a Dicom stack, and gray level to Young's"
                  "modulus relationships. Specifically designed to be used with Ansys .cdb meshes.",
      long_description="Injects material (Young's modulus) to each element, based on a Dicom stack, and gray level to Young's"
                       "modulus relationships. Specifically designed to be used with Ansys .cdb meshes.",
      url='https://github.com/MarcG-LBMC-Lyos/QCTMA',
      author='Marc Gardegaront',
      author_email='m.gardegaront@gmail.com',
      license='GNU GPLv3',
      py_modules=['qctma', 'rw_cdb'],
      install_requires=['matplotlib>=2.2.5', 'numpy>=1.19.5', 'pydicom>=2.1.2', 'quadpy>=0.16.5', 'scipy>=1.5.4',
                        'reportlab>=3.5.66'
                        ],
      python_requires=">=3.6")