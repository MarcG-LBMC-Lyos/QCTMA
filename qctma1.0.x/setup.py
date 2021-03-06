from setuptools import setup


setup(name='qctma',
      version='1.0.8',
      description="Injects material (Young's modulus) to each element, based on a Dicom stack, and gray level to Young's"
                  "modulus relationships. Specifically designed to be used with Ansys .cdb meshes.",
      long_description="Injects material (Young's modulus) to each element, based on a Dicom stack, and gray level to Young's"
                       "modulus relationships. Specifically designed to be used with Ansys .cdb meshes.",
      url='https://github.com/MarcG-LBMC-Lyos/QCTMA',
      author='Marc Gardegaront',
      author_email='m.gardegaront@gmail.com',
      license='GNU GPLv3',
      py_modules=['qctma', 'rw_cdb'],
      install_requires=['ansys-dpf-core>=0.2.1', 'ansys-grpc-dpf>=0.2.2', 'ansys-mapdl-reader>=0.50.7',
                        'appdirs>=1.4.4', 'cachetools>=4.2.1', 'certifi>=2020.12.5', 'chardet>=4.0.0', 'cycler>=0.10.0',
                        'google-api-core>=1.26.1', 'google-api-python-client>=2.0.2', 'google-auth>=1.28.0',
                        'google-auth-httplib2>=0.1.0', 'googleapis-common-protos>=1.53.0', 'grpcio>=1.36.1',
                        'httplib2>=0.19.0', 'idna>=2.10', 'imageio>=2.9.0', 'importlib-metadata>=3.7.3',
                        'kiwisolver>=1.3.1', 'matplotlib>=3.3.4', 'meshio>=4.3.11', 'mpmath>=1.2.1', 'ndim>=0.1.4',
                        'numpy>=1.20.1', 'orthopy>=0.8.5', 'packaging>=20.9', 'pexpect>=4.8.0', 'Pillow>=8.1.2',
                        'protobuf>=3.15.6', 'ptyprocess>=0.7.0', 'pyasn1>=0.4.8', 'pyasn1-modules>=0.2.8',
                        'pydicom>=2.1.2', 'pyparsing>=2.4.7', 'python-dateutil>=2.8.1', 'pytz>=2021.1',
                        'pyvista>=0.29.0', 'quadpy>=0.16.6', 'reportlab>=3.5.66', 'requests>=2.25.1', 'rsa>=4.7.2', 'scipy>=1.6.1',
                        'scooby>=0.5.6', 'six>=1.15.0', 'sympy>=1.7.1', 'tqdm>=4.59.0', 'transforms3d>=0.3.1',
                        'typing-extensions>=3.7.4.3', 'uritemplate>=3.0.1', 'urllib3>=1.26.4', 'vtk>=9.0.1',
                        'zipp>=3.4.1'
                        ],
      python_requires=">=3.6")