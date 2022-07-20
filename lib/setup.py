# setup script
from setuptools import setup

setup(name='libspi',
      version='0.1',
      description='The SCRaMbLE Polymer Interaction model (SPI) library',
      url='https://github.com/stracquadaniolab/spi-nf/',
      author='Giovanni Stracquadanio',
      author_email='stracquadaniolab@gmail.com',
      license='GPL',
      packages=['libspi'],
      install_requires=[
          'numpy',
          'scipy',
      ],
      zip_safe=False
)
