#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
setup the package
To do:
    [type the aims and objectives of this script]
Created on Wed Apr  1 15:03:58 2020

@author: Xiaodong Ming
"""
from setuptools import setup, find_packages
from os import path
here = path.abspath(path.dirname(__file__))
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
      name='hipims_io',
      version='0.2.9.1',
      description='To process input and output files of the HiPIMS model',
      url='https://github.com/xiaxilin/hipims/tree/master/hipims_io',
      author='Xiaodong Ming',
      author_email='xiaodong.ming@ncl.ac.uk',
      long_description=long_description,
      long_description_content_type='text/markdown',
      classifiers=['Development Status :: 3 - Alpha',
                   'Programming Language :: Python :: 3',
                   'License :: OSI Approved :: MIT License',
                   'Programming Language :: Python :: 3',
                   'Operating System :: OS Independent',],
      keywords='hipims model IO',
      license='LICENSE.txt',
      packages=find_packages(),
      include_package_data=True,
      package_data={'hipims_io': ['sample/Example_DEM.asc', 
                                  'sample/Example_data.npy'],},
#      install_requires=['gdal', 'scipy', 'pyshp'],
      python_requires='>=3.6')

"""
#Sometimes you’ll want to use packages that are properly arranged with 
#setuptools, but aren’t published to PyPI. In those cases, you can specify a 
#list of one or more dependency_links URLs where the package can be downloaded,
#along with some additional hints, and setuptools will find and install the
#package correctly.
#For example, if a library is published on GitHub, you can specify it like:
setup(
    ...
    dependency_links=['http://github.com/user/repo/tarball/master#egg=package-1.0']
    ...
)
"""