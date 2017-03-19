#!/usr/bin/env python

from distutils.core import setup

setup(name='genomon_complex_mutation',
      version='0.1.0',
      description='Python tools for detecting complex mutations',
      author='Yuichi Shiraishi',
      author_email='friend1ws@gamil.com',
      url='https://github.com/friend1ws/GenomonComplexMutation.git',
      package_dir = {'': 'lib'},
      packages=['genomon_complex_mutation'],
      scripts=['genomon_complex_mutation'],
      license='GPL-3'
     )

