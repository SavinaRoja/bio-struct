"""
This is the setup file to install relevant bio-struct code.
"""

from distutils.core import setup, Extension

#c_modules = Extension('', sources = [''])

setup (name = 'bio-struct',
       version = '0.1',
       description = 'This installs bio-struct code.',
       author='Paul Barton',
       #packages = [''],
       #ext_modules = [c_module]
       )
