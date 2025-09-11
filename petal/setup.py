from distutils.core import setup
from Cython.Build import cythonize
#import os

for module in ['poscollider']:
    # for item in os.listdir():
    #     for ext_types_to_remove in ['.so', '.pyd', '.c', '.html']:
    #         if module in item and os.path.splitext(item)[-1] in ext_types_to_remove:
    #             os.remove(item)
    setup(name=module, ext_modules=cythonize(module + '.pyx', annotate=True))

# General tips for the Cython-uninitiated...
# ------------------------------------------
# sample call to compile at general console:
# python setup.py build_ext --inplace clean --all
#
# sample call to compile in Spyder console (generally, do so within a fresh console):
# runfile('setup.py', args='build_ext --inplace clean --all')
#
# for Windows compiling, need build tools for visual studio:
# https://visualstudio.microsoft.com/downloads/#build-tools-for-visual-studio-2017
# also see: https://wiki.python.org/moin/WindowsCompilers
