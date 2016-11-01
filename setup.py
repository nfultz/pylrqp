
try:
        from setuptools import setup, Extension
except ImportError:
        from distutils.core import setup, Extension
from glob import glob
import os


#shamefully borrowed from cvxopt

# Modifiy this if BLAS and LAPACK libraries are not in /usr/lib.
BLAS_LIB_DIR = '/usr/lib'

# Default names of BLAS and LAPACK libraries
BLAS_LIB = ['blas']
LAPACK_LIB = ['lapack']
BLAS_EXTRA_LINK_ARGS = []

# Set environment variable BLAS_NOUNDERSCORES=1 if your BLAS/LAPACK do
# not use trailing underscores
BLAS_NOUNDERSCORES = False


# No modifications should be needed below this line.

BLAS_NOUNDERSCORES = int(os.environ.get("CVXOPT_BLAS_NOUNDERSCORES",BLAS_NOUNDERSCORES)) == True
BLAS_LIB = os.environ.get("CVXOPT_BLAS_LIB",BLAS_LIB)
LAPACK_LIB = os.environ.get("CVXOPT_LAPACK_LIB",LAPACK_LIB)
BLAS_LIB_DIR = os.environ.get("CVXOPT_BLAS_LIB_DIR",BLAS_LIB_DIR)
BLAS_EXTRA_LINK_ARGS = os.environ.get("CVXOPT_BLAS_EXTRA_LINK_ARGS",BLAS_EXTRA_LINK_ARGS)
if type(BLAS_LIB) is str: BLAS_LIB = BLAS_LIB.strip().split(',')
if type(LAPACK_LIB) is str: LAPACK_LIB = LAPACK_LIB.strip().split(',')
if type(BLAS_EXTRA_LINK_ARGS) is str: BLAS_EXTRA_LINK_ARGS = BLAS_EXTRA_LINK_ARGS.strip().split(',')


# Macros
MACROS = []
if BLAS_NOUNDERSCORES: MACROS.append(('BLAS_NO_UNDERSCORE',''))


base = Extension('_LowRankQP', libraries = ['m'] + LAPACK_LIB + BLAS_LIB,
    library_dirs = [ BLAS_LIB_DIR ],
    define_macros = MACROS,
    extra_link_args = BLAS_EXTRA_LINK_ARGS,
    sources = ['src/C/LowRankQP.c','src/C/LowRankQPmodule.c']) 

extmods = [base]

setup (name = 'pylrqp', 
    description = 'Low Rank Quadratic Programming',
    version = '1.0.2', 
    long_description = '''
This is a port of the orphanded LowRankQP R package.
See https://cran.r-project.org/web/packages/LowRankQP/index.html
''', 
    author = 'N. Fultz',
    author_email = 'nfultz@openmail.co',
    url = 'http://github.com/nfultz/pylrqp',
    license = 'GNU GPL version 3',
    ext_package = "pylrqp",
    ext_modules = extmods,
    package_dir = {"pylrqp": "src/python"},
    packages = ["pylrqp"],
    setup_requires=['numpy'],
    classifiers=[
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: OS Independent',
        'Programming Language :: C',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering',
        ],
    )
