'''
SSEA: Sample Set Enrichment Analysis
'''
import os
from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

import numpy
numpy_inc = numpy.get_include()

__author__ = "Matthew Iyer, Yashar Niknafs"
__copyright__ = "Copyright 2012-2017"
__credits__ = ["Matthew Iyer", "Yashar Niknafs"]
__license__ = "GPL"
__version__ = "0.1.0"
__maintainer__ = "Yashar Niknafs"
__email__ = "yniknafs@umich.edu"
__status__ = "Development"


ext_modules = [Extension('ssea.lib.kernel',
                         sources=['ssea/lib/kernel.pyx', 'ssea/lib/rng.c'],
                         include_dirs=[numpy_inc],
                         libraries=['m'],
                         extra_compile_args=['-std=c99']),
               Extension('ssea.lib.cfisher',
                         sources=['ssea/lib/cfisher.pyx'],
                         include_dirs=[numpy_inc],
                         extra_compile_args=['-std=c99'])]

extensions = [
    Extension('ssea.lib.ckernel',
              sources=['ssea/lib/ckernel.c', 'ssea/lib/rng.c'],
              extra_compile_args=['-w'],
              include_dirs=[numpy_inc]
             )
]

def main():
    setup(name='SSEA',
          version=__version__,
          description='Sample Set Enrichment Analysis',
          author=__author__,
          author_email=__email__,
          requires=['numpy', 'jinja2', 'cython'],
          license=__license__,
          url='https://github.com/yniknafs/ssea',
          ext_modules=cythonize(ext_modules) + extensions,
          packages=['ssea', 'ssea.lib'],
          package_data={'ssea.templates': ['details.html',
                                           'report.html']},
          scripts=['ssea/ssea', 'ssea/ssea_gui'])

if __name__ == '__main__':
    main()
