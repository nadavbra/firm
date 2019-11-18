from __future__ import absolute_import, division, print_function

import sys

from setuptools import setup
from Cython.Build import cythonize
import numpy as np
        
def readme():
    with open('README.rst', 'r') as f:
        return f.read()

setup(
    name = 'firm',
    version = '1.1.3',
    description = 'FIRM (Functional Impact Rating at the Molecular-level) is a machine-learning model for predicting the functional impact ' + \
            'of genetic variants.',
    long_description = readme(),
    long_description_content_type = 'text/markdown',
    url = 'https://github.com/nadavbra/firm',
    author = 'Nadav Brandes',
    author_email = 'nadav.brandes@mail.huji.ac.il',
    license = 'MIT',
    packages = [
        'firm',
        'firm.ml',
    ],
    package_data = {'firm': [
        '_apply_scale.pyx',
        'data/classifier-py2.pkl',
        'data/classifier-py3.pkl',
    ]},
    install_requires = [
        'numpy',
        'pandas',
        'biopython',
        'scikit-learn',
        'geneffect', # https://github.com/nadavbra/geneffect
    ],
    ext_modules = cythonize('firm/_apply_scale.pyx', compiler_directives = dict(language_level = sys.version_info.major)),
    include_dirs = [np.get_include()],
    zip_safe = False,
)
