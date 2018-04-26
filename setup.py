from __future__ import absolute_import, division, print_function

from setuptools import setup
        
def readme():
    with open('README.rst', 'r') as f:
        return f.read()

setup(
    name = 'firm',
    version = '1.0',
    description = 'FIRM (Functional Impact Rating at the Molecular-level) is a machine-learning model for predicting the functional impact ' + \
            'of genetic variants.',
    long_description = readme(),
    url = 'https://github.com/nadavbra/firm',
    author = 'Nadav Brandes',
    author_email  ='nadav.brandes@mail.huji.ac.il',
    license = 'MIT',
    packages = [
        'firm',
        'firm.ml',
    ],
    package_data = {'firm': [
        '_apply_scale.pyx',
        'data/classifier.pkl',
    ]},
    install_requires = [
        'numpy',
        'pandas',
        'biopython',
        'scikit-learn',
        'geneffect==1.1', # https://github.com/nadavbra/geneffect
    ],
    zip_safe = False,
)
