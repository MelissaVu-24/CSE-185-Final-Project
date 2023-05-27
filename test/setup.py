import os
from setuptools import setup, find_packages

curdir = os.path.abspath(os.path.dirname(__file__))
MAJ = 0
MIN = 0
REV = 0
VERSION = '%d.%d.%d' % (MAJ, MIN, REV)
with open(os.path.join(curdir, 'mypileup/version.py'), 'w') as fout:
        fout.write(
            "\n".join(["",
                       "# THIS FILE IS GENERATED FROM SETUP.PY",
                       "version = '{version}'",
                       "__version__ = version"]).format(version=VERSION)
        )


setup(
    name='diffana',
    version='0.0.0',
    description='CSE185 Project: Diffana',
    author='Owin Gong and Melissa Vu',
    author_email='mlvu@ucsd.edu',
    packages=find_packages(),
    entry_points={
        "console_scripts": [
            "diffana=diffana.diffana:main"
        ],
    },
)
