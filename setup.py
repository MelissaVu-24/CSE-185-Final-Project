import os
from setuptools import setup, find_packages

setup(
    name='diffana',
    version=VERSION,
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
