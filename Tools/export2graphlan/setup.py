import setuptools
from setuptools.command.install import install
from io import open
import os

install_requires = ["pandas==0.13.1"]
setuptools.setup(
    name='export2graphlan',
    version='0.22',
    author='Francesco Asnicar',
    author_email='f.asnicar@unitn.it',
    url='http://github.com/segatalab/export2graphlan',
    packages = setuptools.find_packages(),
    scripts=['export2graphlan.py'],
    package_dir = {'export2graphlan' : '' },
    long_description_content_type='text/markdown',
    long_description=open('README.md').read(),
    description='export2graphlan is a conversion software tool for producing both annotation and tree file for GraPhlAn',
    install_requires=install_requires
)
