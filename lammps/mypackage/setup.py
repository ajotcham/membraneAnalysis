# setup.py

from setuptools import setup, find_packages

setup(
    name='mypackage',
    version='0.1',
    packages=find_packages(),
    install_requires=[],  # Add any dependencies here
    author='Your Name',
    author_email='your.email@example.com',
    description='A simple example package',
    url='https://github.com/yourusername/mypackage',  # URL to the package source
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
)