from setuptools import setup, find_packages  # Always prefer setuptools over distutils
from codecs import open  # To use a consistent encoding
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the relevant file
with open(path.join(here, 'README'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='HumanCircRNADataParser',

    # Versions should comply with PEP440.  For a discussion on single-sourcing
    # the version across setup.py and the project code, see
    # http://packaging.python.org/en/latest/tutorial.html#version
    version='1.0',

    description='A Python project that wants to exist.',
    long_description=long_description,  #this is the

    # The project's main homepage.
    url='https://google.com',

    # Author details
    author='Richard Clayton',
    author_email='rich29@live.unc.edu',

    # Choose your license
    license='UNC',

    # See https://PyPI.python.org/PyPI?%3Aaction=list_classifiers
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 5 - Production/Stable',

        # Indicate who your project is intended for
        'Intended Audience :: Researchers',
        'Topic :: Genome Analysis :: Compare and Parse Tools',

        # Pick your license as you wish (should match "license" above)
        'License :: OSI Approved :: UNC License',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 2.7',
    ],

    # What does your project relate to?
    keywords='circRNA Parser ',

    packages=["hcrdp"],

)