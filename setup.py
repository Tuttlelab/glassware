from setuptools import setup

from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
    name='glassware',
    version='0.0.1',
    description='Tuttlelab/Glassware',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/Tuttlelab/glassware',
    author='Swanson, Alexander van Teijlingen',
    author_email='alexander.van-teijlingen@strath.ac.uk',
    license='BSD 2-clause',
    packages=['hpctools', 'peptideutils', 'rmsd_tools'],
    install_requires=['ase',
                      'pandas',
                      'numpy',
                      'mdtraj',
                      'matplotlib',
                      'argparse'
                      ],

    classifiers=[
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3.10',
    ],
)
