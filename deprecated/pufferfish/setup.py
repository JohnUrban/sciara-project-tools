import os
import sys
from setuptools import setup

version_py = os.path.join(os.path.dirname(__file__), 'pufferfish', 'version.py')
version = open(version_py).read().strip().split('=')[-1].replace('"','')
long_description = """
``pufferfish`` HMM-based approach(es) to finding and analyzing developmentally regulated amplicons (genomic sites that are programmed to increase in copy number over time).'
"""

with open("requirements.txt", "r") as f:
    install_requires = [x.strip() for x in f.readlines()]

setup(
        name="pufferfish",
        version=version,
        install_requires=install_requires,
        requires = ['python (>=2.7, <3.0)'],
        packages=['pufferfish',
                  'pufferfish.scripts'],
        author="John Urban",
        description='HMM-based approach(es) to finding and analyzing developmentally regulated amplicons (genomic sites that are programmed to increase in copy number over time).',
        long_description=long_description,
        url="https://github.com/JohnUrban/poreminion",
        package_dir = {'pufferfish': "pufferfish"},
        package_data = {'pufferfish': []},
        zip_safe = False,
        include_package_data=True,
        #scripts = ['pufferfish/scripts/pufferfish'],
        entry_points = {
            'console_scripts' : [
                 'pufferfish = pufferfish.pufferfish_main:main', 
            ],
        },  
        author_email="mr.john.urban@gmail.com",
        classifiers=[
            'Development Status :: 4 - Beta',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: GNU General Public License (GPL)',
            'Topic :: Scientific/Engineering :: Bio-Informatics']
    )
