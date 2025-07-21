import setuptools
import re

def get_version():
    with open("amrscan/amrscan.py", "r") as f:
        for line in f:
            match = re.search(r'__version__ = "(.*?)"', line)
            if match:
                return match.group(1)
    raise RuntimeError("Version could not be found.")

setuptools.setup(
    name="amrscan-pipeline",
    version=get_version(),
    author="H. Soon Gweon",
    author_email="h.s.gweon@reading.ac.uk",
    description="A Bioinformatics Tool for AMR gene and variant detection from metagenomic data.",
    packages=setuptools.find_packages(),
    
    include_package_data=True,
    package_data={
        'amrscan': ['databases/**/*']
    },

    entry_points={
        'console_scripts': [
            'amrscan = amrscan.amrscan:main',
        ],
    },
    
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires='>=3.7',
)