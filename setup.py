"""Setup.py

NOTE: This is kept to enable `pip install -e`. See PEP 660.

"""
import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pyDO3SE4",
    version="4.40.1",
    author="SEI-York",
    author_email="sam.bland@york.ac.uk",
    description="DO3SE model python API",
    install_requires=[
        'scipy',
        'numpy',
        'click',
        'pandas',
        'matplotlib',
        'deprecated',
        'proflow',
        'data_helpers',
        'thermal_time @ git+ssh://git@github.com/SEI-DO3SE/thermal_time@RELEASE',
        'do3se_phenology @ git+ssh://git@github.com/SEI-DO3SE/do3se_phenology@RELEASE',
        'do3se_met @ git+ssh://git@github.com/SEI-DO3SE/do3se_met@RELEASE',
    ],
    setup_requires=[
        'pytest-cov',
        'pytest-runner',
        'snapshottest'
    ],
    tests_require=['pytest', 'numpy', 'pandas', 'matplotlib', 'pytest-benchmark'],
    extras_require={'test': ['pytest', 'numpy', 'pandas']},
    packages=setuptools.find_packages(),
    package_dir={'pyDO3SE': 'pyDO3SE', 'pyDO3SE_cli': 'pyDO3SE_cli'},
    classifiers=[
        "Programming Language :: Python :: 3.7",
        "Operating System :: OS Independent",
    ],
)
