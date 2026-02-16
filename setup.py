import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="do3se_phenology",
    version="0.17.5",
    author="SEI-York",
    author_email="sam.bland@york.ac.uk",
    description="DO3SE model python API",
    install_requires=[
        'wheel',
        'scipy',
        'numpy',
        'click',
        'pandas',
        'matplotlib',
        'deprecated',
        'data-helpers',
        'thermal_time @ git+ssh://git@github.com/SEI-DO3SE/thermal_time@RELEASE',
        'do3se_met @ git+ssh://git@github.com/SEI-DO3SE/do3se_met@RELEASE',
    ],
    setup_requires=[
        'wheel',
        'pytest-cov',
        'pytest-runner',
        'snapshottest'
    ],
    tests_require=['pytest', 'numpy', 'pandas', 'matplotlib', 'pytest-benchmark'],
    extras_require={'test': ['pytest', 'numpy', 'pandas']},
    packages=setuptools.find_packages(),
    package_dir={'do3se_phenology': 'do3se_phenology'},
    classifiers=[
        "Programming Language :: Python :: 3.7",
        "Operating System :: OS Independent",
    ],
)
