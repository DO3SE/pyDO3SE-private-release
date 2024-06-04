import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="thermal_time",
    version="0.0.15",
    author="SEI-York",
    author_email="sam.bland@york.ac.uk",
    description="Thermal time module",
    long_description=long_description,
    install_requires=[
        'numpy',
        'pandas',
        'deprecated',
        # PROFLOW FROM VENDERS DIR
        # 'proflow'
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
    package_dir={'thermal_time': 'thermal_time'},
    classifiers=[
        "Programming Language :: Python :: 3.7",
        "Operating System :: OS Independent",
    ],
)
