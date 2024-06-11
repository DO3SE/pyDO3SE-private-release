import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="do3se_met",
    version="0.3.4",
    author="SEI-York",
    author_email="sam.bland@york.ac.uk",
    description="DO3SE model python API",
    install_requires=[
        'scipy',
        'numpy',
        'click',
        'deprecated',
        'scipy',
    ],
    setup_requires=[
        'wheel',
        'pytest-cov',
        'pytest-runner',
        'snapshottest',
        'deprecated'
    ],
    tests_require=['pytest', 'numpy','pytest-benchmark'],
    extras_require={'test': ['pytest', 'numpy']},
    packages=setuptools.find_packages(),
    package_dir={'do3se_met': 'do3se_met'},
    classifiers=[
        "Programming Language :: Python :: 3.7",
        "Operating System :: OS Independent",
    ],
)
