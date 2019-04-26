import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pyMolNetEnhancer",
    version="0.0.9",
    author="Madeleine Ernst, Ming Wang, Ricardo Silva",
    author_email="mernst@ucsd.edu",
    description="A python implementation of MolNetEnhancer",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/madeleineernst/pyMolNetEnhancer",
    packages=setuptools.find_packages(),
    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ),
    setup_requires=["pytest-runner"],
    tests_require=["pytest"],
)

