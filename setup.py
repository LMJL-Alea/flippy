import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="flippy", # Replace with your own username
    version="0.0.0.9000",
    author="Aymeric Stamm",
    author_email="aymeric.stamm@math.cnrs.fr",
    description="A flexible permutation framework for making inference such as point estimation, confidence intervals or hypothesis testing, on any kind of data, be it univariate, multivariate, or more complex such as network-valued data, topological data, functional data or density-valued data.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://astamm.github.io/flippy/",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: CC - BY NC",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
