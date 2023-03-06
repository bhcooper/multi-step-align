from setuptools import setup, find_packages
import os

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="multi-step-align tools",
    version="1.0",
    author="Brendon Cooper",
    author_email="bhcooper@usc.edu",
    description="A set of scripts used to align full length SELEX-seq reads based on the protocol described in Cooper et al 2023 (in revision). The package also includes several scripts to validate and extend the results as described in the manuscript.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/bhcooper/multi-step-align",
    packages=find_packages(),
    scripts = [
        f"src/{script}" for script in os.listdir("src")
        if script.endswith(".py") or script.endswith(".R")
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
    install_requires=[
        "numpy",
        "pandas",
        "regex",
        "sharedmem",
        "matplotlib",
        "PyYAML",
        "scipy",
        "scikit-learn",
        "biopython",
        "logomaker",
        "h5py",
        "tensorflow",
        "TopDownCrawl"
    ],
)