from setuptools import setup, find_packages
from setuptools.command.install import install

setup(
    name="oyFlow",
    version="0.0.2",
    description="flow analysis code, OY lab @HMS",
    author="Alon Oyler-Yaniv",
    url="https://github.com/alonyan/oyFlow",
    packages=find_packages(include=["oyFlow", "oyFlow.*"]),
    python_requires=">=3.8",
    dependency_links=[
    ],
    install_requires=[
        "PyYAML",
        "PyQt5",
        "dill==0.3.4",
        "ipython==8.4.0",
        "ipywidgets==7.6.5",
        "matplotlib==3.8.0",
        "magicgui==0.7.3",
        "numpy==1.23.1",
        "pandas==1.4.2",
        "flowcytometrytools==0.5.1",
        "scipy==1.8.1",
        "setuptools==62.3.4",
        "anytree==2.8.0",
        "natsort==8.4.0",
        "QtAwesome==1.2.3"
    ],
    extras_require={
    },
)
