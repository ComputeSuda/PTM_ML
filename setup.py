"""Setup and Install Script."""


from setuptools import setup


def get_readme():
    """Load README.rst for display on PyPI."""
    with open('README.md') as fhandle:
        return fhandle.read()


setup(
    name="PTM-ML",
    version="0.0.1b1",
    description="PTM-ML Package",
    long_description=get_readme(),
    long_description_content_type='text/markdown',
    # url="http://github.com/theochem/procrustes",
    license="GNU (Version 3)",
    author="QC-Devs Community",
    author_email="qcdevs@gmail.com",
    # package_dir={"PTMFun": "features"},
    package_dir={"PTM-ML_Prediction": "PTM-ML"},
    packages=["PTMFun"],
    install_requires=["numpy>=1.18.5", "scipy>=1.5.0", "pytest>=5.4.3", "sphinx>=2.3.0"],
)
