import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="multimds",
    version="0.0.5",
    author="Lila Rieber",
    author_email="lr65358@gmail.com",
    description="Structural inference and alignment of Hi-C datasets",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/seqcode/multimds",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=[
          'scipy',
          'numpy',
          'pymp-pypi',
#          'mayavi',
          'scikit-learn',
          'h5py',
          'pandas',
          'seaborn'
    ]
)
