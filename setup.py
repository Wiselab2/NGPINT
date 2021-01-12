import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="NGPINT", # Replace with your own username
    version="1.0.0",
    author="Sagnik Banerjee",
    author_email="sagnik@iastate.edu",
    description="NGPINT: A next-generation protein-protein interaction software",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Wiselab2/NGPINT/archive/NGPINTv1.0.0.tar.gz",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.7',
)
