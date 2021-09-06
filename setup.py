import setuptools

with open("README.md", "r", encoding="utf-8") as reader:
    long_description = reader.read()

setuptools.setup(
    name="DEPhT",
    version="0.0.1",
    author="Christian Gauthier",
    author_email="christian.gauthier@pitt.edu",
    description="Discovery and Extraction of Phages Tool",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/chg60/DEPhT.git",
    classifiers=["Programming Langauge :: Python :: 3"],
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    python_requires=">=3.6",)
