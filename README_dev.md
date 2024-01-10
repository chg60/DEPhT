# pip build
```
python3 -m build .
twine upload dist/depht-{version}.*
```
# conda build
```
conda-build -c bioconda -c conda-forge -c laa89 .
```
# conda upload
```
# assuming anaconda-client is installed (e.g., conda install anaconda-client)...
anaconda login  # login with anaconda.org account credentials
anaconda upload /path/to/miniconda3/conda-bld/noarch/depht-x.y.z-py_0.tar.bz2
```
