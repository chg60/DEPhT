# pip build
python3 -m build .
twine upload dist/depht-{version}.*

# conda build
conda-build -c bioconda -c conda-forge -c laa89 .
