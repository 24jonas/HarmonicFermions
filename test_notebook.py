import nbformat
from nbconvert.preprocessors import ExecutePreprocessor

with open('Jupyter/BaxterMeth.ipynb', 'r') as f:
    nb = nbformat.read(f, as_version=4)

# Execute only up to cell 6
nb.cells = nb.cells[:7]

ep = ExecutePreprocessor(timeout=600, kernel_name='python3')
try:
    ep.preprocess(nb, {'metadata': {'path': 'Jupyter/'}})
    print("Execution successful!")
except Exception as e:
    print(f"Error executing notebook: {e}")
