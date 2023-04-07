# qiskit-xyz2pdb

This program is designed primarily to convert XYZ files from the results of [Qiskit's protein folding algorithm](https://qiskit-research.github.io/qiskit-research/protein_folding/protein_folding.html), but may also work with other similar XYZ files for protein coordinates. There are two options for the output: 1) a hetero atom (HETATM) format with connect records (CONECT), and 2) a conventional alpha carbon trace protein pdb format. In both outputs, each set of coordinates corresponds to the alpha carbon of each residue in the sequence. 

## Install

The `qiskit-xyz2pdb` package is available on the Python Package Index and Conda and can be installed with the following command line:

```
# Install the package through PyPI
$ pip install qiskit-xyz2pdb

# Install the package through the Conda package manage
$ conda install qiskit-xyz2pdb
```

## Features

The hetero atom format is most useful for quick visualization, using a program such as NGL viewer. Note: this option assumes that the XYZ coordinates are all consecutively bound forming a continous single chain; it assumes that there are no side chains, and the CONECT records reflect this. How to run with this option:

```
$ qiskit-xyz2pdb --in-xyz inputfile.xyz --out-pdb outputfile.pdb --hetero-atoms
```

The conventional alpha carbon trace format is more versatile and can be used for further research beyond visualization - such as reconstructing the full backbone and sidechains, adding a force field, and then performing molecular dynamics simulations. How to run with this option:

```
$ qiskit-xyz2pdb --in-xyz inputfile.xyz --out-pdb outputfile.pdb --alpha-c-traces
```
