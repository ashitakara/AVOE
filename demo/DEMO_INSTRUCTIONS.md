# Demo Instructions for AVOE

This file provides minimal instructions for running the reviewer demo of AVOE using public PDB-derived example structures.

## 1. Demo contents

### `input/`

The demo input directory contains two ligand-series directories:

- `input/ADN/`
  - `4ckj_cleaned_chainAB.pdb`
  - `4ckj_cleaned_chainAB_G810S.pdb`

- `input/Vandetanib/`
  - `2ivu_cleaned_chainAB.pdb`
  - `2ivu_cleaned_chainAB_G810S.pdb`

These PDB files were prepared from public PDB-derived structures.
Unnecessary components were removed using PyMOL, and chain identifiers were reassigned so that:

- receptor chain = `A`
- ligand chain = `B`

The file `input/ligand_list.txt` lists the two ligand-series directories:

```text
demo/input/ADN/
demo/input/Vandetanib/
```

### `expected_output/`

- `expected_output/vdw.txt`
  - provided van der Waals radii file for the demo
- `expected_output/overlap_volume.csv`
  - reviewer-friendly reference output for overlap-volume results
- `expected_output/summary.txt`
  - short summary of the expected demo result

### `output/`

- `output/`
  - directory for writing demo outputs

## 2. System requirements

Tested environment:

- OS:
  - Ubuntu 18.04 LTS
  - Ubuntu 20.04 LTS
- Python:
  - 3.7.12
  - 3.9.19
- PyMOL:
  - 2.5.0
- Biopython:
  - 1.79, 1.80
- NumPy:
  - 1.21.6, 1.26.4
- pandas:
  - 1.3.4, 1.3.5

No non-standard hardware is required.

## 3. Installation guide

For a fresh Ubuntu environment, it is recommended to install Miniforge and create a dedicated conda environment for the demo.

### 3.1 Install Miniforge

```bash
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh
source ~/.bashrc
```

If shell initialization was skipped during installation, initialize conda manually:

```bash
~/miniforge3/bin/conda init bash
source ~/.bashrc
```

### 3.2 Create and activate the demo environment

```bash
conda create -n avoe-demo python=3.9 -y
conda activate avoe-demo
```

### 3.3 Install required packages

Install PyMOL open-source and the required Python packages into the same environment:

```bash
mamba install -c conda-forge pymol-open-source biopython numpy pandas -y
```

If `mamba` is not available, use:

```bash
conda install -c conda-forge pymol-open-source biopython numpy pandas -y
```

## 4. Demo: generate the vdW radii file using `create_vdw_file.py`

From the repository root, run:

```bash
python create_vdw_file.py \
  demo/input/Vandetanib/2ivu_cleaned_chainAB.pdb \
  demo/output/vdw_generated.txt
```

## 5. Demo: run AVOE using `ligand_list.txt`

From the repository root, run:

```bash
python avoe.py \
  --ligand_file demo/input/ligand_list.txt \
  --vdw_file demo/expected_output/vdw.txt \
  --output_dir demo/output/ \
  --receptor_chain A \
  --ligand_chain B
```

In this demo, `--ligand_file` specifies a text file listing the ligand-series directories to be processed.
Each listed directory contains protein-ligand complex PDB files, and the directory name is treated as the ligand name by AVOE.

## 6. Expected output

Reference files are provided in `demo/expected_output/`.

- `vdw.txt`
  - provided van der Waals radii file used for the demo
- `overlap_volume.tsv`
  - overlap-volume results

The original AVOE program should write its native output files to `demo/output/`.

## 7. Expected runtime

On a standard desktop computer, the demo is expected to complete quickly.
Typical runtime should be on the order of seconds to a few tens of seconds, depending on the environment and PyMOL startup overhead.

## 8. Instructions for use on user data

To run AVOE on your own dataset using the `--ligand_file` option:

1. prepare PDB files in which each file contains both receptor and ligand
2. organize the files into one or more directories as needed
3. prepare a text file listing the target PDB files or directories, one path per line
4. provide a van der Waals radii file
5. run `avoe.py` with the appropriate receptor and ligand chain IDs

## 9. Preparation note for the demo input

The demo input files were prepared from public PDB-derived structures.

For this reviewer demo:
- unnecessary components and artifacts were removed using PyMOL
- chain identifiers were reassigned for a simple demo configuration
- the resulting cleaned structures are provided only as minimal demonstration inputs

## 10. Citation for the demo input

### PDB entry and original structure paper for Vandetanib-series input

`PDB 2IVU: Crystal structure of phosphorylated RET tyrosine kinase domain complexed with the inhibitor ZD6474. DOI: 10.2210/pdb2IVU/pdb`

`Knowles, P.P., Murray-Rust, J., Kjaer, S., Scott, R.P., Hanrahan, S., Santoro, M., Ibanez, C.F., McDonald, N.Q. (2006) Structure and chemical inhibition of the RET tyrosine kinase domain. Journal of Biological Chemistry 281: 33577-33587. PubMed: 16928683. DOI: 10.1074/jbc.M605604200`

### PDB entry and original structure paper for ADN-series input

`PDB 4CKJ: Crystal structure of RET tyrosine kinase domain bound to adenosine. DOI: 10.2210/pdb4CKJ/pdb`

`Plaza-Menacho, I., Barnouin, K., Goodman, K., Martinez-Torres, R.J., Borg, A., Murray-Rust, J., Mouilleron, S., Knowles, P., McDonald, N.Q. (2014) Oncogenic RET kinase domain mutations perturb the autophosphorylation trajectory by enhancing substrate presentation in trans. Molecular Cell 53: 738-751. PubMed: 24560924. DOI: 10.1016/j.molcel.2014.01.015`

## 11. License

AVOE is distributed under the MIT License.
See `LICENSE.md` in the repository root.
