# Demo Instructions for AVOE

This file provides minimal instructions for running the reviewer demo of AVOE using a public PDB-derived example.

## 1. Demo contents

- `input/2ivu_cleaned_chainAB.pdb`
  - demo input derived from public PDB entry 2IVU
  - unnecessary components were removed using PyMOL
  - chain identifiers were reassigned so that:
    - receptor chain = `A`
    - ligand chain = `B`

- `expected_output/vdw.txt`
  - expected van der Waals radii file generated from the cleaned demo PDB using `create_vdw_file.py`

- `expected_output/overlap_volume.csv`
  - reviewer-friendly reference output for overlap-volume results

- `expected_output/summary.txt`
  - short summary of the expected demo result

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

Clone or download this repository and prepare a Python environment with the required dependencies.

Example:

```bash
pip install biopython numpy pandas
```

PyMOL must also be available in the Python environment, because the original AVOE workflow uses PyMOL-based processing.

Typical installation time on a standard desktop computer:
- if Python and PyMOL are already available: very short (typically a few minutes or less)

## 4. Demo: generate the vdW radii file

From the repository root, run:

```bash
python create_vdw_file.py \
  demo/2ivu_cleaned_chainAB.pdb \
  demo/output/vdw_generated.txt
```


## 5. Demo: run AVOE

From the repository root, run:

```bash
python avoe.py \
  --ligand demo/input/2ivu_cleaned_chainAB.pdb \
  --vdw_file demo/vdw_generated.txt \
  --output_dir demo/output/ \
  --receptor_chain A \
  --ligand_chain B
```

## 6. Expected output

Reference files are provided in `demo/expected_output/`.

- `vdw.txt`
  - expected output from `create_vdw_file.py`
- `overlap_volume.csv`
  - reviewer-friendly reference table for overlap-volume results
- `summary.txt`
  - short summary of the expected result

The original AVOE program should write its native output files to `demo/output/`.

## 7. Expected runtime

On a standard desktop computer, the demo is expected to complete quickly.
Typical runtime should be on the order of seconds, depending on the environment and PyMOL startup overhead.

## 8. Instructions for use on user data

To run AVOE on another PDB structure:

1. prepare a PDB file in a format suitable for the AVOE workflow
2. generate a van der Waals radii file using `create_vdw_file.py`
3. run `avoe.py` with the appropriate receptor and ligand chain IDs

Example:

```bash
python create_vdw_file.py your_input.pdb your_vdw.txt

python avoe.py \
  --ligand your_input.pdb \
  --vdw_file your_vdw.txt \
  --output_dir your_output_dir/ \
  --receptor_chain A \
  --ligand_chain B
```

## 9. Preparation note for the demo input

The file `2ivu_cleaned_chainAB.pdb` was prepared from public PDB entry 2IVU.

For this reviewer demo:
- unnecessary components and artifacts were removed using PyMOL
- chain identifiers were reassigned for a simple demo configuration
- the resulting cleaned structure is provided only as a minimal demonstration input

## 10. Citation for the demo input

### PDB entry

`PDB 2IVU: Crystal structure of phosphorylated RET tyrosine kinase domain complexed with the inhibitor ZD6474. DOI: 10.2210/pdb2IVU/pdb`

### Original structure paper

`Knowles, P.P., Murray-Rust, J., Kjaer, S., Scott, R.P., Hanrahan, S., Santoro, M., Ibanez, C.F., McDonald, N.Q. (2006) Structure and chemical inhibition of the RET tyrosine kinase domain. Journal of Biological Chemistry 281: 33577-33587. PubMed: 16928683. DOI: 10.1074/jbc.M605604200`
