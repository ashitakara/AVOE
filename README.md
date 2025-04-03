
## Associated Publication

The code in this repository is associated with a scientific publication currently in preparation.

If you use this code or find it helpful in your research, please cite the publication once it becomes available. We will update this README with the full citation details upon publication.

Thank you!

# Atomic Volume Overlap Estimator (AVOE)

## 1. Overview

Atomic Volume Overlap Estimator (AVOE) is a tool designed to estimate the overlapped volume between receptor and ligand in protein-ligand complex structures (PDB files) using the Monte Carlo method. It calculates the volume where atoms of both receptor and ligand are within their van der Waals radii of each other. AVOE is valuable in fields like molecular design, interaction analysis, and drug discovery for evaluating the steric complementarity between proteins and ligands.

## 2. Key Features

- **Monte Carlo Estimation of Overlapped Volume:** Statistically evaluates the overlap of van der Waals volumes between receptor and ligand using a large number of randomly generated pseudo-particles.
- **Flexible Calculation Area Definition:** Allows calculation of overlapped volume in areas around the ligand, specified residues, or on a per-residue basis.
- **Van der Waals Radii Retrieval Options:**
    - **Automatic Retrieval via PyMOL:** AVOE can directly utilize PyMOL to retrieve van der Waals radii for atoms in PDB files.
    - **Pre-calculated Radii File:** Users can provide a pre-calculated Van der Waals radii file.  While `create_vdw_file.py` is provided as a utility script to generate such files, **any file adhering to the specified format can be used**, offering significantly faster execution, especially when processing many PDB files.
- **Parallel Processing for Enhanced Speed:** Supports multi-core parallel processing to efficiently analyze numerous PDB files.
- **Detailed Output:** Provides overlapped volume, Monte Carlo simulation statistics, and computation time.
- **Command-Line Interface:** Offers a flexible command-line interface with extensive options for customization.

## 3. Installation

### 3.1 Dependencies

AVOE requires the following Python libraries:

- **Biopython:** For parsing PDB files.
- **PyMOL:** For retrieving van der Waals radii (required if not using a pre-calculated radii file or if you choose dynamic retrieval).
  - Installation instructions for PyMOL can be found on the [PyMOL website](https://pymol.org/).
  - Ensure `pymol` is executable from your command line.
- **NumPy:** For numerical computations.
- **pandas:** For data frame handling and output of results.

### 3.2 Installation Steps

1. **Download AVOE Scripts:**
   Download `avoe.py` and `create_vdw_file.py` and save them to a suitable directory.
2. **Install Dependencies:** Install the required libraries.
3. **(Optional but Recommended) Prepare Van der Waals Radii File:** For improved performance, especially when processing multiple PDB files, it is highly recommended to use a pre-calculated Van der Waals radii file. While `create_vdw_file.py` is provided as a convenient script to generate these files, you can use any method to create a file that conforms to the specified format (see Section 5.2 "Van der Waals Radii File"). Using a pre-calculated file significantly reduces computational cost during the main AVOE execution.

## 4. Usage

### 4.1 Preparation: Creating Van der Waals Radii File (Using `create_vdw_file.py` Utility)

To enhance the efficiency of AVOE, it's recommended to use a pre-calculated Van der Waals radii file.  `create_vdw_file.py` is provided as a utility script to easily generate such a file from a PDB structure.

**Command (using `create_vdw_file.py`):**

```bash
python create_vdw_file.py <input_pdb_file> <output_file> [--chain_id <chain_id>] [--verbose]
```

**Options (for `create_vdw_file.py`):**

- `<input_pdb_file>`: Path to the input PDB file from which to extract van der Waals radii.
- `<output_file>`: Path to the output file for saving van der Waals radii.
- `--chain_id <chain_id>`: Chain ID to extract radii from (optional). If not specified, radii from all chains will be extracted.
- `--verbose`: Enable verbose output (optional).

**Example (using `create_vdw_file.py`):**

```bash
python create_vdw_file.py receptor.pdb vdw_radii.txt --chain_id A
```

This command uses `create_vdw_file.py` to extract Van der Waals radii from `receptor.pdb` (chain A only) and saves them to `vdw_radii.txt`. You can then use `vdw_radii.txt` with the `--vdw_file` option in `avoe.py` for faster calculations.

### 4.2 Running AVOE (`avoe.py`)

Execute AVOE using the `avoe.py` script.

**Command:**

```bash
python avoe.py [options]
```

**Key Options:**

- `--ligand <ligand_pdb_file or directory>`: Path to a single ligand PDB file or a directory containing ligand PDB files. If a directory is specified, all PDB files within it will be processed. The directory name will be used as the ligand name.
- `--ligand_file <ligand_list_file>`: Path to a text file listing ligand PDB files or directories, one path per line. Cannot be used with `--ligand`. If neither `--ligand` nor `--ligand_file` is specified, AVOE will search for PDB files in the current directory and its subdirectories.
- `--vdw_file <van_der_waals_radii_file>`: Path to a pre-calculated Van der Waals radii file. **Using a pre-calculated file, generated either by `create_vdw_file.py` or any other method that produces a file in the correct format, is highly recommended for reducing computational cost and speeding up calculations, especially when processing multiple PDB files.** If provided, AVOE will use this file instead of dynamically retrieving radii (which may involve calling PyMOL).
- `--output_dir <output_directory>`: Directory to save output files. Default is `./output/`.
- `--pseudo_particles <number_of_pseudo_particles>`: Number of pseudo-particles for Monte Carlo simulation. Default is `10000`. Increasing this number improves accuracy but also increases computation time.
- `--repeat <number_of_repeats>`: Number of times to repeat the Monte Carlo calculation. Default is `1`. Repeating multiple times allows for assessing the statistical reliability of the results.
- `--num_cores <number_of_cores>`: Number of CPU cores to use for parallel processing. Default is `1`. Using more cores can significantly speed up the calculation.
- `--receptor_chain <receptor_chain_id>`: Chain ID for the receptor protein. Default is `A`.
- `--ligand_chain <ligand_chain_id>`: Chain ID for the ligand molecule. Default is `B`.
- `--add_lim <bounding_box_expansion_limit>`: Expansion limit for the bounding box around the calculation area in Ångströms. Default is `2.0` Å.
- `--calculation_area <calculation_area>`: Defines the area for overlap volume calculation. Options are:
    - `ligand` (default): Calculate around the ligand vicinity.
    - `residues`: Calculate around the vicinity of residues specified in `--residue_list_file`.
    - `per_residue`: Calculate overlap volume for each residue specified in `--residue_list_file`.
- `--residue_list_file <residue_list_file>`: Path to a text file containing a list of residue numbers (one residue number per line) when `--calculation_area` is set to `residues` or `per_residue`.
- `--calculate_around_mutant_residue`: Detects mutation name from PDB filename (format `[A-Z]\d+[A-Z]`) and calculates overlapped volume specifically around the mutant residue(s). Skips calculation for "WT".
- `--verbose`: Enable verbose output during calculation (optional).

**For a complete list of options, run `python avoe.py -h` to see the help message.**

**Example:**

```bash
python avoe.py --ligand ligand_dir --vdw_file vdw_radii.txt --output_dir output_avoe --pseudo_particles 5000 --repeat 3 --num_cores 4 --verbose
```

This command calculates the overlapped volume for PDB files in `ligand_dir`, utilizing the Van der Waals radii from `vdw_radii.txt`.  This file could have been generated by `create_vdw_file.py` or created manually, as long as it adheres to the correct format. The calculation is performed with 5000 pseudo-particles, 3 repeats, utilizing 4 CPU cores, and saves results to `output_avoe` directory.

## 5. Input File Formats

### 5.1 PDB Files

AVOE accepts standard PDB (Protein Data Bank) format files as input. It reads atom coordinates and element information from ATOM records.

### 5.2 Van der Waals Radii File

The Van der Waals radii file is a plain text file that specifies the Van der Waals radius for each atom type. The file must be formatted such that each line contains an atom name (e.g., C, N, O) and its corresponding Van der Waals radius (in Ångströms), separated by a tab. While `create_vdw_file.py` is provided to help generate files in this format, you can create this file using any method, as long as it adheres to this structure.

**Example (vdw_radii.txt):**

```
C       1.7
N       1.55
O       1.52
S       1.8
H       1.2
...
```


### 5.3 Residue List File

The residue list file is a text file required when `--calculation_area` is set to `residues` or `per_residue`. Each line should contain a residue number for calculation.

**Example (residue_list.txt):**

```
274
310
312
```

### 5.4 Ligand List File

The ligand list file is used with the `--ligand_file` option. It's a text file with each line specifying a path to a ligand PDB file or a directory containing ligand PDB files.

**Example (ligand_list.txt):**

```
ligand_pdbs/ligand1.pdb
ligand_pdbs/ligand2_dir
ligand_pdbs/ligand3.pdb
```

## 6. Output File Formats

### 6.1 Output Files

AVOE outputs results in tab-separated values (TSV) files (`.tsv`).

**Output Files:**

- `<ligand_name>.tsv`: Contains overlapped volume calculation results for each ligand.
- `<repeat_count>_<pseudo_particle_count>.tsv`: Merged results for all ligands.

### 6.2 Output File Content

The TSV files include the following columns (columns may vary based on calculation area settings):

- `ligand`: Ligand name
- `mutation`: Mutation name extracted from the PDB filename. AVOE attempts to extract a mutation name from the input PDB filename using a regular expression to identify patterns like `[InitialAminoAcid][ResidueNumber][MutatedAminoAcid]` (e.g., `H274Y`). If a mutation name is successfully extracted, it is listed here. If no mutation pattern is found in the filename, or if the extracted name is "WT", then "WT" (Wild Type) is recorded.
- `residue`: Residue number (only when calculation_area="per_residue")
- `overlapping_volume`: Calculated overlapped volume (Å³)

## 7. Important Notes

- **PyMOL Installation:** For dynamic retrieval of van der Waals radii, PyMOL must be installed and executable as `pymol` from the command line. If PyMOL is not available, use a pre-calculated radii file with the `--vdw_file` option.
- **PDB File Format:** Ensure input PDB files are in standard format, with ATOM records correctly specifying atom coordinates and element information.
- **Chain ID Specification:** Correctly specify receptor chain ID (`--receptor_chain`) and ligand chain ID (`--ligand_chain`) according to your PDB files.
- **Calculation Area Settings:** Configure `--calculation_area` and related options (`--residue_list_file`) appropriately for your analysis goals.
- **Parameter Tuning:** Adjust `--pseudo_particles` and `--repeat` to balance calculation accuracy and time. Higher pseudo-particle counts improve accuracy but increase computation time.
- **Error Message Handling:** If error messages appear during execution, review them carefully and check input files, option settings, and dependency installations.

## 8. License

This project is licensed under the MIT License - see the LICENSE file for details.

## 9. References

- [PyMOL Documentation](https://pymolwiki.org/index.php/Main_Page)
- [Biopython Documentation](https://biopython.org/)

## 10. Author

- Kosuke Maruyama
- Affiliation: Division of Genome Biology, National Cancer Center Japan/Graduate School of Medicine, The University of Tokyo
- Contact: komaruya@ncc.go.jp

---
