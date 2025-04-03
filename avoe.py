# Copyright (c) 2025 Kosuke Maruyama
# This source code is licensed under the MIT license found in the LICENSE file in the root directory of this source tree.
import sys
from Bio.PDB.PDBParser import PDBParser
import Bio.PDB.Residue
import Bio.PDB.Atom
import Bio.PDB.Selection
import pymol
import itertools
import random
import time
import os
import numpy as np
import gc
import pandas as pd
import glob
import re
from multiprocessing import Pool
import argparse
from pathlib import Path # Import Pathlib

# Define global variable for Van der Waals radii dictionary
vdw_radius_dict_global = {}

# Function to validate rectangular box dimensions
def validate_rectangular_box(vectors, tolerance=1e-6):
    """
    Validate that the vectors defining a rectangular box are orthogonal and non-zero length.

    Args:
        vectors (list of numpy.ndarray): List of three vectors representing the box axes.
        tolerance (float): Tolerance for orthogonality check (default: 1e-6).

    Returns:
        None: Exits the program if validation fails.
    """
    # Validate orthogonality and non-zero length of box vectors
    vector_combinations = itertools.combinations(vectors, 2)
    for v_c in vector_combinations:
        if abs(np.dot(v_c[0], v_c[1])) > tolerance:
            print(f"Error: Vectors are not orthogonal. Dot product: {np.dot(v_c[0], v_c[1])}")
            sys.exit(1)
    for v in vectors:
        if np.linalg.norm(v) == 0:
            print("Error: Zero length vector detected.")
            sys.exit(1)


# Function to calculate overlapped volume using Monte Carlo method
def calculate_overlapped_area(pdb, mutant, add_lim, pseudo_particles, receptor_chain_id, ligand_chain_id,
                              calculation_area, residue_list_file,
                              calculate_around_mutant_residue, residue,
                              verbose):
    """
    Calculate the overlapped volume between receptor and ligand in a PDB structure using Monte Carlo method,
    allowing for calculation around ligand, specified residues, or per-residue calculation.

    Args:
        pdb (str): Path to the PDB file.
        mutant (str): Mutation name extracted from PDB file name (e.g., "WT", "H274Y").
        add_lim (float): Additional expansion limit for the bounding box around the calculation area (Angstroms).
        pseudo_particles (int): Number of pseudo particles (Monte Carlo points) to use.
        receptor_chain_id (str): Chain ID for the receptor protein (default: "A").
        ligand_chain_id (str): Chain ID for the ligand molecule (default: "B").
        calculation_area (str): Area to calculate overlap.
            - "ligand": Calculate overlap around the ligand.
            - "residues": Calculate overlap around specified residues listed in residue_list_file.
            - "per_residue": Calculate overlap for each residue listed in residue_list_file.
            Defaults to "ligand".
        residue_list_file (str, optional): Path to file with residue list for --calculation_area="residues" or "per_residue". Defaults to None.
        calculate_around_mutant_residue (bool): Calculate overlapped volume specifically around mutant residue(s) if mutant is not "WT". Defaults to False.
        residue (int, optional): Residue number for per-residue calculation when calculation_area="per_residue". Defaults to None.
        verbose (bool): Enable verbose output (default: False).

    Returns:
        tuple: Overlapped volume and elapsed time for the calculation.
    """
    start = time.time()

    pdb_file = pdb
    name = os.path.splitext(os.path.basename(pdb_file))[0]
    if verbose:
        print("open pdb: ", name)

    pdb_parser = PDBParser()
    try:
        with open(pdb_file, "r") as handle:
            struc = pdb_parser.get_structure(name, handle)
    except FileNotFoundError:
        print(f"Error: PDB file not found: {pdb_file}")
        sys.exit(1)
    except Exception as e:
        print(f"Error: Failed to parse PDB file {pdb_file}: {e}")
        sys.exit(1)

    # Check if chain IDs exist in the PDB structure
    if receptor_chain_id not in struc[0]:
        print(f"Error: Receptor chain ID '{receptor_chain_id}' not found in PDB file {pdb_file}.")
        sys.exit(1)
    if ligand_chain_id not in struc[0]:
        print(f"Error: Ligand chain ID '{ligand_chain_id}' not found in PDB file {pdb_file}.")
        sys.exit(1)

    if calculation_area == "per_residue":
        all_recep_atoms_pdb = Bio.PDB.Selection.unfold_entities(struc[0][receptor_chain_id][residue], "A")
    else:
        all_recep_atoms_pdb = Bio.PDB.Selection.unfold_entities(struc[0][receptor_chain_id], "A")
    if verbose:
        print(f"all_recep_atoms (chain {receptor_chain_id}): {len(all_recep_atoms_pdb)}")

    all_ligand_atoms_pdb = Bio.PDB.Selection.unfold_entities(struc[0][ligand_chain_id], "A")
    if verbose:
        print(f"ligand_atoms (chain {ligand_chain_id}): {len(all_ligand_atoms_pdb)}")

    if calculate_around_mutant_residue and mutant == "WT": # Calculate around mutant residue
        if verbose:
            print("Skipping volume calculation for WT mutant.")
        elapsed_time = time.time() - start
        return np.nan, elapsed_time

    elif calculate_around_mutant_residue and mutant != "WT": # Calculate around mutant residue
        try:
            mutant_residue_num = int(mutant[1:-1]) # Extract residue number from mutant name
            mutant_residue = struc[0][receptor_chain_id][mutant_residue_num] # Get mutant residue object
            area_atoms = Bio.PDB.Selection.unfold_entities(mutant_residue, "A") # Get atoms of mutant residue
        except KeyError:
            print(f"Error: Mutant residue {mutant[1:-1]} not found in chain {receptor_chain_id}. Please check PDB file and mutant naming.")
            sys.exit(1)
        except ValueError:
            print(f"Error: Could not parse residue number from mutant name '{mutant}'. Please check PDB filename format.")
            sys.exit(1)

        if not area_atoms: # Error handling if no atoms found for mutant residue
            print(f"Error: No atoms found for mutant residue {mutant[1:-1]} in chain {receptor_chain_id}.")
            sys.exit(1)

        x = [atom.get_coord()[0] for atom in area_atoms]
        y = [atom.get_coord()[1] for atom in area_atoms]
        z = [atom.get_coord()[2] for atom in area_atoms]

    elif calculation_area == "ligand":
        x = [atom.get_coord()[0] for atom in all_ligand_atoms_pdb]
        y = [atom.get_coord()[1] for atom in all_ligand_atoms_pdb]
        z = [atom.get_coord()[2] for atom in all_ligand_atoms_pdb]

    elif calculation_area == "residues":
        # Read residue number from a residue file
        if residue_list_file is None:
            print("Error: --residue_list_file must be specified when --calculation_area='residues'.")
            sys.exit(1)
        try:
            with open(residue_list_file, 'r') as f:
                residue_numbers = [int(line.strip()) for line in f]
        except FileNotFoundError:
            print(f"Error: Residue list file not found: {residue_list_file}")
            sys.exit(1)

        area_atoms = []
        for residue_num in residue_numbers:
            try:
                area_atoms += Bio.PDB.Selection.unfold_entities(struc[0][receptor_chain_id][residue_num], "A")
            except KeyError:
                print(f"Warning: Residue {residue_num} not found in chain {receptor_chain_id}. Skipping.") # Residue not found warning
                continue # Skip to the next residue
        if not area_atoms: # Error handling if no atoms found for specified residues
            print(f"Error: No atoms found for specified residues in chain {receptor_chain_id}.")
            sys.exit(1)

        x = [atom.get_coord()[0] for atom in area_atoms]
        y = [atom.get_coord()[1] for atom in area_atoms]
        z = [atom.get_coord()[2] for atom in area_atoms]

    elif calculation_area == "per_residue":
        if residue is None:
            print("Error: Residue number must be specified when --calculation_area='per_residue'.")
            sys.exit(1)

        area_atoms = all_recep_atoms_pdb

        x = [atom.get_coord()[0] for atom in area_atoms]
        y = [atom.get_coord()[1] for atom in area_atoms]
        z = [atom.get_coord()[2] for atom in area_atoms]

    else:
        print(f"Error: Invalid calculation area: {calculation_area}. Choose 'ligand' or 'residues'.")
        sys.exit(1)

    x_limits = np.array([min(x), max(x)])
    y_limits = np.array([min(y), max(y)])
    z_limits = np.array([min(z), max(z)])
    expansion_limit = np.array([-add_lim, add_lim])
    x_limits += expansion_limit
    y_limits += expansion_limit
    z_limits += expansion_limit

    if calculation_area == "ligand":
        recep_atom_pdb = [
            atom
            for atom in all_recep_atoms_pdb
            if (x_limits[0] -add_lim < atom.get_coord()[0] < x_limits[1] + add_lim) and
               (y_limits[0] -add_lim < atom.get_coord()[1] < y_limits[1] + add_lim) and
               (z_limits[0] -add_lim < atom.get_coord()[2] < z_limits[1] + add_lim)]
    else:
        recep_atom_pdb = [atom for atom in all_recep_atoms_pdb]
    receptor_atoms = np.array([[atom.get_coord()[0], atom.get_coord()[1], atom.get_coord()[2]] for atom in recep_atom_pdb])
    receptor_radii = []

    if not recep_atom_pdb:
        elapsed_time = time.time() - start
        return 0, elapsed_time

    pymol_pdb_loaded = False # PyMOL PDB load flag

    for atom in recep_atom_pdb:
        a_name = atom.get_name()
        a_element = atom.element
        if a_element == None:
            print(f"Error: Failed to get the element of {a_name}")
            sys.exit(1)
        if a_element in vdw_radius_dict_global:
            radius = vdw_radius_dict_global[a_element]
        else:
            if not pymol_pdb_loaded: # Load PDB to PyMOL only when needed
                if verbose: # Verbose output before PyMOL load
                    print(f"Loading PDB '{pdb_file}' to PyMOL for VDW radius retrieval...")
                pymol.pymol_argv = ["pymol", "-c"]
                try:
                    pymol.cmd.load(pdb_file)
                    pymol_pdb_loaded = True # Set flag after successful load
                except Exception as e:
                    print(f"Error: PyMOL failed to load PDB file {pdb_file}: {e}")
                    print("Make sure PyMOL is installed and PDB file is valid.")
                    sys.exit(1)
            resi = atom.get_parent().get_full_id()[3][1]
            myspace = {"van_der_w": []}
            try:
                pymol.cmd.iterate(f"/{name}//{receptor_chain_id}/{resi}/{a_name}", "van_der_w.append(vdw)", space=myspace)
                if verbose: # Verbose output after PyMOL iterate
                    print(f"Retrieved VDW radius for atom '{a_name}' from PyMOL.")
            except Exception as e:
                print(f"Error: PyMOL iterate command failed for atom {a_name}: {e}")
                print("Please check PyMOL installation and PDB file format.")
                sys.exit(1)
            if not myspace["van_der_w"]: # Check if vdw_radius is retrieved from PyMOL
                print(f"Error: Could not retrieve Van der Waals radius for atom {a_name} using PyMOL.")
                print("Please check atom name and PDB file format.")
                sys.exit(1)
            radius = myspace["van_der_w"][0]
            vdw_radius_dict_global[a_element] = radius
        receptor_radii.append(radius)

    if calculation_area == "ligand":
        ligand_atoms_pdb = [atom for atom in all_ligand_atoms_pdb]
        ligand_atoms = np.array([[atom.get_coord()[0], atom.get_coord()[1], atom.get_coord()[2]] for atom in all_ligand_atoms_pdb])
    else:
        ligand_atoms_pdb = [
            atom
            for atom in all_ligand_atoms_pdb
            if (x_limits[0] -add_lim < atom.get_coord()[0] < x_limits[1] + add_lim) and
               (y_limits[0] -add_lim < atom.get_coord()[1] < y_limits[1] + add_lim) and
               (z_limits[0] -add_lim < atom.get_coord()[2] < z_limits[1] + add_lim)]
        ligand_atoms = np.array([[atom.get_coord()[0], atom.get_coord()[1], atom.get_coord()[2]] for atom in ligand_atoms_pdb])
    if verbose:
        print(f"Type of ligand_atoms_pdb: {type(ligand_atoms_pdb)}") # debug
    if not ligand_atoms_pdb:
        elapsed_time = time.time() - start
        return 0, elapsed_time

    ligand_radii = []

    for atom in ligand_atoms_pdb:
        a_name = atom.get_name()
        a_element = atom.element
        if a_element == None:
            print(f"Error: Failed to get the element of {a_name}")
            sys.exit(1)
        if a_element in vdw_radius_dict_global:
            radius = vdw_radius_dict_global[a_element]
        else:
            if not pymol_pdb_loaded: # Load PDB to PyMOL only when needed
                if verbose: # Verbose output before PyMOL load
                    print(f"Loading PDB '{pdb_file}' to PyMOL for VDW radius retrieval...")
                pymol.pymol_argv = ["pymol", "-c"]
                try:
                    pymol.cmd.load(pdb_file)
                    pymol_pdb_loaded = True # Set flag after successful load
                except Exception as e:
                    print(f"Error: PyMOL failed to load PDB file {pdb_file}: {e}")
                    print("Make sure PyMOL is installed and PDB file is valid.")
                    sys.exit(1)
            resi = atom.get_parent().get_full_id()[3][1]
            myspace = {"van_der_w": []}
            try:
                pymol.cmd.iterate(f"/{name}//{ligand_chain_id}/{resi}/{a_name}", "van_der_w.append(vdw)", space=myspace)
                if verbose: # Verbose output after PyMOL iterate
                    print(f"Retrieved VDW radius for atom '{a_name}' from PyMOL.")
            except Exception as e:
                print(f"Error: PyMOL iterate command failed for atom {a_name}: {e}")
                print("Please check PyMOL installation and PDB file format.")
                sys.exit(1)
            if not myspace["van_der_w"]: # Check if vdw_radius is retrieved from PyMOL
                print(f"Error: Could not retrieve Van der Waals radius for atom {a_name} using PyMOL.")
                print("Please check atom name and PDB file format.")
                sys.exit(1)
            radius = myspace["van_der_w"][0]
            vdw_radius_dict_global[a_element] = radius
        ligand_radii.append(radius)

    if pymol_pdb_loaded:
        pymol.cmd.do("delete all") # Conditionally delete all objects in PyMOL

    a = np.array([x_limits[0], y_limits[0], z_limits[0]])
    b = np.array([x_limits[1], y_limits[0], z_limits[0]])
    c = np.array([x_limits[0], y_limits[1], z_limits[0]])
    d = np.array([x_limits[1], y_limits[1], z_limits[0]])
    e = np.array([x_limits[0], y_limits[0], z_limits[1]])
    f = np.array([x_limits[1], y_limits[0], z_limits[1]])
    g = np.array([x_limits[0], y_limits[1], z_limits[1]])
    h = np.array([x_limits[1], y_limits[1], z_limits[1]])

    vectors = [b-a, c-a, e-a]

    validate_rectangular_box(vectors)

    x_length = np.linalg.norm(vectors[0])
    y_length = np.linalg.norm(vectors[1])
    z_length = np.linalg.norm(vectors[2])

    box_volume = x_length * y_length * z_length

    mc_points = pseudo_particles
    p = [[[random.uniform(x_limits[0], x_limits[1]), random.uniform(y_limits[0], y_limits[1]), random.uniform(z_limits[0], z_limits[1])]] for i in range(0, mc_points, 1)]

    # Calculate distances between MC points and receptor/ligand atoms, and check for overlaps
    expanded_receptor_atoms = np.broadcast_to(receptor_atoms, (len(p), len(receptor_atoms), 3))
    receptor_atom_vectors_to_mc_points = expanded_receptor_atoms - p
    receptor_atom_distances = [[np.linalg.norm(j) for j in i] for i in receptor_atom_vectors_to_mc_points]
    receptor_overlap_grid = np.where(np.array(receptor_atom_distances) < receptor_radii, 1, 0)

    expanded_ligand_atoms = np.broadcast_to(ligand_atoms, (len(p), len(ligand_atoms), 3))
    ligand_atom_vectors_to_mc_points = expanded_ligand_atoms - p
    ligand_atom_distances = [[np.linalg.norm(j) for j in i] for i in ligand_atom_vectors_to_mc_points]
    ligand_overlap_grid = np.where(np.array(ligand_atom_distances) < ligand_radii, 1, 0)

    # Pad result arrays to ensure same size for consistent processing
    pad_size = max(len(receptor_overlap_grid[0]), len(ligand_overlap_grid[0])) - min(len(receptor_overlap_grid[0]), len(ligand_overlap_grid[0]))
    if len(receptor_overlap_grid[0]) > len(ligand_overlap_grid[0]):
        ligand_overlap_grid = np.pad(ligand_overlap_grid, [(0, 0), (0, pad_size)], "constant")
    elif len(receptor_overlap_grid[0]) < len(ligand_overlap_grid[0]):
        receptor_overlap_grid = np.pad(receptor_overlap_grid, [(0, 0), (0, pad_size)], "constant")

    result = np.array([receptor_overlap_grid, ligand_overlap_grid])
    mc_point_overlap_results = result.transpose(1, 0, 2)
    total_pseudo_particles = len(mc_point_overlap_results)
    overlapped_pseudo_particles = np.sum([detect_overlap(x) for x in mc_point_overlap_results])

    # Delete variables
    del receptor_overlap_grid, ligand_overlap_grid, result, mc_point_overlap_results, receptor_atoms, ligand_atoms, all_ligand_atoms_pdb, ligand_atoms_pdb, all_recep_atoms_pdb

    if verbose:
        print("ratio: " + str(overlapped_pseudo_particles/total_pseudo_particles))
    overlapped_volume = box_volume * (overlapped_pseudo_particles/total_pseudo_particles)
    if verbose:
        print("Volume: {}.".format(str(overlapped_volume)))
    elapsed_time = time.time() - start
    if verbose:
        print("elapsed_time: {} sec.".format(str(elapsed_time)))

    return overlapped_volume, elapsed_time

# Function to detect overlap for a single pseudo particle
def detect_overlap(two_d_array):
    """
    Detect if a pseudo particle overlaps with both receptor and ligand.

    Args:
        two_d_array (numpy.ndarray): 2D array where two_d_array[0] is overlap grid for receptor and two_d_array[1] for ligand.

    Returns:
        int: 1 if overlap is detected, 0 otherwise.
    """
    if np.sum(two_d_array[0]) > 0 and np.sum(two_d_array[1]) > 0:
        return 1
    else:
        return 0

# Function to repeat volume calculation and aggregate results
def repeat_calculate(pdb, add_lim, pseudo_particles, repeat, ligand, receptor_chain_id, ligand_chain_id,
                      calculation_area, residue_list_file,
                      calculate_around_mutant_residue, verbose):
    """
    Repeat the overlapped volume calculation multiple times and aggregate the results.

    Args:
        pdb (str): Path to the PDB file.
        add_lim (float): Additional expansion limit for the bounding box.
        pseudo_particles (int): Number of pseudo particles for Monte Carlo.
        repeat (int): Number of times to repeat the calculation.
        ligand (str): Ligand name or identifier.
        receptor_chain_id (str): Receptor chain ID.
        ligand_chain_id (str): Ligand chain ID.
        calculation_area (str): Calculation area ("ligand", "residues", or "per_residue").
        residue_list_file (str, optional): Path to residue list file (required for "residues" and "per_residue"). Defaults to None.
        calculate_around_mutant_residue (bool): Calculate around mutant residue.
        verbose (bool): Enable verbose output.

    Returns:
        pandas.DataFrame: DataFrame containing results for each repeat.
    """
    v_list = []
    time_list = []
    if calculation_area=="per_residue":
        results_dict = {"ligand": [], "mutation": [], "residue": [], "overlapping_volume": []}
    else:
        results_dict = {"ligand": [], "mutation": [], "overlapping_volume": []}
    mutant_match = re.findall(r"^.+([A-Z]\d+[A-Z]).*.pdb$", pdb)
    mutant = "WT" # Default mutant name
    if mutant_match and len(mutant_match) > 0:
        mutant = mutant_match[0]
    else:
        if verbose:
            print(f"Warning: Could not extract mutation name from PDB filename '{pdb}'. Using default 'WT'.") # Warning message if regex fails

    for i in range(repeat):
        if calculation_area == "per_residue":
            if residue_list_file is None:
                print("Error: --residue_list_file must be specified when --calculation_area='per_residue'.")
                sys.exit(1)
            try:
                with open(residue_list_file, 'r') as f:
                    residue_numbers = [int(line.strip()) for line in f]
            except FileNotFoundError:
                print(f"Error: Residue list file not found: {residue_list_file}")
                sys.exit(1)
            for residue in residue_numbers:
                v, t = calculate_overlapped_area(
                    pdb, mutant, add_lim, pseudo_particles, receptor_chain_id, ligand_chain_id,
                    calculation_area, residue_list_file,
                    calculate_around_mutant_residue, residue, verbose
                )
                v_list.append(v)
                time_list.append(t)
                gc.collect()
                results_dict["ligand"].append(ligand)
                results_dict["mutation"].append(mutant)
                results_dict["residue"].append(residue)
                results_dict["overlapping_volume"].append(v)
                if verbose and calculation_area != "per_residue":
                    print("Residue: {}".format(str(residue)), "Volume_mean: {}".format(np.mean(np.array(v_list))), "Var: {}".format(np.var(np.array(v_list))))
                del v, t

        else:
            residue = None
            v, t = calculate_overlapped_area(
                pdb, mutant, add_lim, pseudo_particles, receptor_chain_id, ligand_chain_id,
                calculation_area, residue_list_file,
                calculate_around_mutant_residue, residue, verbose
            )
            v_list.append(v)
            time_list.append(t)
            gc.collect()
            results_dict["ligand"].append(ligand)
            results_dict["mutation"].append(mutant)
            results_dict["overlapping_volume"].append(v)
            del v, t

    print(v_list)
    print("Volume_mean: {}".format(np.mean(np.array(v_list))), "Var: {}".format(np.var(np.array(v_list))))
    print("Time_mean: {} sec.".format(np.mean(np.array(time_list))))
    if verbose:
        print(results_dict)
    return pd.DataFrame(results_dict)

# Wrapper function for parallel processing
def calc_wrapper(args):
    """
    Wrapper function to execute repeat_calculate with arguments for parallel processing.

    Args:
        args (tuple): Tuple of arguments to be passed to repeat_calculate.

    Returns:
        pandas.DataFrame: DataFrame returned by repeat_calculate.
    """
    (pdb, add_lim, pseudo_particles, repeat, ligand, receptor_chain_id, ligand_chain_id,
     calculation_area, residue_list_file,
     calculate_around_mutant_residue, verbose) = args
    return repeat_calculate(
        pdb, add_lim, pseudo_particles, repeat, ligand, receptor_chain_id, ligand_chain_id,
        calculation_area, residue_list_file,
        calculate_around_mutant_residue, verbose
    )

def get_ligand_pdb_path_and_name(ligand_input, verbose):
    """
    Get ligand PDB file paths and ligand name based on input.

    Args:
        ligand_input (str): Path to ligand directory or ligand name, or PDB file.
        verbose (bool): Enable verbose output.

    Returns:
        tuple: Ligand name and list of ligand PDB file paths.
    """
    ligand_pdb_path_list = []
    ligand_path = Path(ligand_input)

    if ligand_path.is_dir(): # Ligand input is directory
        ligand_name = ligand_path.name
        ligand_pdb_path_list = glob.glob(str(ligand_path / "**/*.pdb"), recursive=True)
        if verbose:
            print(f"Processing ligands from directory: {ligand_input}")
    elif ligand_path.is_file(): # Ligand input is a file
        if ligand_path.suffix == ".pdb":
            ligand_name = ligand_path.stem
            ligand_pdb_path_list = [ligand_input]
            if verbose:
                print(f"Processing ligand PDB file: {ligand_input}")
        else: # Input file is not a PDB file
            print(f"Error: Input file '{ligand_input}' is not a PDB file. Please provide a PDB file or a directory.")
            sys.exit(1)
    else: # Ligand input is ligand name (assuming PDBs in current directory)
        ligand_name = ligand_input
        ligand_pdb_path_list = glob.glob(f"./{ligand_input}/**/*.pdb")
        if verbose:
            print(f"Processing ligands with name: {ligand_input} from current directory structure")

    return ligand_name, ligand_pdb_path_list

# Main execution block
if __name__ == "__main__":
    t_start = time.time()

    parser = argparse.ArgumentParser(description='Calculate overlapped volume using Monte Carlo method.')
    parser.add_argument('--verbose', action='store_true', help='Enable verbose output during calculation.')
    parser.add_argument('--receptor_chain', type=str, default="A", help='Chain ID for receptor protein. Defaults to "A".')
    parser.add_argument('--ligand_chain', type=str, default="B", help='Chain ID for ligand molecule. Defaults to "B".')
    parser.add_argument('--add_lim', type=float, default=2.0, help='Additional limit for bounding box expansion (Angstroms). Defaults to 2.0.')
    parser.add_argument('--vdw_file', type=str, default=None, help='Path to file containing Van der Waals radii.')
    parser.add_argument('--output_dir', type=str, default="./output/", help='Output directory. Defaults to "./output/".')
    parser.add_argument('--pseudo_particles', type=int, default=10000, help='Number of pseudo particles for Monte Carlo. Defaults to 10000.')
    parser.add_argument('--num_cores', type=int, default=1, help='Number of cores for parallel processing. Defaults to 1.')
    parser.add_argument('--ligand', type=str, default=None, help='Path to a single ligand PDB file or ligand directory/name.')
    parser.add_argument('--ligand_file', type=str, default=None, help='Path to a text file listing ligand PDB files or directories/names.')
    parser.add_argument('--repeat', type=int, default=1, help='Number of times to repeat Monte Carlo calculation. Defaults to 1.')
    parser.add_argument('--calculation_area', type=str, default="ligand", choices=['ligand', 'residues', 'per_residue'],
                        help='Area for overlap calculation: "ligand" (ligand vicinity), "residues" (vicinity of specified residues), or "per_residue" (per residue overlap). Defaults to "ligand".')
    parser.add_argument('--residue_list_file', type=str, default=None,
                        help='Path to a text file containing a list of residue numbers (one residue per line) for calculation area. Required if --calculation_area="residues".')
    parser.add_argument('--calculate_around_mutant_residue', action='store_true', default=False,
                        help='Calculate overlapped volume specifically around mutant residue(s).')


    args = parser.parse_args()
    if args.ligand and args.ligand_file: # Check for mutual exclusion of --ligand and --ligand_file
        parser.error("--ligand and --ligand_file cannot be specified together. Choose one or none.")
    if args.calculation_area == "residues" and args.residue_list_file is None:
        parser.error("--residue_list_file must be specified when --calculation_area='residues'.")
    if args.calculation_area == "per_residue" and args.residue_list_file is None: # Add error check for per_residue
        parser.error("--residue_list_file must be specified when --calculation_area='per_residue'.")


    verbose = args.verbose
    receptor_chain = args.receptor_chain
    ligand_chain = args.ligand_chain
    add_lim = args.add_lim
    vdw_file = args.vdw_file
    output_dir = args.output_dir
    pseudo_particles = args.pseudo_particles
    num_cores = args.num_cores
    ligand_pdb_path = args.ligand
    ligand_file_path = args.ligand_file
    repeat = args.repeat

    calculation_area = args.calculation_area
    residue_list_file = args.residue_list_file
    calculate_around_mutant_residue = args.calculate_around_mutant_residue


    # Command line argument validation
    if add_lim < 0:
        print("Error: --add_lim must be a non-negative value.")
        sys.exit(1)
    if pseudo_particles <= 0:
        print("Error: --pseudo_particles must be a positive integer.")
        sys.exit(1)
    if num_cores <= 0:
        print("Error: --num_cores must be a positive integer.")
        sys.exit(1)
    if repeat <= 0:
        print("Error: --repeat must be a positive integer.")
        sys.exit(1)

    os.makedirs(output_dir, exist_ok=True)

    # Load Van der Waals radii to global dictionary if vdw_file is specified
    if vdw_file:
        if verbose:
            print(f"Loading Van der Waals radii from: {vdw_file} to global dictionary")
        try:
            with open(vdw_file, 'r') as f:
                for line in f:
                    atom_name, radius = line.strip().split()
                    vdw_radius_dict_global[atom_name] = float(radius)
            if verbose:
                print("Van der Waals radii loaded successfully to global dictionary from file.")
        except FileNotFoundError:
            print(f"Error: Van der Waals radii file not found: {vdw_file}")
            sys.exit(1)
        except ValueError:
            print(f"Error: Invalid format in Van der Waals radii file {vdw_file}. Ensure each line contains 'atom_name radius'.")
            sys.exit(1)
        except Exception as e:
            print(f"Error: Failed to load Van der Waals radii from {vdw_file}: {e}")
            sys.exit(1)


    df_list = []
    ligand_list = []
    ligand_name = "ligand" # Default ligand name for default discovery case

    if ligand_pdb_path: # Ligand input from --ligand
        ligand_list = [ligand_pdb_path]
    elif ligand_file_path: # Ligand input from --ligand_file
        try:
            with open(ligand_file_path, 'r') as f:
                ligand_list = [line.strip() for line in f]
        except FileNotFoundError:
            print(f"Error: Ligand list file not found: {ligand_file_path}")
            sys.exit(1)
    else: # Default ligand discovery (no --ligand or --ligand_file provided)
        ligand_list = ["."]


    for ligand in ligand_list:
        ligand_pdb_path_list = []
        if not ligand_pdb_path and not ligand_file_path: # Default ligand discovery
            ligand_pdb_path_list = glob.glob(ligand + "/**/*.pdb")
        else: # Ligand input from --ligand or --ligand_file
            ligand_name, ligand_pdb_path_list = get_ligand_pdb_path_and_name(ligand, verbose) # Process ligand input to get ligand name and PDB paths
        if not ligand_pdb_path_list: # Check if ligand_pdb_path_list is empty
            print(f"Error: No PDB files found for ligand '{ligand}'. Please check ligand input and directory structure.")
            sys.exit(1) # Exit if no PDB files are found

        values = [(ligand_pdb, add_lim, pseudo_particles, repeat, ligand_name, receptor_chain, ligand_chain,
                   calculation_area, residue_list_file,
                   calculate_around_mutant_residue, verbose)
                  for ligand_pdb in ligand_pdb_path_list] # Prepare arguments for parallel processing

        if verbose:
            print("Values for calc_wrapper:", values) # Print values for debugging


        with Pool(num_cores) as p:
            df = pd.concat(p.map(calc_wrapper, values))
        df_list.append(df)
        df.to_csv(os.path.join(output_dir, f"{ligand_name}.tsv"), sep="\t", index=False) # Save result for each ligand

    df_merge = pd.concat(df_list)

    df_merge.to_csv(os.path.join(output_dir, f"{repeat}_{pseudo_particles}.tsv"), sep="\t", index=False) # Save merged results

    t_end = time.time()
    if verbose:
        print("elapsed_time: " + str(t_end - t_start))
