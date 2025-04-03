# Copyright (c) 2025 Kosuke Maruyama
# This source code is licensed under the MIT license found in the LICENSE file in the root directory of this source tree.
import sys
import os
import pymol
import argparse
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Selection import unfold_entities
import time # Import time module

# Function to creates a Van der Waals radii file from a PDB file
def create_vdw_file(pdb_file, output_file, chain_id=None, verbose=False):
    """
    Create a Van der Waals radii file from a PDB file using PyMOL.

    Args:
        pdb_file (str): Path to the input PDB file.
        output_file (str): Path to the output Van der Waals radii file.
        chain_id (str, optional): Chain ID to extract atoms from. Defaults to None (all chains).
        verbose (bool, optional): Enable verbose output. Defaults to False.

    Returns:
        None: Writes Van der Waals radii to the output file.
    """

    start_time = time.time() # Record start time
    vdw_radii = {}

    if verbose:
        print(f"Loading PDB file: {pdb_file}")
    name = os.path.splitext(os.path.basename(pdb_file))[0]

    pdb_parser = PDBParser()
    try:
        with open(pdb_file, "r") as handle:
            structure = pdb_parser.get_structure(name, handle)
    except FileNotFoundError:
        print(f"Error: PDB file not found: {pdb_file}")
        sys.exit(1)
    except Exception as e:
        print(f"Error: Failed to parse PDB file {pdb_file}: {e}")
        sys.exit(1)

    if chain_id is None:
        if verbose:
            print("Extracting atoms from all chains...")
        chains = [chain.id for chain in structure[0].get_chains()]
    else:
        if verbose:
            print(f"Extracting atoms from chain: {chain_id}")
        chains = [chain_id]

    pymol_pdb_loaded = False # PyMOL PDB load flag

    for chain in chains:
        atom_list = unfold_entities(structure[0][chain], 'A')
        for atom in atom_list:
            a_name = atom.get_name()
            a_element = atom.element
            if a_element == None:
                print(f"Error: Failed to get the element of {a_name}")
                sys.exit(1)

            if verbose:
                print(f"Atom name: {a_name}")
            if a_element in vdw_radii.keys():
                pass
            else:
                if not pymol_pdb_loaded: # Load PDB to PyMOL only when needed
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
                    if verbose:
                        print(f"PyMol iterate: /{name}//{chain}/{resi}/{a_name}")
                    pymol.cmd.iterate(f"/{name}//{chain}/{resi}/{a_name}", "van_der_w.append(vdw)", space=myspace)
                except Exception as e:
                    print(f"Error: PyMOL iterate command failed for atom {a_name}: {e}")
                    print("Please check PyMOL installation and PDB file format.")
                    sys.exit(1)
                if not myspace["van_der_w"]: # Check if vdw_radius is retrieved from PyMOL
                    print(f"Error: Could not retrieve Van der Waals radius for atom {a_name} using PyMOL.")
                    print("Please check atom name and PDB file format.")
                    sys.exit(1)
                radius = myspace["van_der_w"][0]
                vdw_radii[a_element] = radius
    if pymol_pdb_loaded:
        pymol.cmd.do("delete all") # Conditionally delete all objects in PyMOL

    if verbose:
        print(f"Writing Van der Waals radii to file: {output_file}")

    try:
        with open(output_file, 'w') as f:
            for atom_name, radius in vdw_radii.items():
                f.write(f"{atom_name}\t{radius}\n")
    except Exception as e:
        print(f"Error: Failed to write Van der Waals radii to file {output_file}: {e}")
        sys.exit(1)

    end_time = time.time()
    elapsed_time = end_time - start_time

    if verbose:
        print(f"Van der Waals radii file created successfully: {output_file}")
        print(f"Elapsed time: {elapsed_time:.2f} seconds")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create Van der Waals radii file from PDB.')
    parser.add_argument('pdb_file', type=str, help='Path to the input PDB file.')
    parser.add_argument('output_file', type=str, help='Path to the output Van der Waals radii file.')
    parser.add_argument('--chain_id', type=str, default=None, help='Chain ID to extract atoms from. Defaults to None (all chains).')
    parser.add_argument('--verbose', action='store_true', help='Enable verbose output.')

    args = parser.parse_args()

    create_vdw_file(args.pdb_file, args.output_file, args.chain_id, args.verbose)
