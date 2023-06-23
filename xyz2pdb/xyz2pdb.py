#!/usr/bin/env python3
"""
Convert a XYZ file from Qiskit output into a PDB file
"""

__authors__ = ("Bryan Raubenolt (raubenb@ccf.org)",
               "Fabio Cumbo (cumbof@ccf.org)",
               "Jayadev Joshi (joshij@ccf.org)",
               "Daniel Blankenberg (blanked2@ccf.org)")

__version__ = "0.1.2"
__date__ = "Jun 23, 2023"

import argparse as ap
import errno
import os
from functools import partial
from typing import List, Union, Optional

import numpy as np

TOOL_ID = "qiskit-xyz2pdb"

# Residue names
# Used in case of --alpha-c-traces only
RESIDUES = {
    "A": "ALA",
    "R": "ARG",
    "N": "ASN",
    "D": "ASP",
    "C": "CYS",
    "E": "GLU",
    "Q": "GLN",
    "G": "GLY",
    "H": "HIS",
    "I": "ILE",
    "L": "LEU",
    "K": "LYS",
    "M": "MET",
    "F": "PHE",
    "P": "PRO",
    "S": "SER",
    "T": "THR",
    "W": "TRP",
    "Y": "TYR",
    "V": "VAL"
}


def read_params():
    """
    Read and test input arguments

    :return:     The ArgumentParser object
    """

    p = ap.ArgumentParser(
        prog=TOOL_ID,
        description="Convert a XYZ file from Qiskit output into a PDB file",
        formatter_class=ap.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "--in-xyz",
        type=os.path.abspath,
        required=True,
        dest="in_xyz",
        help="Path to the XYZ input file from Qiskit output",
    )
    p.add_argument(
        "--out-name",
        type=str,
        dest="out_name",
        help="Name of the output PDB file",
    )
    p.add_argument(
        "--out-folder",
        type=os.path.abspath,
        dest="out_folder",
        help="Path to the output folder",
    )
    p.add_argument(
        "--alpha-c-traces",
        action="store_true",
        default=False,
        dest="alpha_c_trace",
        help="Add C(alpha) traces",
    )
    p.add_argument(
        "--hetero-atoms",
        action="store_true",
        default=False,
        dest="hetero_atoms",
        help="Add hetero atoms",
    )
    p.add_argument(
        "-v",
        "--version",
        action="version",
        version="{} version {} ({})".format(TOOL_ID, __version__, __date__),
        help="Print the {} version and exit".format(TOOL_ID),
    )
    return p.parse_args()


def contains_coordinates(vector: List[str]) -> bool:
    """
    Check whether the input vector contains three orthogonal coordinates X, Y, and Z

    :param vector:  Input list
    :return:        True if the input list contains the three coordinates. False, otherwise
    """

    if len(vector) != 3:
        return False

    for value in vector:
        try:
            float(value)
        except:
            return False

    return True


def format_line(
    record_name: str,
    atom_serial_name: int,
    atom_name: str,
    alternate_location_indicator: str,
    residue_name: str,
    chain_identifier: str,
    residue_sequence_number: int,
    code_for_insertion_of_residues: str,
    orthogonal_coordinates_for_x: float,
    orthogonal_coordinates_for_y: float,
    orthogonal_coordinates_for_z: float,
    occupancy: float,
    temperature_factor: float,
    element_symbol: str,
    charge_on_the_atom: str
) -> str:
    """
    Format a PDB file
    https://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html

    :param record_name:                     Record name
    :param atom_serial_name:                Atom serial number
    :param atom_name:                       Atom name
    :param alternate_location_indicator:    Alternate location indicator
    :param residue_name:                    Residue name
    :param chain_identifier:                Chain identifier
    :param residue_sequence_number:         Residue sequence name
    :param code_for_insertion_of_residues:  Code for insertion of residues
    :param orthogonal_coordinates_for_x:    Orthogonal coordinates for X in Angstroms
    :param orthogonal_coordinates_for_y:    Orthogonal coordinates for Y in Angstroms
    :param orthogonal_coordinates_for_z:    Orthogonal coordinates for Z in Angstroms
    :param occupancy:                       Occupancy
    :param temperature_factor:              Temperature factor
    :param element_symbol:                  Element symbol
    :param charge_on_the_atom:              Charge on the atom
    :return:                                A PDB line
    """

    return "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}".format(
        record_name,
        atom_serial_name,
        atom_name,
        alternate_location_indicator,
        residue_name,
        chain_identifier,
        residue_sequence_number,
        code_for_insertion_of_residues,
        orthogonal_coordinates_for_x,
        orthogonal_coordinates_for_y,
        orthogonal_coordinates_for_z,
        occupancy,
        temperature_factor,
        element_symbol,
        charge_on_the_atom
    )


def load_xyz_data(xyz_filepath: os.path.abspath) -> np.ndarray:
    """
    Load a XYZ file as defined by qiskit

    :param xyz_filepath:    Path to the XYZ file
    :return:                A list of lists as numpy.ndarray with the content of the XYZ file
    """

    if not os.path.isfile(xyz_filepath):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), xyz_filepath)

    return np.array([line.strip().split(" ") for line in open(xyz_filepath).readlines() if line.strip() and len(line.strip().split(" ")) == 4])


def build_pdb(
    xyz_list: Union[List[List[Union[str, int, float]]], np.ndarray],
    out_pdb_name: Optional[str]=None,
    out_pdb_folder: Optional[os.path.abspath]=None,
    alpha_c_trace: bool=False,
    hetero_atoms: bool=False,
    replace: bool=False
) -> None:
    """
    Convert an XYZ list to a PDB file.
    The output PDB file name is optional. In case of None, it is automatically defined based on the series on aminoacids in xyz_list.
    The output folder path is also optional. In case of None, the output file is saved into the current working directory.

    :param xyz_list:        List of lists with the content of the XYZ file (see load_xyz_data)
    :param out_pdb_name:    Name of the output PDB file. Optional
    :param out_pdb_folder:  Path to the folder of the output PDB file. Optional
    :param alpha_c_trace:   Add C(alpha) traces
    :param hetero_atoms:    Add hetero atoms
    :param replace:         Overwrite the output file if it exists
    """

    if isinstance(xyz_list, np.ndarray):
        xyz_list = xyz_list.tolist()

    if not xyz_list:
        raise ValueError("The input XYZ is empty!")

    if not out_pdb_name:
        out_pdb_name = "".join([str(values[0]) for values in xyz_list])

    if not out_pdb_folder:
        out_pdb_folder = os.getcwd()

    out_pdb = os.path.join(out_pdb_folder, "{}.pdb".format(out_pdb_name))

    if os.path.isfile(out_pdb) and not replace:
        raise Exception("The output PDB file already exists")

    if (alpha_c_trace and hetero_atoms) or (not alpha_c_trace and not hetero_atoms):
        raise Exception(
            ("Please use one of the following flags: [--alpha-c-traces; --hetero-atoms]\n"
             "Use --help for additional details")
        )

    with open(out_pdb, "w+") as outfile:
        # Standard headers for the PDB file
        # https://www.wwpdb.org/documentation/file-format-content/format33/sect2.html#COMPND
        outfile.write(
            "{:6s} {:3s}{:70s}\n".format(
                "COMPND",
                "",
                "My Protein Sequence"
            )
        )
        # https://www.wwpdb.org/documentation/file-format-content/format33/sect2.html#AUTHOR
        outfile.write(
            "{:6s}  {:2s}{:69s}\n".format(
                "AUTHOR",
                "",
                "Galaxy User"
            )
        )

        # Count the number of residues
        num_res = 1
        
        for xyz_line in xyz_list:
            if xyz_line:
                if len(xyz_line) != 4 or (len(xyz_line) == 4 and not contains_coordinates(xyz_line[1:])):
                    # Raise an exception in case the input line does not contain 4 columns
                    # or in case of invalid coordinates
                    raise Exception("Malformed input XYZ file")

                format_line_partial = partial(
                    format_line,
                    atom_serial_name=num_res,
                    alternate_location_indicator="",
                    chain_identifier="",
                    residue_sequence_number=num_res,
                    code_for_insertion_of_residues="",
                    orthogonal_coordinates_for_x=float(xyz_line[1]),
                    orthogonal_coordinates_for_y=float(xyz_line[2]),
                    orthogonal_coordinates_for_z=float(xyz_line[3]),
                    occupancy=1.0,
                    temperature_factor=0.0,
                    element_symbol=xyz_line[0],
                    charge_on_the_atom=""
                )

                if alpha_c_trace:
                    # https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
                    outfile.write(
                        "{}\n".format(
                            format_line_partial(
                                record_name="ATOM",
                                atom_name="CA",
                                residue_name=RESIDUES[xyz_line[0]]
                            )
                        )
                    )
                
                elif hetero_atoms:
                    # https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#HETATM
                    outfile.write(
                        "{}\n".format(
                            format_line_partial(
                                record_name="HETATM",
                                atom_name=xyz_line[0],
                                residue_name="PEP"
                            )
                        )
                    )

                num_res += 1

        if hetero_atoms:
            # Add CONECT lines
            # https://www.wwpdb.org/documentation/file-format-content/format33/sect10.html#CONECT
            for i in range(1, num_res - 1):
                outfile.write(
                    "{:6s}{:5d}{:5d}{:5d}{:5d}{:5d}\n".format(
                        "CONECT",   # Record name
                        i,          # Atom serial number
                        i + 1,      # Serial number of bonded atom
                        0,          # Serial number of bonded atom
                        0,          # Serial number of bonded atom
                        0,          # Serial number of bonded atom
                    )
                )

            # Add MASTER line
            # https://www.wwpdb.org/documentation/file-format-content/format33/sect11.html#MASTER
            outfile.write(
                "{:6s}    {:5d}{:5d}{:5d}{:5d}{:5d}{:5d}{:5d}{:5d}{:5d}{:5d}{:5d}{:5d}\n".format(
                    "MASTER",       # Record name
                    0,              # Number of REMARK records
                    0,
                    0,              # Number of HET records
                    0,              # Number of HELIX records
                    0,              # Number of SHEET records
                    0,              # **Deprecated**
                    0,              # Number of SITE records
                    0,              # Number of coordinate transformation records (ORIGX+SCALE+MTRIX)
                    num_res - 1,    # Number of atomic coordinate records (ATOM+HETATM)
                    0,              # Number of TER records
                    num_res - 1,    # Number of CONECT records
                    0               # Number of SEQRES records
                )
            )

        outfile.write("END")


def main() -> None:
    args = read_params()

    # Load the XYZ file into a list of lists
    xyz_list = load_xyz_data(args.in_xyz)

    # Convert the XYZ list to PDB
    build_pdb(
        xyz_list,
        out_pdb_name=args.out_name,
        out_pdb_folder=args.out_folder,
        alpha_c_trace=args.alpha_c_trace,
        hetero_atoms=args.hetero_atoms
    )


if __name__ == "__main__":
    main()
