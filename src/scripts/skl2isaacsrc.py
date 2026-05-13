#!/usr/bin/env python3

# Resource control file for skl2isaacs.py.
#
# This file defines the default parameter values for the
# skl2isaacs script.  A copy of this file is installed to
# $IMAGO_RC during the build process.  Users may place a
# modified copy in their current working directory to
# override the installed defaults on a per-project basis.
#
# To customize: copy this file into your working directory
# and edit the values below.  Command line arguments will
# still override any values set here.


def parameters_and_defaults():
    param_dict = {
        # String: Name of the input imago skeleton file.
        "input_file": "imago.skl",

        # String: Root name for the ISAACS output files.
        #   Two files are produced: <root>.ipf (ISAACS XML
        #   project file) and <root>.chem3d (Chem3D
        #   coordinate file).
        "output_root": "imago",
    }
    return param_dict


if __name__ == '__main__':
    print(parameters_and_defaults())
