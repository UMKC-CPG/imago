#!/usr/bin/env python3
## SPDX-License-Identifier: ECL-2.0
## Copyright (c) 2026 Paul Rulis


import os


def parameters_and_defaults():
    param_dict = {
            # File name parameters. Empty string means use the
            #   automatically constructed default name.
            "dat_file"  : "",  # Str = setup .dat file name.
            "out_file"  : "",  # Str = SYBD .out file name.
            "raw_file"  : "",  # Str = SYBD .raw file name.
            "plot_file" : "",  # Str = Output .plot file name.

            # Naming components used to construct default
            #   file names when the above are empty.
            "basis" : "fb",  # Str = Basis tag (fb, mb, eb).
            "edge"  : "gs",  # Str = Edge tag (gs, es, etc).
            }
    return param_dict


if __name__ == '__main__':
    print(parameters_and_defaults())
