#!/usr/bin/env python3
## SPDX-License-Identifier: ECL-2.0
## Copyright (c) 2026 Paul Rulis


import os


def parameters_and_defaults():
    param_list = [
            "gs_scfV-fb.dat", # Default imago SCF potential file.
            "structure.dat", # Default imago structure file. 
            False # Assume that we are not doing a spin-down calculation.
            ]
    return param_list


if __name__ == '__main__':
    print(parameters_and_defaults())
