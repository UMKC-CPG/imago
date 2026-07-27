#!/usr/bin/env python3
## SPDX-License-Identifier: ECL-2.0
## Copyright (c) 2026 Paul Rulis


import os


def parameters_and_defaults():
    param_dict = {
            "bz" : 1,  # -bz : Which BZ to show
            "real_alpha" : 0.1,  # -rla : Set the real cell alpha.
            "recip_alpha" : 0.4,  # -rca : Set the recip cell alpha.
            "bz_alpha" : 0.4,  # -bza : Set the bz cell alpha.
            "mesh_kp_alpha" : 0.4,  # -mkpa : Set the mesh kpoints alpha.
            "fold_kp_alpha" : 0.4,  # -fkpa : Set the folded kpoints alpha.
            "path_kp_alpha" : 0.4  # -pkpa : Set the path kpoints alpha.
            }
    return param_dict


if __name__ == '__main__':
    print(parameters_and_defaults())
