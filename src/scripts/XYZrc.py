#!/usr/bin/env python3
## SPDX-License-Identifier: ECL-2.0
## Copyright (c) 2026 Paul Rulis


import os


def parameters_and_defaults():
    param_dict = {
            "some" : "Foo", # String
            "thing" : "Bar", # String
            "more" : False, # T (F) = Property is (isn't) true.
            "other" : True, # T (F) = Property is (isn't) true.
            "things" : False, # T (F) = Property is (isn't) true.
            "final_thing" : [1.1, 2.2]  # Real list = Set of numbers.
            }
    return param_dict


if __name__ == '__main__':
    print(parameters_and_defaults())
