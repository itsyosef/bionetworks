#!/usr/env/python3

import os
import sys
from pathlib import Path

path = Path(sys.argv[1]) / "input/wigs/"

if len(os.listdir(path / "tex_notex")) > 0:
    path = path / "tex_notex"
else:
    path = path / "fragment"

lib_names = []

for file in os.listdir(path):
    new_name = file
    number = file.split(".")[-1]
    if number == "wig": number = "1"
    if "pos" in file:
        new_name += f":tex:{number}:a"
    else:
        new_name += f":notex:{number}:a"

    if "forward" in file:
        new_name += ":+"
    else:
        new_name += ":-"

    lib_names.append(str(path / new_name))
    
print(" ".join(lib_names))