#!/usr/bin/env python3

import sys, os
from pathlib import Path

base_path = Path(sys.argv[1])

wigs_path = base_path / "input" / "wigs" / "fragment"

samples = list(set([i.split("_")[0] for i in os.listdir(wigs_path)]))

format_str = ""

for ix, sample in enumerate(samples):
    ix = ix + 1
    format_str += f" {wigs_path / sample}_forward.wig:frag:{ix}:a:+ {wigs_path / sample}_reverse.wig:frag:{ix}:a:-"
print(format_str)