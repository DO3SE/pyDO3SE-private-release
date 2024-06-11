
import os
from pathlib import Path
from distutils.dir_util import copy_tree

import re
from typing import List
import pypandoc
from warnings import warn


EXAMPLE_DIR_OUT = "docs/examples"
EXAMPLE_DIR_IN = "examples"
setups = [
    ["soil_moisture"],
]

if __name__ == "__main__":
    print('Copying example files')

    os.makedirs(EXAMPLE_DIR_OUT)
    for id in setups:
        copy_tree(f"{EXAMPLE_DIR_IN}/{id}", f"{EXAMPLE_DIR_OUT}/{id}")


