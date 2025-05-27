import snapatac2 as snap
import scanpy as sc
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from argparse import ArgumentParser
from pathlib import Path
from snapatac2_peak_calling import peak_calling
import sys


if __name__ == '__main__':
    assert os.path.exists(Path(sys.argv[1])/"fragments.tsv.gz")


