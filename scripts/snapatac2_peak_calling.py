import snapatac2 as snap
import scanpy as sc
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from argparse import ArgumentParser
from pathlib import Path


def concat_atac(snapatac_dir, concated_filename='concated.h5ads'):
    adatas = []
    for filename in os.listdir(snapatac_dir):
        if filename.endswith('h5ad') and 'genematrix' not in filename:
            adatas.append(
                [filename.replace('.h5ad', ''),
                 snap.read(snapatac_dir/filename, backed='r')]
                )
    data = snap.AnnDataSet(
        adatas=[(name, ad) for name, ad in adatas],
        filename=f"{concated_filename}"
    )
    # clustering on binned results
    snap.pp.select_features(data, n_features=50000)
    snap.tl.spectral(data)
    snap.pp.mnc_correct(data, batch="sample")
    snap.pp.harmony(data, batch="sample", max_iter_harmony=20)
    snap.tl.umap(data, use_rep='X_spectral_mnn', key_added='umap_mnn')
    snap.tl.umap(data, use_rep='X_spectral_harmony', key_added='umap_harmony')
    snap.pp.knn(data, use_rep='X_spectral_harmony')
    snap.tl.leiden(data)
    data.close()


def peak_calling(concated_filepath, peakmat_savepath):
    data = snap.read_dataset(concated_filepath)
    # peak calling grouped by leiden clusters
    snap.tl.macs3(data, groupby='leiden')
    merged_peaks = snap.tl.merge_peaks(data.uns['macs3'], chrom_sizes=snap.genome.hg38)
    peak_mat = snap.pp.make_peak_matrix(data, use_rep=merged_peaks['Peaks'])
    peak_mat.write(peakmat_savepath)
    data.close()


def parse():
    parser = ArgumentParser()
    parser.add_argument("--snapatac_dir", dest="snapatac_dir")
    parser.add_argument("--concated_filename", dest="concated_filename")
    parser.add_argument("--peakmat_savepath", dest="peakmat_savepath")
    return parser.parse_args()


def main():
    args = parse()
    concat_atac(
        snapatac_dir=Path(args.snapatac_dir), 
        concated_filename=args.concated_filename
        )
    concated_filepath = f"./{args.concated_filename}"
    peakmat_filepath = f"./{args.peakmat_savepath}"
    peak_calling(
        concated_filepath=concated_filepath, 
        peakmat_savepath=peakmat_filepath
        )


if __name__ =='__main__':
    main()
