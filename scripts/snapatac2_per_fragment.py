import snapatac2 as snap
import sys
from pathlib import Path
import os
import scanpy as sc
from argparse import ArgumentParser

def parse():
    parser = ArgumentParser()
    parser.add_argument("--fragments_dir", dest="fragments_dir")
    parser.add_argument("--name", dest="name")
    return parser.parse_args()


def process_one_fragmentfile(name, fragment_dir):
    # create temp folder
    if not os.path.exists(fragment_dir/f"temp/{name}"):
        os.makedirs(fragment_dir/f"temp/{name}")
    # read in fragment file
    samplepath = Path(fragment_dir/f'{name}/fragments.tsv.gz')
    assert os.path.exists(samplepath)
    # snapatac2 preprocess
    data = snap.pp.import_data(
        samplepath,
        file=f'{name}.h5ad',
        chrom_sizes=snap.genome.hg38,
        min_num_fragments=1000,
        sorted_by_barcode=False,
        tempdir=fragment_dir/f'temp/{name}'
    )
    snap.metrics.tsse(data, snap.genome.hg38)
    snap.pp.filter_cells(data, min_counts=5000, min_tsse=10, max_counts=100000)
    snap.pp.add_tile_matrix(data)
    snap.pp.select_features(data, n_features=250000)
    snap.pp.scrublet(data)
    snap.pp.filter_doublets(data)
    # create gene matrix
    gene_matrix = snap.pp.make_gene_matrix(data, snap.genome.hg38)
    sc.pp.filter_genes(gene_matrix, min_cells=5)
    sc.pp.normalize_total(gene_matrix)
    sc.pp.log1p(gene_matrix)
    sc.external.pp.magic(gene_matrix, solver="approximate")
    data.close()
    # save gene matrix
    gene_matrix.write(f"{name}.genematrix.h5ad", compression='gzip')


def main():
    args = parse()
    fragment_dir = Path(args.fragments_dir)
    process_one_fragmentfile(
        name=args.name.strip(), 
        fragment_dir=fragment_dir
    )

if __name__ == '__main__':
    main()
