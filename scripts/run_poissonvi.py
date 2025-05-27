import os
import scanpy as sc
import anndata
import scvi
from argparse import ArgumentParser
from pathlib import Path
import yaml
scvi.settings.seed = 0


def parse():
    parser = ArgumentParser()
    parser.add_argument("--data_path", dest="data_path")
    parser.add_argument("--model_name", dest="model_name")
    parser.add_argument("--params_path", dest="params_path", default=None)
    return parser.parse_args()


def load_data(peakmat_path, filter_ratio):
    atac_data = sc.read_h5ad(peakmat_path)
    sc.pp.filter_genes(atac_data, min_cells=int(atac_data.shape[0] * filter_ratio))
    return atac_data


class INTEGRATION:
    def __init__(self, datapath, model_name, params_dic):
        self.model_name = model_name
        self.params_dic = params_dic
        self.datapath = datapath
        self.ad = None
        self.model = None

    def load_data(self):
        self.ad = sc.read_h5ad(self.datapath)
        self.ad.obs_names_make_unique()

    def run_poissonvi(self):
        scvi.external.POISSONVI.setup_anndata(
            self.ad, 
            batch_key=self.params_dic['batch_key'], 
            categorical_covariate_keys=self.params_dic['categorical_covariate_keys']
        )
            
        self.model = scvi.external.POISSONVI(self.ad)
        self.model.train(max_epochs=250)
        self.model.save("mod.pt", overwrite=True)
        self.ad.obsm["poissonvi"] = self.model.get_latent_representation()
        sc.pp.neighbors(self.ad, use_rep="poissonvi")
        sc.tl.umap(self.ad, min_dist=0.2)
        self.ad.write(f"adat.h5ad")


    def run_multivi(self):
        scvi.model.MULTIVI.setup_anndata(self.ad, batch_key="modality")
        model = scvi.model.MULTIVI(
            self.ad,
            n_genes=(self.ad.var["modality"] == "gene").sum(),
            n_regions=(self.ad.var["modality"] == "peak").sum(),
        )
        # model.view_anndata_setup()
        self.model.train(max_epochs=250)
        self.model.save("mod.pt", overwrite=True)
        # model = scvi.external.POISSONMULTIVI.load(model_dir, adata=adata_mvi)
        self.ad.obsm["multivi"] = self.model.get_latent_representation()
        sc.pp.neighbors(self.ad, use_rep="multivi")
        sc.tl.umap(self.ad, min_dist=0.2)
        sc.tl.leiden(self.ad, resolution=0.5, key_added='leiden0.5')
        self.ad.write(f'adat.h5ad')

    def run_multimap(self):
        pass

    def run_model(self):
        if 'poissonvi' in self.model_name:
            self.run_poissonvi()
        elif 'multivi' in self.model_name:
            self.run_multivi()
        elif 'multimap' in self.model_name:
            self.run_multimap()
        else:
            raise IOError(f"Not recognized model_name f{self.model_name}")



def main():
    args = parse()
    if not os.path.exists(args.params_path):
        params_dic = {'batch_key': 'kit', 
                      'categorical_covariate_keys': ['sample'],
                      'filter_ratio': 0.1}
    else:
        params_dic = yaml.safe_load(open(args.params_path, "r"))
    data_path = args.data_path
    model_name = args.model_name
    Integration = INTEGRATION(
        data_path, 
        model_name, 
        params_dic)
    Integration.load_data(filter_ratio=args.params_dic['filter_ratio'])
    Integration.run_model()


if __name__ == '__main__':
    main()