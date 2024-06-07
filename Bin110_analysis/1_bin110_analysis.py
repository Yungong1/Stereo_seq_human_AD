import stereo as st
import os
import pandas as pd
import gc
import utils

######### read the data as ann h5ad files
control_list = ["B01806B5", "B01806B6", "B01809A3", "B01809A4", "D02175A4", "D02175A6"]
moderate_list = ["B02008D2", "B02009F6"]
advanced_list = ["B01809C2", "C02248B5"]
severe_list = ["A02092E1", "B02008C6"]

###### Control
# import the stereo data
for num,files in enumerate(control_list):
    data_path = "../../processed_data/{}/GeneExpMatrix/{}.tissue.gef".format(files, files)
    data = st.io.read_gef(file_path=data_path, bin_size=110)
    data.cells["diagnosis"] = "control"
    data.cells["sample"] = files
    globals()[files] = data
    
##### moderate
# import the stereo data
for num,files in enumerate(moderate_list):
    data_path = "../../processed_data/{}/GeneExpMatrix/{}.tissue.gef".format(files, files)
    data = st.io.read_gef(file_path=data_path, bin_size=110)
    data.cells["diagnosis"] = "moderate"
    data.cells["sample"] = files
    globals()[files] = data

#### advanced
for num,files in enumerate(advanced_list):
    data_path = "../../processed_data/{}/GeneExpMatrix/{}.tissue.gef".format(files, files)
    data = st.io.read_gef(file_path=data_path, bin_size=110)
    data.cells["diagnosis"] = "advanced"
    data.cells["sample"] = files
    globals()[files] = data

#### severe_list
for num,files in enumerate(severe_list):
    data_path = "../../processed_data/{}/GeneExpMatrix/{}.tissue.gef".format(files, files)
    data = st.io.read_gef(file_path=data_path, bin_size=110)
    data.cells["diagnosis"] = "severe"
    data.cells["sample"] = files
    globals()[files] = data

# merge the data
merge = st.utils.data_helper.merge(B01806B5, B01806B6, B01809A3, B01809A4, D02175A4, D02175A6,
                                  A02092E1, B02009F6, C02248B5, B02008C6, B02008D2, B01809C2)

utils.cell_t0_res_key(merge)
merge.tl.raw_checkpoint()

merge.tl.normalize_total()
merge.tl.log1p()
merge.tl.pca(use_highly_genes=False, n_pcs=50, res_key='pca')

merge.tl.batches_integrate(pca_res_key='pca', res_key='pca_integrated')

merge.tl.neighbors(
        pca_res_key='pca_integrated',
        n_pcs=30,
        res_key='neighbors'
        )

merge.tl.spatial_neighbors(
        neighbors_res_key='neighbors',
        res_key='spatial_neighbors'
        )

merge.tl.umap(pca_res_key='pca', neighbors_res_key='spatial_neighbors', res_key='umap')

merge.tl.leiden(neighbors_res_key='spatial_neighbors', res_key='spatial_leiden', resolution=1)
merge.plt.cluster_scatter(res_key='spatial_leiden', reorganize_coordinate = 2)

merge.plt.umap(cluster_key = "spatial_leiden")

# save the multiple file
st.io.write_h5ad(
        merge,
        use_raw=True,
        use_result=True,
        key_record=None,
        output='../Result/Stereo/Mult_file/integration.h5ad'
        )

# save one file
merge.cells.to_df().to_csv("../Result/Stereo/One_file/110_integrate_meta.csv")
st.io.write_h5ad(
        merge,
        use_raw=True,
        use_result=True,
        key_record=None,
        output='../Result/Stereo/One_file/integration.h5ad',
        split_batches=False
    
        )