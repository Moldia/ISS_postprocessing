from skimage import (
     color, feature, filters, measure, morphology, segmentation, util
)
import os
import pandas as pd
import numpy as np
import skimage.color
from scipy.sparse import coo_matrix
from scipy.sparse import load_npz, save_npz
from skimage.measure import label, regionprops
import scanpy as sc 
import matplotlib.pyplot as plt
from skimage.segmentation import watershed, expand_labels
from matplotlib.pyplot import rc_context
import matplotlib as mpl

import scanpy as sc
import numpy as np
from matplotlib.pyplot import rc_context
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import os
#import tangram as tg
import scanpy as sc
import pandas as pd
#import squidpy as sq
from math import ceil
# Show plots as part of the notebook

import io
import seaborn as sns

import warnings
warnings.filterwarnings("ignore")
import numpy as np

import numpy as np
import scanpy as sc
import pandas as pd
from sklearn.metrics.pairwise import euclidean_distances

def get_object_info(segementation_result):
    region_proporties = measure.regionprops(segementation_result)
    area = []
    index = []
    x = []
    y = []
    for i in range(len(measure.regionprops(segementation_result))):         
        centroid_intermediate = region_proporties[i].centroid
        centroid_intermediate = list(centroid_intermediate)
        area_intermediate = region_proporties[i].area    
        x.append(centroid_intermediate[1])
        y.append(centroid_intermediate[0])
        area.append(area_intermediate)
        index.append(i)    # create dataframe      
    cell_info = pd.DataFrame(index)
    cell_info['x'] = x
    cell_info['y'] = y
    cell_info['area'] = area
    return cell_info

def assign_spots_to_cells(segmentation_labels, spots):
    from scipy import ndimage as ndi    
    spots1 = spots[["y", "x"]]    
    cell_labels = ndi.map_coordinates(
          segmentation_labels,
          spots1.T,  # assuming spot coords has shape (n, 2)
          order=0,
          )
    spots["cell"] = cell_labels
    return spots

def Diff(li1, li2):
    return list(set(li1) - set(li2)) + list(set(li2) - set(li1))


def create_anndata_obj(spots_file, 
    segmentation_mask, 
    output_file,
    filter_data=True, 
    metric = 'distance', 
    write_h5ad = True,
    value= 1.2,
    convert_coords = True, 
    conversion_factor = 0.1625): 
    print('reading spots file')

    spots = pd.read_csv(spots_file)
    
    if filter_data==True:
        spots_filtered = spots[spots[metric] < value]
    else: 
        spots_filtered = spots
    
    spots_filtered = spots_filtered[['target','xc','yc']]
    
    if convert_coords == True: 
        spots_filtered['x'] = spots_filtered['xc']/conversion_factor
        spots_filtered['y'] = spots_filtered['yc']/conversion_factor
    else: 
        spots_filtered['x'] = spots_filtered['xc']
        spots_filtered['y'] = spots_filtered['yc']
        
    spots_filtered = spots_filtered.rename(columns = {'target':'Gene'})
    spots_filtered = spots_filtered[['Gene','x','y']]
    spots_filtered = spots_filtered.dropna()
    
    coo = load_npz(segmentation_mask)
    print('load coo file')
    assinged = assign_spots_to_cells(coo.toarray(), spots_filtered)
    cells = get_object_info(coo.toarray())
    
    cells[0] = cells[0]+1
    print('assign spots to cells')
    assigned_filt = assinged[assinged.cell != 0]
    hm = assinged.groupby(['Gene','cell']).size().unstack(fill_value=0)
    hm = hm.drop(columns = 0)

    an_sp = sc.AnnData(X=hm.T)
    cells[0]  = cells[0].astype(str)
    cells_filt = cells[cells[0].isin(list(an_sp.obs.index))]
    an_sp.obs = cells_filt
    an_sp.obs = an_sp.obs.drop(columns = 0)
    cells[0].astype(int)-1
    if write_h5ad == True:
        print('write h5ad')
        an_sp.write_h5ad(output_file)
    else: 
        print('not writing')
    
    return an_sp




def recluster_specific_cluster(anndata, 
                               to_cluster, 
                               rerun_umap = False, 
                               resolutions = [0.1,0.2,0.3,0.5]):
    import scanpy as sc
    
    to_cluster = [to_cluster]
    anndata_int = anndata[anndata.obs.cell_type.isin(to_cluster)]

    sc.pp.neighbors(anndata_int, n_neighbors=30, n_pcs=30)
    
    if rerun_umap == True:
        sc.tl.umap(anndata_int, min_dist=1)

    for i in resolutions:
        print('clustering at resolution: '+str(i))
        plt.rcdefaults()
        sc.tl.leiden(anndata_int, resolution =i, key_added = ("cell_type_" + str(i)))

        plt.rcdefaults()
        with rc_context({'figure.figsize': (10, 10), 'figure.dpi': 50}):
            sc.pl.umap(anndata_int, color = ("cell_type_"+str(i)),s=30,legend_loc='on data',legend_fontsize=20,legend_fontoutline=10)
            
    return anndata_int


def plot_umap(anndata,
              color = 'cell_type',
              compute_umap=False, 
              n_neighbors=30, 
              n_pcs=30,
              min_dist=1,
              fig_size = (10, 10),
              fig_dpi = 50, 
              s=20,
              legend_loc='on data',
              legend_fontsize=15,
              legend_fontoutline=10
             ): 
    
    if compute_umap == True: 
        sc.pp.neighbors(anndata, n_neighbors=n_neighbors, n_pcs=n_pcs)
        sc.tl.umap(anndata, min_dist=min_dist)

    plt.rcdefaults()
    with rc_context({'figure.figsize': fig_size, 'figure.dpi': fig_dpi}):
        sc.pl.umap(anndata, color = (color),s=s,legend_loc=legend_loc,legend_fontsize=legend_fontsize,legend_fontoutline=legend_fontoutline, frameon = False, title = ' ')

def plot_marker_genes(anndata, 
                    cluster_label, 
                    method='t-test', 
                    key_added = "t-test",
                    n_genes=25, sharey=False, key = "t-test"):
    sc.tl.rank_genes_groups(anndata, cluster_label, method=method, key_added = key_added)
    plt.rcdefaults()
    sc.pl.rank_genes_groups(anndata, n_genes=n_genes, sharey=sharey, key = key_added)

def plot_clusters(anndata, 
                    clusters_to_map, 
                    broad_cluster,
                    key='t-test', 
                    size = 0.5,
                    number_of_marker_genes = 10, 
                    sample_id_column = 'sample_id', 
                    dim_subplots = [3,3]
                 ): 
    
    mpl.rcParams['text.color'] = 'w'
    plt.style.use('dark_background')

    clusters_to_map = clusters_to_map
    cluster_class = {}
    marker_genes = {}
    for broad in sorted(list(anndata.obs[broad_cluster].unique())): 
        anndata_broad = anndata[anndata.obs[broad_cluster] == broad]
        print('  ')
        print(broad)
        for cluster in sorted(list(anndata_broad.obs[clusters_to_map].unique())):
            print(cluster)
            genes = list(sc.get.rank_genes_groups_df(anndata_broad,group=str(cluster), key=key)['names'].head(number_of_marker_genes))
            print(*list(genes), sep = " ")
            spatial_int = anndata_broad[anndata_broad.obs[clusters_to_map] == str(cluster)]
            fig, axs = plt.subplots(dim_subplots[0],ceil(len(anndata_broad.obs['sample'].unique())/dim_subplots[1]), figsize=(20, 10))
            fig.subplots_adjust(hspace = .5, wspace=.001)
            fig.suptitle('Cluster: '+str(cluster))
            axs = axs.ravel()


            for q, j in enumerate(sorted(list(anndata_broad.obs[sample_id_column].unique()))):


                spatial_celltypes_tag_ = spatial_int[spatial_int.obs[sample_id_column]==j]
                axs[q].plot((anndata[anndata.obs[sample_id_column] == j].obs.y), (anndata[anndata.obs[sample_id_column] == j].obs.x), marker='s', linestyle='', ms=size, color = 'grey', alpha = 0.2)
                axs[q].plot(spatial_celltypes_tag_.obs.x, spatial_celltypes_tag_.obs.y, marker='s', linestyle='', ms=size, color = 'yellow')#spatial_int.uns['leiden_0.4_colors'][0])


            plt.show()

def spatial_neigbourshoods(anndata,
                           cluster_label = 'leiden_0.5',
                           max_distance_allowed = 300, 
                           umap_dist = 1,
                           leiden_resolution = 0.2
                           ):
    
    
    distances_input=np.array([anndata.obs['x'],anndata.obs['y']])
    din=distances_input.transpose()
    distances=euclidean_distances(din, din)
    dist_df=pd.DataFrame(distances)
    max_distance_allowed=max_distance_allowed
    dist_binary=((dist_df<max_distance_allowed)*1)*((dist_df!=0)*1)
    np.sum(np.sum(dist_binary))
    dist_binary['name']=list(anndata.obs[cluster_label])
    distbinsum=dist_binary.groupby('name').sum()
    adata=sc.AnnData(distbinsum.transpose())
    adata.obs=anndata.obs
    sc.tl.umap(adata,min_dist=umap_dist)
    sc.tl.leiden(adata,resolution=leiden_resolution, key_added = 'local_neighborhood')
    return adata
