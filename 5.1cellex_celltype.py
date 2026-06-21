###############################################################################
# Step 5.1: Run cellex
###############################################################################
import numpy as np 
import pandas as pd

import cellex
import scanpy as sc

data = pd.read_csv("scn_combine_filter_harmony_soybean_count.csv", index_col=0 engine="pyarrow",memory_map=True,low_memory=False)
metadata = pd.read_csv("scn_combine_filter_harmony_soybean_metadata.csv", index_col=0)

eso = cellex.ESObject(data=data, annotation=metadata, verbose=True)
eso.compute(verbose=True)
eso.results["esmu"].to_csv("cellex_celltype.esmu.csv.gz")