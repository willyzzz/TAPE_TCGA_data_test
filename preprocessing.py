import pandas as pd
import anndata as ad
import numpy as np

# 读取bulk数据
bulk_data = pd.read_csv('./or_bulk_nature_2012.txt', delimiter='\t', index_col=0)
bulk_data = bulk_data.T
bulk_data = bulk_data.fillna(0.0)

# 处理负值和重复列
min_value = bulk_data.min().min()
offset = -min_value + 0.001
bulk_data = bulk_data + offset
bulk_data = bulk_data.loc[:, ~bulk_data.columns.duplicated()]

# 读取单细胞数据
adata = ad.read_h5ad('single_cell.h5ad')
cell_type_series = adata.obs['celltype_major']
data_matrix = pd.DataFrame(adata.X.toarray(), index=cell_type_series, columns=adata.var['feature_name'], dtype=np.float32)
sc_ref = data_matrix

# 转换数据为数值类型并处理缺失值
sc_ref = sc_ref.apply(pd.to_numeric, errors='coerce')
sc_ref.dropna(axis=1, inplace=True)

# 筛选共同的列
common_columns = bulk_data.columns.intersection(sc_ref.columns)
sc_ref = sc_ref[common_columns]
bulk_data = bulk_data[common_columns]


# 计算每种细胞类型的样本数量
cell_type_counts = sc_ref.index.value_counts()
total_cells = 10000
sample_counts = (cell_type_counts / cell_type_counts.sum() * total_cells).astype(int)

# 采样所需数量的每种细胞类型的细胞
sampled_sc_ref = pd.concat([sc_ref.loc[cell_type].sample(n=count, random_state=42, replace=True)
                            for cell_type, count in sample_counts.items()])

# 保存采样后的单细胞参考数据和bulk数据到CSV
sampled_sc_ref.to_csv('./example_data/example_counts.txt')
bulk_data.T.to_csv('./example_data/bulk_data.txt')

print(f"Original cell type distribution:\n{cell_type_counts}")
print(f"Sampled cell type distribution:\n{sampled_sc_ref.index.value_counts()}")
