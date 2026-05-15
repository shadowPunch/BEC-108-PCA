import gzip
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA


#1 = ER+, 0 = ER-
class_df = pd.read_csv(
    "data/class.tsv",
    sep="\t",
    header=None,
    names=["label"]
)
labels = class_df["label"].values


expr_raw = pd.read_csv(
    "data/filtered.tsv.gz",
    sep="\t",
    compression="gzip",
    index_col=0
)

print("Original matrix shape:", expr_raw.shape)
expr_df = expr_raw.copy()

#gene IDs to clean strings
expr_df.columns = expr_df.columns.map(lambda x: str(x).strip())

print("Expression matrix shape:", expr_df.shape)

print("\nFirst 10 gene IDs:")
print(expr_df.columns[:10])

col_df = pd.read_csv(
    "data/columns.tsv.gz",
    sep="\t",
    comment="#",
    compression="gzip",
    low_memory=False
)


col_df = col_df[["ID", "GeneSymbol"]]

# Remove missing names
col_df = col_df.dropna(subset=["GeneSymbol"])

col_df["ID"] = col_df["ID"].astype(str)
col_df["GeneSymbol"] = col_df["GeneSymbol"].str.strip()

# Build mapping
id_to_name = dict(zip(col_df["ID"], col_df["GeneSymbol"]))
name_to_id = {v: k for k, v in id_to_name.items()}

#Extract XBP1 and GATA3
xbp1_id = name_to_id.get("XBP1")
gata3_id = name_to_id.get("GATA3")

print("XBP1 ID :", xbp1_id)
print("GATA3 ID:", gata3_id)

if xbp1_id is None or gata3_id is None:
    raise ValueError("Could not find XBP1 or GATA3 IDs")

# Directly use gene IDs as column names
xbp1_expr = expr_df[str(xbp1_id)].values.astype(float)
gata3_expr = expr_df[str(gata3_id)].values.astype(float)


#PCA on FULL expression matrix
# rows = samples, cols = genes
X_full = expr_df.values.astype(float)

# Standardize each gene
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X_full)

# PCA
pca = PCA(n_components=2)
X_pca = pca.fit_transform(X_scaled)

# PC scores
pc1_scores = X_pca[:, 0]
pc2_scores = X_pca[:, 1]

print("\nExplained variance ratio:")
print("PC1:", pca.explained_variance_ratio_[0])
print("PC2:", pca.explained_variance_ratio_[1])


# 6. PCA on 2D XBP1/GATA3 data
X2 = np.column_stack([gata3_expr, xbp1_expr])

scaler2 = StandardScaler()
X2_scaled = scaler2.fit_transform(X2)

pca2 = PCA(n_components=2)
pca2.fit(X2_scaled)

v1 = pca2.components_[0]
v2 = pca2.components_[1]

fig, axes = plt.subplots(1, 3, figsize=(16, 5))

ax = axes[0]

mask_pos = labels == 1
mask_neg = labels == 0

ax.scatter(
    gata3_expr[mask_neg],
    xbp1_expr[mask_neg],
    c="black",
    s=20,
    label="ER-"
)

ax.scatter(
    gata3_expr[mask_pos],
    xbp1_expr[mask_pos],
    c="red",
    s=20,
    label="ER+"
)

ax.set_xlabel("GATA3", fontsize=16, fontstyle="italic")
ax.set_ylabel("XBP1", fontsize=16, fontstyle="italic")

ax.set_title("Figure 1a")

ax = axes[1]

ax.scatter(
    gata3_expr[mask_neg],
    xbp1_expr[mask_neg],
    c="black",
    s=20
)

ax.scatter(
    gata3_expr[mask_pos],
    xbp1_expr[mask_pos],
    c="red",
    s=20
)

# Draw PC directions
scale = 4

ax.arrow(
    0, 0,
    v1[0] * scale,
    v1[1] * scale,
    color="black",
    width=0.02,
    head_width=0.15,
    length_includes_head=True
)

ax.arrow(
    0, 0,
    v2[0] * scale,
    v2[1] * scale,
    color="black",
    width=0.02,
    head_width=0.15,
    length_includes_head=True
)

ax.text(
    v1[0] * scale * 1.1,
    v1[1] * scale * 1.1,
    "PC1",
    fontsize=14
)

ax.text(
    v2[0] * scale * 1.1,
    v2[1] * scale * 1.1,
    "PC2",
    fontsize=14
)

ax.set_xlabel("GATA3", fontsize=16, fontstyle="italic")
ax.set_ylabel("XBP1", fontsize=16, fontstyle="italic")

ax.set_title("Figure 1b")

ax = axes[2]

# y-levels
y_all = np.ones(len(pc1_scores)) * 2
y_neg = np.ones(np.sum(mask_neg)) * 1
y_pos = np.ones(np.sum(mask_pos)) * 0

# All samples
ax.scatter(
    pc1_scores[mask_neg],
    y_all[mask_neg],
    c="black",
    s=20
)

ax.scatter(
    pc1_scores[mask_pos],
    y_all[mask_pos],
    c="red",
    s=20
)

# ER-
ax.scatter(
    pc1_scores[mask_neg],
    y_neg,
    c="black",
    s=20
)

# ER+
ax.scatter(
    pc1_scores[mask_pos],
    y_pos,
    c="red",
    s=20
)

ax.set_yticks([2, 1, 0])
ax.set_yticklabels(["All", "ER-", "ER+"])

ax.set_xlabel("Projection onto PC1", fontsize=14)

ax.set_title("Figure 1c")

plt.tight_layout()

plt.savefig(
    "figure1_reproduction.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()

print("\nSaved figure to: figure1_reproduction.png")
