# make_plots.py (SIMPLE & COMPATIBLE)
import os, sys
import scanpy as sc
import matplotlib.pyplot as plt
import hdf5plugin  # ΑΠΑΡΑΙΤΗΤΟ πριν από αναγνώσεις h5ad με plugins

# ---- Ρυθμίσεις εισόδου/εξόδου ----
infile = sys.argv[1] if len(sys.argv) > 1 else "data/filtered.h5ad"
outdir = "figs"
os.makedirs(outdir, exist_ok=True)

print(f"[i] Διαβάζω: {infile}")
adata = sc.read(infile)

# ---- QC metrics αν λείπουν ----
need_qc = not {"total_counts", "n_genes_by_counts"}.issubset(adata.obs.columns) or "pct_counts_mt" not in adata.obs.columns
if need_qc:
    adata.var["mt"] = adata.var_names.str.upper().str.startswith(("MT-", "MT."))
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

# ---- Ανάλυση αν λείπει ----
if "X_pca" not in adata.obsm_keys():
    print("[i] Τρέχω normalize/log1p/HVGs/PCA...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, subset=True)
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver="arpack")

if "neighbors" not in adata.uns:
    sc.pp.neighbors(adata, n_pcs=30)

if "X_umap" not in adata.obsm_keys():
    sc.tl.umap(adata)

if "leiden" not in adata.obs.columns:
    print("[i] Τρέχω Leiden (χρειάζονται python-igraph & leidenalg)...")
    sc.tl.leiden(adata, resolution=1.0)

# Κατηγορική στήλη κλάστερ
if adata.obs["leiden"].dtype.name != "category":
    adata.obs["leiden"] = adata.obs["leiden"].astype("category")

# ---- 1) QC violin ----
print("[i] Φτιάχνω QC violin...")
sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    groupby="leiden",
    multi_panel=True,
    show=False,      # ΜΟΝΟ show=False
)
fig = plt.gcf()      # παίρνουμε το τρέχον figure
fig.savefig(os.path.join(outdir, "qc_violin.png"), dpi=200, bbox_inches="tight")
plt.close(fig)

# ---- 2) UMAP ----
print("[i] Φτιάχνω UMAP...")
sc.pl.umap(
    adata,
    color=["leiden"],
    # legend_loc παραλείπεται για μέγιστη συμβατότητα
    show=False,
)
fig = plt.gcf()
fig.savefig(os.path.join(outdir, "umap.png"), dpi=200, bbox_inches="tight")
plt.close(fig)

# ---- 3) Heatmap με top markers ανά κλάστερ ----
print("[i] Υπολογίζω markers & heatmap...")
sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon")

# Προσπαθούμε με το helper (νέες εκδόσεις), αλλιώς fallback
marker_list = []
top_n = 3
try:
    df = sc.get.rank_genes_groups_df(adata, None)
    for g in adata.obs["leiden"].cat.categories:
        genes = df[df["group"] == g].head(top_n)["names"].tolist()
        for gene in genes:
            if gene not in marker_list:
                marker_list.append(gene)
except Exception:
    names = adata.uns["rank_genes_groups"]["names"]
    groups = list(adata.obs["leiden"].cat.categories)
    for idx, g in enumerate(groups):
        # names είναι structured array ή 2D list-like
        arr = names[idx] if hasattr(names, "__getitem__") else []
        for i in range(min(top_n, len(arr))):
            gene = arr[i]
            if gene not in marker_list:
                marker_list.append(gene)

sc.pl.heatmap(
    adata,
    var_names=marker_list,
    groupby="leiden",
    swap_axes=True,
    cmap="viridis",
    show=False,
)
fig = plt.gcf()
fig.savefig(os.path.join(outdir, "markers_heatmap.png"), dpi=200, bbox_inches="tight")
plt.close(fig)

# (προαιρετικό) σώσε και το επεξεργασμένο h5ad
os.makedirs("data", exist_ok=True)
adata.write("data/analysis_with_umap_leiden.h5ad")

print(f"[✓] Έτοιμα:\n - {os.path.join(outdir,'qc_violin.png')}\n - {os.path.join(outdir,'umap.png')}\n - {os.path.join(outdir,'markers_heatmap.png')}")
