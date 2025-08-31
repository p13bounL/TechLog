from __future__ import annotations
from typing import List, Optional, Union
from pathlib import Path
import os, sys, platform, importlib, shutil

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import hdf5plugin

# -------- Intro helpers --------
def check_versions() -> dict:
    def v(mod: str):
        try: return importlib.import_module(mod).__version__
        except Exception: return "not installed"
    return {
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "numpy": v("numpy"),
        "pandas": v("pandas"),
        "matplotlib": v("matplotlib"),
        "scanpy": v("scanpy"),
        "anndata": v("anndata"),
        "hdf5plugin": v("hdf5plugin"),
        "decoupler": v("decoupler"),
        "scikit-learn": v("sklearn"),
    }

def list_dir(path: str = ".") -> list[str]:
    return sorted(os.listdir(path))

def disk_free(path: str = ".") -> float:
    total, used, free = shutil.disk_usage(path)
    return round(free / (1024**3), 2)

# -------- Raw → H5AD --------
def _looks_like_10x_dir(p: Path) -> bool:
    names = {x.name.lower() for x in p.glob("*")}
    has_mtx = any(n.startswith("matrix.mtx") for n in names)
    has_barc = any(n.startswith("barcodes.tsv") for n in names)
    has_feat = any(n.startswith("features.tsv") or n.startswith("genes.tsv") for n in names)
    return has_mtx and has_barc and has_feat

from pathlib import Path
import shutil, tempfile
import anndata as ad
import scanpy as sc

def read_10x_mtx_folder(folder: str, var_names: str = "gene_symbols", condition: str | None = None) -> ad.AnnData:
    p = Path(folder)
    if not p.exists():
        raise FileNotFoundError(f"Ο φάκελος δεν υπάρχει: {p}")

    # Βρες το matrix (με ή χωρίς .gz)
    matrix_files = list(p.glob("*matrix.mtx.gz")) + list(p.glob("*matrix.mtx"))
    if not matrix_files:
        raise FileNotFoundError("Δεν βρέθηκε αρχείο matrix.mtx ή matrix.mtx.gz μέσα στον φάκελο.")

    mtx = sorted(matrix_files, key=lambda x: len(x.name))[0]
    name = mtx.name
    gz = name.endswith(".mtx.gz")

    # Προσδιόρισε prefix (ό,τι υπάρχει πριν από το 'matrix.mtx')
    pivot = "matrix.mtx"
    idx = name.find(pivot)
    prefix = name[:idx] if idx >= 0 else ""
    suffix = ".gz" if gz else ""

    # Έλεγξε ότι υπάρχουν barcodes & features/genes με το ΙΔΙΟ prefix
    def has_file(basename: str) -> bool:
        return (p / f"{prefix}{basename}").exists() or (p / f"{prefix}{basename}.gz").exists()

    if not has_file("barcodes.tsv"):
        raise FileNotFoundError("Λείπει barcodes.tsv(.gz) με το ίδιο prefix.")
    if not (has_file("features.tsv") or has_file("genes.tsv")):
        raise FileNotFoundError("Λείπει features.tsv(.gz) ή genes.tsv(.gz) με το ίδιο prefix.")

    # 1η προσπάθεια: νέα Scanpy (δέχεται prefix/suffix)
    try:
        adata = sc.read_10x_mtx(p, var_names=var_names, make_unique=True, prefix=prefix, suffix=suffix)
    except TypeError:
        # Παλαιότερη Scanpy: φτιάξε προσωρινό φάκελο με canonical ονόματα και διάβασέ τον
        tmpdir = Path(tempfile.mkdtemp(prefix="tenx_"))
        shutil.copyfile(mtx, tmpdir / f"matrix.mtx{suffix}")

        # Αντιγραφή σχετικών αρχείων με canonical ονόματα
        def copy_first(candidates: list[Path], dst: Path) -> None:
            for f in candidates:
                if f.exists():
                    shutil.copyfile(f, dst)
                    return
            raise FileNotFoundError(f"Δεν βρέθηκε: {dst.name}")

        # barcodes
        copy_first(
            [p / f"{prefix}barcodes.tsv{suffix}", p / f"{prefix}barcodes.tsv"],
            tmpdir / f"barcodes.tsv{suffix}"
        )
        # features ή genes
        feats = [p / f"{prefix}features.tsv{suffix}", p / f"{prefix}features.tsv",
                 p / f"{prefix}genes.tsv{suffix}",    p / f"{prefix}genes.tsv"]
        copy_first(feats, tmpdir / f"features.tsv{suffix}")

        adata = sc.read_10x_mtx(tmpdir, var_names=var_names, make_unique=True)

    if condition is not None:
        adata.obs["condition"] = condition
    return adata


def read_10x_h5_file(h5_path: str, condition: str | None = None) -> ad.AnnData:
    p = Path(h5_path)
    if not p.is_file():
        raise FileNotFoundError(f"Δεν βρέθηκε αρχείο: {h5_path}")
    adata = sc.read_10x_h5(str(p))
    if condition is not None:
        adata.obs["condition"] = condition
    return adata

def read_10x_auto(path: str, var_names: str = "gene_symbols", condition: str | None = None) -> ad.AnnData:
    p = Path(path)
    if p.is_dir():
        return read_10x_mtx_folder(str(p), var_names=var_names, condition=condition)
    if p.is_file() and p.name.lower().startswith("matrix.mtx"):
        return read_10x_mtx_folder(str(p.parent), var_names=var_names, condition=condition)
    if p.suffix.lower() == ".h5ad" and p.is_file():
        adata = sc.read_h5ad(str(p))
        if condition is not None:
            adata.obs["condition"] = condition
        return adata
    if p.suffix.lower() in {".h5", ".hdf5"}:
        return read_10x_h5_file(str(p), condition=condition)
    raise ValueError("Δώσε φάκελο 10x mtx ή αρχείο .h5 (CellRanger) ή .h5ad.")

def save_h5ad(adata: ad.AnnData, out_path: str) -> str:
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out_path)
    return out_path

def raw2h5ad(mtx_folder: str, out_path: str, condition: str = "disease", var_names: str = "gene_symbols") -> str:
    adata = read_10x_mtx_folder(mtx_folder, var_names=var_names, condition=condition)
    return save_h5ad(adata, out_path)

# -------- Preprocessing (batch) --------
def _ensure_dir(p: str | Path) -> None:
    Path(p).mkdir(parents=True, exist_ok=True)

def _compression_arg(name: str | None):
    if name is None: return None
    name = str(name).lower()
    if name == "zstd": return hdf5plugin.FILTERS["zstd"]
    if name in {"gzip", "lzf"}: return name
    return hdf5plugin.FILTERS["zstd"]

def preprocess_h5ad_file(
    input_path: str | Path,
    output_path: str | Path,
    n_genes_min: int = 1000,
    n_genes_max: int = 10000,
    compression: str = "zstd",
) -> str:
    try:
        import adata_preprocessor as ap
        adata_filtered = ap.adata_preprocessor(str(input_path), n_genes_min=n_genes_min, n_genes_max=n_genes_max)
    except Exception:
        adata = sc.read_h5ad(str(input_path))
        sc.pp.calculate_qc_metrics(adata, inplace=True)
        mask = (adata.obs["n_genes_by_counts"] >= n_genes_min) & (adata.obs["n_genes_by_counts"] <= n_genes_max)
        adata_filtered = adata[mask].copy()
    _ensure_dir(Path(output_path).parent)
    comp = _compression_arg(compression)
    adata_filtered.write_h5ad(str(output_path), compression=comp)
    return str(output_path)

def batch_preprocess_h5ad_dir(
    input_dir: str | Path,
    output_dir: str | Path,
    n_genes_min: int = 1000,
    n_genes_max: int = 10000,
    compression: str = "zstd",
) -> list[str]:
    input_dir = Path(input_dir); output_dir = Path(output_dir); _ensure_dir(output_dir)
    written: list[str] = []
    for file in sorted(os.listdir(input_dir)):
        if file.endswith(".h5ad"):
            in_path = input_dir / file
            out_path = output_dir / f"filtered_{file}"
            saved = preprocess_h5ad_file(in_path, out_path, n_genes_min=n_genes_min, n_genes_max=n_genes_max, compression=compression)
            written.append(saved)
    return written

def preprocess_pipeline(
    adata: ad.AnnData,
    min_counts_cells: int = 1000,
    min_cells_genes: int = 3,
    n_top_genes: int = 2000,
    n_neighbors: int = 15,
    n_pcs: int = 40,
) -> ad.AnnData:
    adata = adata.copy()
    sc.pp.filter_cells(adata, min_counts=min_counts_cells)
    sc.pp.filter_genes(adata, min_cells=min_cells_genes)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, subset=True)
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver="arpack")
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    sc.tl.umap(adata)
    return adata

# -------- Concatenation --------
def concatenate_adatas(adatas: List[ad.AnnData]) -> ad.AnnData:
    return ad.concat(adatas, join="outer", label="batch", fill_value=0)

# -------- Integration (placeholder) --------
def integrate_scanorama(adata: ad.AnnData) -> ad.AnnData:
    return adata  # βάλε εδώ τον πραγματικό σου κώδικα

# -------- UMAP helpers --------
def ensure_umap(adata: ad.AnnData, n_neighbors: int = 15, n_pcs: int = 40) -> ad.AnnData:
    adata = adata.copy()
    if "X_pca" not in adata.obsm: sc.tl.pca(adata)
    if "connectivities" not in adata.obsp: sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    if "X_umap" not in adata.obsm: sc.tl.umap(adata)
    return adata

def umap_figure(adata: ad.AnnData, color: Optional[Union[str, list[str]]] = None, point_size: int = 10):
    sc.pl.umap(adata, color=color, size=point_size, show=False)
    return plt.gcf()

# -------- Annotation --------
def annotate_cells(adata: ad.AnnData, resolution: float = 1.0) -> ad.AnnData:
    adata = adata.copy()
    candidates = ["cell_type", "celltype", "CellType", "annotation", "labels", "celltype_major", "leiden", "louvain", "cluster"]
    for col in candidates:
        if col in adata.obs:
            series = adata.obs[col].astype(str)
            if series.nunique() > 1 and not series.str.lower().eq("unknown").all():
                adata.obs["cell_type"] = series
                return adata
    if "X_pca" not in adata.obsm: sc.tl.pca(adata)
    if "connectivities" not in adata.obsp: sc.pp.neighbors(adata, n_neighbors=15)
    sc.tl.leiden(adata, key_added="cluster", resolution=resolution)
    adata.obs["cell_type"] = adata.obs["cluster"].astype(str).map(lambda x: f"cluster_{x}")
    return adata

# -------- DEG 2 groups --------
def deg_two_groups(
    adata: ad.AnnData,
    groupby: str,
    group1,
    group2,
    method: str = "wilcoxon",
) -> pd.DataFrame:
    """
    Διαφορική έκφραση group1 vs group2 (robust):
    - Χωρίζει ΔΕΝ βασίζεται σε stringification των τιμών.
    - Φιλτράρει σε δύο ομάδες, ρίχνει NaN, θέτει κατηγορική στήλη σωστά.
    - Επιστρέφει columns: gene, log2FC, pval, padj (ό,τι λείπει -> NaN).
    """
    if group1 == group2:
        raise ValueError("Οι δύο ομάδες πρέπει να είναι διαφορετικές.")

    # 1) Drop NaN στη στήλη groupby
    if groupby not in adata.obs:
        raise KeyError(f"Το '{groupby}' δεν υπάρχει στο obs.")
    ad_ok = adata[adata.obs[groupby].notna()].copy()
    if ad_ok.n_obs < 3:
        raise ValueError("Πολύ λίγα κελιά μετά το φιλτράρισμα NaN.")

    # 2) Κράτα ΜΟΝΟ τις δύο ομάδες
    mask_2 = ad_ok.obs[groupby].isin([group1, group2])
    ad2 = ad_ok[mask_2].copy()
    if ad2.obs[groupby].nunique() < 2:
        raise ValueError("Μία από τις ομάδες δεν βρέθηκε στο dataset.")

    # 3) Κάνε την groupby κατηγορική με ακριβείς τιμές (όχι strings)
    if not pd.api.types.is_categorical_dtype(ad2.obs[groupby]):
        ad2.obs[groupby] = pd.Categorical(ad2.obs[groupby])
    # αφαίρεσε αχρησιμοποίητες κατηγορίες και βεβαιώσου ότι περιλαμβάνει μόνο τις 2 ζητούμενες
    ad2.obs[groupby] = ad2.obs[groupby].cat.remove_unused_categories()

    # 4) Ελάχιστο μέγεθος ομάδας
    n1 = int((ad2.obs[groupby] == group1).sum())
    n2 = int((ad2.obs[groupby] == group2).sum())
    if n1 < 2 or n2 < 2:
        raise ValueError(f"Πολύ μικρές ομάδες για DEG: A={n1}, B={n2}.")

    # 5) Τρέξε Scanpy DEG
    sc.tl.rank_genes_groups(
        ad2,
        groupby=groupby,
        groups=[group1],
        reference=group2,
        method=method,
    )

    # 6) Πάρε αποτελέσματα. Προσπαθούμε πρώτα το επίσημο API
    try:
        df_all = sc.get.rank_genes_groups_df(ad2)  # όλα τα groups
        df = df_all[df_all["group"].astype(object) == group1].copy()
    except Exception:
        # Fallback: διαχείριση dict/recarray/ndarray από ad2.uns['rank_genes_groups']
        rg = ad2.uns["rank_genes_groups"]

        def _field(rg_obj, key):
            try:
                if isinstance(rg_obj, dict):
                    return rg_obj.get(key, None)
                # structured/recarray: index με όνομα πεδίου
                return rg_obj[key]
            except Exception:
                return None

        def _col(arr, gval):
            """Επιστρέφει Series με τιμές για το συγκεκριμένο group, από structured/2D/1D array."""
            if arr is None:
                return None
            a = np.asarray(arr)
            # Structured array με ονομαστικά πεδία (π.χ. ονόματα groups)
            if a.dtype.names:
                # Αν υπάρχει πεδίο με ΑΚΡΙΒΗ όνομα του group (όχι str(group))
                if gval in a.dtype.names:
                    return pd.Series(a[gval])
                # αλλιώς δοκίμασε stringified
                if str(gval) in a.dtype.names:
                    return pd.Series(a[str(gval)])
                # αλλιώς πάρε την πρώτη στήλη ως fallback
                return pd.Series(a[a.dtype.names[0]])
            # Κανονικό ndarray
            if a.ndim == 2:
                # Αν υπάρχει λίστα column labels; διαφορετικά επέλεξε τη στήλη που αντιστοιχεί στο group
                # Θα προσπαθήσουμε να βρούμε ποια στήλη ταιριάζει συγκρίνοντας τα "names"
                return pd.Series(a[:, 0])
            if a.ndim == 1:
                return pd.Series(a)
            return None

        names_col = _col(_field(rg, "names"), group1)
        lfc_col   = _col(_field(rg, "logfoldchanges"), group1)
        p_col     = _col(_field(rg, "pvals"), group1)
        padj_col  = _col(_field(rg, "pvals_adj"), group1)

    # 7) Κανονικοποίηση/renames αν ήρθαν διαφορετικά
    if 'df' not in locals():
        df = pd.DataFrame()
        if names_col is not None: df["gene"] = names_col.astype(str)
        else:
            n = 0
            for c in (lfc_col, p_col, padj_col):
                if c is not None:
                    n = len(c)
                    break
            df["gene"] = pd.RangeIndex(n).astype(str)

        df["log2FC"] = pd.to_numeric(lfc_col, errors="coerce") if lfc_col is not None else np.nan
        df["pval"]   = pd.to_numeric(p_col,   errors="coerce") if p_col   is not None else np.nan
        df["padj"]   = pd.to_numeric(padj_col,errors="coerce") if padj_col is not None else np.nan
    else:
        # Αν το πήραμε από sc.get.rank_genes_groups_df
        rename_map = {}
        if "names" in df.columns and "gene" not in df.columns: rename_map["names"] = "gene"
        if "logfoldchanges" in df.columns and "log2FC" not in df.columns: rename_map["logfoldchanges"] = "log2FC"
        if "pvals" in df.columns and "pval" not in df.columns: rename_map["pvals"] = "pval"
        if "pvals_adj" in df.columns and "padj" not in df.columns: rename_map["pvals_adj"] = "padj"
        if rename_map: df = df.rename(columns=rename_map)
        for c in ("log2FC", "pval", "padj"):
            if c not in df.columns: df[c] = np.nan
        df["gene"] = df["gene"].astype(str)
        for c in ("log2FC", "pval", "padj"):
            df[c] = pd.to_numeric(df[c], errors="coerce")

    return df[["gene", "log2FC", "pval", "padj"]]





# -------- Volcano --------
def volcano_figure(
    df: pd.DataFrame, log2fc_col: str = "log2FC", p_col: str = "padj",
    fc_thr: float = 1.0, p_thr: float = 0.05, annotate_top: int = 0,
):
    d = df.copy()
    d[log2fc_col] = pd.to_numeric(d[log2fc_col], errors="coerce")
    d[p_col]      = pd.to_numeric(d[p_col], errors="coerce")
    d = d.dropna(subset=[log2fc_col, p_col])
    d["neglog10p"] = -np.log10(np.clip(d[p_col].values, 1e-300, 1.0))
    x, y = d[log2fc_col].values, d["neglog10p"].values
    sig_mask  = (np.abs(x) >= fc_thr) & (d[p_col].values <= p_thr)
    up_mask   = (x >=  fc_thr) & (d[p_col].values <= p_thr)
    down_mask = (x <= -fc_thr) & (d[p_col].values <= p_thr)

    fig, ax = plt.subplots()
    ax.scatter(x[~sig_mask], y[~sig_mask], s=8,  alpha=0.3, label="μη-σημαντικά")
    ax.scatter(x[up_mask],   y[up_mask],   s=10, alpha=0.9, label="σημαντικά (up)")
    ax.scatter(x[down_mask], y[down_mask], s=10, alpha=0.9, label="σημαντικά (down)")
    p_thr_safe = max(float(p_thr), 1e-300)
    ax.axvline(fc_thr, linestyle="--"); ax.axvline(-fc_thr, linestyle="--"); ax.axhline(-np.log10(p_thr_safe), linestyle="--")
    ax.set_xlabel("log2 Fold Change"); ax.set_ylabel("-log10(padj)")
    ax.set_title(f"Volcano — σημαντικά: {int(sig_mask.sum())} (↑{int(up_mask.sum())}, ↓{int(down_mask.sum())})")
    ax.legend(loc="best")

    if annotate_top and sig_mask.any() and "gene" in d.columns:
        sig_df = d.loc[sig_mask].sort_values(p_col).head(int(annotate_top))
        for _, row in sig_df.iterrows():
            ax.annotate(str(row.get("gene","")), (row[log2fc_col], row["neglog10p"]),
                        xytext=(3,3), textcoords="offset points", fontsize=8)
    return fig

# -------- Expression plots --------
def violin_expression_figure(adata: ad.AnnData, gene: str, groupby: Optional[str] = None):
    sc.pl.violin(adata, keys=[gene], groupby=groupby, show=False)
    return plt.gcf()
