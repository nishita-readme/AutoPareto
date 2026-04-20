import scanpy as sc
import pandas as pd
import numpy as np
import anndata as ad
import partipy as pt
import matplotlib.pyplot as plt
import scanpy.external as sce
import scipy.sparse as sp
from scipy.spatial.distance import cdist
import gseapy as gp
from IPython.display import display
import plotly.io as pio
from itertools import permutations
from itertools import combinations

def check_raw_integers_in_adataX(adata, tol = 1e-6):
    """Check if adata.X contains only integers (or close to integers for sparse matrices).
    
    Parameters:
        adata : AnnData object
        tol : float, tolerance for floating point rounding

    Returns:
        str : message indicating whether adata.X is integer"""
    
    X = adata.X
    
    if X is None:
            return "adata.X is not present"
        
    if sp.issparse(X):
        # Find first non-integer in the stored data
        for idx, val in enumerate(X.data):
            if not val.is_integer():
                row, col = X.nonzero()[0][idx], X.nonzero()[1][idx]
                return f"❌ non-integer found at ({row}, {col}) = {val}"
        return "adata.X contains only integers ✅"
    else:
        # Dense matrix
        rows, cols = X.shape
        for i in range(rows):
            for j in range(cols):
                if not float(X[i, j]).is_integer():
                    return f"❌ Non-integer found at ({i}, {j}) = {X[i, j]}"
        return "adata.X contains only integers ✅"
    
    
def preprocess_adata(adata, exclude_quality_genes=False, custom_exclude_genes=None, n_pcs=50, pca_seed=123,
    apply_harmony=False, batch_key="Sample"):
    """
    Preprocess an AnnData object in-place.
    
    Steps:
    -----------
    - Normalize total counts per cell
    - log-transform
    - Select highly variable genes
    - Optionally exclude quality-associated genes or custom list
    - Scaling
    - Run PCA
    - Optionally run Harmony batch correction
    

    Parameters:
    -----------
    adata : AnnData
        AnnData object to preprocess
    exclude_quality_genes : bool
        Whether to exclude predefined quality-associated genes (default False)
    custom_exclude_genes : list of str
        Custom list of gene names to exclude (overrides quality genes)
    n_pcs : int
        Number of principal components to compute
    pca_seed : int
        add seed to pca function (default 123)
    apply_harmony : bool
        Whether to perform Harmony batch correction (default False)
    batch_key : str
        Column in adata.obs containing batch info (default "Sample")
        

    Notes:
    ------
    - Harmony requires the batch column to exist in adata.obs. The default batch_key is "Sample"
    """

    # Normalization
    sc.pp.normalize_total(adata)
    print("adata.X normalized")
    
    # log1p transformation 
    sc.pp.log1p(adata)
    print("adata.X log transformed")
    
    # Highly variable genes
    sc.pp.highly_variable_genes(adata)
    print("Highly variable genes retrieved")
    
    np.random.seed(0)

    # Exclude quality genes
    if exclude_quality_genes or (custom_exclude_genes is not None):
        if custom_exclude_genes is not None:
            genes_to_exclude = custom_exclude_genes
        else:
            genes_to_exclude = set()  # Replace with actual lab list

        if genes_to_exclude:
            adata.var['highly_variable'] = adata.var['highly_variable'] & (~adata.var_names.isin(genes_to_exclude))
            print("QC-associated genes excluded from highly variable gene list")
    
    
    # Scaling
    adata.layers["z_scaled"] = sc.pp.scale(adata.X, copy=True)
    print("adata scaled and stored in layer z_scaled")
    
    # PCA
    sc.pp.pca(adata, mask_var="highly_variable")
    print("PCA ✅")

    # Harmony (optional)
    if apply_harmony:
        if batch_key not in adata.obs:
            raise ValueError(f"Batch key '{batch_key}' not found in .obs slot")
        sce.pp.harmony_integrate(adata, key=batch_key)
        print("Batch correction performed using harmony")

    print("Preprocessing complete. AnnData object updated in-place.")
    
def get_partipy_selection_metrics(adata, n_dims, k_min=2, k_max=10):
    """
    Compute archetype selection metrics over a range of k values.

    Sets the PCA embedding in partipy and computes selection metrics
    (e.g., RSS, variance explained) for each k in [k_min, k_max].
    Results are stored in adata and used downstream for archetype
    number selection.

    Parameters
    ----------
    adata : AnnData
        Annotated data object containing PCA coordinates in obsm['X_pca'].
    n_dims : int
        Number of PCA dimensions to use for archetype fitting.
    k_min : int, optional
        Minimum number of archetypes to evaluate. Default is 2.
    k_max : int, optional
        Maximum number of archetypes to evaluate. Default is 10.

    Returns
    -------
    adata : AnnData
        Updated AnnData object with selection metrics stored in adata.uns.

    Example
    -------
    >>> adata = get_partipy_selection_metrics(adata, n_dims=20, k_min=3, k_max=7)
    >>> plots_for_n_archetypes_selection(adata, n_archetype_range=range(3, 8))
    """
    pt.set_obsm(adata=adata, obsm_key="X_pca", n_dimensions=n_dims)
    pt.compute_selection_metrics(
        adata=adata, 
        n_archetypes_list=list(range(k_min, k_max + 1))
    )
    return adata
    
def plots_for_n_archetypes_selection(adata, n_archetype_range=range(3, 7), color=None):
    """
    Generate diagnostic plots to aid selection of the optimal number of archetypes.

    Produces the following outputs:
        1. Variance explained plot across archetype solutions
        2. Information criterion (IC) plot across archetype solutions
        3. Bootstrap variance summary across all n (stability plot)
        4. Bar plot of t-ratio & RSS significance p-values by number of archetypes
        5. 2D archetype scatter plots for each n in n_archetype_range
        6. Per-n bootstrap 2D plots for each n in n_archetype_range

    Parameters
    ----------
    adata : anndata.AnnData
        Annotated data object. Must have PCA in adata.obsm and AA_selection_metrics
        in adata.uns (run pt.compute_selection_metrics first).
    n_archetype_range : range or list of int, default range(3, 7)
        The number of archetypes to evaluate. Must start at 3+ (t-ratio requires > 2).
    color : str or None
        Column in adata.obs to color cells by in 2D scatter plots (e.g. "CellType").
    """

    n_archetype_range = list(n_archetype_range)
    assert min(n_archetype_range) >= 3, "n_archetype_range must start at 3 (t-ratio requires > 2)"

    # ── 1. Variance explained + IC ───────────────────────────────────────────
    print("=== Variance Explained ===")
    display(pt.plot_var_explained(adata).draw())

    print("=== Information Criterion ===")
    display(pt.plot_IC(adata).draw())

    # ── 2. Compute archetypes ────────────────────────────────────────────────
    print("Computing archetypes...")
    for i in n_archetype_range:
        print(f"  n_archetypes={i}")
        pt.compute_archetypes(adata, n_archetypes=i, delta=0.5, archetypes_only=False)

    # ── 3. Bootstrap variance ────────────────────────────────────────────────
    print("Computing bootstrap variance...")
    pt.compute_bootstrap_variance(
        adata=adata,
        n_bootstrap=50,
        n_archetypes_list=n_archetype_range,
        delta=0.5
    )

    # ── 4. Bootstrap variance summary plot (across all n) ───────────────────
    print("=== Bootstrap Variance Summary (stability across n) ===")
    display(pt.plot_bootstrap_variance(adata, summary_method="median").draw())

    # ── 5. T-ratio & RSS bar plot ────────────────────────────────────────────
    print("Computing t-ratio significance...")
    p_val_dict = {}
    for i in n_archetype_range:
        sig = pt.t_ratio_significance(adata, result_filters={"n_archetypes": i, "delta": 0.5})
        p_val_dict[i] = sig

    x       = list(p_val_dict.keys())
    x_pos   = list(range(len(x)))
    width   = 0.35
    t_pvals = [p_val_dict[i]['t_ratio_p_value'] for i in x]
    r_pvals = [p_val_dict[i]['rss_p_value']     for i in x]

    print("=== T-ratio & RSS Significance ===")
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.bar([p - width/2 for p in x_pos], t_pvals, width, label='T-ratio p-value', color='steelblue')
    ax.bar([p + width/2 for p in x_pos], r_pvals, width, label='RSS p-value',     color='coral')
    ax.axhline(y=0.05, color='red', linestyle='--', label='p = 0.05')
    ax.set_xlabel("Number of Archetypes")
    ax.set_ylabel("P-value")
    ax.set_title("T-ratio & RSS Significance by Number of Archetypes")
    ax.set_xticks(x_pos)
    ax.set_xticklabels(x)
    ax.legend()
    plt.tight_layout()
    plt.show()

    # ── 6. Per-n 2D scatter + bootstrap 2D ──────────────────────────────────
    for i in n_archetype_range:
        filters = {"n_archetypes": i, "delta": 0.5}
        print(f"=== 2D Plots: n_archetypes={i} ===")
        display(pt.plot_archetypes_2D(adata, color=color, result_filters=filters).draw())


    
def plot_archetypes_3D_range(adata, n_archetype_range=range(3, 7), color=None):
    """
    Plot 3D archetype scatter plots for a range of archetype solutions.

    Parameters
    ----------
    adata : anndata.AnnData
        Annotated data object with AA results in adata.uns["AA_results"].
    n_archetype_range : range or list of int, default range(3, 7)
        The number of archetypes to plot.
    color : str or None
        Column in adata.obs to color cells by (e.g. "CellType").
    """
    pio.renderers.default = "notebook"  # force render in Jupyter

    for i in list(n_archetype_range):
        print(f"n_archetypes={i}")
        fig = pt.plot_archetypes_3D(
            adata,
            color=color,
            result_filters={"n_archetypes": i, "delta": 0.5}
        )
        fig.show()
        
        
def get_top_cells_per_archetype(adata, n_archetypes, top_n=200, n_dims=8):
    """
    Annotate the top_n closest cells to each archetype in PC space.

    Parameters
    ----------
    adata : anndata.AnnData
        Must have PCA in adata.obsm["X_pca"] and AA results in adata.uns["AA_results"].
    n_archetypes : int
        Number of archetypes to use.
    top_n : int, default 200
        Number of top cells to annotate per archetype.
    n_dims : int, default 8
        Number of PC dimensions to use when computing distances.

    Returns
    -------
    adata : anndata.AnnData
        Same object with "archetype" column added to adata.obs.
        Values are 1-indexed archetype number for top cells, 0 for all other cells.
    """
    aa_result = pt.get_aa_result(adata, n_archetypes=n_archetypes, delta=0.5)
    Z = aa_result["Z"]  # (n_archetypes, n_pca_dims)

    X = adata.obsm["X_pca"][:, :Z.shape[1]]  # match dims of Z

    dists = cdist(X, Z, metric="euclidean")   # (n_cells, n_archetypes)

    archetype_labels = np.zeros(adata.n_obs, dtype=int)
    for arch_idx in range(n_archetypes):
        top_idx = np.argsort(dists[:, arch_idx])[:top_n]
        archetype_labels[top_idx] = arch_idx + 1

    adata.obs["archetype"] = archetype_labels
    return adata

def plot_top_cells_per_archetype(adata, color_col="archetype", dims=(0, 1)):
    """
    Plot top cells per archetype on a 2D PCA scatter.
    Non-assigned cells (archetype=0) are shown in grey in the background.

    Parameters
    ----------
    adata : anndata.AnnData
        Must have PCA in adata.obsm["X_pca"] and "archetype" in adata.obs
        (run get_top_cells_per_archetype first).
    color_col : str, default "archetype"
        Column in adata.obs to color by.
    dims : tuple of int, default (0, 1)
        Which PC dimensions to plot (0-indexed).
    """


    pc_x, pc_y = dims
    df = pd.DataFrame({
        f"PC{pc_x + 1}": adata.obsm["X_pca"][:, pc_x],
        f"PC{pc_y + 1}": adata.obsm["X_pca"][:, pc_y],
        "archetype":     adata.obs["archetype"].values
    })

    n_archetypes = df["archetype"].max()
    colors = plt.cm.tab10.colors  # up to 10 archetypes

    fig, ax = plt.subplots(figsize=(8, 6))

    # background: unassigned cells in grey
    bg = df[df["archetype"] == 0]
    ax.scatter(bg[f"PC{pc_x + 1}"], bg[f"PC{pc_y + 1}"],
               c="lightgrey", s=8, alpha=0.4, label="unassigned", zorder=1)

    # foreground: top cells per archetype
    for arch in range(1, n_archetypes + 1):
        sub = df[df["archetype"] == arch]
        ax.scatter(sub[f"PC{pc_x + 1}"], sub[f"PC{pc_y + 1}"],
                   c=[colors[(arch - 1) % len(colors)]],
                   s=15, alpha=0.8, label=f"Archetype {arch}", zorder=2)

    ax.set_xlabel(f"PC{pc_x + 1}")
    ax.set_ylabel(f"PC{pc_y + 1}")
    ax.set_title(f"Top cells per archetype (PC{pc_x + 1} vs PC{pc_y + 1})")
    ax.legend(markerscale=2, bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()
    plt.show()

def run_deg_per_archetype(adata, lfc_threshold=1.0, pval_threshold=0.05):
    """
    Subset to top archetype-assigned cells and run DEG (1 vs rest) for each archetype.

    Parameters
    ----------
    adata : anndata.AnnData
        Must have "archetype" column in adata.obs (run get_top_cells_per_archetype first).
        Values of 0 are unassigned and excluded.
    lfc_threshold : float, default 1.0
        Minimum absolute log fold change to keep.
    pval_threshold : float, default 0.05
        Maximum adjusted p-value to keep.

    Returns
    -------
    dict : {archetype (int): pd.DataFrame of significant DEGs sorted by score}
    """

    # subset to only assigned cells
    adata_sub = adata[adata.obs["archetype"] != 0].copy()
    adata_sub.obs["archetype"] = adata_sub.obs["archetype"].astype(str)

    print(f"Cells after subsetting: {adata_sub.n_obs}")
    print(adata_sub.obs["archetype"].value_counts())

    # run wilcoxon 1 vs rest across all genes
    sc.tl.rank_genes_groups(
        adata_sub,
        groupby="archetype",
        method="wilcoxon",
        reference="rest",
        n_genes=adata_sub.n_vars,  # test all genes
        key_added="rank_genes_groups"
    )

    # extract + filter per archetype
    results = {}
    groups = adata_sub.uns["rank_genes_groups"]["names"].dtype.names

    for group in groups:
        df = sc.get.rank_genes_groups_df(
            adata_sub,
            group=group,
            key="rank_genes_groups"
        )

        df = df[
            (df["logfoldchanges"] > lfc_threshold) &
            (df["pvals_adj"] < pval_threshold)
        ].sort_values("scores", ascending=False)

        print(f"Archetype {group}: {len(df)} significant DEGs")
        results[int(group)] = df

    return results

def run_pairwise_deg_per_archetype(adata, lfc_threshold=1.0, pval_threshold=0.05):
    """
    Run pairwise DEG for each archetype pair (A vs B for all A != B).

    Parameters
    ----------
    adata : anndata.AnnData
        Must have "archetype" column in adata.obs with 0 = unassigned.
    lfc_threshold : float, default 1.0
        Minimum log fold change to keep.
    pval_threshold : float, default 0.05
        Maximum adjusted p-value to keep.

    Returns
    -------
    dict : {(i, j): pd.DataFrame} — genes upregulated in archetype i vs archetype j
    """

    adata_sub = adata[adata.obs["archetype"] != 0].copy()
    adata_sub.obs["archetype"] = adata_sub.obs["archetype"].astype(str)

    archetypes = sorted(adata_sub.obs["archetype"].unique())
    results = {}

    for ref, tgt in permutations(archetypes, 2):
        pair_adata = adata_sub[adata_sub.obs["archetype"].isin([ref, tgt])].copy()

        sc.tl.rank_genes_groups(
            pair_adata,
            groupby="archetype",
            groups=[ref],
            reference=tgt,
            method="wilcoxon",
            n_genes=pair_adata.n_vars,
            key_added="rank_genes_groups"
        )

        df = sc.get.rank_genes_groups_df(
            pair_adata,
            group=ref,
            key="rank_genes_groups"
        )

        df = df[
            (df["logfoldchanges"] > lfc_threshold) &
            (df["pvals_adj"] < pval_threshold)
        ].sort_values("scores", ascending=False)

        print(f"Archetype {ref} vs {tgt}: {len(df)} significant DEGs")
        results[(int(ref), int(tgt))] = df

    return results


def get_strict_archetype_genes(deg_dict, pairwise_deg_dict):
    """
    Get strict archetype marker genes by intersecting 1-vs-rest and all pairwise DEG results.

    A gene is a strict marker for archetype i if it is:
        - significant in 1-vs-rest DEG for archetype i
        - significant in archetype i vs every other archetype in pairwise DEG

    Parameters
    ----------
    deg_dict : dict
        Output of run_deg_per_archetype. {archetype (int): pd.DataFrame}
    pairwise_deg_dict : dict
        Output of run_pairwise_deg_per_archetype. {(i, j): pd.DataFrame}

    Returns
    -------
    pd.DataFrame
        Flattened dataframe with columns: archetype, gene, plus stats from 1-vs-rest.
    """

    archetypes = sorted(deg_dict.keys())
    strict_results = []

    for arch in archetypes:
        # genes from 1 vs rest
        one_vs_rest_genes = set(deg_dict[arch]["names"].values
                                if "names" in deg_dict[arch].columns
                                else deg_dict[arch].index)

        # intersect with genes significant in arch vs every other archetype
        other_archetypes = [a for a in archetypes if a != arch]
        pairwise_gene_sets = []

        for other in other_archetypes:
            key = (arch, other)
            if key in pairwise_deg_dict and len(pairwise_deg_dict[key]) > 0:
                genes = set(pairwise_deg_dict[key]["names"].values
                            if "names" in pairwise_deg_dict[key].columns
                            else pairwise_deg_dict[key].index)
                pairwise_gene_sets.append(genes)
            else:
                pairwise_gene_sets.append(set())

        # strict = in 1vrest AND in all pairwise comparisons
        if pairwise_gene_sets:
            strict_genes = one_vs_rest_genes.intersection(*pairwise_gene_sets)
        else:
            strict_genes = one_vs_rest_genes

        print(f"Archetype {arch}: {len(strict_genes)} strict marker genes")

        # pull stats from 1-vs-rest df
        df = deg_dict[arch]
        name_col = "names" if "names" in df.columns else df.index.name
        if name_col and name_col in df.columns:
            subset = df[df["names"].isin(strict_genes)].copy()
        else:
            subset = df[df.index.isin(strict_genes)].copy()
            subset = subset.reset_index()

        subset["archetype"] = arch
        strict_results.append(subset)

    return pd.concat(strict_results, ignore_index=True).sort_values(
        ["archetype", "scores"], ascending=[True, False]
    )

def run_go_analysis(deg_dict, adata, organism="mouse", n_top_genes=200):
    """
    Run GO enrichment analysis on DEG results per archetype using gseapy.

    Parameters
    ----------
    deg_dict : dict
        Output of run_deg_per_archetype. {archetype (int): pd.DataFrame}
    adata : anndata.AnnData
        Used to extract background gene list (all genes in the dataset).
    organism : str, default "mouse"
        Organism for GO analysis. "mouse" for mouse, "human" for human.
    n_top_genes : int, default 200
        Number of top genes (by score) to use per archetype for enrichment.

    Returns
    -------
    dict : {archetype (int): pd.DataFrame of GO enrichment results}
    """

    gene_sets  = [
        "GO_Biological_Process_2023",
        "GO_Molecular_Function_2023",
        "GO_Cellular_Component_2023"
    ]

    background = adata.var_names.tolist()
    print(f"Background: {len(background)} genes")

    results = {}

    for arch, df in deg_dict.items():
        print(f"\nRunning GO for archetype {arch}...")

        name_col = "names" if "names" in df.columns else None
        genes = (df.sort_values("scores", ascending=False)
                   .head(n_top_genes)[name_col]
                   .tolist()) if name_col else (
                 df.sort_values("scores", ascending=False)
                   .head(n_top_genes).index.tolist())

        if len(genes) == 0:
            print(f"  No genes for archetype {arch}, skipping.")
            results[arch] = pd.DataFrame()
            continue

        print(f"  {len(genes)} genes submitted")

        try:
            enr = gp.enrichr(
                gene_list       = genes,
                gene_sets       = gene_sets,
                organism        = organism,
                background      = background,
                outdir          = None,
                cutoff          = 0.05
            )

            res_df  = enr.results
            adj_col = [c for c in res_df.columns if "adjusted" in c.lower()][0]
            res_df  = res_df[res_df[adj_col] < 0.05].sort_values(adj_col)
            print(f"  {len(res_df)} significant GO terms")
            results[arch] = res_df

        except Exception as e:
            print(f"  Failed for archetype {arch}: {e}")
            results[arch] = pd.DataFrame()

    return results

def run_strict_go_analysis(strict_genes_df, adata, organism="mouse"):
    """
    Run GO enrichment on strict archetype marker genes.

    Parameters
    ----------
    strict_genes_df : pd.DataFrame
        Output of get_strict_archetype_genes.
        Must have "archetype" column and a gene name column ("names" or index).
    adata : anndata.AnnData
        Used to extract background gene list (all genes in dataset).
    organism : str, default "mouse"

    Returns
    -------
    dict : {archetype (int): pd.DataFrame of GO enrichment results}
    """
    gene_sets  = [
        "GO_Biological_Process_2023",
        "GO_Molecular_Function_2023",
        "GO_Cellular_Component_2023"
    ]

    background = adata.var_names.tolist()
    archetypes = sorted(strict_genes_df["archetype"].unique())
    results    = {}

    print(f"Background: {len(background)} genes\n")

    for arch in archetypes:
        print(f"Archetype {arch}...")

        arch_df  = strict_genes_df[strict_genes_df["archetype"] == arch]
        name_col = "names" if "names" in arch_df.columns else arch_df.columns[0]
        genes    = arch_df[name_col].tolist()

        print(f"  {len(genes)} strict marker genes")

        if len(genes) == 0:
            print(f"  No genes — skipping.")
            results[arch] = pd.DataFrame()
            continue

        try:
            enr = gp.enrichr(
                gene_list  = genes,
                gene_sets  = gene_sets,
                organism   = organism,
                background = background,
                outdir     = None,
                cutoff     = 0.05
            )

            res_df  = enr.results
            adj_col = [c for c in res_df.columns if "adjusted" in c.lower()][0]
            res_df  = res_df[res_df[adj_col] < 0.05].sort_values(adj_col)
            print(f"  {len(res_df)} significant GO terms")
            results[arch] = res_df

        except Exception as e:
            print(f"  Failed: {e}")
            results[arch] = pd.DataFrame()

    return results

