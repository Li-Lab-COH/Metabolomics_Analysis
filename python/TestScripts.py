import re
import numpy as np
import pandas as pd
from scipy.stats import hypergeom, chi2

# ---------- helpers ----------
def parse_kegg_maps(cell):
    """Return a list of KEGG pathway IDs normalized to the 5-digit core (e.g., '04211').
    Accepts 'ko04211', 'map04211', 'hsa04211' and multiple separated by , ; | space."""
    if pd.isna(cell):
        return []
    tokens = re.split(r'[,\s;|]+', str(cell))
    out = []
    for t in tokens:
        t = t.strip()
        if not t:
            continue
        m = re.match(r'(?:ko|map|hsa)?(\d{5})$', t)
        if m:
            out.append(m.group(1))
    # unique + sorted for stability
    return sorted(set(out))

def add_parsed_maps(df, kegg_col='kegg_maps', out_col='kegg_ids'):
    return df.assign(**{out_col: df[kegg_col].map(parse_kegg_maps)})

def benjamini_hochberg(pvals):
    p = np.asarray(pvals, float)
    n = p.size
    order = np.argsort(p)
    ranks = np.empty(n, int); ranks[order] = np.arange(1, n+1)
    adj = p * n / ranks
    # monotone: from largest to smallest, then reverse back
    adj_sorted = np.minimum.accumulate(adj[order[::-1]])[::-1]
    out = np.empty(n, float); out[order] = np.clip(adj_sorted, 0, 1)
    return out

# ---------- core ORA ----------
def enrich_kegg_maps(
    df,
    feature_col,                 # unique feature id per metabolite (e.g., 'Index' or metabolite name)
    kegg_col='kegg_maps',
    pval_col=None,               # 'P.Value' if no adj.P given
    adjp_col=None,               # 'adj.P.Val' preferred if present
    logfc_col=None,              # optional for direction filters
    alpha=0.05,
    direction='any',             # 'any' | 'up' | 'down'
    min_set_size=3,
    max_set_size=5000
):
    tmp = add_parsed_maps(df, kegg_col=kegg_col, out_col='kegg_ids')

    # background = metabolites with at least one KEGG pathway mapping
    universe = tmp[tmp['kegg_ids'].str.len() > 0].copy()
    if universe.empty:
        raise ValueError(f"No features with KEGG pathway IDs found in column '{kegg_col}'.")

    # choose significance
    if adjp_col and adjp_col in tmp.columns:
        tmp['adjP'] = tmp[adjp_col].astype(float)
    elif pval_col and pval_col in tmp.columns:
        tmp['adjP'] = benjamini_hochberg(tmp[pval_col].astype(float).values)
    else:
        raise ValueError("Provide either 'adjp_col' or 'pval_col' present in df.")

    # significant set
    sig = universe[universe['adjP'] <= alpha].copy()
    if logfc_col and direction in ('up', 'down'):
        sig = sig[sig[logfc_col] > 0] if direction == 'up' else sig[sig[logfc_col] < 0]

    # explode to long form
    uni_long = universe[[feature_col, 'kegg_ids']].explode('kegg_ids')
    sig_long = sig[[feature_col, 'kegg_ids']].explode('kegg_ids')

    # pathway membership
    members_by_pw = uni_long.groupby('kegg_ids')[feature_col].apply(lambda s: set(s)).to_dict()
    sig_by_pw     = sig_long.groupby('kegg_ids')[feature_col].apply(lambda s: set(s)).to_dict()

    # sizes
    M = universe[feature_col].nunique()
    n = sig[feature_col].nunique()

    rows = []
    for pw, members in members_by_pw.items():
        K = len(members)
        if K < min_set_size or K > max_set_size:
            continue
        k = len(sig_by_pw.get(pw, set()) & members)
        pval = hypergeom.sf(k - 1, M, K, n) if k > 0 else 1.0
        rows.append({
            'pathway_id': pw,                # 5-digit KEGG code, e.g., '04211'
            'M_universe': M,
            'K_pathway': K,
            'n_sig': n,
            'k_overlap': k,
            'overlap_ratio': (k / n) if n else np.nan,   # fraction of sig landing in this pathway
            'pathway_coverage': (k / K) if K else np.nan,# fraction of pathway hit by your sig
            'p_overrep': pval
        })

    res = pd.DataFrame(rows)
    if res.empty:
        # nothing enriched; still return a table with universe/sig counts for transparency
        return pd.DataFrame({'pathway_id': [], 'M_universe':[M], 'n_sig':[n]})[:0]

    res['fdr'] = benjamini_hochberg(res['p_overrep'].values)
    res = res.sort_values(['fdr', 'p_overrep', 'pathway_id']).reset_index(drop=True)
    return res

# ---------- (later) combine metabolite + proteomic pathway p-values ----------
def fishers_method(pvals):
    pvals = [p for p in pvals if (p is not None and np.isfinite(p) and p > 0)]
    if not pvals:
        return np.nan
    stat = -2 * np.sum(np.log(pvals))
    df = 2 * len(pvals)
    return chi2.sf(stat, df)

def combine_pathway_pvalues(meta_res, prot_res,
                            meta_pcol='p_overrep', prot_pcol='p_overrep'):
    # both inputs must have 'pathway_id' using the same 5-digit normalization
    merged = pd.merge(
        meta_res[['pathway_id', meta_pcol]],
        prot_res[['pathway_id', prot_pcol]],
        on='pathway_id', how='outer'
    )
    merged[meta_pcol] = merged[meta_pcol].fillna(1.0)
    merged[prot_pcol] = merged[prot_pcol].fillna(1.0)
    merged['p_fisher'] = merged.apply(
        lambda r: fishers_method([r[meta_pcol], r[prot_pcol]]), axis=1
    )
    merged['fdr_fisher'] = benjamini_hochberg(merged['p_fisher'].values)
    return merged.sort_values('fdr_fisher')
