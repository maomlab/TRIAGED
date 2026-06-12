import pandas as pd
import json
import os
import argparse
from datetime import datetime


##########################
# Feature Block Definitions
##########################

def get_feature_blocks(df):
    """
    Maps block names to column lists based on column prefixes/patterns.
    """
    cols = df.columns.tolist()

    blocks = {
        'tgcpl_covalent_ifp_res':       [c for c in cols if c.startswith('tgcpl_covalent') and len(c.split('_')) > 3],
        'hscpl_covalent_ifp_res':       [c for c in cols if c.startswith('hscpl_covalent') and len(c.split('_')) > 3],
        'tgcpl_noncovalent_ifp_res':    [c for c in cols if c.startswith('tgcpl_noncovalent') and len(c.split('_')) > 3],
        'hscpl_noncovalent_ifp_res':    [c for c in cols if c.startswith('hscpl_noncovalent') and len(c.split('_')) > 3],

        'tgcpl_covalent_ifp_counts':    [c for c in cols if c.startswith('tgcpl_covalent') and len(c.split('_')) == 3],
        'hscpl_covalent_ifp_counts':    [c for c in cols if c.startswith('hscpl_covalent') and len(c.split('_')) == 3],
        'tgcpl_noncovalent_ifp_counts': [c for c in cols if c.startswith('tgcpl_noncovalent') and len(c.split('_')) == 3],
        'hscpl_noncovalent_ifp_counts': [c for c in cols if c.startswith('hscpl_noncovalent') and len(c.split('_')) == 3],

        'tgcpl_covalent_preds':         [c for c in cols if c.startswith('tgcpl_covalent') and any(x in c for x in ['pred', 'binding'])],
        'hscpl_covalent_preds':         [c for c in cols if c.startswith('hscpl_covalent') and any(x in c for x in ['pred', 'binding'])],
        'tgcpl_noncovalent_preds':      [c for c in cols if c.startswith('tgcpl_noncovalent') and any(x in c for x in ['pred', 'binding'])],
        'hscpl_noncovalent_preds':      [c for c in cols if c.startswith('hscpl_noncovalent') and any(x in c for x in ['pred', 'binding'])],

        'tgcpl_covalent_conf':          [c for c in cols if c.startswith('tgcpl_covalent') and any(x in c for x in ['confidence', 'ptm', 'iptm', 'plddt', 'pde'])],
        'hscpl_covalent_conf':          [c for c in cols if c.startswith('hscpl_covalent') and any(x in c for x in ['confidence', 'ptm', 'iptm', 'plddt', 'pde'])],
        'tgcpl_noncovalent_conf':       [c for c in cols if c.startswith('tgcpl_noncovalent') and any(x in c for x in ['confidence', 'ptm', 'iptm', 'plddt', 'pde'])],
        'hscpl_noncovalent_conf':       [c for c in cols if c.startswith('hscpl_noncovalent') and any(x in c for x in ['confidence', 'ptm', 'iptm', 'plddt', 'pde'])],

        'metadata':                     [c for c in cols if c in [
                                            'molecular_weight', 'log_p', 'log_d', 'log_s',
                                            'num_aromatic_rings', 'num_h_bond_donors',
                                            'num_h_bond_acceptors', 'num_rule_of_5_violations',
                                            'pka', 'pka_basic', 'heavy_atom_count', 'tpsa',
                                            'num_rotatable_bonds', 'cns_mpo_score', 'bbb2_score',
                                            'fsp3', 'pka_acidic']],
        'rdkit_metadata':               [c for c in cols if c in [
                                            'molwt', 'numheavy', 'numrotbonds',
                                            'numrings', 'hydrophobicatoms', 'hbondacceptors']],
        'rgroup':                       [c for c in cols if c.startswith('R1_')
                                            or c.startswith('R3a_')
                                            or c.startswith('R3b_')],
    }

    for name, block_cols in blocks.items():
        if not block_cols:
            print(f"Warning: block '{name}' has no columns — check prefixes")

    return blocks


##########################
# Apply Block Filters
##########################

def apply_block_filters(blocks, block_filters):
    """
    Apply include/exclude filters to block columns from JSON BLOCK_FILTERS.
    
    Supports per block:
      - include_prefixes : only keep columns starting with these prefixes
      - include_columns  : whitelist specific columns
      - exclude_columns  : blacklist specific columns
    """
    filtered = {}
    for block_name, cols in blocks.items():
        if block_name not in block_filters:
            filtered[block_name] = cols
            continue

        f = block_filters[block_name]

        # apply include_prefixes
        if 'include_prefixes' in f:
            prefixes = f['include_prefixes']
            cols = [c for c in cols if any(c.startswith(p) for p in prefixes)]

        # apply include_columns (whitelist)
        if 'include_columns' in f:
            whitelist = set(f['include_columns'])
            cols = [c for c in cols if c in whitelist]

        # apply exclude_columns (blacklist)
        if 'exclude_columns' in f:
            blacklist = set(f['exclude_columns'])
            cols = [c for c in cols if c not in blacklist]

        if not cols:
            print(f"Warning: block '{block_name}' has no columns after filtering")

        filtered[block_name] = cols

    return filtered


##########################
# Column Selection
##########################

def select_columns(df, blocks, requested):
    selected_cols = []
    resolved = {}
    unknown = []

    for item in requested:
        if isinstance(item, str):
            if item in blocks:
                resolved[item] = blocks[item]
                selected_cols.extend(blocks[item])
            elif item in df.columns:
                resolved[item] = [item]
                selected_cols.append(item)
            else:
                unknown.append(item)

        elif isinstance(item, dict):
            block_name     = item.get('block')
            cols_requested = item.get('columns', [])

            if block_name not in blocks:
                unknown.append(block_name)
                continue

            block_cols = blocks[block_name]
            valid      = [c for c in cols_requested if c in block_cols]
            invalid    = [c for c in cols_requested if c not in block_cols]

            if invalid:
                print(f"Warning: these columns not found in block '{block_name}': {invalid}")

            key = f"{block_name}[subset]"
            resolved[key] = valid
            selected_cols.extend(valid)

    if unknown:
        print(f"Warning: unknown blocks or columns: {unknown}")

    # deduplicate preserving order
    seen   = set()
    deduped = []
    for c in selected_cols:
        if c not in seen:
            seen.add(c)
            deduped.append(c)

    return deduped, resolved


##########################
# Main
##########################

def main(args):
    with open(args.json_file) as f:
        cfg = json.load(f)

    feature_matrix_path = os.path.expandvars(cfg['FEATURE_MATRIX'])
    output_base_dir     = os.path.expandvars(cfg['OUTPUT_DIR'])
    requested           = cfg['INCLUDE']
    target_cols         = cfg.get('TARGETS', [])
    block_filters       = cfg.get('BLOCK_FILTERS', {})

    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    run_dir   = os.path.join(output_base_dir, timestamp)
    os.makedirs(run_dir, exist_ok=True)

    df     = pd.read_csv(feature_matrix_path, index_col='substance_id')
    # --- compute selectivity and add to feature matrix ---
    tgcpl_col = 'mean_tgcpl_log_ic50 (uM)'
    hscpl_col = 'mean_hscpl_log_ic50 (uM)'

    if tgcpl_col in df.columns and hscpl_col in df.columns:
        # selectivity = log(HsCPL IC50 / TgCPL IC50) = HsCPL - TgCPL in log space
        # positive value = more selective for TgCPL (higher HsCPL IC50 = harder to inhibit human)
        df['selectivity'] = df[hscpl_col] - df[tgcpl_col]
        n_valid = df['selectivity'].notna().sum()
        print(f"Computed selectivity for {n_valid} compounds "
            f"(range: {df['selectivity'].min():.2f} to {df['selectivity'].max():.2f})")
    else:
        print(f"Warning: could not compute selectivity — "
            f"missing {tgcpl_col!r} or {hscpl_col!r}")

    blocks = get_feature_blocks(df)
    blocks = get_feature_blocks(df)

    # apply JSON block filters
    if block_filters:
        blocks = apply_block_filters(blocks, block_filters)
        print(f"Applied BLOCK_FILTERS for: {list(block_filters.keys())}")

    selected_cols, resolved = select_columns(df, blocks, requested)

    if not selected_cols:
        raise ValueError("No valid columns selected — check your INCLUDE list.")

    keep_cols = selected_cols.copy()
    for t in target_cols:
        if t in df.columns and t not in keep_cols:
            keep_cols.append(t)

    X = df[keep_cols]

    fname   = f"{timestamp}.csv"
    outpath = os.path.join(run_dir, fname)
    X.to_csv(outpath)

    print(f"\nRun: {timestamp}")
    print(f"Blocks/columns requested: {requested}")
    if block_filters:
        print(f"Block filters applied: {json.dumps(block_filters, indent=2)}")
    print(f"\nResolved columns per block:")
    for item, cols in resolved.items():
        print(f"  {item}: {len(cols)} columns")
    print(f"\nTotal features selected: {len(selected_cols)}")
    print(f"Total columns in output (inc. targets): {X.shape[1]}")
    print(f"Saved to: {outpath}")

    manifest = {
        'timestamp':    timestamp,
        'requested':    requested,
        'block_filters': block_filters,
        'targets':      target_cols,
        'resolved':     {k: v for k, v in resolved.items()},
        'output_file':  outpath,
        'n_features':   len(selected_cols),
        'n_samples':    X.shape[0],
    }
    manifest_path = os.path.join(run_dir, f"{timestamp}_manifest.json")
    with open(manifest_path, 'w') as f:
        json.dump(manifest, f, indent=2)
    print(f"Manifest saved to: {manifest_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Feature selector: retains specified blocks/columns from feature matrix."
    )
    parser.add_argument("--json_file", required=True, help="Path to JSON config.")
    args = parser.parse_args()
    main(args)