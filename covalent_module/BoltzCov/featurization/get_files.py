import json
import os

#####################
# JSON / Config
#####################

def read_json_args(json_file):
    """
    Reads the JSON input file. All keys are required.
    Always loads both covalent and noncovalent data — no mechanism flag needed.
    """
    with open(json_file, 'r') as jf:
        a = json.load(jf)

    def ep(key):
        """Path helper — expands env vars, raises if missing."""
        val = a.get(key)
        if val is None:
            raise ValueError(f"Missing required JSON key: '{key}'")
        return os.path.expandvars(val) if isinstance(val, str) else val

    def lp(key):
        """List helper — raises if missing or not a list."""
        val = a.get(key)
        if val is None:
            raise ValueError(f"Missing required JSON key: '{key}'")
        if not isinstance(val, list):
            raise ValueError(f"JSON key '{key}' must be a list, got: {type(val)}")
        return val

    return dict(
        record_path            = ep("RECORD_PATH"),
        metadata_feats         = lp("METADATA"),
        boltz_pred_feats       = lp("BOLTZ_PRED"),
        boltz_confidence_feats = lp("BOLTZ_CONF"),
        sdf_file               = ep("SDF_FILE"),
        avg_preds_tgcpl        = ep("TGCPL_PRED"),
        avg_preds_hscpl        = ep("HSCPL_PRED"),
        ifp_tgcpl_res          = ep("TGCPL_IFP_RES"),
        ifp_hscpl_res          = ep("HSCPL_IFP_RES"),
        ifp_tgcpl_types        = ep("TGCPL_IFP_TYPE"),
        ifp_hscpl_types        = ep("HSCPL_IFP_TYPE"),
        avg_preds_tgcpl_noncov = ep("TGCPL_PRED_NONCOV"),
        avg_preds_hscpl_noncov = ep("HSCPL_PRED_NONCOV"),
        ifp_tgcpl_res_noncov   = ep("TGCPL_IFP_RES_NONCOV"),
        ifp_hscpl_res_noncov   = ep("HSCPL_IFP_RES_NONCOV"),
        ifp_tgcpl_types_noncov = ep("TGCPL_IFP_TYPE_NONCOV"),
        ifp_hscpl_types_noncov = ep("HSCPL_IFP_TYPE_NONCOV"),
    )


def get_mechanism_files(cfg):
    """
    Always returns both covalent and noncovalent file paths.
    Validates all paths exist before returning.
    """
    files = {
        'covalent': {
            'tgcpl':          cfg['avg_preds_tgcpl'],
            'hscpl':          cfg['avg_preds_hscpl'],
            'tgcpl_ifp_res':  cfg['ifp_tgcpl_res'],
            'hscpl_ifp_res':  cfg['ifp_hscpl_res'],
            'tgcpl_ifp_type': cfg['ifp_tgcpl_types'],
            'hscpl_ifp_type': cfg['ifp_hscpl_types'],
        },
        'noncovalent': {
            'tgcpl':          cfg['avg_preds_tgcpl_noncov'],
            'hscpl':          cfg['avg_preds_hscpl_noncov'],
            'tgcpl_ifp_res':  cfg['ifp_tgcpl_res_noncov'],
            'hscpl_ifp_res':  cfg['ifp_hscpl_res_noncov'],
            'tgcpl_ifp_type': cfg['ifp_tgcpl_types_noncov'],
            'hscpl_ifp_type': cfg['ifp_hscpl_types_noncov'],
        },
    }

    for mech, paths in files.items():
        missing = [k for k, v in paths.items() if v is None]
        if missing:
            raise ValueError(f"[{mech}] Missing paths for: {missing}")
        for k, path in paths.items():
            if not os.path.exists(path):
                raise FileNotFoundError(f"[{mech}] File not found for '{k}': {path}")

    return files
