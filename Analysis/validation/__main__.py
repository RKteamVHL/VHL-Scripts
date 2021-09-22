import pandas as pd

from .core import *

if __name__ == "__main__":
    masterlist_ref_df = pd.read_csv("all_refs.csv")
    masterlist_ref_df = masterlist_ref_df.rename(columns={MASTERLIST_PMID_COL: "MASTER_PMID"})

    # umd_df = get_umd_variants()
    umd_df = pd.read_csv("umd_df.csv", dtype={'UMD_PMID': str})
    # umd_df.to_csv("umd_df.csv")

    vhldb_df = get_vhldb_df("vhldbMuts_germline.csv")

    vhldb_umd_pmid_comp_df = compare_pmids([vhldb_df, umd_df, masterlist_ref_df],
                                           ["VHLDB_PMID", "UMD_PMID", "MASTER_PMID"])

    vhldb_umd_pmid_comp_df.to_csv("vhldb_umd_pmid_summary.csv")
    vhldb_variant_comp_df = compare_vhldb_variants("vhldbMuts_germline.csv", "all_variants.csv")
    vhldb_variant_comp_df.to_csv("vhldb_variant_summary.csv")
