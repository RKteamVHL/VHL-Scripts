import pandas as pd
import requests
from lxml import html

cdna_xpath = '/html/body/div/table[1]/tr[2]/td[1]'
aa_xpath = '/html/body/div/table[1]/tr[2]/td[2]'

pubmed_xpath = '//td/a'

base_link = 'http://www.umd.be/VHL/4DACTION'

MASTERLIST_PMID_COL = "PMID"
MASTERLIST_HGVS_COL = "HGVS_transcript"
VHLDB_PMID_COL = "PubMed ID"
VHLDB_VARIANT_COL = "Variant"
UMD_PMID_COL = "PMID"
UMD_VARIANT_COL = "Mutation Event c.DNA."


def get_vhldb_df(filename):
    df = pd.read_csv(filename, delimiter="\t")
    vhldb_df = df.copy()
    vhldb_df.loc[:, VHLDB_PMID_COL] = vhldb_df[VHLDB_PMID_COL].str.replace('[\"\[\]]', '', regex=True).str.split(",")
    vhldb_df = vhldb_df.explode(column=VHLDB_PMID_COL)
    vhldb_df = vhldb_df.rename(columns={VHLDB_PMID_COL: "VHLDB_PMID"})
    return vhldb_df

def compare_vhldb_variants(vhldb_filename, masterlist_filename):
    vhldb_df = get_vhldb_df(vhldb_filename)
    masterlist_ref_df = pd.read_csv(masterlist_filename)

    vhldb_refs = set(vhldb_df[VHLDB_VARIANT_COL])
    masterlist_refs = set(masterlist_ref_df[MASTERLIST_HGVS_COL])

    union = vhldb_refs.union(masterlist_refs)
    summary_df = pd.DataFrame(columns=["In VHLdb", "In Masterlist"], index=union)

    summary_df.loc[:, "In VHLdb"] = [ref in vhldb_refs for ref in summary_df.index]
    summary_df.loc[:, "In Masterlist"] = [ref in masterlist_refs for ref in summary_df.index]


    return summary_df

def compare_vhldb_pmids(vhldb_df, masterlist_df):

    vhldb_refs = set(vhldb_df[VHLDB_PMID_COL])
    masterlist_refs = set(masterlist_df[MASTERLIST_PMID_COL])

    union = vhldb_refs.union(masterlist_refs)
    summary_df = pd.DataFrame(columns=["In VHLdb", "In Masterlist"], index=union)

    summary_df.loc[:, "In VHLdb"] = [ref in vhldb_refs for ref in summary_df.index]
    summary_df.loc[:, "In Masterlist"] = [ref in masterlist_refs for ref in summary_df.index]


    # intersect = vhldb_refs.intersection(masterlist_refs)
    # missing = vhldb_refs.difference(masterlist_refs)
    # extra = masterlist_refs.difference(vhldb_refs)

    return summary_df

def compare_pmids(df_list, df_cols):
    pmid_refs = []
    ref_union = set()
    for i in range(len(df_list)):
        refs = set(df_list[i][df_cols[i]].astype(str))
        pmid_refs.append(df_list[i][df_cols[i]])
        ref_union = ref_union.union(refs)

    summary_df = pd.DataFrame(columns=df_cols, index=ref_union)

    for i in range(len(df_list)):
        summary_df.loc[:, df_cols[i]] = [ref in set(pmid_refs[i]) for ref in summary_df.index]

    return summary_df



def get_umd_variants():

    page = requests.get('http://www.umd.be/VHL/4DACTION/W_DMDT1/1')
    tree = html.fromstring(page.content)
    all_links = tree.xpath("//a")
    variant_link_set = set()
    for link in all_links:
        href = link.attrib["href"]
        if href.startswith("../../4DACTION/WV/"):
            abslink = href.replace("../../4DACTION", base_link)
            variant_link_set.add(abslink)

    variants = []

    for variant_link in variant_link_set:
        variant = {}

        v_page = requests.get(variant_link)
        v_tree = html.fromstring(v_page.content)

        cdna_ele = v_tree.xpath(cdna_xpath)
        aa_ele = v_tree.xpath(aa_xpath)

        if cdna_ele:
            variant["Mutation Event c.DNA."] = cdna_ele[0].text

        if aa_ele:
            variant["Predicted Consequence Protein Change"] = aa_ele[0].text

        pubmed_ele = v_tree.xpath(pubmed_xpath)
        if pubmed_ele and pubmed_ele[0].text.isdecimal():
            variant["UMD_PMID"] = pubmed_ele[0].text


        if variant:
            variants.append(variant)

    variant_df = pd.DataFrame.from_dict(variants)

    return variant_df
