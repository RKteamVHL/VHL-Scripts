import requests
import pandas as pd
from lxml import html

cdna_xpath = '/html/body/div/table[1]/tr[2]/td[1]'
aa_xpath = '/html/body/div/table[1]/tr[2]/td[2]'

pubmed_xpath = '//td/a'

base_link = 'http://www.umd.be/VHL/4DACTION'


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
            variant["PMID"] = pubmed_ele[0].text


        if variant:
            variants.append(variant)

    variant_df = pd.DataFrame.from_dict(variants)

    return variant_df
