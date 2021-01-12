"""
The function `fetch_from_api` searches the PubMed database for publications about clinical trials related to a substance and a disease. It returns a DataFrame of all found publications and additional information about them (publication year, publication type(s), NCT id's of clinical trials associated with the publication).
"""

import pandas as pd
from lxml import etree
from Bio import Entrez
from urllib.error import HTTPError


def fetch_from_api(df, retmax=1000, mindate='2000/01', maxdate='2020/08', email='vanboefer@gmail.com'):
    """
    Take a DataFrame with the columns 'active_substance' and 'disease_name', and search the PubMed database for publications about clinical trials related to the substance and the disease.
    Return a DataFrame of all found publications and additional information about them (publication year, publication type(s), NCT id's of clinical trials associated with the publication).
    """

    Entrez.email = email

    # get pmids
    pmids = df.apply(get_pmids, axis=1, retmax=retmax, mindate=mindate, maxdate=maxdate)

    # create pmid_df
    pmid_df = pd.DataFrame(pmids.explode().dropna(), columns=['pmid'])
    pmid_df.index.name = 'idx'

    # get additional info
    info_pmids = pmid_df.pmid.apply(get_pub_info, retmax=retmax).apply(pd.Series, index=['info_pub_year', 'info_ncts', 'info_pub_types'])
    pmid_df = pmid_df.assign(
        info_pub_year = info_pmids.info_pub_year,
        info_ncts = info_pmids.info_ncts,
        info_pub_types = info_pmids.info_pub_types,
)

    return pmid_df

def compose_query(drug, disease):
    return f"""
        {drug}[SUBS]
        AND {disease}[MESH]
        AND English[LANG]
        AND (
            Clinical Study[PTYP] OR
            Clinical Trial[PTYP] OR
            Randomized Controlled Trial[PTYP] OR
            Controlled Clinical Trial[PTYP] OR
            Clinical Trial, Phase I[PTYP] OR
            Clinical Trial, Phase II[PTYP] OR
            Clinical Trial, Phase III[PTYP] OR
            Clinical Trial, Phase IV[PTYP]
    )
    """

def get_pmids(row, **kwargs):
    """
    Take a row from a DataFrame with the columns 'active_substance' and 'disease_name' and return a list of publications about clinical trials related to the substance and the disease.
    """
    # compose query
    query = compose_query(row.active_substance, row.disease_name)

    # get pmid's
    searchhandle = Entrez.esearch(db="pubmed", term=query, usehistory="y", **kwargs)
    searchresults = Entrez.read(searchhandle)
    searchhandle.close()
    return [str(item) for item in searchresults['IdList']]

def get_pub_info(pmid, **kwargs):
    """
    For one PMID, fetch info from Entrez.

    Returns
    -------
    pub_year: {None, str}
        Publication year
    nct: {None, set}
        NCT id's of clinical trials associated with the publication
    pubtype: {None, set}
        Publication types
    """

    pub_year = None
    nct = None
    pubtype = None

    try:
        fetchhandle = Entrez.efetch(db="pubmed", id=str(pmid), retmode="xml", **kwargs)
        data = fetchhandle.read()
        fetchhandle.close()
        root = etree.fromstring(data)
    except (HTTPError):
        print(f"PMID {pmid} could not be fetched: bad request")
        return pub_year, nct, pubtype

    # retrieve nct's
    databank_list = root.find('PubmedArticle/MedlineCitation/Article/DataBankList')
    if databank_list is not None:
        nct = set()
        for databank in databank_list:
            if databank.find('DataBankName').text == 'ClinicalTrials.gov':
                for item in databank.find('AccessionNumberList'):
                    if hasattr(item, 'text'):
                        nct.add(item.text)

    # retrieve publication year
    node_pub_year = root.find('PubmedArticle/MedlineCitation/Article/Journal/JournalIssue/PubDate/Year')
    if node_pub_year is not None:
        pub_year = node_pub_year.text

    # retrieve publication types
    pubtype_list = root.find('PubmedArticle/MedlineCitation/Article/PublicationTypeList')
    if pubtype_list is not None:
        pubtype = set()
        for item in pubtype_list:
            pubtype.add(item.text)

    return pub_year, nct, pubtype


if __name__ == '__main__':

    # load ema data
    ema = pd.read_pickle('../data/ema_cuis.pkl')

    # fetch api data and pickle df
    df = fetch_from_api(ema)
    df.to_pickle('../data/results/pubmed_api.pkl')
