PubMed search
======================
The script [`api_query.py`](code/api_query.py) searches the PubMed database for publications about clinical trials related to a substance and a disease. It returns a DataFrame of all found publications and additional information about them: publication year, publication type(s), NCT id's of clinical trials associated with the publication.

The [data](data/) folder contains a sample of 72 EMA drug authorizations whose substance and disease is used in the query.

For more details on the experiment and the results, see [here](https://github.com/vanboefer/link_ct_results).
