import os
import urllib
import requests

SEARCH_API_URL = "https://search.rcsb.org/rcsbsearch/v2/query"


def protein_query(gene):
    gene = gene.lower()
    return {
        "query": {
                "type": "group",
                "logical_operator": "and",
                "nodes": [
                    {
                        "type": "terminal",
                        "service": "text",
                        "parameters": {
                            "attribute": "struct.title",
                            "operator": "contains_phrase",
                            "value": gene
                        }
                    },
                    {
                        "type": "terminal",
                        "service": "text",
                        "parameters": {
                            "attribute": "struct.title",
                            "operator": "contains_phrase",
                            "value": "apo"
                        }
                    }
                ]

        },
        "return_type": "entry",
        "request_options": {
            "results_content_type": [
                "experimental"
            ],
            "sort": [
                {
                    "sort_by": "em_3d_reconstruction.resolution",
                    "direction": "asc"
                },
                {
                    "sort_by": "score",
                    "direction": "desc"
                }
            ],
            "scoring_strategy": "combined"
        }
    }


def search_protein(search_term):
    query_data = protein_query(search_term)
    response = requests.post(SEARCH_API_URL, json=query_data)
    results = response.json()
    if "total_count" in results and results["total_count"] > 0:
        return results["result_set"][0]["identifier"]
    else:
        return None


def write_pdb(pdb_id, path=None):
    structure_url = 'https://files.rcsb.org/download/'
    file_name = pdb_id + '.pdb'

    if path:
        os.makedirs(path, exist_ok=True)
        output_path = os.path.join(path, file_name)
    else:
        output_path = file_name

    urllib.request.urlretrieve(structure_url + pdb_id + '.pdb', output_path)
    return output_path
