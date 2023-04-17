import requests
from pathlib import Path

from util.util import clean_directory


def __protein_query(gene):
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
    query_data = __protein_query(search_term)
    response = requests.post("https://search.rcsb.org/rcsbsearch/v2/query", json=query_data)
    results = response.json()
    if "total_count" in results and results["total_count"] > 0:
        return results["result_set"][0]["identifier"]
    else:
        return None


def write_pdb(pdb_id, path='../pdbs'):
    path = Path(path)
    output_path = path/(pdb_id + '.pdb')

    with open(output_path, 'wb') as o:
        o.write(requests.get('https://files.rcsb.org/download/' + pdb_id + '.pdb').content)
    return output_path
