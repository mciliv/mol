import requests
from pathlib import Path

import util


def local_data_dir():
    return util.local_data_dir(__file__)


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
                        "value": gene,
                    },
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "struct.title",
                        "operator": "contains_phrase",
                        "value": "apo",
                    },
                },
            ],
        },
        "return_type": "entry",
        "request_options": {
            "results_content_type": ["experimental"],
            "sort": [
                {"sort_by": "em_3d_reconstruction.resolution", "direction": "asc"},
                {"sort_by": "score", "direction": "desc"},
            ],
            "scoring_strategy": "combined",
        },
    }


def search_protein(search_term):
    query_data = __protein_query(search_term)
    response = requests.post(
        "https://search.rcsb.org/rcsbsearch/v2/query", json=query_data
    )
    results = response.json()
    if "total_count" in results and results["total_count"] > 0:
        return results["result_set"][0]["identifier"]
    else:
        return None


def write_pdb(pdb_id, path=local_data_dir(), prepended_label=""):
    path = Path(path)
    if prepended_label is not "":
        prepended_label += "_"
    output_path = (path / (prepended_label + pdb_id)).with_suffix(".pdb")

    with open(output_path, "wb") as o:
        o.write(
            requests.get("https://files.rcsb.org/download/" + pdb_id + ".pdb").content
        )
    return output_path


def apoprotein(names, destination=local_data_dir() / "apoproteins"):
    util.mkdirs(destination, clean=True)
    for name in names:
        write_pdb(search_protein(name), path=destination, prepended_label=name + "_apo")
    return destination


if __name__ == "__main__":
    write_pdb("2nnq")
