from __future__ import print_function

import argparse

import pandas as pd
from sklearn import cluster
import libscol

COMPOUND_START_ROW, COMPOUND_END_ROW = 176, 276
GLUCOSE_10_COLUMN = "D"
NANOTEMPER_START_COLUMN, NANOTEMPER_END_COLUMN = "J", "M"

GLUCOSE_RANGE = (
    GLUCOSE_10_COLUMN
    + str(COMPOUND_START_ROW)
    + ":"
    + GLUCOSE_10_COLUMN
    + str(COMPOUND_END_ROW)
)
NANOTEMPER_RANGE = (
    NANOTEMPER_START_COLUMN
    + str(COMPOUND_START_ROW)
    + ":"
    + NANOTEMPER_END_COLUMN
    + str(COMPOUND_END_ROW)
)

RELEVANT_ROWS = slice(COMPOUND_START_ROW - 1, COMPOUND_END_ROW)
RELEVANT_COLS = [
    libscol.scol2int(GLUCOSE_10_COLUMN),
    *range(
        libscol.scol2int(NANOTEMPER_START_COLUMN),
        libscol.scol2int(NANOTEMPER_END_COLUMN) + 1,
    ),
]


def get_args():
    parser = argparse.ArgumentParser(description="Run K-means.")
    parser.add_argument("-K", nargs="*", type=int)
    parser.add_argument("-U", "--update-spreadsheet")
    parser.add_argument("--random_state", type=int, default=0)
    return parser.parse_args()


def get_kmeans_results(data, K, random_state=0):
    X = data.dropna()
    kmeans = (
        cluster.KMeans(random_state=random_state)
        if K == None
        else cluster.KMeans(n_clusters=K, random_state=random_state)
    )
    kmeans_result = kmeans.fit(X)
    centers = pd.DataFrame(kmeans.cluster_centers_.round(2))
    centers.columns = ["Glu", "Na.4.5", "Na.4.25", "Na.5.5", "Na.5.25"]
    centers.index.name = "Center label"
    labeled_data = pd.concat(
        [
            data,
            pd.Series(
                kmeans.labels_, index=X.index, dtype="Int64", name="kmeans cluster"
            ),
        ],
        axis=1,
    )
    labeled_data.index.name = "Spreadsheet row"
    return centers, labeled_data


def get_results_for_spreadsheet(spreadsheet, centers, labeled_data, K, random_state):
    last_used_column = libscol.int2scol(len(spreadsheet.columns))
    labels_range = (
        last_used_column
        + str(COMPOUND_START_ROW)
        + ":"
        + last_used_column
        + str(COMPOUND_END_ROW)
    )
    range_names = [
        last_used_column + "1",
        last_used_column + "2",
        last_used_column + "3",
        labels_range,
    ]
    values = [
        [["K = " + str(K)]],
        [["Random state = " + str(random_state)]],
        [[centers.to_string()]],
        [list(labeled_data["kmeans cluster"].astype(str))],
    ]
    return list(map(spreadsheet.create_value_range, range_names, values))


def display_to_terminal(data, X, centers, results):
    print("Shape without NaN removal ", data.shape)
    print("Shape after NaN removal ", X.shape)
    print("\n", centers)
    print("\n", results)


if __name__ == "__main__":
    args = get_args()
    for K in args.K:
        spreadsheet = spreadsheet.get_spreadsheet()
        relevant = spreadsheet.iloc[RELEVANT_ROWS, RELEVANT_COLS]
        centers, labeled_data = get_kmeans_results(relevant, K)
        results_for_spreadsheet = get_results_for_spreadsheet(
            spreadsheet, centers, labeled_data, K, args.random_state
        )
        if args.update_spreadsheet:
            spreadsheet.batch_update_values(
                spreadsheet.SPREADSHEET_ID, "USER_ENTERED", results_for_spreadsheet
            )
        else:
            display_to_terminal(
                relevant, labeled_data, centers, results_for_spreadsheet
            )
