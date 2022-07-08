from __future__ import print_function

import os
import sys
import pdb
import argparse

import google.auth
from google.auth.transport.requests import Request
from google.oauth2.credentials import Credentials
from google_auth_oauthlib.flow import InstalledAppFlow
from googleapiclient.discovery import build
from googleapiclient.errors import HttpError
import numpy as np
from sklearn import cluster 

# If modifying these scopes, delete the file token.json.
SCOPES = ['https://www.googleapis.com/auth/spreadsheets.readonly']

SPREADSHEET_ID = '1VRl8dghA5t1JiLEa47OWW-hr-Qbvn9yksB4Hrga19Ck'
GLUCOSE_10_COLUMN = "D"
COMPOUND_START_COLUMN = 176
COMPOUND_END_COLUMN = 276
NANOTEMPER_START_COLUMN = "J"
NANOTEMPER_END_COLUMN = "M"

GLUCOSE_RANGE = GLUCOSE_10_COLUMN + str(COMPOUND_START_COLUMN) + ":" + GLUCOSE_10_COLUMN + str(COMPOUND_END_COLUMN)
NANOTEMPER_RANGE = NANOTEMPER_START_COLUMN + str(COMPOUND_START_COLUMN) + ":" + NANOTEMPER_END_COLUMN + str(COMPOUND_END_COLUMN)

def get_creds():
    creds = None
    # The file token.json stores the user's access and refresh tokens, and is
    # created automatically when the authorization flow completes for the first
    # time.
    if os.path.exists('token.json'):
        creds = Credentials.from_authorized_user_file('token.json', SCOPES)
    # If there are no (valid) credentials available, let the user log in.
    if not creds or not creds.valid:
        if creds and creds.expired and creds.refresh_token:
            creds.refresh(Request())
        else:
            flow = InstalledAppFlow.from_client_secrets_file(
                'credentials.json', SCOPES)
            creds = flow.run_local_server(port=0)
        # Save the credentials for the next run
        with open('token.json', 'w') as token:
            token.write(creds.to_json())
    return creds

def get_values(spreadsheet_id, range_names):
    """
    Creates the batch_update the user has access to.
    Load pre-authorized user credentials from the environment.
    TODO(developer) - See https://developers.google.com/identity
    for guides on implementing OAuth2 for the application.\n"
        """
    # creds, _ = google.auth.default()
    # pylint: disable=maybe-no-member
    try:
        service = build('sheets', 'v4', credentials=get_creds())
        if type(range_names) is str:
            result = service.spreadsheets().values().get(
                spreadsheetId=spreadsheet_id, range=range_names, valueRenderOption='UNFORMATTED_VALUE').execute()
            rows = result.get('values', [])
            print(f"{len(rows)} rows retrieved")
        elif type(range_names) is list:
            result = service.spreadsheets().values().batchGet(
                spreadsheetId=spreadsheet_id, ranges=range_names, valueRenderOption='UNFORMATTED_VALUE').execute()
        else:
            raise Error("'range_names' should be of type str or list")
        return result
    except HttpError as error:
        print(f"An error occurred: {error}")
        return error


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Run K-means.")
    parser.add_argument('-K', type=int)
    args = parser.parse_args()

    # Pass: spreadsheet_id, and range_name
    tables = get_values(SPREADSHEET_ID, [GLUCOSE_RANGE, NANOTEMPER_RANGE])
    if 'range' in tables.keys():
#        X = np.array(tables['range'])
        X = pd.DataFrame(tables['range'])
        print(X)
        exit()
    elif 'valueRanges' in tables.keys():
        valueRanges = tables['valueRanges']
        arrays = []
        for valueRange in valueRanges:
            values = valueRange['values']
            array = np.array(values)
            arrays.append(array)
        X = np.concatenate(arrays, axis=1) 
    else:
        raise Error("The key 'range' or 'valueRanges' isn't found")

    X[X == '#N/A (No matches are found in FILTER evaluation.)'] = np.nan
    X = X.astype('float')
    print("Shape without Nan removal ", X.shape)
    X = X[~np.isnan(X).any(axis=1)]
    print("Shape after NaN removal ", X.shape)
    kmeans = cluster.KMeans(n_clusters=args.K, random_state=0).fit(X)
    print("kmeans.labels_  \n", kmeans.labels_)
    print("kmeans.cluster_centers_  \n", kmeans.cluster_centers_)

