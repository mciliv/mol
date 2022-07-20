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
import pandas as pd
from sklearn import cluster 
import libscol

import config

# If modifying these scopes, delete the file token.json.
SCOPES = ['https://www.googleapis.com/auth/spreadsheets.readonly']

SPREADSHEET_ID = '1VRl8dghA5t1JiLEa47OWW-hr-Qbvn9yksB4Hrga19Ck'
GLUCOSE_10_COLUMN = "D"
COMPOUND_START_ROW = 176
COMPOUND_END_ROW = 276
NANOTEMPER_START_COLUMN = "J"
NANOTEMPER_END_COLUMN = "M"

GLUCOSE_RANGE = GLUCOSE_10_COLUMN + str(COMPOUND_START_ROW) + ":" + GLUCOSE_10_COLUMN + str(COMPOUND_END_ROW)
NANOTEMPER_RANGE = NANOTEMPER_START_COLUMN + str(COMPOUND_START_ROW) + ":" + NANOTEMPER_END_COLUMN + str(COMPOUND_END_ROW)

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
                config.CREDENTIALS, SCOPES)
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
        restResource = service.spreadsheets().values()
        if type(range_names) is str:
            result = service.spreadsheets().values().get(
                spreadsheetId=spreadsheet_id, range=range_names, majorDimension='COLUMNS', valueRenderOption='UNFORMATTED_VALUE').execute()
        elif type(range_names) is list: result = service.spreadsheets().values().batchGet(
                spreadsheetId=spreadsheet_id, ranges=range_names, majorDimension='COLUMNS', valueRenderOption='UNFORMATTED_VALUE').execute()
        else:
            raise Error("'range_names' should be of type str or list")
        return result
    except HttpError as error:
        print(f"An error occurred: {error}")
        return error

def update_values(spreadsheet_id, range_names, value_input_option,
                  values):
    """
    Creates the batch_update the user has access to.
    Load pre-authorized user credentials from the environment.
    TODO(developer) - See https://developers.google.com/identity
    for guides on implementing OAuth2 for the application.\n"
        """
    # pylint: disable=maybe-no-member
    try:
        service = build('sheets', 'v4', credentials=get_creds())
        body = {
            'values': values
        }
        result = service.spreadsheets().values().update(
            spreadsheetId=spreadsheet_id, range=range_name,
            valueInputOption=value_input_option, body=body).execute()
        print(f"{(result.get('updates').get('updatedCells'))} cells updated.")
        return result
    except HttpError as error:
        print(f"An error occurred: {error}")
        return error

def batch_update_values(spreadsheet_id, range_names,
                        value_input_option, data):
    """
        Creates the batch_update the user has access to.
        Load pre-authorized user credentials from the environment.
        TODO(developer) - See https://developers.google.com/identity
        for guides on implementing OAuth2 for the application.\n"
            """
    # pylint: disable=maybe-no-member
    try:
        service = build('sheets', 'v4', credentials=get_creds())
        body = {
            'valueInputOption': value_input_option,
            'data': data
        }
        result = service.spreadsheets().values().batchUpdate(
            spreadsheetId=spreadsheet_id, body=body).execute()
        print(f"{(result.get('totalUpdatedCells'))} cells updated.")
        return result
    except HttpError as error:
        print(f"An error occurred: {error}")
        return error

def set_pandas_to_max_view():
    pd.set_option('display.max_rows', None)
    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', None)  # more options can be specified also

def spreadsheet_to_dataframe():
    set_pandas_to_max_view()
    
    # Pass: spreadsheet_id, and range_name
    tables = get_values(SPREADSHEET_ID, 'Sheet1') 
    if 'values' in tables.keys():
        spreadsheet = pd.concat(map(pd.Series, tables['values']), axis=1)
    elif 'valueRanges' in tables.keys():
        valueRanges = tables['valueRanges']
        ranges = []
        for valueRange in valueRanges:
            values = valueRange['values']
            range = pd.concat(map(pd.Series, values), axis=1)
            ranges.append(range)
        spreadsheet = pd.concatenate(map(pd.Series, ranges), axis=1)
    else:
        raise Error("The key 'range' or 'valueRanges' isn't found")
    spreadsheet[spreadsheet == '#N/A (No matches are found in FILTER evaluation.)'] = np.nan
    return spreadsheet

def display_terminal(relevant, X, centers, results):
    print("Shape without NaN removal ", relevant.shape)
    print("Shape after NaN removal ", X.shape)
    print("\n", centers)
    print("\n", results)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Run K-means.")
    parser.add_argument('-K', type=int)
    args = parser.parse_args()
    rows = slice(COMPOUND_START_ROW - 1, COMPOUND_END_ROW)
    cols = [libscol.scol2int(GLUCOSE_10_COLUMN)] + [*range(libscol.scol2int(NANOTEMPER_START_COLUMN), libscol.scol2int(NANOTEMPER_END_COLUMN) + 1)]
    spreadsheet = spreadsheet_to_dataframe()
    relevant = spreadsheet.iloc[rows, cols]
    X = relevant.dropna()
    kmeans = cluster.KMeans(random_state=0) if args.K == None else cluster.KMeans(n_clusters=args.K, random_state=0)
    kmeans_result = kmeans.fit(X)
    centers = pd.DataFrame(kmeans.cluster_centers_.round(2))
    centers.index.name = 'Center label'
    labeled_data = pd.concat([relevant, pd.Series(kmeans.labels_, index=X.index, dtype='Int64', name='kmeans cluster')], axis=1)
    labeled_data.index.name = "Spreadsheet row"
    last_used_column = libscol.int2scol(len(spreadsheet.columns))
    range_name = last_used_column + str(COMPOUND_START_ROW) + ':' + last_used_column + str(COMPOUND_END_ROW)
    if input(f'Do you want to write to the range: {range_name} and its title range? yes?') == 'yes':
        range_names = [
        batch_update_values(SPREADSHEET_ID,
                      range_names, "USER_ENTERED",
                      labeled_data['kmeans cluster'])

    
