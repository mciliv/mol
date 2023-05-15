import os

from google.auth.transport.requests import Request
from google.oauth2.credentials import Credentials
from google_auth_oauthlib.flow import InstalledAppFlow
from googleapiclient.discovery import build
from googleapiclient.errors import HttpError
import pandas as pd
import numpy as np

# If modifying these scopes, delete token.json
SCOPES = ['https://www.googleapis.com/auth/spreadsheets']
CREDENTIALS = './credentials.json'
SPREADSHEET_ID = '1VRl8dghA5t1JiLEa47OWW-hr-Qbvn9yksB4Hrga19Ck'


def get_creds():
    creds = None
    # The file token.json stores the user's access and refresh tokens, and is
    # created automatically when the authorization flow completes for the first
    # time.
    if os.path.exists('../token.json'):
        creds = Credentials.from_authorized_user_file('../token.json', SCOPES)
    # If there are no (valid) credentials available, let the user log in.
    if not creds or not creds.valid:
        if creds and creds.expired and creds.refresh_token:
            creds.refresh(Request())
        else:
            flow = InstalledAppFlow.from_client_secrets_file(CREDENTIALS, SCOPES)
            creds = flow.run_local_server(port=0)
        # Save the credentials for the next run
        with open('../token.json', 'w') as token:
            token.write(creds.to_json())
    return creds


def get_values(spreadsheet_id, range_names):
    # creds, _ = google.auth.default()
    # pylint: disable=maybe-no-member
    try:
        service = build('sheets', 'v4', credentials=get_creds())
        restResource = service.spreadsheets().values()
        if type(range_names) is str:
            result = service.spreadsheets().values().get(
                spreadsheetId=spreadsheet_id, range=range_names, majorDimension='COLUMNS',
                valueRenderOption='UNFORMATTED_VALUE').execute()
        elif type(range_names) is list:
            result = service.spreadsheets().values().batchGet(
                spreadsheetId=spreadsheet_id, ranges=range_names, majorDimension='COLUMNS',
                valueRenderOption='UNFORMATTED_VALUE').execute()
        else:
            raise Error("'range_names' should be of type str or list")
        return result
    except HttpError as error:
        print(f"An error occurred: {error}")
        return error


def get_spreadsheet() -> pd.DataFrame:
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
        raise Exception("The key 'range' or 'valueRanges' isn't found")
    spreadsheet[spreadsheet == '#N/A (No matches are found in FILTER evaluation.)'] = np.nan
    return spreadsheet


def create_value_range(range, values):
    return {
        'range': range,
        'majorDimension': 'COLUMNS',
        'values': values
    }


def update_values(spreadsheet_id, range_names, value_input_option,
                  values):
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


def batch_update_values(spreadsheet_id,
                        value_input_option, data):
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
