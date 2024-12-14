'''

This module tests the app_script.py module. Since that module is using mainly the function 
generate_assembly_library from main.py (../../../src), those tests can be found in test_main.py
Tests you can find here: 
    - test_optimize_no_file: Test if there is no file uploaded
    - test_optimize_wrong_file_type: Test if the file is uploaded, but it is the wrong file type
    
'''

import pytest
from unittest.mock import patch
import io
import pandas as pd
import sys
import os

# Add the GUI directory to the Python path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../")))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../test_data')))

from app_script import app  # Now it should find app_script.py

#Configure the app for testing mocking a client to send requests
@pytest.fixture
def client():
    app.config['TESTING'] = True
    with app.test_client() as client:
        yield client

# Test if there is no file uploaded
def test_optimize_no_file(client):
    # Create the data dictionary
    data = {
        'enzyme_name': 'BbsI',
        'min_oligo_size': 20,
        'max_oligo_size': 100
    }
    # Send the POST request without file data
    response = client.post('/optimize', data=data)
    # Verify the response
    assert response.status_code == 400
    assert b"Error: Both gene file and mutations file are required" in response.data

# Test if the file is uploaded, but the wrong file type
def test_optimize_wrong_file_type(client):
    # Open the mock files
    # Create the data dictionary
    with open('test_data/test_seq_single.csv', 'rb') as gene_file, \
         open('test_data/mixed_base_mut_list.csv', 'rb') as mutations_file:
        data = {
            'enzyme_name': 'BbsI',
            'min_oligo_size': 20,
            'max_oligo_size': 100,
            'gene_file': (gene_file, 'mock_data.csv'),
            'mutations_file': (mutations_file, 'mock_mutation_data.csv'),
        }

        # Send the POST request with real file data
        response = client.post('/optimize', data=data, content_type='multipart/form-data')
    
    # Verify the response
    assert response.status_code == 500
    assert b"An error occurred" in response.data

