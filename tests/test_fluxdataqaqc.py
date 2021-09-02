# -*- coding: utf-8 -*-

import pkg_resources
import pytest
from pathlib import Path
from shutil import move, copy, rmtree

import numpy as np
import pandas as pd

from fluxdataqaqc import Data
from fluxdataqaqc import QaQc
from fluxdataqaqc import Plot
from fluxdataqaqc import Convert
from fluxdataqaqc import util

@pytest.fixture(scope="session")
def data(request):
    """Prepare input data used for tests"""
    d = {}
    package_root_dir = Path(
        pkg_resources.resource_filename(
            'fluxdataqaqc', '.'
        )
    ).parent

    d['package_root_dir'] = package_root_dir      
    # use example data from tutorial for most tests
    examples_dir = package_root_dir / 'examples'
    test_data_files = [
        f for f in examples_dir.rglob('*') if f.suffix in [
            '.csv','.ini','.xlsx'] and not 'gridMET_data' in str(f)
    ]
    
    temp_data_dir = Path('tests') / 'test_data'
    if not temp_data_dir.is_dir():
        temp_data_dir.mkdir(parents=True, exist_ok=True)

    print(temp_data_dir, Path())

    for f in test_data_files:
        copy(f, temp_data_dir)




    def teardown():
        rmtree(temp_data_dir)

    request.addfinalizer(teardown)

    return d
    
    
class TestData(object):

#    def setup_method(self, data):
#        self.root_dir = data['package_root_dir']
        
    def test_ACSE_refET(self, data):
        config = data['package_root_dir']/'examples'/'Basic_usage'/'US-Tw3_config.ini'
        data_obj = Data(config)
        ts = data_obj.hourly_ASCE_refET()
        assert len(ts) == 47544
        
