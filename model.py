import sys
sys.path.append('.')

import pandas as pd
import numpy as np
import pickle as pkl

from rdkit import Chem
from mordred import Calculator, descriptors

from src import models
from src import experimental_setup

def get_ld50(smiles_list):
    if smiles_list == '':
        return 'unknown'
    experimental_setup.path_prefix = '/../'

    _benchmarks = {
        'rf_mordred': {'model': models.RF, 'encoding': 'mordred'},
    }
    
    # benchmarks to train/validate, check _benchmark_dict for options
    run_benchmarks = ['rf_mordred']
    
    # `random` or `stratified`
    sampling_type = 'random'
    
    kfold = experimental_setup.CrossValidator(
        splits = 5, # dont change without re-running data preprocessing
        sampling_type = sampling_type,
    )
    
    converter = experimental_setup.LD50UnitConverter()
    
    fn = 'rf_mordred' + str(0) + '_' + sampling_type
    mordred_rf = _benchmarks['rf_mordred']['model']()
    mordred_rf.load_weights('../../data/benchmark-models/chkpts/%s.chkpt' % fn)
    
    with open('./src/mordred_cols.pkl', 'rb') as f:
        m_new = pkl.load(f)
    
    calc_2d = Calculator(descriptors, ignore_3D=False)
    calc_3d = Calculator(descriptors, ignore_3D=True)
    
    mols = [Chem.MolFromSmiles(smi) for smi in smiles_list]
    
    df = calc_2d.pandas(mols)
    df_3d = calc_3d.pandas(mols)
    y_hat = experimental_setup.scaler.inverse_transform(mordred_rf.predict(df[m_new].to_numpy()))
    return y_hat[0][0]