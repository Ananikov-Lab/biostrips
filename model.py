import sys
sys.path.append('.')

import pandas as pd
import numpy as np
import pickle as pkl

from rdkit import Chem
from mordred import Calculator, descriptors
from rdkit.Chem import Descriptors

from src import models
from src import experimental_setup

def get_ld50(smiles_list, mesure_variable):
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
    
    fn = 'rf_mordred' + str(4) + '_' + sampling_type
    mordred_rf = _benchmarks['rf_mordred']['model']()
    mordred_rf.load_weights('src/%s.chkpt' % fn)
    
    with open('./src/mordred_cols.pkl', 'rb') as f:
        m_new = pkl.load(f)
    
    calc_2d = Calculator(descriptors, ignore_3D=False)
    calc_3d = Calculator(descriptors, ignore_3D=True)
    
    if '.' in smiles_list[0]:
        max_wt = 0
        heavy_cmpd = ''
    for cmpd in smiles_list[0].split('.'):
        wt = Chem.Descriptors.MolWt(Chem.MolFromSmiles(cmpd))
        if max_wt < wt:
            max_wt = wt
            heavy_cmpd = cmpd
            
    smiles_list = [heavy_cmpd]

    mols = [Chem.MolFromSmiles(Chem.MolToSmiles(Chem.MolFromSmiles(smi))) for smi in smiles_list]
    if mols[0] is None:
        return 'unknown'

    df = calc_2d.pandas(mols)
    df_3d = calc_3d.pandas(mols)

    with open('./src/scaler.pkl', 'rb') as f:
        y_train = pkl.load(f)

    y_train = experimental_setup.scaler.fit_transform(y_train)

    y_hat = experimental_setup.scaler.inverse_transform(mordred_rf.predict(df[m_new].to_numpy()))
    if mesure_variable == 'mol/kg or mol/L':
        return str(10**(-1*round(y_hat[0][0], 2)))
    else:
        molwt = Descriptors.MolWt(Chem.MolFromSmiles(smiles_list[0]))
    return  str((10**(-1*round(y_hat[0][0], 2)))*molwt)
