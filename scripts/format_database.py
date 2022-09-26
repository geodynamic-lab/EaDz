#!/usr/bin/env python3
#
# =======================================================================
#
#               Python Scripts for EaDz web-based database
#                  -----------------------------------
#
#                                Authors:
#                 Bo Zhang, Shaofeng Liu, Chenxi Zhang
#             China university of Geoscience(Beijing) 2022
#                          ALL RIGHTS RESERVED
#
#
# =====================================================================
#
#
#  Last update: 24th September 2022 by Bo Zhang


# import necessary Libs

import os
import pathlib
from pathlib import Path
import time
import math
import glob
from collections import defaultdict
import sys
import shutil

import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import scipy
from scipy.signal import savgol_filter
import detritalpy
import detritalpy.detritalFuncs as dFunc
import pygplates


def get_paleo_position(row, model, static):

    try:
        anchor_plate = 601
        #
        points = [(row['longitude'], row['latitude'])]
        point_features = []
        for lon, lat in points:
            point_feature = pygplates.Feature()
            point_feature.set_geometry(pygplates.PointOnSphere(lat, lon))
            point_features.append(point_feature)

        point_features_cc = pygplates.partition_into_plates(static, model,
                                                            point_features,
                                                            properties_to_copy=[
                                                                pygplates.PartitionProperty.reconstruction_plate_id])
        try:
            reconstruction_time = row['age(ma)']
        except:
            reconstruction_time = row['tda']

        # for anchor plate 601
        reconstructed_feature_geometries = []
        pygplates.reconstruct(point_features_cc, model, reconstructed_feature_geometries, reconstruction_time,
                              anchor_plate)

        paleo_lat, paleo_lon = reconstructed_feature_geometries[0].get_reconstructed_geometry().to_lat_lon()
        row['paleo_lat'] = paleo_lat
        row['paleo_lon'] = paleo_lon

        # for anchor plate 0
        reconstructed_feature_geometries = []
        pygplates.reconstruct(point_features_cc, model, reconstructed_feature_geometries, reconstruction_time)

        paleo_lat, paleo_lon = reconstructed_feature_geometries[0].get_reconstructed_geometry().to_lat_lon()
        row['paleo_lat_mantle'] = paleo_lat
        row['paleo_lon_mantle'] = paleo_lon
    except:
        pass
    print(row[['paleo_lat', 'paleo_lat']])

    return row


def get_tda_of_detrital_zircon_samples(data):
    # convert detrital mda to mda
    data.loc[~data['YC1S WM'].isna(), 'mda_tag'] = 'Yc1sWm'
    data.loc[~data['YC1S WM'].isna(), 'mda'] = data.loc[~data['YC1S WM'].isna(), 'YC1S WM']
    data.loc[~data['YC1S WM'].isna(), 'mda_1se'] = data.loc[~data['YC1S WM'].isna(), 'YC1S WM 1serr']
    data.loc[~data['YC1S WM'].isna(), 'mda_MSWD'] = data.loc[~data['YC1S WM'].isna(), 'YC1S WM MSWD']

    data.loc[data['mda'].isna(), 'mda_tag'] = 'Yc2sWm'
    data.loc[data['mda'].isna(), 'mda'] = data.loc[data['mda'].isna(), 'YC2S WM']
    data.loc[data['mda'].isna(), 'mda_1se'] = data.loc[data['mda'].isna(), 'YC2S WM 1serr']
    data.loc[data['mda'].isna(), 'mda_MSWD'] = data.loc[data['mda'].isna(), 'YC2S WM MSWD']

    data.loc[data['mda'].isna(), 'mda_tag'] = 'YSG'
    data.loc[data['mda'].isna(), 'mda'] = data.loc[data['mda'].isna(), 'YSG']
    data.loc[data['mda'].isna(), 'mda_1se'] = data.loc[data['mda'].isna(), 'YSG_err1s']
    data.loc[data['mda'].isna(), 'mda_MSWD'] = data.loc[data['mda'].isna(), 'YC1S WM']

    data.loc[data['mda'].isna(), 'mda_tag'] = 'NotGet'

    # determine the true depositional age (tda)
    # mda represent the tda
    data.loc[(data['mda'] < data['rgt_max']) &
             (data['mda'] > data['rgt_min']), 'tda_tag'] = 'mda'

    data.loc[(data['mda'] < data['rgt_max']) &
             (data['mda'] > data['rgt_min']), 'tda'] = \
        data.loc[(data['mda'] < data['rgt_max']) &
                 (data['mda'] > data['rgt_min']), 'mda']

    data.loc[data['tda'].isna() &
             ~data['mda'].isna() &
             (data['mda'] < data['rgt_mean'])
    , 'tda_tag'] = 'mda'

    data.loc[data['tda'].isna() &
             ~data['mda'].isna() &
             (data['mda'] < data['rgt_mean'])
    , 'tda'] = \
        data.loc[data['tda'].isna() &
                 ~data['mda'].isna() &
                 (data['mda'] < data['rgt_mean'])
        , 'mda']

    # rgt_mean represent the tda
    data.loc[data['tda'].isna() &
             ~data['rgt_mean'].isna(), 'tda_tag'] = 'rgt_mean'

    data.loc[data['tda'].isna() &
             ~data['rgt_mean'].isna(), 'tda'] = \
        data.loc[data['tda'].isna() &
                 ~data['rgt_mean'].isna(), 'rgt_mean']

    # mad represent the tda, but no record limit
    data.loc[data['tda'].isna() &
             ~data['mda'].isna(), 'tda_tag'] = 'mda_no_constrain'

    data.loc[data['tda'].isna() &
             ~data['mda'].isna(), 'tda'] = \
        data.loc[data['tda'].isna() &
                 ~data['mda'].isna(), 'mda']

    # add the 1se error and MSWD for tda which based on MDA
    data.loc[data['tda_tag'] == 'mda', 'tda_1se'] = \
        data.loc[data['tda_tag'] == 'mda', 'mda_1se']

    data.loc[data['tda_tag'] == 'mda_no_constrain', 'tda_MSWD'] = \
        data.loc[data['tda_tag'] == 'mda_no_constrain', 'mda_MSWD']

    data.loc[data['tda_tag'] == 'mda_no_constrain', 'tda_1se'] = \
        data.loc[data['tda_tag'] == 'mda_no_constrain', 'mda_1se']

    data.loc[data['tda_tag'] == 'mda_no_constrain', 'tda_MSWD'] = \
        data.loc[data['tda_tag'] == 'mda_no_constrain', 'mda_MSWD']

    data = data.round({'tda': 2, 'mda': 2, 'tda_1se': 2, 'mda_1se': 2, 'tda_MSWD': 4, 'mda_MSWD': 4})

    return data


def convert_2se_into_1se(trim_data):
    for col in trim_data.columns:
        if '2se' in col:
            if col[:-3] + '1se' in trim_data.columns:
                trim_data.loc[trim_data[col[:-3] + '1se'].isna(), col[:-3] + '1se'] = trim_data[col] / 2
            else:
                trim_data[col[:-3] + '1se'] = trim_data[col] / 2

    drop_list = list()
    for col in trim_data.columns:
        if '2se' in col:
            drop_list.append(col)
    trim_data.drop(columns=drop_list, inplace=True)

    return trim_data


def cal_the_best_age(trim_data, threshold=1500):
    # the default threshold is 1500 according to Spencer et al., 2016, 
    # if age > 1500 use the age_207Pb/206Pb, else use the age_206Pb/238U
    # note the age_207Pb/206pb the pb is small case

    # process the age based on threshold

    trim_data.loc[trim_data['age_206Pb/238U'] <= threshold, 'cal_best_age'] = \
        trim_data.loc[trim_data['age_206Pb/238U'] <= threshold, 'age_206Pb/238U']

    trim_data.loc[trim_data['age_206Pb/238U'] <= threshold, 'cal_best_age_1se'] = \
        trim_data.loc[trim_data['age_206Pb/238U'] <= threshold, 'age_206Pb/238U_1se']

    trim_data.loc[trim_data['age_206Pb/238U'] > threshold, 'cal_best_age'] = \
        trim_data.loc[trim_data['age_206Pb/238U'] > threshold, 'age_207Pb/206Pb']

    trim_data.loc[trim_data['age_206Pb/238U'] > threshold, 'cal_best_age_1se'] = \
        trim_data.loc[trim_data['age_206Pb/238U'] > threshold, 'age_207Pb/206Pb_1se']

    return trim_data


def convert_tdm_to_hf_tdm(trim_data):
    col_list = ['TDM1', 'TDM1_1se', 'TDM2', 'TDM2_1se']
    for col in col_list:
        if col in trim_data.columns:
            if 'Hf_' + col in trim_data.columns:
                trim_data.loc[trim_data['Hf_' + col].isna(), 'Hf_' + col] = \
                    trim_data.loc[trim_data['Hf_' + col].isna(), col]
            else:
                trim_data['Hf_' + col] = trim_data[col]
            trim_data.drop(columns=col, inplace=True)
    return trim_data


def format_hf_tdm(trim_data):
    col_list = ['Hf_TDM1', 'Hf_TDM1_1se', 'Hf_TDM2', 'Hf_TDM2_1se']
    for col in col_list:
        if col in trim_data.columns:
            if "se" in col:
                threshold = 1
            else:
                threshold = 5

            trim_data.loc[trim_data[col] <= threshold, col] = trim_data.loc[trim_data[col] <= threshold, col] * 1000

    return (trim_data)


def cal_lu_hf_isotopic(trim_data):
    # define the constant
    lam = 0.00001867  # 176Lu decay rate (per myr), ref: 1.867 (Soderlund et al., 2004);   1.865 (Scherer et al., 2001)
    c_hf_chur_0 = 0.282772  # (present day) 176Hf/177Hf CHUR, ref: 0.282772 (BlichertToft, J & Albarede, F,
    # 1997; Griffin et al., 2000),   0.282793 (Iizuka et al., 2015)
    c_hf_dm = 0.28325  # 176Hf/177Hf DM, ref:0.283250 (Amelin et al., 1999; BlichertToft, J & Albarede, F,
    # 1997; Griffin et al., 2000)
    c_lu_chur_0 = 0.0332  # (present day) 176Lu/177Hf CHUR, ref: 0.332 (BlichertToft, J & Albarede, F, 1997; Griffin
    # et al., 2000),  0.0338 (Iizuka et al., 2015)
    c_lu_dm = 0.0384  # 176Lu/177Hf DM
    flu_hf = 0.0150  # 176Lu/177Hf (crustal avg. for fcc), ref: 0.0150 (Amelin et al. 1999)
    fcc = -0.5482  # fcc    (for TDM2 calc), ref: Duan, C et al. (2018). Minerals 8: 547.
    # https://doi.org/10.3390/min8120547
    fdm = 0.1566  # Duan, C et al. (2018). Minerals 8: 547. https://doi.org/10.3390/min8120547

    # calculation    
    trim_data['cal_176Hf/177Hf(t)'] = trim_data['176Hf/177Hf'] - trim_data['176Lu/177Hf'] * (
                np.exp(lam * trim_data['cal_best_age']) - 1)
    trim_data['cal_Hf(t)_Denominator'] = c_hf_chur_0 - c_lu_chur_0 * (np.exp(lam * trim_data['cal_best_age']) - 1)
    trim_data['cal_epsilon_Hf(t)'] = ((trim_data['cal_176Hf/177Hf(t)'] / trim_data[
        'cal_Hf(t)_Denominator']) - 1) * 10000
    trim_data['cal_lu/hf'] = (trim_data['176Lu/177Hf'] / c_lu_chur_0) - 1
    trim_data['cal_Hf_TDM1(Ma)'] = (1 / lam) * np.log(
        1 + ((trim_data['176Hf/177Hf'] - c_hf_dm) / (trim_data['176Lu/177Hf'] - c_lu_dm)))
    trim_data['cal_Hf_TDM2(Ma)'] = trim_data['cal_Hf_TDM1(Ma)'] - (
                trim_data['cal_Hf_TDM1(Ma)'] - trim_data['cal_best_age']) * (
                                               (fcc - trim_data['cal_lu/hf']) / (fcc - fdm))
    trim_data['disc_with_pub_of_epsilon_hf(t)'] = (trim_data['cal_epsilon_Hf(t)'] / trim_data[
        'epsilon_Hf(t)'] - 1) * 100
    return trim_data


def convert_u_th_to_th_u(trim_data):
    if 'U/Th' in trim_data.columns:
        trim_data.loc[trim_data['Th/U'].isna(), 'Th/U'] = trim_data.loc[trim_data['Th/U'].isna(), 'U/Th']
        trim_data.drop(columns=['U/Th'], inplace=True)
    if 'U/Th_1se' in trim_data.columns:
        trim_data.loc[trim_data['Th/U_1se'].isna(), 'Th/U_1se'] = trim_data.loc[trim_data['Th/U'].isna(), 'U/Th_1se']
        trim_data.drop(columns=['U/Th_1se'], inplace=True)
    return trim_data


def format_error_for_isotopic(trim_data,
                              col_list=['176Hf/177Hf_1se', '176Lu/177Hf_1se', '176Yb/177Hf_1se', '176Hf/177Hf(t)_1se'],
                              threshold=1, ratio=1000000):
    for col in col_list:
        if col in trim_data.columns:
            trim_data.loc[trim_data[col] > threshold, col] = trim_data.loc[trim_data[col] > threshold, col] / ratio
    return trim_data


def format_error_for_u_pb_isotopic(trim_data, threshold=1):
    col_list = ['ele_206Pb/238U_1se', 'ele_207Pb/206Pb_1se',
                'ele_207Pb/235U_1se', 'ele_208Pb/232Th_1se']
    for col in col_list:
        if col in trim_data.columns:
            trim_data.loc[trim_data[col] >= threshold, col] = trim_data.loc[trim_data[col] >= threshold, col] * \
                                                              trim_data.loc[trim_data[col] >= threshold, col[:-4]] / 100
    return trim_data


def fill_miss_spot(trim_data):
    """
    fill the spot holes, call the pd.series.values 
    """
    samples = list(set(trim_data.loc[trim_data['spot'].isna(), 'sample']))
    for sample in samples:
        spot_num = len(trim_data.loc[trim_data['sample'] == sample, 'spot'])

        trim_data.loc[trim_data['sample'] == sample, 'spot'] = 'spot_' + pd.Series(np.arange(spot_num)).astype(
            str).values
    return trim_data


def generate_streamlit_sample_spots(trim_data):
    trim_data['st_sample'] = trim_data['sample'].copy()
    trim_data['st_spot'] = trim_data['spot'].copy()

    for col in ['st_sample', 'st_spot']:
        trim_data.loc[trim_data[col].astype(str).str.isdecimal(), col] = 'sam-' + \
                                                                         trim_data.loc[trim_data[col].astype(
                                                                             str).str.isdecimal(), col].astype(str)

    return trim_data


def format_data_set(trim_data):
    trim_data = convert_2se_into_1se(trim_data)
    trim_data = cal_the_best_age(trim_data, threshold=1500)
    trim_data = convert_tdm_to_hf_tdm(trim_data)
    trim_data = format_hf_tdm(trim_data)
    trim_data = cal_lu_hf_isotopic(trim_data)
    trim_data = convert_u_th_to_th_u(trim_data)
    trim_data = format_error_for_isotopic(trim_data)
    trim_data = format_error_for_u_pb_isotopic(trim_data, threshold=1)
    trim_data = fill_miss_spot(trim_data)
    trim_data = generate_streamlit_sample_spots(trim_data)

    return trim_data


# parameters for gplates
work_dir = Path('<replace with your data directory>')
sample_file = '<replace with your sample information file>'
data_file = '<replace with you data file>'
sample_inf = pd.read_excel(os.path.join(work_dir, sample_file), sheet_name='Sheet1')
zircon_data = pd.read_excel(os.path.join(work_dir, data_file), sheet_name='Sheet1')

# calculate the paleoposition of samples
gplates_data_dir = Path('<replace with your data directory>')

# here we use the model presented in Young et al., 2019
rot_file = os.path.join(gplates_data_dir, 'FeatureCollections/Rotations',
                        'Muller2019-Young2019-Cao2020_CombinedRotations.rot')
static_polygon_file = os.path.join(gplates_data_dir, 'FeatureCollections/StaticPolygons',
                                   'Global_EarthByte_GPlates_PresentDay_StaticPlatePolygons.gpmlz'
                                   )
rotation_model_M19 = pygplates.RotationModel(rot_file)
sample_inf.apply(get_paleo_position, model=rotation_model_M19, static=static_polygon_file, axis=1)

# format data
zircon_data = format_data_set(zircon_data)
zircon_data = cal_the_best_age(zircon_data)

# calculate the mda based on detrialpy
detripy_data = zircon_data[['unique_sample_id', 'spot', 'cal_best_age', 'cal_best_age_1se']]
rename_dict = {'unique_sample_id': 'Sample_ID',
               'spot': 'Analysis_ID',
               'cal_best_age': 'BestAge',
               'cal_best_age_1se': 'BestAge_err'
               }
detripy_data.rename(columns=rename_dict, inplace=True)
samples = pd.DataFrame(detripy_data['Sample_ID'].unique(), columns=['Sample_ID'])
with pd.ExcelWriter('detritalpy.xlsx') as writer:
    samples.to_excel(writer, sheet_name='Samples', index=None)
    detripy_data.to_excel(writer, sheet_name='ZrUPb', index=None)

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
print('detritalPy version: ', detritalpy.__version__)

dataToLoad = ['detritalpy.xlsx']
main_df, main_byid_df, samples_df, analyses_df = dFunc.loadDataExcel(dataToLoad)
sampleList = main_df.Sample_ID.to_list()
ages, errors, numGrains, labels = dFunc.sampleToData(sampleList, main_byid_df, sigma='1sigma')
fileName = 'EAGC_MDA.csv'

makePlot = False
# Specify how grains are sorted: by best age, best age + 1 sigma error, or best age + 2 sigma error
sortBy = '1sigma'

# Specify plot parameters
plotWidth = 8
plotHeight = 5
barWidth = 0.25
ageColors = ['blue', 'red', 'green']
fillMDACalcs = True
alpha = 0.25

figMDA = dFunc.MDAtoCSV(sampleList,
                        ages, errors, numGrains, labels, fileName,
                        sortBy, barWidth, plotWidth, plotHeight,
                        ageColors, alpha, makePlot, fillMDACalcs)
