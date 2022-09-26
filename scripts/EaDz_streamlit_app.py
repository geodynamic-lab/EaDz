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
#  Last update: 25rd September 2022 by Bo Zhang


# import necessary Libs
import os
from pathlib import Path
import streamlit as st
import pandas as pd
import numpy as np
import base64
import matplotlib.pyplot as plt
import seaborn as sns
import detritalpy
import detritalpy.detritalFuncs as dFunc
import pathlib
import geopandas as gpd
import pygplates
import cartopy.crs as ccrs
from shapely.geometry.polygon import LinearRing, Point
import leafmap.foliumap as leafmap
from math import atan2, degrees
import math
from collections import defaultdict
import warnings
import mas_function
import socket
import timeit
import io

warnings.filterwarnings('ignore')


# function set
# load data
@st.cache(allow_output_mutation=True)
def load_data(sample, data):
    data_pd = pd.read_pickle(data)
    sample_pd = pd.read_pickle(sample)
    count_data = data_pd.groupby(['unique_sample_id'])['cal_best_age'].count()
    drop_list = list(count_data[count_data < 60].index)
    data_pd = data_pd.set_index('unique_sample_id')
    data_pd = data_pd.drop(drop_list, axis=0)
    data_pd = data_pd.reset_index()

    data_pd.dropna(subset=['cal_best_age'], axis=0, inplace=True)
    return sample_pd, data_pd


def uploaded_file_to_gdf(data):
    import tempfile
    import os
    import uuid

    _, file_extension = os.path.splitext(data.name)
    file_id = str(uuid.uuid4())
    file_path = os.path.join(tempfile.gettempdir(), f"{file_id}{file_extension}")

    with open(file_path, "wb") as file:
        file.write(data.getbuffer())

    if file_path.lower().endswith(".kml"):
        gpd.io.file.fiona.drvsupport.supported_drivers["KML"] = "rw"
        gdf = gpd.read_file(file_path, driver="KML")
    else:
        gdf = gpd.read_file(file_path)

    return file_path, gdf


def download_data_as_csv(data_inf, sample_inf, dl_col, dl_re_dict, file_name):
    @st.cache
    def prepare_download_file(data_df, download_cols, download_data_rename_dicts, sample_df):
        if 'custom_group_name' in data_df.columns:
            download_cols.append('custom_group_name')
        data_df = data_df[download_cols]
        data_df.rename(columns=download_data_rename_dicts, inplace=True)
        merge_sample = sample_df[['unique_sample_id', 'formation', 'group', 'reference', 'DOI']]

        download_data = pd.merge(download_data,
                                 merge_sample,
                                 left_on='unique_sample_id',
                                 right_on='unique_sample_id',
                                 how='left')

        download_data.sort_values(by='unique_sample_id', inplace=True)
        return download_data.to_csv(index=None).encode('utf-8')

    download_file = prepare_download_file(data_inf, dl_col, dl_re_dict, sample_inf)
    st.download_button(
        label="Download data as CSV",
        data=download_file,
        file_name=file_name,
        mime='text/csv',
    )


def download_fig_as_pdf(fig, fig_name, key='st'):
    fn = fig_name
    fig.savefig(fn)
    with open(fn, "rb") as img:
        btn = st.download_button(
            label="Download image as pdf",
            data=img,
            file_name=fn,
            mime="image/pdf",
            key=key
        )


def add_leafmap(data, lat_col, lon_col, center=[36, 120], zoom=4, user_geo=None):
    """
    add map based on the leafmap
    :param data: Pd.DataFrame
    :param lat_col: str, latitude column name
    :param lon_col: str, longitude column name
    :param center: list, [latitude, longitude] define the plot center
    :param zoom: int, define the map zoom lever, default: 4
    :return:
    """
    plot_data = data[['unique_sample_id', lat_col, lon_col, 'rock_type', 'formation', 'MDA(Ma)', 'reference', 'DOI']]
    plot_data.dropna(subset=[lon_col, lat_col], inplace=True)
    if center == 'cal':
        center = [plot_data[lat_col].mean(), plot_data[lon_col].mean()]
    map_sample = leafmap.Map(center=center, zoom=zoom, height=400, width=600,
                             draw_export=True,
                             basemap="HYBRID",
                             plugin_Draw=True,
                             locate_control=True,
                             plugin_LatLngPopup=False,
                             )
    map_sample.add_points_from_xy(plot_data,
                                  x=lon_col,
                                  y=lat_col,
                                  layer_name='detrital sample')
    # map_sample.add_basemap("Esri.NatGeoWorldMap")
    map_sample.add_basemap("HYBRID")
    map_sample.add_basemap("Esri.WorldStreetMap")

    # add china geology map
    map_sample.add_geojson('data/geo3al/major_basins.geojson', layer_name='Major basin')
    if user_geo:
        map_sample.add_geojson(user_geo, layer_name='User define')

    # map_sample.add_shp()
    st.session_state['map_num'] += 1
    return map_sample.to_streamlit()


def prepare_visualization_data(show_samples_inf, sample_list,
                               st_sel_data, st_sel_samples, data_suffix, key_prefix='samples'):
    """
    prepare data for visualization
    :param key_prefix:
    :param show_samples_inf: essential information of the sample
    :param sample_list: user selected samples
    :param st_sel_data: the data based on space and time selection criteria
    :param st_sel_samples: the sample data based on space and time selection criteria
    :param data_suffix: the name of the data subset
    :return: samples and associated data selected by user
    """
    vis_sample_data = show_samples_inf[show_samples_inf['unique_sample_id'].isin(sample_list)]
    vis_dp_data = st_sel_data[st_sel_data['unique_sample_id'].isin(sample_list)]

    st.dataframe(show_samples_inf[show_samples_inf['unique_sample_id'].isin(sample_list)])

    download_data_as_csv(vis_dp_data, vis_sample_data, download_col, download_data_rename_dict,
                         '{} Visualization Data From EaDz database {}.csv'.format(key_prefix, data_suffix))
    if st.button('show samples in map', key=key_prefix):
        add_leafmap(st_sel_samples[st_sel_samples['unique_sample_id'].isin(sample_list)],
                    lat_col, lon_col, center='cal', zoom=7)
    dp_sample = st_sel_samples[st_sel_samples['unique_sample_id'].isin(sample_list)]
    dp_data = st_sel_data[st_sel_data['unique_sample_id'].isin(sample_list)]
    dp_data = trim_data_for_detritalPy(dp_data)

    return dp_sample, dp_data


@st.cache
def trim_data_for_detritalPy(in_data):
    in_data.loc[in_data['cal_best_age_1se'].isna(), 'cal_best_age_1se'] = \
        in_data.loc[in_data['cal_best_age_1se'].isna(), 'cal_best_age'] * 0.01
    in_data.loc[in_data['cal_discordance'].isna(), 'cal_discordance'] = 5
    return in_data


@st.cache(allow_output_mutation=True)
def load_data_into_detritalpy(sel_sample, st_sel_data, sample_id='unique_sample_id'):
    if sample_id == 'unique_sample_id':
        samples = sel_sample[
            ['sample', sample_id, 'formation', 'group', 'baisn', 'System', 'longitude', 'latitude', 'reference']]
    else:
        samples = sel_sample[[sample_id, 'longitude', 'latitude']]

    sample_col_rename_dict = {sample_id: 'Sample_ID',
                              'formation': 'Unit',
                              'group': 'Group',
                              'baisn': 'Basin',
                              'System': 'Age',
                              'longitude': 'Longitude',
                              'latitude': 'Latitude'}
    samples.rename(columns=sample_col_rename_dict, inplace=True)
    # samples['Age'] = samples['Age'].replace(0, np.nan)

    zircon_data = st_sel_data[
        [sample_id, 'sample_spots_id', 'cal_best_age', 'cal_best_age_1se', 'cal_disc']]
    zircon_data_col_rename_dict = {sample_id: 'Sample_ID',
                                   'sample_spots_id': 'Analysis_ID',
                                   'cal_best_age': 'BestAge',
                                   'cal_best_age_1se': 'BestAge_err',
                                   'cal_disc': 'Disc',
                                   'Th/U': 'Th_U'}
    zircon_data.rename(columns=zircon_data_col_rename_dict, inplace=True)

    id_col = 'Sample_ID'
    main_df = samples.reset_index()
    samples_df = main_df.copy()
    analyses_df = zircon_data.reset_index()

    for sample_ind in range(main_df.shape[0]):  # loop through entries in main_df
        active_sample_id = main_df.loc[sample_ind, id_col]
        active_UPb_data = analyses_df.loc[analyses_df[id_col].isin([active_sample_id]), :]

        for col_name in active_UPb_data:
            if col_name not in [id_col]:  # Skip if the indexing column
                # Check column naming overlap with the Samples table (having the same column name will otherwise
                # result in an error)
                if col_name in samples_df.columns:
                    # New name for col_name
                    col_name_adj = col_name + '_' + 'ZrUPb'
                    if col_name_adj not in main_df.columns:
                        main_df[col_name_adj] = (np.nan * np.empty(shape=(len(main_df), 1))).tolist()
                        main_df[col_name_adj] = np.asarray(main_df[col_name_adj])
                    main_df.at[sample_ind, col_name_adj] = active_UPb_data[col_name].values
                else:
                    if col_name not in main_df.columns:  # Make col_name with revised name if already in samples table
                        main_df[col_name] = (np.nan * np.empty(shape=(len(main_df), 1))).tolist()
                        main_df[col_name] = np.asarray(main_df[col_name])
                    main_df.at[sample_ind, col_name] = active_UPb_data[col_name].values

    # Make a copy of the dataset and set the sample ID as index
    main_byid_df = main_df.copy()
    main_byid_df.set_index(id_col, inplace=True, drop=False)

    return main_df, main_byid_df, samples_df, analyses_df, samples, zircon_data


@st.cache(allow_output_mutation=True)
def sampleToData(sampleList, main_byid_df, sampleLabel='Sample_ID', bestAge='BestAge', bestAgeErr='BestAge_err',
                 sigma='1sigma', ID_col='Sample_ID', verify=False):
    """
    after detritalpy
    :param sampleList:
    :param main_byid_df:
    :param sampleLabel:
    :param bestAge:
    :param bestAgeErr:
    :param sigma:
    :param ID_col:
    :param verify:
    :return:
    """

    N = len(sampleList)
    ages = []
    errors = []
    numGrains = []
    labels = []
    stop = False

    if type(sampleList[0]) == tuple:
        for i in range(N):
            samples = sampleList[i][0]
            # Verify that all samples are in the database
            if verify:
                if not all(sample in list(main_byid_df[ID_col]) for sample in sampleList[i][0]):
                    print('These samples are not in the database - check for typos!')
                    print(list(np.setdiff1d(sampleList[i][0], list(main_byid_df[ID_col]))))
                    print('Function stopped')
                    stop = True
                    break
            sampleAges = []
            sampleErrors = []
            for sample in samples:
                sampleAges = np.append(sampleAges, main_byid_df.loc[sample, bestAge])
                if sigma == '2sigma':
                    sampleErrors = np.append(sampleErrors, main_byid_df.loc[sample, bestAgeErr] / 2)
                else:
                    sampleErrors = np.append(sampleErrors, main_byid_df.loc[sample, bestAgeErr])
            ages.append(sampleAges)
            errors.append(sampleErrors)
            numGrains.append(len(sampleAges))
            labels.append(sampleList[i][1])

    else:
        for sample in sampleList:
            # Verify that all samples are in the database
            if verify:
                if not all(sample in list(main_byid_df[ID_col]) for sample in sampleList):
                    print('These samples are not in the database - check for typos!')
                    print(list(np.setdiff1d(sampleList, list(main_byid_df[ID_col]))))
                    print('Function stopped')
                    stop = True
                    break
            ages.append(main_byid_df.loc[sample, bestAge])
            if sigma == '2sigma':
                errors.append(main_byid_df.loc[sample, bestAgeErr] / 2.)
            else:
                errors.append(main_byid_df.loc[sample, bestAgeErr])
            numGrains.append(len(main_byid_df.loc[sample, bestAge]))
            labels.append(main_byid_df.loc[sample, sampleLabel])

    if verify:
        if not stop:  # Only check for missing data if function not terminated beforehand
            if np.min([len(x) for x in ages]) == 0:  # Check whether any of the samples returned no data
                samples_no_data = np.asarray(sampleList)[
                    [len(x) == 0 for x in ages]]  # Return a list of the samples with no data
                print('Warning! These samples have no data: ', samples_no_data)
                print('Please check for consistency in sample naming')
                print('and/or verify that all data have not been filtered')

    return ages, errors, numGrains, labels


def plot_age_distribution(sampleList, main_byid_df, data_subset_name, key_prefix='rc', sample_submit=False):
    st.markdown("""---""")
    st.write("#### Zircon age distribution")
    col1, col2 = st.columns([1, 2])
    ages, errors, numGrains, labels = sampleToData(sampleList, main_byid_df, sigma='1sigma')
    with col1:
        with st.form('age_dis_form_{}'.format(key_prefix)):
            whatToPlot = st.selectbox('Chose the figure to plot',
                                      ('both', 'cumulative', 'relative'),
                                      key=key_prefix
                                      )

            # Specify the age range (Myr) that you want to plot
            plot_age_range = st.slider('Specify the age range (Myr) that you want to plot',
                                       0, 4500,
                                       (0, 3000),
                                       key=key_prefix
                                       )

            b = st.slider('Specify histogram bin size (Myr)',
                          5, 100, 30,
                          key=key_prefix
                          )  # Specify the histogram bin size (Myr)

            with st.expander('Options for plot elements'):
                # Specify the plot dimensions
                w = st.slider('Specify width of figure', 5, 15, 10,
                              key=key_prefix
                              )  # width of the plot
                c = st.slider('Specify height of CDF panel', 2, 10, 4,
                              key=key_prefix
                              )  # height of CDF panel

                h = st.slider('Specify height of relative panel', 4, 16, 6,
                              key=key_prefix
                              )  # height of CDF panel

                plotKDE = st.selectbox('Whether plot KDE?',
                                       (True, False),
                                       key=key_prefix
                                       )  # True # Set to True if want to plot KDE
                colorKDE = False  # Will color KDE according to same coloration as used in CDF plotting
                bw = 'optimizedFixed'
                # True Will color KDE according to age populations if set to True
                colorKDEbyAge = st.selectbox('Whether color KDE according to age populations?',
                                             (True, False),
                                             key=key_prefix
                                             )

                plotPDP = st.selectbox('Whether plot PDP?',
                                       (False, True),
                                       key=key_prefix
                                       )
                colorPDP = False

                colorPDPbyAge = st.selectbox('Whether color PDP according to age populations?',
                                             (True, False),
                                             key=key_prefix
                                             )

                plotHist = st.selectbox('Whether plot histogram?  \n  only available when separateSubplots is True',
                                        (True, False),
                                        key=key_prefix
                                        )

                plotPIE = st.selectbox('Whether plot pie diagram?  \n   only available when separateSubplots is True',
                                       (True, False),
                                       key=key_prefix
                                       )

            age_bins = st.text_input('Define the age categories as following format:',
                                     '0, 100, 140, 180, 220, 260, 500, 600, 900, 1500, 2000, 2200, 2600, 2900, 4500',
                                     key=key_prefix
                                     )

            age_bins_color = st.text_input("Define the age categories associated color as following format:",
                                           "#ffffff,  #ffffbf, #f46d43, #92c5de, "
                                           "#c2a5cf, #a6dba0,  #fddbc7, #018571,"
                                           "#f46d43, #dfc27d, #f1b6da, #3288bd,"
                                           "#e0f3f8, #fee08b",
                                           key=key_prefix
                                           )

            agePeakOptions = ['KDE', 0.05, 5, 2, True]

            submitted = st.form_submit_button("Submit for calculation")

        if len(age_bins) > 10:
            try:
                age_bins = [int(x.strip()) for x in age_bins.split(',')]
            except:
                st.write('please check the input text')
        else:
            age_bins = [0, 80, 160, 200, 600, 1500, 2000, 2500, 4500]

        if len(age_bins_color) > 10:
            try:
                age_bins_color = [str(x.strip()) for x in age_bins_color.split(',')]
            except:
                st.write('please check the input text', key=key_prefix)

            if len(age_bins) - len(age_bins_color) != 1:
                st.write('please check the age categories should be one more than color categories',
                         key=key_prefix
                         )
        else:
            age_bins_color = ['#f0f9e8', '#2b8cbe', '#fdcc8a', '#fc8d59', '#f2f0f7', '#cbc9e2', '#6a51a3', '#ce1256']
        plotColorBar = False
        plotAgePeaks = False

        separateSubplots = True
        x1 = plot_age_range[0]
        x2 = plot_age_range[1]
        plotLog = False

        xdif = 1

        plotCDF = False
        plotCPDP = False
        plotCKDE = True
        plotDKW = False
        normPlots = False

    def plot_age_dis(sampleList, ages, errors, numGrains, labels, whatToPlot, separateSubplots, plotCDF,
                     plotCPDP,
                     plotCKDE,
                     plotDKW, normPlots, plotKDE, colorKDE, colorKDEbyAge, plotPDP, colorPDP, colorPDPbyAge,
                     plotColorBar,
                     plotHist, plotLog, plotPIE, x1, x2, b, bw, xdif, age_bins, age_bins_color, w, c, h,
                     plotAgePeaks,
                     agePeakOptions,
                     CDFlw=3, KDElw=1, PDPlw=1,
                     key_prefix=key_prefix
                     ):

        fig = dFunc.plotAll(sampleList, ages, errors, numGrains, labels, whatToPlot, separateSubplots, plotCDF,
                            plotCPDP,
                            plotCKDE,
                            plotDKW, normPlots, plotKDE, colorKDE, colorKDEbyAge, plotPDP, colorPDP, colorPDPbyAge,
                            plotColorBar,
                            plotHist, plotLog, plotPIE, x1, x2, b, bw, xdif, age_bins, age_bins_color, w, c, h,
                            plotAgePeaks,
                            agePeakOptions,
                            CDFlw=CDFlw, KDElw=KDElw, PDPlw=PDPlw)

        st.session_state['age_dis_{}'.format(key_prefix)] = fig
        return fig

    with col2:
        if 'age_dis_{}'.format(key_prefix) not in st.session_state:
            st.session_state['age_dis_{}'.format(key_prefix)] = 0

        if submitted or \
                st.session_state['age_dis_{}'.format(key_prefix)] == 0 or sample_submit:
            fig = plot_age_dis(sampleList, ages, errors, numGrains, labels, whatToPlot, separateSubplots, plotCDF,
                               plotCPDP,
                               plotCKDE,
                               plotDKW, normPlots, plotKDE, colorKDE, colorKDEbyAge, plotPDP, colorPDP, colorPDPbyAge,
                               plotColorBar,
                               plotHist, plotLog, plotPIE, x1, x2, b, bw, xdif, age_bins, age_bins_color, w, c, h,
                               plotAgePeaks,
                               agePeakOptions,
                               CDFlw=3, KDElw=1, PDPlw=1)
        else:
            fig = st.session_state['age_dis_{}'.format(key_prefix)]
        st.pyplot(fig)
        download_fig_as_pdf(fig, 'Age distribution of  {}.pdf'.format(data_subset_name),
                            key="{}_age_dis".format(key_prefix))


def plot_maximum_depositional_age(sampleList, main_byid_df, data_subset_name, key_prefix='rc', sample_submit=False):
    st.markdown("""---""")
    st.write("#### Maximum Depositional Age")
    col1, col2 = st.columns([1, 2])

    ages, errors, numGrains, labels = dFunc.sampleToData(sampleList, main_byid_df, sigma='1sigma')

    fileName = 'MDA_{}.csv'.format(data_subset_name)
    pathlib.Path('Output').mkdir(parents=True, exist_ok=True)

    with col1:

        makePlot = True

        sortBy = '1sigma'
        with st.form('mda_form_{}'.format(key_prefix)):
            plotWidth = st.slider('Specify width of MDA figure',
                                  5, 15, 8,
                                  key=key_prefix + 'MDA'
                                  )

            plotHeight = st.slider('Specify width of MDA figure',
                                   3, 10, 5,
                                   key=key_prefix + 'MDA'
                                   )
            barWidth = 0.25

            ysg_c = st.color_picker('Pick a color for horizontal bars for YSG',
                                    '#67a9cf',
                                    key=key_prefix + 'MDA'
                                    )
            ysg_1s = st.color_picker('Pick a color for horizontal bars for YC1S',
                                     '#fc8d59',
                                     key=key_prefix + 'MDA'
                                     )
            ysg_2s = st.color_picker('Pick a color for horizontal bars for YC2s',
                                     '#74c476',
                                     key=key_prefix + 'MDA'
                                     )

            alpha = st.slider('Specify the transparency of MDA plots colors',
                              0.1, 1.0, 0.25,
                              key=key_prefix
                              )

            submitted = st.form_submit_button("Submit for calculation")

        ageColors = [ysg_c, ysg_1s, ysg_2s]

        fillMDACalcs = True

    def plot_mda(sampleList, ages, errors, numGrains, labels, fileName, sortBy,
                 barWidth, plotWidth, plotHeight, ageColors, alpha, makePlot, fillMDACalcs, key_prefix=key_prefix):
        figMDA = dFunc.MDAtoCSV(sampleList,
                                ages, errors,
                                numGrains, labels,
                                fileName, sortBy,
                                barWidth, plotWidth,
                                plotHeight,
                                ageColors, alpha,
                                makePlot, fillMDACalcs)
        st.session_state['mda_fig_{}'.format(key_prefix)] = figMDA
        return figMDA

    with col2:
        mda_tag = 'waiting for calculation'
        if 'mda_fig_{}'.format(key_prefix) not in st.session_state:
            st.session_state['mda_fig_{}'.format(key_prefix)] = 0

        if submitted or sample_submit or st.session_state['mda_fig_{}'.format(key_prefix)] == 0:
            figMDA = plot_mda(sampleList,
                              ages, errors,
                              numGrains, labels,
                              fileName, sortBy,
                              barWidth, plotWidth,
                              plotHeight,
                              ageColors, alpha,
                              makePlot, fillMDACalcs)
        else:
            figMDA = st.session_state['mda_fig_{}'.format(key_prefix)]

        st.pyplot(figMDA)
        download_fig_as_pdf(figMDA, 'MDA of  {}.pdf'.format(data_subset_name),
                            key="{}_mda".format(key_prefix)
                            )
        mda_tag = 'completed'

    with col1:
        if mda_tag == 'completed':
            if os.path.exists('Output/{}'.format(fileName)):
                mda_pd = pd.read_csv('Output/{}'.format(fileName))

                st.dataframe(mda_pd)
                download_file = mda_pd.to_csv(index=None).encode('utf-8')

                st.download_button(
                    label="Download data as CSV",
                    data=download_file,
                    file_name=fileName,
                    mime='text/csv',
                    key=key_prefix
                )


def plot_multi_dimensional_scaling(sampleList, main_byid_df, data_subset_name, key_prefix='rc', sample_submit=False):
    st.markdown("""---""")
    st.write("#### Multi-dimensional scaling")
    col1, col2, col3 = st.columns(3)
    with col1:
        with st.form('mds_form_{}'.format(key_prefix)):

            metric = st.selectbox('Whether use the metric MDS?',
                                  (True, False),
                                  key=key_prefix + 'MDS'
                                  )

            criteria = st.selectbox('Select basis for comparison',
                                    ['Vmax', 'Dmax',
                                     'R2-PDP', 'R2-KDE',
                                     'similarity-PDP', 'similarity-KDE',
                                     'likeness-PDP', 'likeness-KDE'],
                                    0,
                                    key=key_prefix + 'MDS'
                                    )
            with st.expander('Click to see the detail information'):
                st.write(
                    '**Vmax**: the sum of teh maximum positive and negative differences between the two CDF Plots ')
                st.write('**Dmax**: the maximum difference between two PDP plots')
                st.write('**R2-PDP**: the Cross-correlation coefficient between two PDP plots')
                st.write('**R2-KDE**: the Cross-correlation coefficient between two KDE plots')
                st.write('**similarity-PDP**: the similarity value between two PDP plots')
                st.write('**similarity-KDE**: the similarity value between two KDE plots')
                st.write('**likeness-PDP**: the likeness value between two PDP plots')
                st.write('**likeness-KDE**: the likeness value between two KDE plots')
                st.write('see Saylor et al. (2017) for more detail')

            x_min = st.slider('Define the low limit (Ma) of the age distributions',
                              0, 5400, 0,
                              key=key_prefix + 'MDS'
                              )
            x_max = st.slider('Define the upper limit (Ma) of the age distributions',
                              0, 5400, 4500,
                              key=key_prefix + 'MDS'
                              )

            plot_mds = st.selectbox('Always recalculate age MDS',
                                    (True, False),
                                    key=key_prefix + 'MDS'
                                    )
            submitted = st.form_submit_button("Submit for calculation")

        def qq_plot():
            fig_qq = model.QQplot(figsize=(8, 8),
                                  savePlot=False, fileName='QQplot.pdf',
                                  halfMatrix=True
                                  )
            fig_qq = plt.gcf()
            st.session_state['qq_fig_{}'.format(key_prefix)] = fig_qq
            return fig_qq

        def heat_plot():
            fig_heatmap = model.heatMap(figsize=(8, 8),
                                        savePlot=False, fileName='HeatMapPlot.pdf',
                                        plotValues=True,
                                        plotType='distance', fontsize=10)
            fig_heatmap = plt.gcf()
            st.session_state['heat_fig_{}'.format(key_prefix)] = fig_heatmap
            return fig_heatmap

        def stress_plot():
            fig_stress = model.stressPlot(figsize=(6, 6),
                                          savePlot=False, fileName='stressPlot.pdf',
                                          stressType='sklearn'
                                          )
            fig_stress = plt.gcf()
            st.session_state['stress_fig_{}'.format(key_prefix)] = fig_stress
            return fig_stress

        def shepard_plot():
            fig_shepard = model.shepardPlot(figsize=(6, 6),
                                            savePlot=False, fileName='shepardPlot.pdf',
                                            plotOneToOneLine=True)
            fig_shepard = plt.gcf()
            st.session_state['shepard_fig_{}'.format(key_prefix)] = fig_shepard
            return fig_shepard

        def mds_plot():
            fig_mds = model.MDSplot(figsize=(8, 6),
                                    savePlot=False, fileName='MDSplot.pdf', plotLabels=True,
                                    equalAspect=False,
                                    stressType='sklearn')
            fig_mds = plt.gcf()
            st.session_state['mds_fig_{}'.format(key_prefix)] = fig_mds
            return fig_mds

        if 'qq_fig_{}'.format(key_prefix) not in st.session_state:
            st.session_state['qq_fig_{}'.format(key_prefix)] = 0
        if 'heat_fig_{}'.format(key_prefix) not in st.session_state:
            st.session_state['heat_fig_{}'.format(key_prefix)] = 0
        if 'stress_fig_{}'.format(key_prefix) not in st.session_state:
            st.session_state['stress_fig_{}'.format(key_prefix)] = 0
        if 'shepard_fig_{}'.format(key_prefix) not in st.session_state:
            st.session_state['shepard_fig_{}'.format(key_prefix)] = 0

        if 'mds_fig_{}'.format(key_prefix) not in st.session_state:
            st.session_state['mds_fig_{}'.format(key_prefix)] = 0

        if st.session_state['mds_fig_{}'.format(key_prefix)] == 0 or submitted or sample_submit:
            ages, errors, numGrains, labels = dFunc.sampleToData(sampleList, main_byid_df, sigma='1sigma')
            model = dFunc.MDS_class(ages, errors, labels, sampleList,
                                    metric=metric, criteria=criteria, bw='optimizedFixed',
                                    n_init='metric',
                                    max_iter=1000, x1=x_min, x2=x_max, xdif=1,
                                    min_dim=1, max_dim=3, dim=2
                                    )
            fig_qq = qq_plot()
            fig_heatmap = heat_plot()
            fig_stress = stress_plot()
            fig_shepard = shepard_plot()
            fig_mds = mds_plot()
        else:
            fig_qq = st.session_state['qq_fig_{}'.format(key_prefix)]
            fig_heatmap = st.session_state['heat_fig_{}'.format(key_prefix)]
            fig_stress = st.session_state['stress_fig_{}'.format(key_prefix)]
            fig_shepard = st.session_state['shepard_fig_{}'.format(key_prefix)]
            fig_mds = st.session_state['mds_fig_{}'.format(key_prefix)]

    with col2:

        st.pyplot(fig_qq)
        download_fig_as_pdf(fig_qq, 'QQplot of {}.pdf'.format(data_subset_name))

        st.pyplot(fig_heatmap)
        download_fig_as_pdf(fig_heatmap, 'heatMap of {}.pdf'.format(data_subset_name))
    with col3:

        st.pyplot(fig_stress)
        download_fig_as_pdf(fig_stress, 'stressPlot of {}.pdf'.format(data_subset_name))

        st.pyplot(fig_shepard)
        download_fig_as_pdf(fig_shepard, 'shepardPlot of {}.pdf'.format(data_subset_name))

    col1, col2 = st.columns([1, 2])
    with col2:
        st.pyplot(fig_mds)
        download_fig_as_pdf(fig_mds, 'MdaMap of  {}.pdf'.format(data_subset_name))


@st.cache
def prepare_data_for_lu_hf_ploter(data, sample_id='unique_sample_id', custom_col='custom'):
    """
    prepare the data for lu-hf visualization
    :param sample_id:
    :param custom_col:
    :param data: the data with the columns define by cols
    :return: the samples which contain lu-hf isotopic measurements
    """
    if custom_col not in data.columns:
        data[custom_col] = None
    cols = ['unique_sample_id', 'sample_spots_id', 'cal_best_age', 'cal_best_age_1se', 'cal_disc', 'cal_epsilon_Hf(t)',
            'cal_Hf_TDM1(Ma)', 'cal_Hf_TDM2(Ma)', custom_col]

    col_rename_dict = {'unique_sample_id': 'unique_sample_id',
                       'cal_epsilon_Hf(t)': 'epsilon_Hf(t)',
                       'cal_Hf_TDM1(Ma)': 'Hf_TDM1(Ma)',
                       'cal_Hf_TDM2(Ma)': 'Hf_TDM2(Ma)',
                       'sample_spots_id': 'Spot',
                       'cal_best_age': 'BestAge',
                       'cal_best_age_1se': 'BestAge_err',
                       'cal_disc': 'Disc',
                       'Th/U': 'Th/U'}
    count_data = data.groupby([sample_id])['cal_epsilon_Hf(t)'].count()
    drop_list = list(count_data[count_data < 20].index)
    data = data.set_index(sample_id)
    data = data.drop(drop_list, axis=0)
    data = data.reset_index()

    temp_data = data.dropna(subset=['cal_epsilon_Hf(t)', 'cal_best_age'], axis=0)
    sample_list = list(set(temp_data[sample_id]))
    clean_data = data[cols]
    clean_data = clean_data.loc[clean_data[sample_id].isin(sample_list)]
    clean_data.rename(columns=col_rename_dict, inplace=True)
    sample_list = list(set(clean_data[sample_id]))

    return sample_list, clean_data


def plot_kde_lu_hf(subplot_list, data, data_subset_name, sample_id='unique_sample_id',
                   x_col='BestAge', y_col='epsilon_Hf(t)', key_prefix='lu-hf-kde', sample_submit=False):
    """
    plot the  2d-kde of lu-hf isotopic data
    :param key_prefix:
    :param data_subset_name:
    :param subplot_list: which should be a sequence defining the sample list or group list
    :param data: the clean lu-hf data
    :param sample_id: the sample column name
    :param x_col:the x-axis variable name
    :param y_col:the x-axis variable name
    :return: plt
    """
    st.markdown("""---""")
    st.write("#### KDE plots of Lu-Hf isotopic data")
    col1, col2 = st.columns([1, 2])

    data = data
    sample_id = sample_id
    subplot_list = subplot_list
    x_col = x_col
    y_col = y_col
    col = 2
    row = math.ceil(len(subplot_list) / col)
    width = col * 5
    hight = row * 3

    # st options
    with col1:
        with st.form('kde_form_{}'.format(key_prefix)):
            ked_color_map = st.selectbox('Pick a color pattern for 2D-KDE?',
                                         ('RdBu_r', 'vlag', 'rocket_r', 'Blues', 'Greens', 'Reds', 'coolwarm'),
                                         key=key_prefix
                                         )

            plot_scatter = st.selectbox('Whether plot data point on the 2D-KDE?',
                                        (True, False),
                                        key=key_prefix
                                        )
            plot_95 = st.selectbox('Whether plot 95% percent contour?',
                                   (True, False),
                                   key=key_prefix
                                   )
            # lint_95_color = '#fb6a4a'
            lint_95_color = st.color_picker('Pick a color for 95% percent contour line',
                                            '#fb6a4a',
                                            key=key_prefix + 'MDA'
                                            )
            evolution_line = st.selectbox('Whether plot the Evolution Lines ?',
                                          (True, False),
                                          key=key_prefix
                                          )
            submitted = st.form_submit_button("Submit for calculation")

    sharex = False
    sharey = False
    custom_xlim = (0, 4000)
    custom_ylim = (-30, 30)

    # constant value
    Decay_const_176Lu = 0.01867
    # 176Lu decay constant (Scherer et al., 2001) 1.867*10^-11 (same as Soderland et al., 2004)

    DM_176Hf_177Hf = 0.283225
    # Vervoort and Blichert-Toft, 1999

    DM_176Lu_177Hf = 0.0383
    # Vervoort and Blichert-Toft, 1999
    BSE_176Hf_177Hf = 0.282785
    # Bouvier et al. 2008
    BSE_176Lu_177Hf = 0.0336
    # Bouvier et al. 2008

    np.arange(0, 5400, 5)
    hf_evol_df = pd.DataFrame(np.arange(-2000, 5600, 5), columns=['age'])

    hf_evol_df['t_176Hf_177Hf'] = DM_176Hf_177Hf - (
            DM_176Lu_177Hf * (np.exp(Decay_const_176Lu * hf_evol_df['age'] / 1000) - 1))
    hf_evol_df['t_chur'] = BSE_176Hf_177Hf - (
            BSE_176Lu_177Hf * (np.exp(Decay_const_176Lu * hf_evol_df['age'] / 1000) - 1))

    hf_evol_df['epsilon_dm'] = 10000 * ((hf_evol_df['t_176Hf_177Hf'] / hf_evol_df['t_chur']) - 1)
    hf_evol_df['epsilon_chur'] = 0

    def plot_2dkde():
        fig, axs = plt.subplots(row, col, figsize=(width, hight), sharex=sharex, sharey=sharey)
        subplot_num = 0
        # subplot_num = 0
        for i in range(row):
            for j in range(col):

                if row == 1 and col == 1:
                    c_ax = axs
                elif row == 1:
                    c_ax = axs[j]
                elif col == 1:
                    c_ax = axs[i]
                else:
                    c_ax = axs[i, j]
                if subplot_num >= len(subplot_list):
                    fig.delaxes(axs[i, j])
                    break

                plot_data = data.loc[data[sample_id] == subplot_list[subplot_num]]

                sns.kdeplot(data=plot_data,
                            x=x_col,
                            y=y_col,
                            kind="kde",
                            ax=c_ax,
                            fill=True,
                            #                     levels=6,
                            cmap=ked_color_map)  # vlag, RdBu_r, rocket_r, Blues, Greens, Reds,coolwarm
                if plot_95:
                    sns.kdeplot(data=plot_data,
                                x=x_col,
                                y=y_col,
                                kind="kde",
                                ax=c_ax,
                                fill=False,
                                levels=[0.05],
                                color=lint_95_color,
                                linestyle='dashed',
                                linewidth=2
                                )
                if plot_scatter:
                    #             s=10, color='#f47e26'
                    sns.scatterplot(data=plot_data,
                                    x=x_col,
                                    y=y_col,
                                    s=10,
                                    color='white',
                                    ax=c_ax)

                if evolution_line:
                    ax_xlim = c_ax.get_xlim()
                    plot_evol = hf_evol_df.loc[(hf_evol_df['age'] >= ax_xlim[0]) &
                                               (hf_evol_df['age'] <= ax_xlim[1])]
                    c_ax.plot(plot_evol['age'],
                              plot_evol['epsilon_dm'],
                              color='gray',
                              linewidth=2,
                              label='DM'
                              )
                    c_ax.plot(plot_evol['age'],
                              plot_evol['epsilon_chur'],
                              color='gray',
                              linestyle='dashed',
                              linewidth=1,
                              label='CHUR'
                              )

                c_ax.set_xlabel('Age(Ma)')
                c_ax.set_ylabel('$\epsilon$Hf(t)')
                c_ax.set_xlabel('Age(Ma)')

                c_ax.set_title(subplot_list[subplot_num] + '  ',
                               y=0.1, loc='right', fontsize=11, fontweight='bold')

                subplot_num += 1
        plt.subplots_adjust(wspace=0.3,
                            hspace=0.35
                            )
        st.session_state['lu_hf_kde_fig_{}'.format(key_prefix)] = fig
        return fig

    with col2:
        if 'lu_hf_kde_fig_{}'.format(key_prefix) not in st.session_state:
            st.session_state['lu_hf_kde_fig_{}'.format(key_prefix)] = 0
        if submitted or sample_submit or st.session_state['lu_hf_kde_fig_{}'.format(key_prefix)] == 0:
            fig = plot_2dkde()
        else:
            fig = st.session_state['lu_hf_kde_fig_{}'.format(key_prefix)]

        st.pyplot(fig)
        download_fig_as_pdf(fig, 'Lu-Hf kde plot of  {}.pdf'.format(data_subset_name),
                            key="{}_age_dis".format(key_prefix))


def plot_violins(sample_list, data, data_subset_name, sample_id='unique_sample_id',
                 x_col='BestAge', y_col='epsilon_Hf(t)', key_prefix='lu-hf', sample_submit=False):
    """
    plot violins plot of the lu-hf isotopic data
    :param data_subset_name:
    :param key_prefix:
    :param sample_list: which should be a sequence defining the sample list or group list
    :param data: the clean lu-hf data
    :param sample_id: the sample column name
    :param x_col: x-axis variable name
    :param y_col: y-axis variable name
    :return: plt
    """
    st.markdown("""---""")
    st.write("#### Violins plots of Lu-Hf isotopic data")
    col1, col2 = st.columns([1, 2])

    data = data
    sample_id = sample_id
    sample_list = sample_list
    x_col = x_col
    y_col = y_col
    with col1:
        with st.form('vl_form_{}'.format(key_prefix)):
            # st option
            # y_lim = [-30, 30]
            inner = st.selectbox('Pick up the style of the data points in the violin interior?',
                                 ('quartile', 'box', 'point', 'stick'),
                                 key=key_prefix
                                 )

            y_lim = st.slider('Specify ɛHf(t) range  that you want to plot',
                              -60, 60,
                              (-30, 30),
                              key=key_prefix
                              )
            # violins_num_of_each_subplot = 6
            violins_num_of_each_subplot = st.slider('Specify sample number in each panel',
                                                    4, 8, 6,
                                                    key=key_prefix
                                                    )
            submitted = st.form_submit_button("Submit for calculation")

    sharey = False
    data = data.loc[(data[y_col] <= y_lim[1]) & (data[y_col] >= y_lim[0])]
    sub_plot_num = math.ceil(len(sample_list) / violins_num_of_each_subplot)
    if len(sample_list) <= violins_num_of_each_subplot:
        col = 1
    else:
        col = 1
    row = math.ceil(sub_plot_num / col)

    width = col * 8
    height = row * 3

    def plot_vl():
        fig, axs = plt.subplots(row, col, figsize=(width, height), sharey=sharey)

        subplot_num = 0
        for i in range(row):
            for j in range(col):

                if row == 1 and col == 1:
                    c_ax = axs
                elif row == 1:
                    c_ax = axs[j]
                elif col == 1:
                    c_ax = axs[i]
                else:
                    c_ax = axs[i, j]
                s = violins_num_of_each_subplot * subplot_num
                e = violins_num_of_each_subplot * (subplot_num + 1)
                if e <= len(sample_list):
                    subplot_sample = sample_list[s:e]
                else:
                    subplot_sample = sample_list[s:]

                if violins_num_of_each_subplot * subplot_num > len(sample_list):
                    break

                plot_data = data.loc[data[sample_id].isin(subplot_sample)]
                plot_data.dropna(subset=[y_col], axis=0, inplace=True)
                #         print(subplot_sample)

                sns.violinplot(x=sample_id,
                               y=y_col,
                               data=plot_data,
                               palette="Set2",
                               inner=inner,
                               ax=c_ax)

                c_ax.set_xlabel("")
                c_ax.set_ylim(y_lim)
                c_ax.set_ylabel('$\epsilon$Hf(t)')
                c_ax.grid(ls='--')
                for tick in c_ax.get_xticklabels():
                    tick.set_rotation(20)

                subplot_num += 1

        plt.subplots_adjust(wspace=0.2,
                            hspace=0.6
                            )
        st.session_state['lu_hf_vl_fig_{}'.format(key_prefix)] = fig
        return fig

    with col2:
        if 'lu_hf_vl_fig_{}'.format(key_prefix) not in st.session_state:
            st.session_state['lu_hf_vl_fig_{}'.format(key_prefix)] = 0
        if submitted or sample_submit or st.session_state['lu_hf_vl_fig_{}'.format(key_prefix)] == 0:
            fig = plot_vl()
        else:
            fig = st.session_state['lu_hf_vl_fig_{}'.format(key_prefix)]

        st.pyplot(fig)
        download_fig_as_pdf(fig, 'Lu-Hf violin plot of  {}.pdf'.format(data_subset_name),
                            key="{}_age_dis".format(key_prefix))


def load_gplates_model(model_name, pb=False):
    """
    load the gplates model
    :param model_name: Gplates model name
    :param pb: weather read the plate boundary file
    :return:
    """
    plate_boundaries = None
    # location of files we want to use
    if model_name == 'Merdith_2022':
        model_dir = 'Gplates_models/{}/'.format(model_name)
        rotation_filename = ['{}1000_0_rotfile_Merdith_et_al.rot'.format(model_dir)]
        static_polygon_file = '{}shapes_static_polygons_Merdith_et_al.gpml'.format(model_dir)
        plate_boundaries_file = ['{}250-0_plate_boundaries_Merdith_et_al.gpml'.format(model_dir),
                                 '{}410-250_plate_boundaries_Merdith_et_al.gpml'.format(model_dir),
                                 '{}TopologyBuildingBlocks_Merdith_et_al.gpml'.format(model_dir)
                                 ]
        coastline_filename = 'Gplates_models/Global_EarthByte_GPlates_PresentDay_Coastlines.gpmlz'
    elif model_name == 'Young_2018':
        model_dir = 'Gplates_models/{}/'.format(model_name)
        rotation_filename = ['{}Global_250-0Ma_Young_et_al.rot'.format(model_dir),
                             '{}Global_410-250Ma_Young_et_al.rot'.format(model_dir)
                             ]

        static_polygon_file = '{}Global_EarthByte_GPlates_PresentDay_StaticPlatePolygons.gpmlz'.format(model_dir)

        plate_boundaries_file = ['{}Global_Mesozoic-Cenozoic_plate_boundaries_Young_et_al.gpml'.format(model_dir),
                                 '{}410-250_plate_boundaries_Merdith_et_al.gpml'.format(model_dir),
                                 '{}TopologyBuildingBlocks_Young_et_al.gpml'.format(model_dir)]
        coastline_filename = 'Gplates_models/Global_EarthByte_GPlates_PresentDay_Coastlines.gpmlz'
    else:
        print('Please check the gplate model name')

    # read into pygplates
    rotation_model = pygplates.RotationModel(rotation_filename)
    static_polygons = pygplates.FeatureCollection(static_polygon_file)
    coastlines = pygplates.FeatureCollection(coastline_filename)
    if pb:
        print(plate_boundaries_file)
        plate_boundaries = pygplates.FeatureCollection()
        for file in plate_boundaries_file:
            print(file)
            plate_boundaries.add(file)

    # add record to the st.session_state
    st.session_state['rt_gplates'] = rotation_model
    st.session_state['sp_gplates'] = static_polygons
    st.session_state['cs_gplates'] = coastlines
    st.session_state['pb_gplates'] = plate_boundaries
    st.session_state['model_name'] = model_name

    return rotation_model, static_polygons, coastlines, plate_boundaries


def data_to_gplates(data, feature_attribute, static_polygons, rotation_model, lat='latitude', lon='longitude'):
    """
    convert dataframe to gplates feature
    :param data: pd.Dataframe
    :param feature_attribute: str, list of str, the attribute will add  to  Gplates feature
    :param static_polygons: gplates feature, used to assign the Plate-ID for the new feature
    :param rotation_model: gplates rotation model
    :param lat: latitude column name
    :param lon: longitude column name
    :return: gplates feature
    """
    data_point_features = []
    data.dropna(subset=['YC1S WM', lat, lon], inplace=True)
    data = data.loc[(data['YC1S WM'] >= 1) & (data[lon] <= 180)]

    for index, row in data.iterrows():
        point = pygplates.PointOnSphere(float(row[lat]),
                                        float(row[lon]))  # use the lat and lng columns as our coordinates
        data_point_feature = pygplates.Feature()
        data_point_feature.set_geometry(point)

        # set other attributes that aren't native to pygplates
        for i in feature_attribute:
            if i in data.columns:
                data_point_feature.set_shapefile_attribute(i, row[i])

        # set time the point is valid for
        data_point_feature.set_valid_time(row['YC1S WM'], 0)  # assign ages for each point
        data_point_features.append(data_point_feature)  # save to pbdb_point_features

    # add plate id to features
    data_point_features = pygplates.partition_into_plates(static_polygons,
                                                          rotation_model,
                                                          data_point_features,
                                                          properties_to_copy=[
                                                              pygplates.PartitionProperty.reconstruction_plate_id])

    return data_point_features


def reconstructed_to_time(rotation_model, feature, reconstruction_time):
    """
    reconstruct Gplates feature to a specific time
    :param rotation_model: str, the model name
    :param feature: gplates object
    :param reconstruction_time: int, the reconstruction time
    :return:
    """
    reconstructed_data = []
    pygplates.reconstruct(feature, rotation_model, reconstructed_data, reconstruction_time)
    return reconstructed_data


def create_geodataframe_gplates(pygplates_recon_geom, feature_attribute, reconstruction_time):
    """
    convert Gplates feature to GeoDataFrame, the input geometry must be a point
    :param pygplates_recon_geom: Gplates feature, which must be a point
    :param feature_attribute: str, list of str, the attributes which convert to  DataFrame columns
    :param reconstruction_time: int, the target reconstruction time
    :return:
    """

    # create new and empy geodataframe
    recon_gpd = gpd.GeoDataFrame()
    arrti_dict = dict()
    for attri in feature_attribute:
        recon_gpd[attri] = None
        arrti_dict[attri] = list()
    # set crs
    arrti_dict['geometry'] = list()
    arrti_dict['reconstruction_time'] = list()
    arrti_dict['plate_id'] = list()
    arrti_dict['start_time'] = list()
    arrti_dict['end_time'] = list()

    recon_gpd['geometry'] = None
    recon_gpd['reconstruction_time'] = None
    recon_gpd = recon_gpd.set_crs(epsg=4326)

    for i in pygplates_recon_geom:
        point = i.get_feature()
        # contains present-day coordinates
        # plat, plon = point.get_geometry().to_lat_lon()

        # reconstructed coordinates
        recon_point = i.get_reconstructed_geometry()
        rlat, rlon = recon_point.to_lat_lon()

        # get feature attributes
        arrti_dict['start_time'].append(point.get_valid_time()[0])
        arrti_dict['end_time'].append(point.get_valid_time()[1])
        arrti_dict['plate_id'].append(point.get_reconstruction_plate_id())

        # get other attributes
        for attri in feature_attribute:
            arrti_dict[attri].append(point.get_shapefile_attribute(attri))

        arrti_dict['reconstruction_time'].append(reconstruction_time)

        # convert point into shapely geometry
        geometry = Point(rlon, rlat)
        arrti_dict['geometry'].append(geometry)

    # print(accepted_names)
    # write to geodataframe

    recon_gpd['geometry'] = arrti_dict['geometry']
    recon_gpd['reconstruction_time'] = arrti_dict['reconstruction_time']
    recon_gpd['plate_id'] = arrti_dict['plate_id']
    recon_gpd['start_time'] = arrti_dict['start_time']
    recon_gpd['end_time'] = arrti_dict['end_time']
    for attri in feature_attribute:
        recon_gpd[attri] = arrti_dict[attri]

    return recon_gpd


def reconstruct_and_plot(coastlines, rotation_model, point_features,
                         gpl_point_att, reconstruction_time, ax,
                         marker_size=30, alpha=0.6, color='red'
                         ):
    """
    plot the reconstructed data points and coastlines
    :param coastlines: Gplates feature
    :param rotation_model: Gplates rotation model
    :param point_features: Gplates feature, the point feature
    :param gpl_point_att: str, str_list
    :param reconstruction_time: int
    :param ax: the matplotlib axis
    :return:
    """
    # reconstruct fossils and coastlines
    reconstruction_time = reconstruction_time

    # reconstruct data point
    reconstructed_data = reconstructed_to_time(rotation_model,
                                               point_features,
                                               reconstruction_time)
    # convert reconstructed gplates points data to geopandas
    reconstructed_gpd = create_geodataframe_gplates(reconstructed_data,
                                                    gpl_point_att,
                                                    reconstruction_time)

    # reconstructed coastline
    reconstructed_coastlines = []
    pygplates.reconstruct(coastlines,
                          rotation_model,
                          reconstructed_coastlines,
                          reconstruction_time,
                          export_wrap_to_dateline=True)

    # ---- plotting parts of the function
    # define colours
    marker_size = marker_size
    alpha = alpha
    color = color

    # ax.set_global()
    # plot present-day coastlines
    ax.coastlines(color='black', linewidth=0.75, resolution='50m', alpha=0.2)

    # plot reconstructed coastlines
    coastlines_for_plotting = []
    date_line_wrapper = pygplates.DateLineWrapper()
    for polygon in reconstructed_coastlines:
        wrapped_polygons = date_line_wrapper.wrap(polygon.get_reconstructed_geometry())
        for poly in wrapped_polygons:
            coastlines_for_plotting.append(
                LinearRing([(p.get_longitude(), p.get_latitude())
                            for p in poly.get_exterior_points()]))
    ax.add_geometries(coastlines_for_plotting, ccrs.PlateCarree(), facecolor='grey', alpha=0.25)

    # plot pbdb points
    reconstructed_gpd.plot(ax=ax, alpha=alpha, edgecolor='w', color=color,
                           markersize=marker_size, transform=ccrs.PlateCarree())
    #     ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)

    ax.set_title('  %s Ma' % reconstruction_time, fontsize=12, loc='left', y=0.90, fontweight='bold')


def sequence_reconstruction_plot(coastlines, rotation_model, point_features,
                                 gpl_point_att, key_prefix='rm'):
    """
    plot  sequence reconstructed points and coastlines
    :param coastlines: Gplates feature, coastline
    :param rotation_model: Gplates object, rotation_model
    :param point_features: Gplates feature, point feature
    :param gpl_point_att: str, str_list
    :param time_sequence: int, int_list, define the reconstructed time point
    :return:
    """
    st.markdown("""---""")
    st.write("#### Plot the Reconstructed Map")
    col1, col2 = st.columns([1, 2])
    with col1:
        with st.form('{}_rmp'.format(key_prefix)):
            time_sequence = st.text_input('Define the reconstruction time sequence as following format:',
                                          '600, 500, 400, 300, 200, 100, 50, 10',
                                          key='rm_seq_{}'.format(key_prefix)
                                          )
            plot_type = st.selectbox('Define which samples will be plotted',
                                     ('total', 'selected'),
                                     key='pt_{}'.format(key_prefix)
                                     )
            marker_size = st.slider('Define points marker size',
                                    10, 80, 30,
                                    key='ms_{}'.format(key_prefix)
                                    )
            map_type = st.selectbox('Pick up map type',
                                    ('global', 'regional'),
                                    key='map_t_{}'.format(key_prefix)
                                    )
            # with st.expander('If regional plot, please click to define the map extend'):
            lon_extend = list(st.slider('Define longitude extend of map',
                                        -180, 180,
                                        (30, 170),
                                        key='long_ex_{}'.format(key_prefix)
                                        ))
            lat_extend = st.slider('Define latitude extend of map',
                                   -90, 90,
                                   (-50, 85),
                                   key='lat_ex_{}'.format(key_prefix)
                                   )
            alpha = st.slider('Define points transparency',
                              0.0, 1.0, 0.6,
                              key='p_alpha_{}'.format(key_prefix)
                              )
            color = st.color_picker('Define the clor of the Point',
                                    '#b30b00',
                                    key='pc_{}'.format(key_prefix)
                                    )
            submitted = st.form_submit_button("Submit for calculation")

    try:
        time_sequence = [int(x.strip()) for x in time_sequence.split(',')]
    except:
        print('Please check the input text in the time sequence box')
        time_sequence = [600, 500, 400, 300, 200, 100, 50, 0]

    if plot_type == 'total':
        point_features = st.session_state['total_sample_feature']

    def plot_rm():

        map_extent = [lon_extend[0], lon_extend[1], lat_extend[0], lat_extend[1]]
        if len(time_sequence) <= 1:
            col = 1
        else:
            col = 2
        row = math.ceil(len(time_sequence) / col)
        # print(row)
        width = col * 6
        height = row * 4

        # map_type = 'regional'
        if map_type == 'regional':
            projection = ccrs.PlateCarree()
            fig, axs = plt.subplots(row, col, figsize=(width, height),
                                    subplot_kw={'projection': projection}
                                    )


        elif map_type == 'global':
            projection = ccrs.Robinson(central_longitude=120)
            fig, axs = plt.subplots(row, col, figsize=(width, height),
                                    subplot_kw={'projection': projection}
                                    )

        subplot_num = 0
        for i in range(row):
            for j in range(col):

                if row == 1 and col == 1:
                    c_ax = axs
                elif row == 1:
                    c_ax = axs[j]
                elif col == 1:
                    c_ax = axs[i]
                else:
                    c_ax = axs[i, j]

                if subplot_num >= len(time_sequence):
                    fig.delaxes(c_ax)
                    break
                if map_type == 'regional':
                    c_ax.set_extent(map_extent, crs=ccrs.PlateCarree())
                elif map_type == 'global':
                    c_ax.set_global()

                reconstruct_and_plot(coastlines,
                                     rotation_model,
                                     point_features,
                                     gpl_point_att,
                                     time_sequence[subplot_num],
                                     c_ax,
                                     marker_size=marker_size, alpha=alpha, color=color)
                subplot_num += 1
        plt.subplots_adjust(wspace=0.1, hspace=0.1)
        st.session_state['rm_map'] = fig
        return fig

    with col2:
        if 'rm_map' not in st.session_state:
            st.session_state['rm_map'] = 0

        if st.session_state['rm_map'] == 0 or submitted:
            fig = plot_rm()
        else:
            fig = st.session_state['rm_map']

        st.pyplot(fig)
        download_fig_as_pdf(fig,
                            'Robinson Map of reconstructed feature {}.pdf'.format(data_subset_name),
                            key="paleo_map")


def get_the_paleolocation_of_point(rotation_model, point_features,
                                   gpl_point_att, time_sequence):
    """
    get the paleolocation of points data
    :param rotation_model: Gplates object, rotation model
    :param point_features: Gplates feature, point feature
    :param gpl_point_att: Gplates feature's attributes
    :param time_sequence: int, int_list
    :return:
    """
    if len(point_features) > 10:
        point_features = point_features[:10]

    out_gpd = gpd.GeoDataFrame()

    for time in time_sequence:

        re_points = reconstructed_to_time(rotation_model, point_features, time)
        re_gpd = create_geodataframe_gplates(re_points, gpl_point_att, time)

        if out_gpd.empty:
            out_gpd = re_gpd
        else:
            out_gpd = out_gpd.append(re_gpd)

    out_gpd['paleo_lon'] = out_gpd['geometry'].x
    out_gpd['paleo_lat'] = out_gpd['geometry'].y
    out_gpd.reset_index(inplace=True)

    return out_gpd


def plot_paleo_lat_lon(rotation_model, point_features, gpl_point_att, key_prefix='rpp'):
    st.markdown("""---""")
    st.write("#### Plot the reconstructed Palaeolongitude and Palaeolatitude")
    col1, col2 = st.columns([1, 2])

    with col1:
        with st.form('{}_rpp'.format(key_prefix)):
            time_range = st.slider('Specify reconstruction time range (Myr)',
                                   0, 1000,
                                   (0, 300),
                                   key='time_range'
                                   )
            time_increment = st.slider('Pick up time increment of reconstruction time sequence',
                                       1, 30, 2,
                                       key='time increment'
                                       )

            lw = st.slider('Define the line width',
                           0.5, 6.0, 3.0,
                           key=key_prefix
                           )
            len_p = st.slider('Define the legend location',
                              0.1, 1.2, 1.02,
                              key=key_prefix
                              )
            draw_style = st.selectbox('Pick draw style of lines',
                                      ('default', 'steps', 'steps-pre', 'steps-mid', 'steps-post'),
                                      key='ds_{}'.format(key_prefix)
                                      )

            submitted = st.form_submit_button("Submit for calculation")

        time_sequence = range(time_range[0], time_range[1], time_increment)
        paleo_location = get_the_paleolocation_of_point(rotation_model,
                                                        point_features,
                                                        gpl_point_att,
                                                        time_sequence
                                                        )

        paleo_location_pd = pd.DataFrame(paleo_location.drop(columns='geometry'))

        def plot_rpp():
            fig, axs = plt.subplots(1, 2, figsize=(10, 8))
            sns.lineplot(data=paleo_location_pd,
                         x='reconstruction_time',
                         y='paleo_lon',
                         hue='sample',
                         ax=axs[0],
                         linewidth=lw,
                         palette="Paired",
                         drawstyle=draw_style
                         )
            axs[0].legend([], [], frameon=False)
            axs[0].grid(linestyle='--')
            axs[0].invert_xaxis()
            axs[0].set_xlabel('Reconstruction Age (Myr)')
            axs[0].set_ylabel('Palaeolongitude')
            sns.lineplot(data=paleo_location_pd,
                         x='reconstruction_time',
                         y='paleo_lat',
                         hue='sample',
                         ax=axs[1],
                         linewidth=lw,
                         palette="Paired",
                         drawstyle=draw_style
                         )
            axs[1].grid(linestyle='--')
            axs[1].invert_xaxis()
            axs[1].legend(bbox_to_anchor=(len_p, 1), loc='upper left', borderaxespad=0)
            axs[1].set_xlabel('Reconstruction Age (Myr)')
            axs[1].set_ylabel('Palaeolatitude')
            st.session_state['rpp_fig_{}'.format(key_prefix)] = fig
            return fig

    with col2:
        if 'rpp_fig_{}'.format(key_prefix) not in st.session_state:
            st.session_state['rpp_fig_{}'.format(key_prefix)] = 0

        if submitted or st.session_state['rpp_fig_{}'.format(key_prefix)] == 0:
            fig = plot_rpp()
        else:
            fig = st.session_state['rpp_fig_{}'.format(key_prefix)]

        st.pyplot(fig)
        download_fig_as_pdf(fig,
                            'Robinson Map of reconstructed feature {}.pdf'.format(data_subset_name),
                            key="paleo_lon_lat")

        download_file = paleo_location_pd.to_csv(index=None).encode('utf-8')
        st.download_button(
            label="Download data as CSV",
            data=download_file,
            file_name='Paleo_location.csv',
            mime='text/csv',
            key='paleo_lon_lon_dl'
        )


def prepare_template_file(template_file):
    """
    get the template file URL
    :param template_file: file path of the template file
    :return:
    """
    sample = pd.read_excel(template_file, sheet_name='sample')
    zircon_data = pd.read_excel(template_file, sheet_name='data')
    output = io.BytesIO()
    with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
        sample.to_excel(writer, index=False, sheet_name='sample')
        zircon_data.to_excel(writer, index=False, sheet_name='data')
    writer.save()

    val = output.getvalue()
    b64 = base64.b64encode(val)
    return f'<a href="data:application/octet-stream;base64,{b64.decode()}" ' \
           f'download="EADZ_template.xlsx">Download file</a>'  # decode b'abc' => abc


def read_user_file(user_data_file, data_base_col=None):
    """
    read the user upload file
    :param data_base_col: the database sample name list
    :param user_data_file:
    :return:
    """
    sample = pd.read_excel(user_data_file, sheet_name='sample')
    zircon_data = pd.read_excel(user_data_file, sheet_name='data')
    data_dict = {'best_age': 'cal_best_age',
                 'best_age_1se': 'cal_best_age_1se',
                 'discordance': 'cal_discordance',
                 'sample_spots_id': 'sample_spots_id',
                 'f(lu/hf)': 'cal_lu/hf',
                 'epsilon_Hf(t)': 'cal_epsilon_Hf(t)',
                 'Hf_TDM1(Ma)': 'cal_Hf_TDM1(Ma)',
                 'Hf_TDM2(Ma)': 'cal_Hf_TDM2(Ma)'
                 }
    zircon_data.rename(columns=data_dict, inplace=True)
    for i in sample['unique_sample_id']:
        if i in data_base_col:
            st.sidebar.write('Sample Name: {} already in the database, Please rename'.format(i))
    return sample, zircon_data


st.set_option('deprecation.showPyplotGlobalUse', False)
st.set_page_config(page_title="Detrital  Zircon Database of East Asia", layout="wide")

# data read sequence
work_dir = Path('<replace with your work directory>')
os.chdir(work_dir)
data_dir = "<replace with your data directory>"
sample_file = "<replace with your sample information file name>"
data_file = "<replace with your data file name>"

columns_list = ['id_order', 'order', 'sample_spots_id', 'record_num', 'sample_id', 'spot', 'sample', 'data_tag',
                'unique_sample_id', 'age(Ma)', 'Age(Ma)_1se', 'age_206Pb/238U', 'age_206Pb/238U_1se', 'age_207Pb/206Pb',
                'age_207Pb/206Pb_1se', 'age_207Pb/235U', 'age_207Pb/235U_1se', 'age_208Pb/232Th', 'age_208Pb/232Th_1se',
                'ele_206Pb/238U', 'ele_206Pb/238U_1se', 'ele_207Pb/206Pb', 'ele_207Pb/206Pb_1se',
                'ele_207Pb/235U', 'ele_207Pb/235U_1se', 'ele_208Pb/232Th', 'ele_208Pb/232Th_1se',
                'ele_238Pb/206U', 'ele_238Pb/206U_1se', '176Hf/177Hf', '176Hf/177Hf(cor)', '176Hf/177Hf(t)',
                '176Hf/177Hf_1se', '176Lu/177Hf', '176Lu/177Hf_1se', '176Yb/177Hf', '176Yb/177Hf_1se',
                'Ce', 'Dy', 'Er', 'Eu', 'Gd', 'Hf', 'Hg', 'Ho', 'La', 'Lu', 'Nb', 'Nd', 'P', 'Pb', 'Pr', 'Sc', 'Sm',
                'Ta', 'Tb', 'Th', 'Th/U', 'Ti', 'Tm', 'Tzr(°C)', 'U', 'Y', 'Yb', 'Zr',
                '176Hf/177Hf(cor)_1se', '176Hf/177Hf(t)_1se', '204Pb/206Pb_1se', '206Pb/204Pb_1se',
                '208Pb/206Pb_1se', 'O18/O16', 'delta_O18(%)_1se', 'cal_best_age', 'cal_best_age_1se', 'best_age_tag',
                'cal_discondance', 'cal_dis_tag', 'cal_176Hf/177Hf(t)', 'cal_epsilon_Hf(t)', 'cal_lu/hf',
                'cal_Hf_TDM1(Ma)', 'cal_Hf_TDM2(Ma)', 'T_zr_Ti', 'st_sample', 'st_spot', 'Tzr_Ti(°C)']

age_cols = ['st_sample', 'sample_spots_id', 'unique_sample_id', 'age_206Pb/238U', 'age_206Pb/238U_1se',
            'age_207Pb/206Pb',
            'age_207Pb/206Pb_1se',
            'age_207Pb/235U', 'age_207Pb/235U_1se', 'age_208Pb/232Th', 'age_208Pb/232Th_1se', 'ele_206Pb/238U',
            'ele_206Pb/238U_1se', 'ele_207Pb/206Pb', 'ele_207Pb/206Pb_1se', 'ele_207Pb/235U', 'ele_207Pb/235U_1se',
            'ele_208Pb/232Th', 'ele_208Pb/232Th_1se']

lu_hf_cols = ['176Hf/177Hf', '176Hf/177Hf_1se', '176Hf/177Hf(t)_1se', '176Lu/177Hf',
              '176Lu/177Hf_1se', '176Yb/177Hf', '176Yb/177Hf_1se', 'cal_176Hf/177Hf(t)',
              'cal_epsilon_Hf(t)', 'cal_lu/hf', 'cal_Hf_TDM1(Ma)', 'cal_Hf_TDM2(Ma)']

trace_elements_cols = ['Nb', 'La', 'Ce', 'Pr', 'Nd', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho',
                       'Er', 'Tm', 'Yb', 'Lu', 'Hf',
                       'Ta', 'Pb', 'Th', 'U', 'Th/U', 'Hg', 'P', 'Y', 'Sc', 'Ti', 'Zr']

show_data_cols = age_cols + lu_hf_cols + trace_elements_cols
download_col = show_data_cols

download_data_rename_dict = {'st_sample': 'sample',
                             'sample_spots_id': 'spot',
                             'cal_epsilon_Hf(t)': 'epsilon_Hf(t)',
                             'cal_lu/hf': 'f(Lu/Hf)',
                             'cal_Hf_TDM1(Ma)': 'Hf_TDM1(Ma)',
                             'cal_Hf_TDM2(Ma)': 'Hf_TDM2(Ma)',
                             'latitude': 'Latitude',
                             'cal_176Hf/177Hf(t)': '176Hf/177Hf(t)'}

mda_col = 'MDA'
lat_col = 'latitude'
lon_col = 'longitude'

sample_info_cols = ['sample', 'unique_sample_id', 'formation', 'group',
                    'rock_type', 'simplified_lithology', 'lithology', mda_col, 'MDA_1se', 'MDA_MSWD', 'MDA_type',
                    'Method',
                    'TDA', 'TDA_1se', 'TDA_MSWD', 'TDA_type', 'reference', 'DOI']

# load data
sample_data, detrital_data = load_data(sample_file, data_file)

# set the basic information
st.title('East Asian Detrital Zircon Database')
st.sidebar.header('User Input Features')
st.sidebar.subheader('Depositional Age and Location Filter ')

# init st.session_state
if 'map_num' not in st.session_state:
    st.session_state['map_num'] = 0

with st.form("st_sample_selected"):
    with st.sidebar:
        data_subset_name = st.text_input('Please define the name of this data subset , '
                                         'which will be used as suffix of download file name')
        data_subset_name = str(data_subset_name.strip())
        tag = st.selectbox('Pick a data set to show in the map',
                           ['total', 'selected'],
                           key='map'
                           )

        mda_range = st.slider('Define depositional age span',
                              0, 2000,
                              (0, 200),
                              key='mda_range'
                              )

        lat_range = st.slider('Define latitude span',
                              10, 65,
                              (25, 35),
                              key='lat_range'
                              )

        lon_range = st.slider('Define longitude span',
                              60, 150,
                              (105, 115),
                              key='lon_range'
                              )
        user_geo = st.file_uploader("Upload a GeoJSON file to select samples within the polygon",
                                    type=["geojson", "kml", "zip"],
                                    )
        st.write('Steps for select samples with polygon'
                 ': Draw a polygon on the map -> Export it as a GeoJSON -> Upload it back to the app')
        submitted = st.form_submit_button("Submit for calculation")
        st.markdown("""---""")

selected_sample = sample_data[(sample_data[mda_col] > mda_range[0]) &
                              (sample_data[mda_col] < mda_range[1]) &
                              (sample_data[lat_col] > float(lat_range[0])) &
                              (sample_data[lat_col] < float(lat_range[1])) &
                              (sample_data[lon_col] > float(lon_range[0])) &
                              (sample_data[lon_col] < float(lon_range[1]))]

# for spatial selection
sample_gdf = gpd.GeoDataFrame(sample_data,
                              geometry=gpd.points_from_xy(sample_data['longitude'], sample_data['latitude'])
                              )
sample_gdf.crs = {"init": 'epsg:4326'}

# for user upload geojoson
if user_geo:
    user_file, user_gpd = uploaded_file_to_gdf(user_geo)
    data_within = sample_gdf[sample_gdf.geometry.within(user_gpd.geometry.unary_union)]
    selected_sample = sample_data.loc[sample_data['unique_sample_id'].isin(data_within['unique_sample_id'])]
    # st.dataframe(data_within['unique_sample_id'])
else:
    user_file = None

# for stratigraphic information selection
st.sidebar.subheader('Stratigraphic information Filter')

with st.form("sti_sample_selected"):
    with st.sidebar:
        formation = st.multiselect('Define the formations',
                                   list(selected_sample['formation'].unique()),
                                   key='formation'
                                   )
        group = st.multiselect('Define the groups',
                               list(selected_sample['group'].unique()),
                               key='group'
                               )
        submitted = st.form_submit_button("Submit for calculation")
        st.markdown("""---""")

if len(formation) >= 1 or len(group) >= 1:
    selected_sample = selected_sample[(selected_sample['formation'].isin(formation)) |
                                      (selected_sample['group'].isin(formation))
                                      ]

space_time_selected_samples = selected_sample['unique_sample_id']
st_selected_data = detrital_data[detrital_data['unique_sample_id'].isin(space_time_selected_samples)]

st.sidebar.subheader('User Upload Data')

with st.sidebar:
    with st.expander('Click for more options'):
        with st.form("user_upload_data"):
            user_data_tag = st.selectbox('Whether use uploading data?',
                                         [True, False],
                                         key='udt'
                                         )

            database_data_tag = st.selectbox('Whether use database?',
                                             [True, False],
                                             key='ddt'
                                             )

            database_contribution_tag = st.selectbox('Whether contribute your data to the database?',
                                                     [True, False],
                                                     key='dct'
                                                     )

            st.markdown("Prepare data file according to this template {}"
                        .format(prepare_template_file('data/template.xlsx')),
                        unsafe_allow_html=True)
            user_data = st.file_uploader("Upload excel file",
                                         type=["xlsx"],
                                         key='user_upload_data'
                                         )
            submitted = st.form_submit_button("Submit for calculation")
        st.markdown("""---""")

if user_data and user_data_tag:
    if 'user_sample' not in st.session_state:
        st.session_state['user_sample'] = 0
    if 'user_data' not in st.session_state:
        st.session_state['user_data'] = 0
    if submitted or st.session_state['user_sample'] == 0:
        user_upload_sample, user_upload_data = read_user_file(user_data,
                                                              selected_sample['unique_sample_id'].to_list())
    else:
        user_upload_sample = st.session_state['user_sample']
        user_upload_data = st.session_state['user_data']

    if not user_upload_sample.empty:
        st.write('combine')
        selected_sample = pd.concat([user_upload_sample, selected_sample],
                                    ignore_index=True,
                                    axis=0
                                    )
        # selected_sample
    if not user_upload_data.empty:
        st_selected_data = pd.concat([user_upload_data, st_selected_data],
                                     ignore_index=True,
                                     axis=0
                                     )
    if not database_data_tag:
        user_upload_sample = user_upload_sample.reindex(selected_sample.columns, axis='columns', fill_value=None)
        user_upload_data = user_upload_data.reindex(st_selected_data.columns, axis='columns', fill_value=None)
        selected_sample = user_upload_sample
        st_selected_data = user_upload_data

if tag == 'total':
    data = sample_data
else:
    data = selected_sample

# convert 
data['MDA(Ma)'] = data['MDA'].astype(str) + " ± " + data['MDA_1se'].astype(str) + "Ma"

# add leafmap
add_leafmap(data,
            lat_col,
            lon_col,
            center=[36, 120],
            zoom=4,
            user_geo=user_file)

data_c1, data_c2 = st.columns([1, 1])
# st.session_state
with data_c1:
    st.write("### Selected samples information")

    show_data_samples = selected_sample[sample_info_cols]
    with st.expander("Click to see samples information table"):
        st.dataframe(show_data_samples)

        st.write("Now, **{}** samples been selected, and zircon age measurements are **{}**, "
                 "zircon trace elements measurements are **{}**, zircon Lu-Hf isotopic measurements are  **{}**."
                 .format(len(space_time_selected_samples),
                         len(st_selected_data['cal_best_age']),
                         len(st_selected_data['Nd'].dropna()),
                         len(st_selected_data['cal_epsilon_Hf(t)'].dropna())))

# show the detrital zircon data and do some analysis

with data_c2:
    st.write("### Measurements data")

    # show the selected data
    show_dz_data = st_selected_data[download_col]
    show_dz_data.rename(columns=download_data_rename_dict, inplace=True)

    with st.expander("Click to see measurements data table"):
        st.dataframe(show_dz_data)

    download_data_as_csv(st_selected_data,
                         selected_sample,
                         download_col,
                         download_data_rename_dict,
                         file_name='Data From East Asia Detrital Zircon Database of {}.csv'.format(data_subset_name)
                         )

st.write("### Visualizations and Applications")

st.sidebar.header('Visualization analysis')
st.sidebar.subheader('Zircon ages visualization')
with st.sidebar:
    with st.form("age_dist"):
        sampleList = st.multiselect('select the samples for analysis, recommend less than 20',
                                    list(set(st_selected_data['unique_sample_id'])),
                                    list(set(st_selected_data['unique_sample_id']))[0:7]
                                    )
        submitted = st.form_submit_button("Submit")

with st.expander("Click to see the information of samples selected in the left sidebar"):
    # st.dataframe(show_dz_data)
    input_dp_samples, input_dp_data = prepare_visualization_data(show_data_samples,
                                                                 sampleList,
                                                                 st_selected_data,
                                                                 selected_sample,
                                                                 data_subset_name
                                                                 )

main_df, main_byid_df, samples_df, analyses_df, dp_samples, dp_zircon_data = \
    load_data_into_detritalpy(input_dp_samples,
                              input_dp_data,
                              sample_id='unique_sample_id')

plot_age_distribution(sampleList,
                      main_byid_df,
                      data_subset_name,
                      key_prefix='rc',
                      sample_submit=submitted
                      )

plot_maximum_depositional_age(sampleList,
                              main_byid_df,
                              data_subset_name,
                              key_prefix='rc',
                              sample_submit=submitted
                              )

plot_multi_dimensional_scaling(sampleList,
                               main_byid_df,
                               data_subset_name,
                               key_prefix='rc',
                               sample_submit=submitted
                               )

st.sidebar.subheader('Lu-Hf isotope visualization')
lu_hf_list, lu_data = prepare_data_for_lu_hf_ploter(st_selected_data)

with st.sidebar:
    with st.form("age_lu_hf"):
        plot_lu_hf_sample_list = st.multiselect('select the samples for lu-hf isotopic visualization,'
                                                ' recommend less than 20',
                                                lu_hf_list,
                                                lu_hf_list[:12],
                                                key='lu_hf'
                                                )
        submitted = st.form_submit_button("Submit")

plot_violins(plot_lu_hf_sample_list,
             lu_data,
             data_subset_name,
             sample_id='unique_sample_id',
             x_col='BestAge',
             y_col='epsilon_Hf(t)',
             key_prefix='lu-hf',
             sample_submit=submitted
             )

plot_kde_lu_hf(plot_lu_hf_sample_list,
               lu_data,
               data_subset_name,
               sample_id='unique_sample_id',
               x_col='BestAge',
               y_col='epsilon_Hf(t)',
               key_prefix='lu-hf',
               sample_submit=submitted
               )

st_key = 'custom_group'
st.markdown("""---""")
st.write("### Visualization based on custom sample groups")
st.sidebar.subheader('For custom sample groups visualization')

# get the user define  group number


with st.sidebar:
    with st.expander('Click to see the options'):
        with st.form("customs_group"):
            group_num = st.slider('Specify the groups number for samples grouping', 2, 20, 6, key='sdg1')
            group_dict = dict()
            group_dict = defaultdict(dict)
            user_define_col_name = 'custom_group_name'
            user_define_group_samples_list = list()
            custom_sample_list = list()

            for i in range(group_num):
                group_dict['sub_{}'.format(i)]['name'] = str(
                    st.text_input('Define the name of group {}'.format(i + 1),
                                  'group_{}'.format(i + 1),
                                  key=st_key + 'gn{}'.format(i + 1)
                                  )).strip()

                group_dict['sub_{}'.format(i)]['sample_list'] = st.multiselect(
                    'Select the samples to label as group {}'.format(i + 1),
                    list(set(st_selected_data['unique_sample_id'])),
                    list(set(st_selected_data['unique_sample_id']))[5 * i:5 * (i + 1)],
                    key=st_key + 'g{}'.format(i + 1)
                )
                # add label to detrital zircon data
                st_selected_data.loc[st_selected_data['unique_sample_id'].isin(
                    group_dict['sub_{}'.format(i)]['sample_list']),
                                     user_define_col_name] = group_dict['sub_{}'.format(i)]['name']

                selected_sample.loc[selected_sample['unique_sample_id'].isin(
                    group_dict['sub_{}'.format(i)]['sample_list']),
                                    user_define_col_name] = group_dict['sub_{}'.format(i)]['name']

                user_define_group_samples_list.extend(group_dict['sub_{}'.format(i)]['sample_list'])
                custom_sample_list.append(group_dict['sub_{}'.format(i)]['name'])

            submitted = st.form_submit_button("Submit")
    st.markdown("""---""")

with st.expander("Click to see the information of samples selected in the left sidebar"):
    sample_info_cols.append(user_define_col_name)
    sample_info = selected_sample[sample_info_cols]
    custom_samples, custom_data = prepare_visualization_data(sample_info,
                                                             user_define_group_samples_list,
                                                             st_selected_data,
                                                             selected_sample,
                                                             data_subset_name,
                                                             key_prefix='Custom',
                                                             )

custom_samples = custom_samples.groupby(user_define_col_name)['longitude', 'latitude'].mean().reset_index()

main_df, main_byid_df, samples_df, analyses_df, dp_samples, dp_zircon_data = \
    load_data_into_detritalpy(custom_samples,
                              custom_data,
                              sample_id=user_define_col_name)

plot_age_distribution(custom_sample_list,
                      main_byid_df,
                      data_subset_name,
                      key_prefix='custom',
                      sample_submit=submitted
                      )

plot_multi_dimensional_scaling(custom_sample_list,
                               main_byid_df,
                               data_subset_name,
                               key_prefix='custom',
                               sample_submit=submitted
                               )

lu_hf_list, lu_data = prepare_data_for_lu_hf_ploter(custom_data,
                                                    sample_id=user_define_col_name,
                                                    custom_col=user_define_col_name
                                                    )
# lu_hf_list = custom_samples[user_define_col_name].to_list()

plot_violins(lu_hf_list,
             lu_data,
             data_subset_name,
             sample_id=user_define_col_name,
             x_col='BestAge',
             y_col='epsilon_Hf(t)',
             key_prefix='lu-hf-custom',
             sample_submit=submitted
             )

plot_kde_lu_hf(lu_hf_list,
               lu_data,
               data_subset_name,
               sample_id=user_define_col_name,
               x_col='BestAge',
               y_col='epsilon_Hf(t)',
               key_prefix='lu-hf-custom',
               sample_submit=submitted
               )

# Connection with Gplates
st.sidebar.header('Paleo-position Reconstruction with Gplates')
with st.sidebar:
    with st.form("gplates"):
        selected_model = st.selectbox('Pick up Gplates Plate reconstruction model',
                                      ('Merdith_2022', 'Young_2018'),
                                      key='gpl_model'
                                      )
        gpl_samples_list = st.multiselect('Select the samples for analysis, recommend less than 20',
                                          list(set(st_selected_data['unique_sample_id'])),
                                          list(set(st_selected_data['unique_sample_id']))[0:10],
                                          key='gpl'
                                          )
        submitted_gpl = st.form_submit_button('Submit')
# load Gplates reconstruction model


if 'rt_gplates' not in st.session_state:
    st.session_state['rt_gplates'] = 0
if 'rt_gplates' not in st.session_state:
    st.session_state['sp_gplates'] = 0
if 'rt_gplates' not in st.session_state:
    st.session_state['cs_gplates'] = 0
if 'rt_gplates' not in st.session_state:
    st.session_state['pb_gplates'] = 0

if 'model_name' not in st.session_state:
    st.session_state['model_name'] = 0

if selected_model == st.session_state['model_name']:
    rotation_model = st.session_state['rt_gplates']
    static_polygons = st.session_state['sp_gplates']
    coastlines = st.session_state['cs_gplates']
    plate_boundaries = st.session_state['pb_gplates']
else:
    rotation_model, static_polygons, coastlines, plate_boundaries = load_gplates_model(selected_model)

# Begin plotting
# for rpp


gpl_samples = selected_sample.loc[selected_sample['unique_sample_id'].isin(gpl_samples_list)]
gpl_point_att = ['sample', 'unique_sample_id', 'formation', 'group']
point_features = data_to_gplates(gpl_samples,
                                 gpl_point_att,
                                 static_polygons,
                                 rotation_model
                                 )
plot_paleo_lat_lon(rotation_model, point_features, gpl_point_att)

# create the total samples Gplates features
if 'total_sample_feature' not in st.session_state:
    st.session_state['total_sample_feature'] = 0
if st.session_state['total_sample_feature'] == 0:
    st.session_state['total_sample_feature'] = data_to_gplates(sample_data,
                                                               gpl_point_att,
                                                               static_polygons,
                                                               rotation_model
                                                               )

sequence_reconstruction_plot(coastlines,
                             rotation_model,
                             point_features,
                             gpl_point_att,
                             )
