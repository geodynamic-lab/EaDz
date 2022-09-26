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
import math
import time
from pathlib import Path
import glob
from collections import defaultdict
import sys
import shutil

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import numpy as np
import scipy
from scipy.signal import savgol_filter
import xlsxwriter


def map_col_name(obs, map_dict, file_name=None, no_record_col=set(), task_tag='mapping'):
    error_records = dict()
    error_records = defaultdict(list)
    no_use_list = list()

    col_index = 0
    for item in obs.iloc[0, :]:
        new_item = item
        # print(item)
        for key in map_dict.keys():
            if item in map_dict[key]:
                new_item = key
                break
        if new_item in map_dict.keys():
            # rename for the no use columns
            if new_item == "NoUse":
                no_use_list.append(item)
                #     # new_item = "NoUse_" + item
                new_item = item

            obs.iloc[0, col_index] = new_item
            # pass
            # print('we cannot find the standard name')
        else:
            if not task_tag == 'prefix':
                print("No Find %s of %s" % (new_item, file_name))
            #             print(len(new_item))
            #             print(file_name)
            if len(file_name) > 3:
                error_records['file_name'].append(str(file_name))
                error_records['col_name'].append(str(new_item))
            # error_records.append({'file_name':file_name,'col_name':new_item},ignore_index = True)
            #             print(error_records)
            no_record_col.update([new_item])
            obs.iloc[0, col_index] = new_item
        col_index += 1

    return obs, no_record_col, error_records, no_use_list


def read_the_columns_map(map_file_name, task_tag):
    if task_tag == "prefix":
        sheet_name = 'Sheet2'
    else:
        sheet_name = 'Sheet1'
    # col_map = pd.read_excel(map_file_name,index=None,header=None,sheet_name=sheet_name)
    col_map = pd.read_excel(map_file_name, index_col=None, header=None, sheet_name=sheet_name)
    # define the col map dict
    col_map_dict = dict()
    col_map_dict = defaultdict(list)
    n_row = col_map.shape[0]
    n_col = col_map.shape[1]
    for row in range(n_row):
        row = list(pd.Series(col_map.iloc[row, :]).dropna())
        key = row[0]
        #         mylist = list(dict.fromkeys(mylist))
        item = list(dict.fromkeys(row))
        col_map_dict[key] = item

    return col_map_dict


def add_prefix_for_TDM(obs, file_name=None):
    """
    :param obs:  panda dataframe
    :param file_name: the filename
    :return: the dataframe with the columns added prefix
    """
    col_index = 0
    for item in obs.iloc[0, :]:
        # print(item)
        if str(item) in ['TDM1', 'TDM1_1se', 'TDM1_2se', 'TDM2', 'TDM2_1se', 'TDM2_2se']:

            for i in range(4):
                check_row = abs(col_index - i)

                if "ND" in str(obs.iloc[0, check_row]).upper():
                    prefix = 'Nd'
                    break
                elif "HF" in str(obs.iloc[0, check_row]).upper():
                    prefix = "Hf"
                    break
                elif "Sr" in str(obs.iloc[0, check_row]).upper():
                    prefix = "Sr"
                    break
                else:
                    #                     pass
                    prefix = "Not_find"

            #             print("get the prefix for TDM is %s" % prefix)
            if prefix == "Not_find":
                print("{0}, not find in the file {1} during adding prefix".format(item, file_name))
                prefix = ""

            obs.iloc[0, col_index] = "%s_%s" % (prefix, obs.iloc[0, col_index])

        else:
            pass

        col_index += 1
    return obs


def add_prefix_for_1se(obs):
    """
    add the prefix for the columns include 1se and 2se
    :param obs:
    :return:
    """
    col_index = 0
    check_status = 'no'
    for item in obs.iloc[0, :]:
        if item.lower().strip() in ['1se', '2se']:

            prefix = obs.iloc[0, abs(col_index - 1)]
            print(prefix)
            obs.iloc[0, col_index] = "%s_%s" % (prefix, obs.iloc[0, col_index])
            check_status = 'yes'
        else:
            check_status = 'no'

        col_index += 1

    return obs, check_status


def traverse_the_data_dir(read_dir, write_dir, back_up_dir, col_dict_file_name, task_tag, task_name=None,
                          cleaned_data_dir=None):
    """

    :param cleaned_data_dir:
    :param task_name: task identified information
    :param read_dir:
    :param write_dir:
    :param back_up_dir:
    :param col_dict_file_name:
    :param task_tag: review, prefix, mapping
    :return:
    """
    no_record_col = set()
    duplicated_item = set()
    error_records = dict()
    error_records = defaultdict(list)
    col_map_dict = read_the_columns_map(col_dict_file_name, task_tag)
    data_check = defaultdict(dict)

    # map the col_name
    for obs_file in glob.glob(read_dir + '/*.xlsx'):
        file_name = obs_file.split('\\')[-1]
        # print(file_name)
        # if os.path.exists(os.path.join(handled_mti_dir,file_name)):
        #     continue

        obs = pd.read_excel(obs_file, sheet_name='Sheet1', header=None)
        obs, return_record_col, return_error_records, drop_list = map_col_name(obs, col_map_dict, file_name,
                                                                               no_record_col=no_record_col,
                                                                               task_tag=task_tag)
        # print(obs)

        if task_tag == 'review':
            pass
            # print('there still are {} columns need update'.format(len(return_record_col)))

        elif task_tag == 'prefix':
            obs = add_prefix_for_TDM(obs, file_name)
            obs, check = add_prefix_for_1se(obs)
            obs.to_excel(os.path.join(write_dir, file_name), header=None, index=None)

        elif task_tag == 'mapping':
            obs.columns = obs.iloc[0]
            obs.to_excel(os.path.join(back_up_dir, file_name), header=None, index=None)

            obs.drop(drop_list, inplace=True, axis=1)

            # add  the sample_spots_id
            if 'spot' in obs.columns:
                obs['sample_spots_id'] = obs['sample'].astype(str) + '_' + obs['spot'].astype(str)
            # add record_num
            record_num = file_name.split('.')[0].split('-')[1]
            obs['record_num'] = str(record_num)
            obs['sample_id'] = obs['record_num'].astype(str) + '_' + obs['sample'].astype(str)

            # deplete the row one
            obs.drop([0], inplace=True)

            obs.to_excel(os.path.join(write_dir, file_name), index=None)

            # data check
            check_list = ['spot', 'sample', 'age(Ma)']
            for i in check_list:
                if i in obs.columns:
                    data_check[file_name][i] = 'Yes'
                else:
                    data_check[file_name][i] = 'No'
            if len(obs.columns) == len(set(obs.columns)):
                data_check[file_name]['duplicated'] = 'No'
            else:
                data_check[file_name]['duplicated'] = 'Yes'
                data_check[file_name]['duplicated_item'] = set(
                    [x for x in list(obs.columns) if list(obs.columns).count(x) > 1])
                duplicated_item.update(data_check[file_name]['duplicated_item'])

            # copy data the single paper dir
            # determine the file type
            # age_data
            if file_name[0:2] in ['a_', 'a-']:
                paper_dir = os.path.join(age_cleaned_dir, record_num)
            # major_element for mti data
            elif file_name[0:2] in ['mt', 'm_', 't_', 'ti', 'm-', 't-']:
                paper_dir = os.path.join(major_element_dir, record_num)
            # tracer_elements for zircon data
            elif file_name[0:2] in ['ts']:
                paper_dir = os.path.join(tracer_element_dir, record_num)
            elif file_name[0:2] in ['i_', 'i-', 'is']:
                paper_dir = os.path.join(isotopic_dir, record_num)
            else:
                print('please check the filename: {}'.format(file_name))

            if not os.path.exists(paper_dir):
                os.mkdir(paper_dir)
            source = os.path.join(write_dir, file_name)
            target = os.path.join(paper_dir, file_name)
            shutil.copy(source, target)

            # write data to the total data dir
            paper_dir = os.path.join(total_dir, record_num)
            if not os.path.exists(paper_dir):
                os.mkdir(paper_dir)
            source = os.path.join(write_dir, file_name)
            target = os.path.join(paper_dir, file_name)
            shutil.copy(source, target)

        else:
            print("please check the task_tag information")

        # write the record information
        no_record_col.update(return_record_col)
        if len(return_error_records['file_name']) > 0:
            error_records['file_name'].extend(return_error_records['file_name'])
            error_records['col_name'].extend(return_error_records['col_name'])

    if not task_tag == 'prefix':
        pd.DataFrame(no_record_col).to_excel(
            os.path.join(work_dir, 'log_files', 'Error01-{0}_{1}-No_record_columns.xlsx'.format(task_name, task_tag)),
            header=None, index=None)
        pd.DataFrame(error_records).to_excel(
            os.path.join(work_dir, 'log_files', 'Error02-{0}_{1}-record_columns.xlsx'.format(task_name, task_tag)),
            index=None)

    if task_tag == 'mapping':
        pd.DataFrame.from_dict(data_check, orient='index').to_excel(
            os.path.join(work_dir, 'log_files', 'CheckFile-{0}_{1}.xlsx'.format(task_name, task_tag)))
        if len(duplicated_item) > 0:
            pd.DataFrame(duplicated_item).to_excel(
                os.path.join(work_dir, 'log_files', 'Error03-{0}_{1}-duplicated.xlsx'.format(task_name, task_tag)),
                header=None, index=None)

    return len(no_record_col)


# function for merged files
def merge_files_into_one_of_single_paper(source_dir, target_dir, file_prefix, index_col='sample', file_tag='age'):
    for sub_dir in os.listdir(source_dir):
        combine_df = pd.DataFrame()
        output_file = os.path.join(target_dir, '{}_merged_{}_{}.xlsx'.format(str(sub_dir), file_prefix, file_tag))

        last_file_prefix = "no"
        for obs_file in glob.glob(os.path.join(source_dir, sub_dir) + '/*.xlsx'):

            file_name = obs_file.split('\\')[-1]
            print(file_name)
            if last_file_prefix == "no":
                last_file_prefix = file_name.split('-')[0]
            obs = pd.read_excel(obs_file, sheet_name='Sheet1', index_col=index_col)

            if combine_df.empty:
                combine_df = obs

            # if the data type same, we use the concat function, otherwise use the merge function
            else:
                if last_file_prefix == file_name.split('-')[0]:
                    combine_df = pd.concat([combine_df, obs], join='outer', ignore_index=False)
                else:
                    combine_df = pd.merge(combine_df, obs, left_index=True, right_index=True, how='outer',
                                          suffixes=['', '_dup'])
            last_file_prefix = file_name.split('-')[0]

        combine_df = combine_df.reset_index()
        combine_df['paper_ref'] = str(sub_dir)
        combine_df['sample_paired_id'] = combine_df['paper_ref'].astype(str) + '_' + combine_df['sample'].astype(str)
        # print(combine_df)
        # combine_df['sample_id'] = combine_df['sample'] +  combine_df['paper_ref']
        combine_df.to_excel(output_file, index=None)


def merge_files_into_one_database(source_dir, target_dir, file_name, file_tag='age'):
    combine_df = pd.DataFrame()
    for obs_file in glob.glob(source_dir + '/*{}*.xlsx'.format(file_tag)):
        obs = pd.read_excel(obs_file, sheet_name='Sheet1', index_col='sample_paired_id')
        print(obs_file.split('/')[-1])
        if combine_df.empty:
            combine_df = obs
        else:
            combine_df = pd.concat([combine_df, obs], join='outer', ignore_index=False)

    # combine_df = combine_df.reset_index()
    # writer = pd.ExcelWriter(os.path.join(target_dir,file_name),engine='xlsxwriter')
    combine_df.to_excel(os.path.join(target_dir, file_name), index=None)
    # combine_df.to_excel(writer,index=False)
    return combine_df


## user config region

work_dir = Path("replace with your data directory")
os.chdir(str(work_dir))

mapping_file_name = 'update_column_map.xlsx'
original_dir_name = 'obs_data'
task_name = 'detrital_zircon'
suffix = 'v1'
new_dir_list = ['prefix_completion', 'column_mapping', 'drop_no_use']

# define the dir
cleaned_data_dir = os.path.join(work_dir, '0_cleaned_data_dir_{}'.format(task_name))
mapping_file = os.path.join(work_dir, mapping_file_name)
origin_data_dir = os.path.join(work_dir, original_dir_name)
prefix_completion_dir = os.path.join(work_dir, '{0}_{1}'.format(task_name, new_dir_list[0]))
column_mapping_dir = os.path.join(work_dir, '{0}_{1}'.format(task_name, new_dir_list[1]))
drop_no_use_dir = os.path.join(work_dir, '{0}_{1}'.format(task_name, new_dir_list[2]))
single_paper_merged_dir = os.path.join(work_dir, '02_single_paper_merged_{}'.format(task_name))
final_result_dir = os.path.join(work_dir, '00_database')

# mkdir
for subdir in new_dir_list:
    if not os.path.exists(os.path.join(work_dir, '{0}_{1}'.format(task_name, subdir))):
        os.mkdir(os.path.join(work_dir, '{0}_{1}'.format(task_name, subdir)))

if not os.path.exists(os.path.join(work_dir, 'log_files')):
    os.mkdir(os.path.join(work_dir, 'log_files'))
if not os.path.exists(cleaned_data_dir):
    os.mkdir(cleaned_data_dir)
if not os.path.exists(single_paper_merged_dir):
    os.mkdir(single_paper_merged_dir)
if not os.path.exists(final_result_dir):
    os.mkdir(final_result_dir)

age_cleaned_dir = os.path.join(cleaned_data_dir, 'age_data')
tracer_element_dir = os.path.join(cleaned_data_dir, 'tracer_element_data')
major_element_dir = os.path.join(cleaned_data_dir, 'major_element_data')
isotopic_dir = os.path.join(cleaned_data_dir, 'isotopic_data')
total_dir = os.path.join(cleaned_data_dir, 'total')

cleaned_data_list = [age_cleaned_dir, tracer_element_dir, major_element_dir, isotopic_dir, total_dir]
for subdir in cleaned_data_list:
    if not os.path.exists(subdir):
        os.mkdir(subdir)

# step 1,review the columns dict
status = traverse_the_data_dir(origin_data_dir,
                               prefix_completion_dir,
                               prefix_completion_dir,
                               mapping_file_name,
                               'review',
                               task_name=task_name)
if not status == 0:
    print("by now we not complete the mapping dict, there still {} columns name need process".format(status))
    sys.exit(0)
# step 2, add the prefix
traverse_the_data_dir(origin_data_dir,
                      prefix_completion_dir,
                      prefix_completion_dir,
                      mapping_file_name,
                      'prefix',
                      task_name=task_name)
# step 3, mapping
traverse_the_data_dir(prefix_completion_dir,
                      drop_no_use_dir,
                      column_mapping_dir,
                      mapping_file_name,
                      'mapping',
                      task_name=task_name,
                      cleaned_data_dir=cleaned_data_dir)
# sample_spots_id for zircon data, sample_id for whole rock data
database_name = '<replace_with_your_database_name>_{}_{}.xlsx'.format(task_name, suffix)
combine_df = merge_files_into_one_database(single_paper_merged_dir,
                                           final_result_dir,
                                           database_name,
                                           'total')
