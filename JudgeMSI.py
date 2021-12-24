#!/usr/bin/env python
import sys 
import os  
from os.path import basename, dirname, exists, join
from glob import glob 
import pandas as pd 
from collections import OrderedDict
import copy
import numpy as np
from datetime import datetime

def readTagDepthFromNGS(path, peak_threshold=0, threshold_as_relative=True, peak_smoothing=False):
    '''
    input: info under MSI folder. 
    '''
#     print(peak_threshold)
    sample_name = basename(dirname(dirname(path)))
    df_tags = pd.read_csv(path, sep='\t', header=None)
    df_tags.columns = ['tag','length','depth']
    dfs_for_each_tag = []
    for tag,df_this_tag in df_tags.groupby('tag'):
        df_this_tag = df_this_tag.sort_values('length')
        df_this_tag.reset_index(inplace=True, drop=True) 
        peak_indices = findPeakLocation(df_this_tag['depth'], peak_threshold=peak_threshold, \
                                        threshold_as_relative=threshold_as_relative, peak_smoothing=peak_smoothing) 
        df_this_tag.loc[:,'is_peak'] = df_this_tag.index.map(lambda x:1 if x in peak_indices else 0)
        dfs_for_each_tag.append(df_this_tag)
    return dfs_for_each_tag
        
def findPeakLocation(list_in, peak_threshold=0, threshold_as_relative=True, peak_smoothing=0):
    '''
    Given a list, output the index of the peaks. 
    Definition of peak:
        larger than (or equal) both left and the right.
    
    peak_threshold: 
        Only those peaks larger than or equal to peak_threshold will be defined as peak. 
    
    If threshold_as_relative:
        use peak_threshold * average depth of the list_in. 
    
    peak_smoothing:
        Int. How many numbers will be used as the left flank or right flank to smooth.  
    '''
    if threshold_as_relative:
        peak_threshold = sum(list_in) / len(list_in) * peak_threshold
    
#     print('using {} as peak_threshold'.format(peak_threshold))
    
    list_in = [0] + list(list_in) + [0] 
    # Merge numbers that are same and sequent. 
    df_ordered_unique = pd.DataFrame(columns=['depth','same'])    
     
    list_smooth = list_in[:peak_smoothing]
    for i in range(peak_smoothing,len(list_in)-peak_smoothing):
        smoothed = sum(list_in[i-peak_smoothing:i+peak_smoothing+1]) / (peak_smoothing * 2 +1)
        list_smooth.append(smoothed)
    list_smooth.extend(list_in[-peak_smoothing:])
    list_in = list_smooth
    
#     if peak_smoothing: 
#         list_smooth = [0]
#         for i in range(1,len(list_in)-1):
#             smoothed = sum(list_in[i-1:i+2]) / 3
#             list_smooth.append(smoothed)
#         list_smooth.append(0)
#         list_in = list_smooth
        
    
#     print(list_in)
    for i in range(0,len(list_in)):
        if i==0 or list_in[i] != list_in[i-1]: 
            df_ordered_unique.loc[i,'depth'] = list_in[i]
            df_ordered_unique.loc[i,'same'] = 1
        elif list_in[i] == list_in[i-1]:
            df_ordered_unique.loc[df_ordered_unique.index[-1],'same'] += 1 
    
#     print(df_ordered_unique) 
        
    list_in = list(df_ordered_unique['depth'])
    # It will shift
    peak_indices = []
    for i in range(1,len(list_in)-1):
        if list_in[i] < peak_threshold:
            continue
        
        # If its in the left end or the right end, more strict. 
        if i in [1, len(list_in)-2]:
            if list_in[i] > list_in[i-1] and list_in[i] > list_in[i+1]: 
                peak_indices.append(i) 
        else:
            if list_in[i] >= list_in[i-1] and list_in[i] >= list_in[i+1]: 
                peak_indices.append(i)
    
    df_ordered_unique.loc[:,'is_peak'] = 0
#     print(peak_indices)
    df_ordered_unique.iloc[peak_indices,2] = 1
#     print(df_ordered_unique)
    peak_indices = []
    
    df_ordered_unique.loc[:,'peak_indices'] = df_ordered_unique.apply(lambda row: row.name + int(row['same'] / 2), axis=1)
    df_peaks = df_ordered_unique[df_ordered_unique['is_peak'] == 1]
    peak_indices = list(df_peaks.loc[:,'peak_indices'])
    peak_indices = [i-1 for i in peak_indices]
    
#     print(peak_indices)
    return peak_indices

def determineStabilility(path_ngs_info_F, path_ngs_info_B, \
                         diff_threshold_for_one_tag=2, peak_threshold=0.2, threshold_as_relative=True, \
                         tags_to_use=['Bat25','Bat26','Mono27','NR21','NR24','NR27'], show_details=False, peak_smoothing=1,
                         raise_exception_if_tag_missing=False):
    dfs_ngs_F = readTagDepthFromNGS(path_ngs_info_F, peak_threshold=peak_threshold, threshold_as_relative=threshold_as_relative, peak_smoothing=peak_smoothing)
    dfs_ngs_B = readTagDepthFromNGS(path_ngs_info_B, peak_threshold=peak_threshold, threshold_as_relative=threshold_as_relative, peak_smoothing=peak_smoothing)
    
    unstable_tags = []
    
    dfs_ngs_F_in_use = []
    dfs_ngs_B_in_use = []
    
    number_tags_should_use = len(tags_to_use)
    tags_available_for_F = map(lambda df:df['tag'][0], dfs_ngs_F)
    tags_available_for_B = map(lambda df:df['tag'][0], dfs_ngs_F)
    tags_to_use = sorted(set(tags_to_use) & set(tags_available_for_F) & set(tags_available_for_B))

    if raise_exception_if_tag_missing and len(tags_to_use) < number_tags_should_use:
        raise Exception('too less tags available')
    
    if show_details:
        print(tags_to_use)
    for df in dfs_ngs_F:
        tag = df['tag'][0]
        if tag in tags_to_use:
            dfs_ngs_F_in_use.append(df)
    for df in dfs_ngs_B:
        tag = df['tag'][0]
        if tag in tags_to_use:
            dfs_ngs_B_in_use.append(df)
    
    for df_tag_F, df_tag_B in zip(dfs_ngs_F_in_use,dfs_ngs_B_in_use):
        tag = df_tag_F['tag'][0] 
        df_tag_F_peak = df_tag_F[df_tag_F['is_peak']==1]
        df_tag_B_peak = df_tag_B[df_tag_B['is_peak']==1]
        if show_details:
            print('same') 
            print('F') 
            print(df_tag_F)
            print('B') 
            print(df_tag_B)
        
        if df_tag_F_peak.shape[0] > df_tag_B_peak.shape[0]:  # != or > 
            unstable_tags.append(tag)
        else:
            if df_tag_F_peak.shape[0] < df_tag_B_peak.shape[0]:
                # Ignore the minor peaks for B. 
                F_num_for_this_tag = df_tag_B_peak.shape[0] - df_tag_F_peak.shape[0]
                print('dropping due to B is more than F')
                df_tag_B_peak.sort_values('depth', inplace=True, ascending=False)
                df_tag_B_peak = df_tag_B_peak.iloc[:F_num_for_this_tag, :]
                df_tag_B_peak.sort_values('length', inplace=True, ascending=True)
            
            abs_gap_between = (df_tag_F_peak['length'].reset_index(drop=True) - df_tag_B_peak['length'].reset_index(drop=True)).map(abs).max()
#             print(abs_gap_between)
            if abs_gap_between >= diff_threshold_for_one_tag:
                unstable_tags.append(tag)
    print(unstable_tags)
    return unstable_tags, tags_to_use
    
if __name__ == '__main__':
    print('''
    ''')
    TAGS_SHOULD_USE = ['Bat25','Bat26','Mono27','NR21','NR24','NR27']
    path_ngs_info_F, path_ngs_info_B, path_csv = sys.argv[1:4]
    unstable_tags, tags_to_use = determineStabilility(path_ngs_info_F, path_ngs_info_B, diff_threshold_for_one_tag=2, peak_threshold=0.2, threshold_as_relative=True, \
                         tags_to_use=TAGS_SHOULD_USE, show_details=False, peak_smoothing=1)
    str_tumor_normal_with_ids = '{}_{}'.format(path_ngs_info_F, path_ngs_info_B).replace('/','')
    
    df_unstable_tag_for_samples = pd.DataFrame()
    df_unstable_tag_for_samples.loc[str_tumor_normal_with_ids,'tags_to_use'] = ','.join(tags_to_use)
    df_unstable_tag_for_samples.loc[str_tumor_normal_with_ids,'unstable'] = ','.join(unstable_tags)

    df_unstable_tag_for_samples.loc[:,'num_unstable'] = len(unstable_tags)
    df_unstable_tag_for_samples.loc[:,'result'] = df_unstable_tag_for_samples['num_unstable'].map(lambda x:'MSI_H' if x > 1 else 'MSI_L')

    if len(tags_to_use) < len(TAGS_SHOULD_USE):
        df_unstable_tag_for_samples.loc[str_tumor_normal_with_ids,'result'] += '_SomeMarkersMissing_AskWangBo'
    df_unstable_tag_for_samples.loc[str_tumor_normal_with_ids, 'bin_for_msi_version'] = datetime.fromtimestamp(os.path.getmtime(sys.argv[0])).strftime('%Y%m%d_%H%M%S')
    
    df_unstable_tag_for_samples.index.name = 'sample_pair'
    # df_unstable_tag_for_samples.loc[:,'sample'] = df_unstable_tag_for_samples.index

    df_unstable_tag_for_samples.to_csv(path_csv)
