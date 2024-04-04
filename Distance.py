import pandas as pd
import numpy as np
import os

def Extract_SVs(celltype):
    
    if celltype == 'K562':
        #data = pd.read_table('/lanec2_home/zhangx/svs/K562/s4_K562_GM12878_hg19.bed',encoding = 'utf-8' )
        data = pd.read_table('/lanec2_home/zhangx/svs/data/K562_sv.bed')
        data['type'] = [ _.replace(_,_.strip(' " ').strip(' “ ').strip(' ” '))  for _ in data['type']]
        region = []
        for i in range(data.shape[0]):
            r = str(data['chr1'].iloc[i]) + ':' + str(data['breakpoint1'].iloc[i]) + '-' + str(data['breakpoint2'].iloc[i])
            region.append(r)
        data['region'] = region
        #type = pd.read_csv('/lanec2_home/zhangx/svs/data/K562_svtype.csv',encoding='utf-8')
        #data['type'] = type['type']
        data['chr1'] = [_.replace(_,'chr'+_) for _ in data['chr1']]
        result = data
        #result = pd.merge(data,s4_cell[['region','type']],on = 'region',how = 'inner')
    else:
        path = '/lanec2_home/zhangx/svs/'
        s4 = pd.read_csv(path + 'sciadv_s4.csv')
        s4 = s4.rename(columns={'Unnamed: 0':'cell', 'chrom1':'chr1', 'chrom2':'chr2'})
        s4= s4[s4['chr1'] == s4['chr2']] # 删除inter-
        s4_cell = s4[s4['cell'] == celltype]
        s4 = s4.rename(columns={'Unnamed: 0':'cell', 'chrom1':'chr1', 'chrom2':'chr2'})
        s4= s4[s4['chr1'] == s4['chr2']] # 删除inter-
        s4_cell = s4[s4['cell'] == celltype]
        #chr = list(set(s4_cell['chr1']))
        s4_cell = s4_cell[s4_cell['chr1'] != 'chr9'] # 删除9号染色体
        # add type to dataframe
        region = []
        for i in range(s4_cell.shape[0]):
            r = str(s4_cell['chr1'].iloc[i][3:]) + ':' + str(s4_cell['breakpoint1'].iloc[i]) + '-' + str(s4_cell['breakpoint2'].iloc[i])
            region.append(r)
        s4_cell['region'] = region

        choose = ['++','+-','-+','--']
        type = []
        for i in range(s4_cell.shape[0]):
            for k in range(5,9):
                if s4_cell.iloc[i,k] > 0.5:
                    type.append(choose[k-5])
        s4_cell['type'] = type
        s4_cell = s4_cell[['cell','chr1','breakpoint1','breakpoint2','region','type']]
        s4_cell = s4_cell.reset_index(drop = True)
        result = s4_cell
        
        #result.loc[result['svtype'] == '++','SVtype'] = 'dup33'
        #result.loc[result['svtype'] == '+-','SVtype'] = 'deletion'
        #result.loc[result['svtype'] == '-+','SVtype'] = 'duplication'
        #result.loc[result['svtype'] == '--','SVtype'] = 'dup55'

    return result


def FindRelatedTADS(chrom,startpoint,endpoint,tads_dataset):
    data = tads_dataset
    test = data[ ((data['chr'] == chrom) & (data['start'] > startpoint) & (data['start'] < endpoint)) 
        |  ((data['chr'] == chrom) & (data['end'] > startpoint) & (data['end'] < endpoint))
        |((data['chr'] == chrom) & (data['start'] < startpoint) & (data['end'] > endpoint))]
    return test


def CreateExcel(condition1,condition2,method_change):
    tads = pd.read_table(f'/lanec2_home/zhangx/DiffCompare/result/{condition1}_vs_{condition2}_10000_f1b5_BH.txt')
    tads_change = pd.read_table(f'/lanec2_home/zhangx/DiffCompare/result/{condition1}_vs_{condition2}_10000_{method_change}_f1b5_BH.txt')
    sv = Extract_SVs(condition2)
    sv_region = sv['region']
    #print(sv.shape[0])
    related_regions = []
    #Pvalue = []
    for _ in sv['region']:
        #print(_)
        chr = _.split(':')[0]
        start = int(_.split('-')[0].split(':')[1])
        end = int(_.split('-')[1])
        related_tads = FindRelatedTADS(chrom = chr,startpoint=start,endpoint=end,tads_dataset=tads)

        if not related_tads.shape[0] == 0:
            rg = ';'.join(related_tads['region'])
            #bh = ';'.join(related_tads['BH'])
            related_regions.append(rg)
            #Pvalue.append(bh)
        else:
            related_regions.append('None')
            #Pvalue.append(None)
    sv['related_tads'] = related_regions

    # 将 relate_tads 列按分号拆分，并在原数据框中添加多行信息
    sv['related_tads'] = sv['related_tads'].str.split(';')
    expanded_data = []
    for idx, row in sv.iterrows():
        for rt in row['related_tads']:
            expanded_data.append([row['chr1'], row['breakpoint1'], row['breakpoint2'],row['region'],row['type'] ,rt])

    expanded_df = pd.DataFrame(expanded_data, columns=['chr1', 'breakpoint1', 'breakpoint2', 'region','type','related_tads'])
    #expanded_df['region'] = expanded_df['region'].where(~expanded_df['region'].duplicated(), '') # 删除重复的region避免数据过于冗余

    # 加入原始方法的pvalue，为pvalue_raw
    tads = tads.rename(columns={'region':'related_tads'})
    merge_df = pd.merge(expanded_df,tads[['related_tads','BH']],on = 'related_tads',how = 'inner')
    #merge_df = expanded_df.join(tads.set_index('region'),on = 'region',lsuffix = '_left',rsuffix = '_right')
    merge_df = merge_df.rename(columns={'BH':'pvalue_raw'})

    # 加入插补处理后的pvalue，为pvalue_to1
    tads_change = tads_change.rename(columns={'region':'related_tads'})
    merge_df = pd.merge(merge_df,tads_change[['related_tads','BH']],on = 'related_tads',how = 'inner')
    merge_df = merge_df.rename(columns={'BH':'pvalue_to1'})

    for r in sv_region:
        pv_raw = merge_df[merge_df['region'] == r]['pvalue_raw']
        pv_to1 =  merge_df[merge_df['region'] == r]['pvalue_to1']
        if np.any(pv_raw <= 0.05):
            merge_df.loc[merge_df['region'] == r,'retads_related_raw'] = 'yes'
        if np.any(pv_to1 <= 0.05):
            merge_df.loc[merge_df['region'] == r,'retads_related_to1'] = 'yes'
        else:
            merge_df.loc[merge_df['region'] == r,'retads_related_raw'] = 'no'
            merge_df.loc[merge_df['region'] == r,'retads_related_to1'] = 'no'
            
    merge_df = merge_df.sort_values(by = 'region')
    # 删除重复的region避免数据过于冗余
    #index_to_delete = merge_df[merge_df.duplicated(subset='region')].index
    #columns_to_delete = ['chr1','breakpoint1','breakpoint2','region','type']
    #merge_df.loc[index_to_delete,columns_to_delete] = None

    return merge_df
