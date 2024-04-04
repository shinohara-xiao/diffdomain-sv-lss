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

#file_path = '/lanec2_home/zhangx/svs/data/Summary_of_sv_tads.xlsx'
method_change = 'to1'

#with pd.ExcelWriter(file_path) as writer:
condition1,condition2 = 'GM12878', 'K562'
data = CreateExcel(condition1,condition2,method_change)
#print(data[data['region'] == '17:26897018-27272018'])
data.to_csv('/lanec2_home/zhangx/svs/data/K562_summary.csv')

condition1,condition2 = 'NHA', 'DIPG007'
data = CreateExcel(condition1,condition2,method_change)
data.to_csv('/lanec2_home/zhangx/svs/data/DIPG007_summary.csv')

condition1,condition2 = 'NHA', 'DIPGXIII'
data = CreateExcel(condition1,condition2,method_change)
data.to_csv('/lanec2_home/zhangx/svs/data/DIPGXIII_summary.csv')

# ------ 计算tads 与 sv 之间是否有overlap ----- #

def generate_region(svregion,tadsregion,svtype,K):
    # svregion and tadsregion format should be : [chr:start-end]
    #sv_start, sv_end = svdataset['start'], svdataset['end']
    # 定义sv的影响区域
    sv_start = int(svregion.split('-')[0].split(':')[1])
    sv_end = int(svregion.split('-')[1])

    tad_start = int(tadsregion.split('-')[0].split(':')[1])
    tad_end = int(tadsregion.split('-')[1])
    #svtype = 

    ltad = tad_end - tad_start
    if svtype == '+-':  #deletion
        k = K
        arm = k*(sv_end-sv_start)
        x1, x2 = sv_start - arm, sv_start
        y1, y2 = sv_end,  sv_end + arm
    elif svtype == '-+':   # duplication
        k = K
        arm = k*(sv_end-sv_start)
        x1, x2= sv_start, sv_start + arm
        y1, y2 = sv_end - arm, sv_end
    elif svtype == '++':  # '++'
        k = K
        arm = k*(sv_end-sv_start)
        x1, x2 = sv_start, sv_start+arm
        y1, y2 = sv_end, sv_end + arm
    else : # "--"
        k = K
        arm = k*(sv_end-sv_start)
        x1, x2 = sv_start -arm,  sv_start
        y1, y2 = sv_end - arm,  sv_end
    sv_row = (x1,x2)
    sv_col = (y1,y2)
    tad_row = (tad_start,tad_start+ltad)
    tad_col = (tad_end-ltad,tad_end)
    return sv_row,sv_col,tad_row,tad_col

def are_intersecting(seg1,seg2):
    x1_seg1, x2_seg1 = seg1
    x1_seg2, x2_seg2 = seg2

    # 确保seg1的左端点在seg2的左侧
    if x1_seg1 > x2_seg1:
        x1_seg1, x2_seg1 = x2_seg1, x1_seg1
    # 确保seg2的左端点在seg1的左侧
    if x1_seg2 > x2_seg2:
        x1_seg2, x2_seg2 = x2_seg2, x1_seg2
    # 判断是否有重叠
    if x2_seg1 >= x1_seg2 and x1_seg1 <= x2_seg2:
        return True
    return False

def Caldistance(celltype,Kvalue):
    #tads = pd.read_table(f'/lanec2_home/zhangx/DiffCompare/result/{condition1}_vs_{condition2}_10000_to1_f1b5_BH.txt')
    #svs = pd.read_table(f'/lanec2_home/zhangx/DiffCompare/result/{condition2}_sv_{svtype}.bed',header = None)
    #svs = svs.rename(columns={0:'chr',1:'start',2:'end'})

    #file_path = '/lanec2_home/zhangx/svs/data/Summary_of_sv_tads.xlsx'
    #data = pd.read_excel(file_path,sheet_name=celltype)

    data = pd.read_csv(f'/lanec2_home/zhangx/svs/data/{celltype}_summary.csv')
    overlap = []
    for i in range(data.shape[0]):
        sg = data['region'][i]
        tg = data['related_tads'][i]
        st = data['type'][i]
        sv_row,sv_col,tad_row,tad_col = generate_region(svregion=sg, tadsregion=tg,svtype = st,K = Kvalue)
        if are_intersecting(sv_row,tad_row) and are_intersecting(sv_col,tad_col):
            overlap.append('yes')
        else:
            overlap.append('no')
    data['overlap'] = overlap
    #print(data[data['region'] == '17:26897018-27272018'])

    # 删除重复的region避免数据过于冗余
    index_to_delete = data[data.duplicated(subset='region')].index
    columns_to_delete = ['chr1','breakpoint1','breakpoint2','region','type','retads_related_raw','retads_related_to1']
    data.loc[index_to_delete,columns_to_delete] = None
    return data


# ---------------- USE Function and CODEs ---------------- #

#distance_path = '/lanec2_home/zhangx/svs/data/Summary_distance.xlsx'
#with pd.ExcelWriter(distance_path) as writer:
data = Caldistance('K562',Kvalue=0.7)
data.to_csv('/lanec2_home/zhangx/svs/data/K562_distance.csv',index = False)

data = Caldistance('DIPG007',Kvalue = 0.7)
data.to_csv('/lanec2_home/zhangx/svs/data/DIPG007_distance.csv',index = False)

data = Caldistance('DIPGXIII',Kvalue = 0.7)
data.to_csv('/lanec2_home/zhangx/svs/data/DIPGXIII_distance.csv',index = False)
