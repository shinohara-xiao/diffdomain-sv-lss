import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator
import seaborn as sns
from scipy import stats
from scipy import integrate
from scipy.stats import norm  #导入正态分布
import matplotlib.mlab as mlab
import math

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

def sv_with_retads(svset,tadset):
    overlap = lying = 0
    for i in range(svset.shape[0]):
        #print(i)
        chrom = svset['chr1'].iloc[i].split('r')[1]
        start = svset['breakpoint1'].iloc[i]
        end = svset['breakpoint2'].iloc[i]
        lying1 = tadset[(tadset['chr'] == chrom) & (tadset['start'] < start) & (tadset['end'] > start)]
        lying2 = tadset[(tadset['chr'] == chrom) & (tadset['start'] < end) & (tadset['end'] > end)]
        overlapp = tadset[(tadset['chr'] == chrom) & (tadset['start'] > start) & (tadset['end'] < end)]
        if lying1.shape[0] > 0:
            lying += 1
        if lying2.shape[0] > 0:
            lying += 1
        if overlapp.shape[0] > 0:
            overlap += 1

    lying_ratio = 100 * lying / (2*svset.shape[0])
    overlap_ratio = 100 * overlap / svset.shape[0]

    return lying_ratio,overlap_ratio    


def Calculate_Ratio(svset,tadset,svtype = None,tadtype = None):
    # Just input the whole svset dataset and also ,all tadsset , they are all dataset that pd.read_table get
    # all sv and all retads
    retads = tadset[(tadset['significant'] == 1) & (tadset['origin'] == 'condition1')]  # 所有condition1 下的reorganized tads

    if (svtype == None) & (tadtype == None):
        svdf = svset
        # 计算有多少breakpoint 在retads里面 (breakpoint lying in tads)
        lying_ratio,overlap_ratio = sv_with_retads(svset = svdf,tadset = retads)

    # svtype 对应的sv 和 all retads
    elif tadtype == None:
        #print('yes')
        svdf = svset[svset['type'] == svtype]
        #print('...')
        #print(svdf)
        lying_ratio,overlap_ratio = sv_with_retads(svset = svdf,tadset = retads)

    
    # svtype 对应的sv 和 subtype of retads
    else:
        svdf = svset[svset['type'] == svtype]

        if tadtype == 'strength' or tadtype == 'zoom':
            retads = retads[retads['subtype'] == tadtype]
            lying_ratio,overlap_ratio = sv_with_retads(svset = svdf,tadset = retads)
        else:
            retads = retads[retads['type'] == tadtype]
            lying_ratio,overlap_ratio = sv_with_retads(svset = svdf,tadset = retads)
    return lying_ratio,overlap_ratio



def Random_Simulation(svset,tadset,svtype = None,tadtype = None,times = 10000):
    tads_condition1 = tadset[tadset['origin'] == 'condition1']
    retadsnum = tadset[(tadset['significant'] == 1) & (tadset['origin'] == 'condition1')].shape[0]
    # 真实的retads 与 sv 之间的关系
    lying_ratio,overlap_ratio = Calculate_Ratio(svset = svset,tadset = tadset,svtype = svtype,tadtype = tadtype)
    
    # simulation  for times and plot  - 随机选择tads 作为retads
    Lying_sim,Overlap_sim = [],[]
    for t in range(times):
        sn = np.random.randint(0,tads_condition1.shape[0],retadsnum)
        random_choose_retads = tads_condition1.iloc[sn]
        lying_ratio_sim,overlap_ratio_sim = Calculate_Ratio(svset = svset,tadset = random_choose_retads,svtype = svtype,tadtype = tadtype)
        Lying_sim.append(lying_ratio_sim)
        Overlap_sim.append(overlap_ratio_sim)

    return lying_ratio,overlap_ratio,Lying_sim,Overlap_sim

condition1,condition2,resolution = 'GM12878','K562','10000'
sv = Extract_SVs(celltype = condition2)
svtype = None
tadtype = None
tads = pd.read_table(f'/lanec2_home/zhangx/DiffCompare/result/{condition1}_vs_{condition2}_{resolution}_to1_classification.txt_types.txt')

lying_ratio,overlap_ratio,Lying_sim,Overlap_sim = Random_Simulation(svset = sv,tadset = tads,times = 100)

miu_lying , cigema_lying = np.mean(Lying_sim), np.std(Lying_sim)
miu_Over , cigema_Over = np.mean(Overlap_sim), np.std(Overlap_sim) 

# break_points lying in retads
x = np.linspace(0,lying_ratio+0.25,10000)
miu,cigema = miu_lying,cigema_lying
fx = 1 / (cigema * (2 * np.pi)**0.5) * np.exp(-(x - miu)**2 / (2 * cigema**2))
plt.plot(x,fx,color='dodgerblue')
x_major_locator=MultipleLocator(10)#以每15显示
ax=plt.gca()
ax.xaxis.set_major_locator(x_major_locator)
plt.axvline(x = lying_ratio, c='r', ls='--', lw=2)  # 垂直于x轴的参考线
plt.title(f'{condition1}_vs_{condition2}_{resolution}_lying')
# 添加pvlaue
pvalue = stats.ttest_1samp(Lying_sim,lying_ratio)[1]
c = 0.7 * (min(x) + max(x))
d = min(fx) + 0.9 * (max(fx) - min(fx))
plt.text(c,d,pvalue,fontsize = 4)
plt.show()
plt.savefig(f'/lanec2_home/zhangx/svs/simulation/{condition2}_{resolution}_lying_{svtype}_{tadtype}.pdf')

# Sv overlap the retads
x = np.linspace(0,overlap_ratio+0.25,10000)
miu,cigema = miu_Over,cigema_Over
fx = 1 / (cigema * (2 * np.pi)**0.5) * np.exp(-(x - miu)**2 / (2 * cigema**2))
plt.plot(x,fx,color='dodgerblue')
x_major_locator=MultipleLocator(10)#以每15显示
ax=plt.gca()
ax.xaxis.set_major_locator(x_major_locator)
plt.axvline(x = overlap_ratio, c='r', ls='--', lw=2)  # 垂直于x轴的参考线
plt.title(f'{condition1}_vs_{condition2}_{resolution}_overlap')

pvalue = stats.ttest_1samp(Lying_sim,lying_ratio)[1]
c = 0.7 * (min(x) + max(x))
d = min(fx) + 0.9 * (max(fx) - min(fx))
plt.text(c,d,pvalue,fontsize = 4)

plt.show()
plt.savefig(f'/lanec2_home/zhangx/svs/simulation/{condition2}_{resolution}_overlap_{svtype}_{tadtype}.pdf')

