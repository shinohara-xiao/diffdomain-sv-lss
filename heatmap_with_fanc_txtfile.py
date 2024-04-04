import pandas as pd
import numpy as np
import os
import seaborn as sns
import matplotlib.pyplot as plt

# ---------- all tads ------------ #
celltype = 'GM12878'
df = pd.read_table(f'/lanec2_home/zhangx/FANC/result/fanc_{celltype}_alltads_25000_r0.5_e_rescale.agg.txt',header = None)

multi_array = None
for i in range(df.shape[0]):
    str_array = list(map(float,df[0][i].split(' ')))
    if multi_array is None:
        multi_array = str_array
    else:
        multi_array = np.vstack((multi_array,str_array))
multi_array

plt.figure(figsize = (10,10))
ax = sns.heatmap(multi_array,cmap ="YlOrBr", annot=False, fmt=".2f", square=True
            ,vmin = 0,vmax = 0.1)
# 不显示x，y坐标
ax.axis('off')
# 对齐颜色条和图形
cbar_ax = ax.figure.axes[-1]
cbar_ax.set_position([ax.get_position().x1 + 0.01, ax.get_position().y0, 0.03, ax.get_position().height])

plt.title(f'{celltype}_25000_alltads')
plt.savefig(f'/lanec2_home/zhangx/FANC/result/snsplot/{celltype}_25000_alltads.pdf')
plt.show()

# ----- subtypes for each type of SV ----- #
#condition1,condition2 = 'GM12878','K562'
#condition1,condition2 = 'GM12878','IMR90'
#condition1,condition2 = 'NHA','DIPG007'
#condition1,condition2 = 'NHA','DIPGXIII'
condition1,condition2 = 'HMEC','NHEK'

types = ['complex','loss','merge','split','strength','zoom']
subtype = types[0]


df = pd.read_table(f'/lanec2_home/zhangx/FANC/result/retads/subtypes_of_celllines/{condition1}_vs_{condition2}_25000_{subtype}_retads_r0.5_e_rescale.agg.txt',header = None)
multi_array = None
for i in range(df.shape[0]):
    str_array = list(map(float,df[0][i].split(' ')))
    if multi_array is None:
        multi_array = str_array
    else:
        multi_array = np.vstack((multi_array,str_array))
retads = multi_array

df = pd.read_table(f'/lanec2_home/zhangx/FANC/result/fanc_{condition1}_alltads_25000_r0.5_e_rescale.agg.txt',header = None)
multi_array = None
for i in range(df.shape[0]):
    str_array = list(map(float,df[0][i].split(' ')))
    if multi_array is None:
        multi_array = str_array
    else:
        multi_array = np.vstack((multi_array,str_array))
multi_array
alltads = multi_array

difftads = np.log2(retads/alltads)

# --- plot
plt.figure(figsize = (10,10))
ax = sns.heatmap(retads,cmap ="YlOrBr", annot=False, fmt=".2f", square=True
            ,vmin = 0,vmax = 0.1)
# 不显示x，y坐标
ax.axis('off')
# 对齐颜色条和图形
cbar_ax = ax.figure.axes[-1]
cbar_ax.set_position([ax.get_position().x1 + 0.01, ax.get_position().y0, 0.03, ax.get_position().height])

plt.title(f'{condition1}_vs_{condition2}_25000_{subtype}_retads')
plt.savefig(f'/lanec2_home/zhangx/FANC/result/snsplot/{condition1}_vs{condition2}_25000_{subtype}_retads.pdf')
plt.show()



plt.figure(figsize = (10,10))
ax = sns.heatmap(difftads,cmap ="coolwarm", annot=False, fmt=".2f", square=True
            ,vmin = -1,vmax = 1)
# 不显示x，y坐标
ax.axis('off')
# 对齐颜色条和图形
cbar_ax = ax.figure.axes[-1]
cbar_ax.set_position([ax.get_position().x1 + 0.01, ax.get_position().y0, 0.03, ax.get_position().height])
plt.title(f'{condition1}_vs_{condition2}_25000_{subtype}_Diff')
plt.savefig(f'/lanec2_home/zhangx/FANC/result/snsplot/{condition1}_vs{condition2}_25000_{subtype}_Diff.pdf')
plt.show()

print(f'scp zhangx@172.25.49.139:/lanec2_home/zhangx/FANC/result/snsplot/{condition1}_vs{condition2}_25000_{subtype}_retads.pdf /Users/shinohara_xiao/Documents/diffdomain/first-review/APA/snsplot/celllines-25kb/Subtypes/')
print(f'scp zhangx@172.25.49.139:/lanec2_home/zhangx/FANC/result/snsplot/{condition1}_vs{condition2}_25000_{subtype}_Diff.pdf /Users/shinohara_xiao/Documents/diffdomain/first-review/APA/snsplot/celllines-25kb/Subtypes/')
