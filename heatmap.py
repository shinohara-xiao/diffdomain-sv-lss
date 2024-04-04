data = np.array([
    [8,8,3,1,0,3],
    [6,2,5,6,0,1],
    [1,11,0,1,1,0],
    [9,8,2,1,4,0],
    [7,1,12,13,0,0],
    [0,1,0,2,2,0],
    [6,1,2,2,0,0],
    [7,2,6,4,1,1],
    [1,9,3,0,1,0],
    [40,10,1,3,4,14],
    [4,0,4,3,1,0,],
    [0,3,0,0,0,0]
])

row_sums = data.sum(axis = 1)
ratios = data / row_sums[:,np.newaxis]
ratios

plt.figure(figsize = (10,10))
ax = sns.heatmap(data,cmap ="YlOrBr", annot=False, fmt=".2f", square=True
            ,vmin = 0,vmax = 10)
# 不显示x，y坐标
ax.axis('off')
# 对齐颜色条和图形
cbar_ax = ax.figure.axes[-1]
cbar_ax.set_position([ax.get_position().x1 + 0.01, ax.get_position().y0, 0.03, ax.get_position().height])


plt.title('Num of Reorganzied TADs of each Subtypes and svtype')


plt.figure(figsize = (10,10))
ax = sns.heatmap(ratios,cmap ="YlOrBr", annot=False, fmt=".2f", square=True
            ,vmin = 0,vmax = 100)
# 不显示x，y坐标
ax.axis('off')
# 对齐颜色条和图形
cbar_ax = ax.figure.axes[-1]
cbar_ax.set_position([ax.get_position().x1 + 0.01, ax.get_position().y0, 0.03, ax.get_position().height])

plt.title('Percent(%) of Reorganzied TADs of each Subtypes and svtype')
