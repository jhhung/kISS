import matplotlib.pyplot as plt
import numpy as np
import csv
import scienceplots
import matplotlib.ticker
import seaborn

algos = ['kiss', 'kiss-new', 'psais', 'sacak', 'pardss']
ks = [2, 4, 8, 16, 32, 64, 128, 0]
fn = './sa_new_20231024.csv'
datat = ["example", "african_clawed_frog", "cat", "chicken", "dna", "dog", "hs37d5", "marmoset", "mouse", "protein", "zebrafish"]
nn = [100000, 2718433805, 2521863845, 1065365425, 104857600, 2312802198, 3117292070, 2897824427, 2728222451, 104857600, 1679203469]
ax1_v_ = ["example", 410, 400, 160, 11, 380, 500, 490, 440, 12, 260]
ax2_l_ = ["example", 40, 30, 10, 1, 30, 40, 40, 30, 1, 20]
ax2_r_ = ["example", 140, 140, 65, 4.5, 130, 170, 160, 150, 5.5, 100]

n = {}
ax1_v = {}
ax2_l = {}
ax2_r = {}

for i, name in enumerate(datat):
    n[name] = nn[i]
    ax1_v[name] = ax1_v_[i]
    ax2_l[name] = ax2_l_[i]
    ax2_r[name] = ax2_r_[i]

data = {}
memory = {}
for d in datat:
    data[d] = {}
    memory[d] = {}
    for algo in algos:
        data[d][algo] = {}
        memory[d][algo] = {}
        for k in ks:
            data[d][algo][k] = []
            memory[d][algo][k] = []

with open(fn) as csvFile:
    csvDictReader = csv.DictReader(csvFile)
    for row in csvDictReader:
        algo = row['algo']
        test = row['test']
        k = int(row['k'])
        t = float(row['time'])
        m = float(row['space'])
        if k in ks:
            data[test][algo][k].append(t)
            memory[test][algo][k].append(m)

# calculate the average for each test, algo, k
for d in datat:
    for algo in algos:
        for k in ks:
            times = data[d][algo][k]
            times = sorted(times)
            data[d][algo][k] = sum(times) / 3

            memories = memory[d][algo][k]
            memories = sorted(memories)
            memory[d][algo][k] = sum(memories) / 3
        new_times = []
        for k in ks:
            new_times.append(data[d][algo][k])
        data[d][algo] = new_times

        # if algo == 'psais' or algo == 'sacak' or algo == 'pardss':
        #     memory[d][algo] = memory[d][algo][0]
        # else:
        #     memory[d][algo] = memory[d][algo][256]


x = np.arange(8)  # the label locations
width = 0.4  # the width of the bars

label = [str(2**i) for i in range(1, 8)] + ['u']

plt.style.use(['science', 'nature'])
plt.rcParams.update({"font.size": 30})

cnt = 0

fig, ax = plt.subplots(6, 2, figsize=(30,40), gridspec_kw={'height_ratios': [1, 3, 1, 3, 1, 3], 'width_ratios': [1, 1]})
ax[-1, -1].axis('off')
ax[-2, -1].axis('off')
color = seaborn.color_palette("Set1")

idx = 0
# for d in ['dna', 'protein', 'hs37d5', 'african_clawed_frog', 'cat', 'dog', 'chicken', 'marmoset', 'mouse', 'zebrafish']:
for d in ['dna', 'protein', 'hs37d5', 'mouse', 'zebrafish']:
    if d == 'example':
        continue
    multiplier = 0

    ax1 = ax[(idx // 2) * 2][idx % 2]
    ax2 = ax[(idx // 2) * 2 + 1][idx % 2]
    idx = idx + 1
    fig.subplots_adjust(hspace=0.3)
    
    for algo, measurement in data[d].items():
        if algo == 'psais' or algo == 'sacak' or algo == 'pardss':
            if algo == 'sacak':
                ax1.axhline(y = measurement[7], linewidth=1, linestyle='dashed', color=color[6], label='SACAK')
                ax2.axhline(y = measurement[7], linewidth=1, linestyle='dashed', color=color[6], label='SACAK')
                # psais
                ax1.axhline(y = measurement[7] / 9.305 * 2.849, linewidth=1, linestyle='dashed', color=color[4], label='pSAIS')
                ax2.axhline(y = measurement[7] / 9.305 * 2.849, linewidth=1, linestyle='dashed', color=color[4], label='pSAIS')
                # psacak+
                ax1.axhline(y = measurement[7] / 12.513 * 2.849, linewidth=1, linestyle='dashed', color=color[2], label='pSACAK+')
                ax2.axhline(y = measurement[7] / 12.513 * 2.849, linewidth=1, linestyle='dashed', color=color[2], label='pSACAK+')
            elif algo == 'pardss': #pardss
                ax1.axhline(y = measurement[7], linewidth=1, linestyle='dashed', color=color[7], label='pDSS')
                ax2.axhline(y = measurement[7], linewidth=1, linestyle='dashed', color=color[7], label='pDSS')
        else:
            if algo == 'kiss-new':
                offset = width * multiplier
                ax1.bar(x + offset, measurement, width, label='kISS-2', color=color[1])
                ax2.bar(x + offset, measurement, width, label='kISS-2', color=color[1])
                multiplier += 1
            else: # kiss
                offset = width * multiplier
                ax1.bar(x + offset, measurement, width, label='kISS-1', color=color[0])
                ax2.bar(x + offset, measurement, width, label='kISS-1', color=color[0])
                multiplier += 1

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax2.set_ylabel('time (s)', fontsize=30)
    ax2.set_xlabel('k', fontsize=30)
    ax1.set_title('CHM13v2.0' if d == 'hs37d5' else ('frog' if d == 'african_clawed_frog' else d))
    r = (10 if d != 'dna' and d != 'protein' else 1)
    ax1.set_ylim(ax1_v[d] - r, ax1_v[d] + r)
    ax2.set_ylim(ax2_l[d], ax2_r[d])
    ax2.set_xticks(x + 0.5 * width, label, fontsize=26)
    ax1.spines.bottom.set_visible(False)
    ax2.spines.top.set_visible(False)
    ax1.xaxis.tick_top()
    ax1.tick_params(labeltop=False)  # don't put tick labels at the top
    ax2.xaxis.tick_bottom()

    d_ = .5  # proportion of vertical to horizontal extent of the slanted line
    kwargs = dict(marker=[(-1, -d_), (1, d_)], markersize=8,
              linestyle="none", color='k', mec='k', mew=1, clip_on=False)
    ax1.plot([0, 1], [0, 0], transform=ax1.transAxes, **kwargs)
    ax2.plot([0, 1], [1, 1], transform=ax2.transAxes, **kwargs)

    ax2.set_yscale('log')
    ax1.yaxis.set_ticks(np.arange(ax1_v[d], ax1_v[d] + 50, 50), np.arange(ax1_v[d], ax1_v[d] + 50, 50), fontsize=26)
    ax1.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax2.yaxis.set_ticks(np.arange(ax2_l[d], ax2_r[d], r), np.arange(ax2_l[d], ax2_r[d], r), fontsize=26)
    ax2.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    if idx == 5:
        ax2.legend(loc='center', ncols=1, bbox_to_anchor=(1.7, 0.8), title = 'Methods', fontsize=40)
plt.savefig('./sa/total_runtime.png', dpi=300)
for d in datat:
    for algo in algos:
        print(d, algo, memory[d][algo] / n[d])