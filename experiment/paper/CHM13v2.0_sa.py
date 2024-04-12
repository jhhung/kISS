import matplotlib.pyplot as plt
import numpy as np
import csv
import scienceplots
import matplotlib.ticker
import seaborn

algos = ['kiss-1', 'kiss-2', 'libsais', 'saca-k', 'pardss', 'libsais-32', 'pardss-32']
ks = [2, 4, 8, 16, 32, 64, 128, 256, -1]
fn = './combined_20240411.csv'
datat = ['example', 'chm13v2.0', 'mouse', 'zebrafish', 'human_protein_5640', 'mouse_protein_589', 'zebrafish_protein_437']
nn = [0, 3117292070, 2728222451, 1679203469, 11415812, 11747732, 15185681]
# datat = ['"example", "african_clawed_frog", "cat", "chicken", "dna", "dog", "hs37d5", "marmoset", "mouse", "protein", "zebrafish"']
# nn = [100000, 2718433805, 2521863845, 1065365425, 104857600, 2312802198, 3117292070, 2897824427, 2728222451, 104857600, 1679203469]
ax1_v_ = ["example", 480, 400, 160, 11, 380, 500, 490, 440, 12, 260]
ax2_l_ = ["example", 20, 30, 10, 1, 30, 40, 40, 30, 1, 20]
ax2_r_ = ["example", 145, 140, 65, 4.5, 130, 170, 160, 150, 5.5, 100]

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
        if row['k'] == '':
            k = -1
        else:
            k = int(row['k'])
        num_threads = int(row['num-threads'])
        t = float(row['time'])
        m = float(row['space'])
        if (num_threads == 32 and k in ks) or (algo == 'saca-k'):
            data[test][algo][k].append(t)
            memory[test][algo][k].append(m)

# calculate the average for each test, algo, k
for d in datat:
    for algo in algos:
        for k in ks:
            times = data[d][algo][k]            
            if len(times) == 0:
                continue

            if len(times) > 3:
                times = times[-3:]
            # print(d, algo, k, times)
            assert(len(times) == 3)
            # for i, p in enumerate(times):
            #     p[1] = i + 1
            # print(algo, k, times)
            times = sorted(times)
            data[d][algo][k] = sum(times) / len(times) 

            memories = memory[d][algo][k]
            memories = sorted(memories)
            memory[d][algo][k] = sum(memories) / len(memories)
        new_times = []
        for k in ks:
            new_times.append(data[d][algo][k])
        data[d][algo] = new_times

        if algo == 'kiss-1' or algo == 'kiss-2':
            memory[d][algo] = memory[d][algo][256]
        else:
            memory[d][algo] = memory[d][algo][-1]
        # print(d, algo, memory[d][algo])


x = np.arange(9)  # the label locations
width = 0.4  # the width of the bars

label = [str(2**i) for i in range(1, 9)] + ['unbounded']

plt.style.use(['science', 'nature'])
plt.rcParams.update({"font.size": 30})

cnt = 0

fig, ax = plt.subplots(2, 1, figsize=(17,13), gridspec_kw={'height_ratios': [1, 7], 'width_ratios': [1]})
color = seaborn.color_palette("Set1")

idx = 0
for d in ['chm13v2.0']:
    multiplier = 0

    ax1 = ax[(idx // 2) * 2]
    ax2 = ax[(idx // 2) * 2 + 1]
    idx = idx + 1
    fig.subplots_adjust(hspace=0.3)
    
    for algo, measurement in data[d].items():
        if algo == 'libsais' or algo == 'saca-k' or algo == 'pardss':
            if algo == 'saca-k':
                # print(measurement)
                ax1.axhline(y = measurement[8], linewidth=1, linestyle='dashed', color=color[6], label='SACAK')
                ax2.axhline(y = measurement[8], linewidth=1, linestyle='dashed', color=color[6], label='SACAK')
            elif algo == 'pardss': #pardss
                ax1.axhline(y = measurement[8], linewidth=1, linestyle='dashed', color=color[7], label='pDSS')
                ax2.axhline(y = measurement[8], linewidth=1, linestyle='dashed', color=color[7], label='pDSS')
            elif algo == 'libsais': # libsais
                # print(measurement)
                ax1.axhline(y = measurement[8], linewidth=1, linestyle='dashed', color=color[2], label='pSACAK+')
                ax2.axhline(y = measurement[8], linewidth=1, linestyle='dashed', color=color[2], label='pSACAK+')
        else:
            if algo == 'kiss-2':
                offset = width * multiplier
                ax1.bar(x + offset, measurement, width, label='kISS-2', color=color[1])
                ax2.bar(x + offset, measurement, width, label='kISS-2', color=color[1])
                multiplier += 1
            elif algo == 'kiss-1': # kiss
                offset = width * multiplier
                ax1.bar(x + offset, measurement, width, label='kISS-1', color=color[0])
                ax2.bar(x + offset, measurement, width, label='kISS-1', color=color[0])
                multiplier += 1

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax2.set_ylabel('time (seconds)', fontsize=30)
    ax2.set_xlabel('k', fontsize=30)
    ax1.set_title('CHM13v2.0' if d == 'chm13v2.0' else ('frog' if d == 'african_clawed_frog' else d))
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
    if idx == 1:
        ax2.legend(loc='lower center', ncols=3, bbox_to_anchor=(0.5, -0.4), title = 'Methods', fontsize=30)
# plt.savefig('./human_runtime_20240411.png', dpi=300)
for d in datat[1:]:
    for algo in algos:
        if type(memory[d][algo]) == float:
            print(d, algo, memory[d][algo] / 1024 / 1024 / 1024, memory[d][algo] / n[d])