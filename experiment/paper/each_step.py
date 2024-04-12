import numpy as np
import matplotlib.pyplot as plt
import csv
import seaborn
import scienceplots

items = ['classifyLSType', 'GetLMSPositions', 'parallelkOrderedLMSSuffixSort', 'getCharacterBucket','putLMSSuffixes','InduceSort','Misc']

def write_data(file, m, algo):
    with open(file) as f:
        lines = f.readlines()
        for line in lines:
            line_list = line.split()
            step = line_list[4]
            if ((step not in items) and (step != algo)):
                continue
            spend_time = float(line_list[6])
            if step not in m:
                m[step] = spend_time
            else:
                m[step] += spend_time
    total_time = 0
    for x in m:
        if x != algo:
            total_time += m[x]
    remain_time = m[algo] - total_time
    m.pop(algo)
    m['Misc'] = remain_time

data = {}
each_step = {}
datat = ["example", "african_clawed_frog", "cat", "chicken", "dna", "dog", "hs37d5", "marmoset", "mouse", "protein", "zebrafish"]
r_ = ["example", 80, 70, 40, 4, 70, 110, 80, 80, 5, 70]
algos = ['kiss', 'kiss-new']
ks = [2, 4, 8, 16, 32, 64, 128, 0]

r = {}
for i, d in enumerate(datat):
    r[d] = r_[i]

for d in datat:
    data[d] = {}
    each_step[d] = {}
    for algo in algos:
        data[d][algo] = {}
        each_step[d][algo] = {}
        for k in ks:
            data[d][algo][k] = []
            each_step[d][algo][k] = {}
            for item in items:
                each_step[d][algo][k][item] = 0.0

fn = './sa_new.csv'
with open(fn) as csvFile:
    csvDictReader = csv.DictReader(csvFile)
    for row in csvDictReader:
        algo = row['algo']
        if algo != 'kiss' and algo != 'kiss-new':
            continue
        test = row['test']
        k = int(row['k'])
        t = float(row['time'])
        if k in ks:
            times = len(data[test][algo][k]) + 1
            data[test][algo][k].append((t, times))

# calculate the average for each test, algo, k
for d in datat:
    for algo in algos:
        for k in ks:
            times = data[d][algo][k]
            times = sorted(times)
            for _, x in times:
                algo_ = ('kiss_new' if algo == 'kiss-new' else 'kiss')
                file_name = ('./sa_/log_%s_%s_%d_%d.txt' % (algo_, d, k, x))
                m1 = {}
                write_data(file_name, m1, algo)
                for step in m1:
                    each_step[d][algo][k][step] += m1[step] / 3

x = np.arange(8)  # the label locations
width = 0.4  # the width of the bars

label = [str(2**i) for i in range(1, 8)] + ['npos']
multiplier = 0

color = seaborn.color_palette("Set3")
f = {}
for i, step in enumerate(items):
    f[step] = color[i]

cnt = 0

plt.style.use(['science', 'nature'])
plt.rcParams.update({"font.size":30})

idx = 0
fig, ax_ = plt.subplots(3, 2, figsize=(30,40))
ax_[-1, -1].axis('off')
fig.subplots_adjust(hspace=0.3)
# for d in ['dna', 'protein', 'hs37d5', 'african_clawed_frog', 'cat', 'dog', 'chicken', 'marmoset', 'mouse', 'zebrafish']:
for d in ['dna', 'protein', 'hs37d5', 'mouse', 'zebrafish']:
    if d == 'example':
        continue
    ax = ax_[idx // 2][idx % 2]
    multiplier = 0
    idx = idx + 1
    for algo in algos:
        bottom = np.zeros(8)
        offset = width * multiplier
        for step in items:
            measurement = np.array([each_step[d][algo][k][step] for k in ks])
            ax.bar(x + offset, measurement, width, label=('InducedSort' if step == 'InduceSort' else step), bottom=bottom, color=f[step])
            bottom += measurement
        if idx == 5 and multiplier == 0:
            ax.legend(loc='center', ncols=1, bbox_to_anchor=(1.7, 0.6), title = 'Steps', fontsize=40)
        multiplier += 1
    ax.set_ylabel('time (s)', fontsize=30)
    ax.set_xlabel('k', fontsize=30)
    ax.set_title('CHM13v2.0' if d == 'hs37d5' else ('frog' if d == 'african_clawed_frog' else d))
    ax.set_xticks(x + 0.5 * width, label, fontsize=28)
    ax.set_yticks(np.arange(0, r[d], (10 if d != 'dna' and d != 'protein' else 1)),
     np.arange(0, r[d], (10 if d != 'dna' and d != 'protein' else 1)), fontsize=28)
plt.savefig('./sa/total_each_step.png', dpi=300)