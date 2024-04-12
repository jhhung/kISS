import numpy as np
import matplotlib.pyplot as plt
import csv
import seaborn
import scienceplots

items = ['LMS Positions Finding', 'LMS suffixes k-ordered Sorting', 'Induced Sort', 'Misc']

def write_data(file, m, algo):
    with open(file) as f:
        lines = f.readlines()
        for line in lines:
            line_list = line.split()
            step = line_list[4]

            if line_list[5] != 'elapsed':
                continue

            # print(step, algo, step == algo)
            
            if step == algo:
                step_ = algo
            elif step == 'get_lms':
                step_ = 'LMS Positions Finding'
            elif step == 'lms_suffix_direct_sort' or step == 'prefix_doubling':
                step_ = 'LMS suffixes k-ordered Sorting'
            elif step == 'induced_sort' or step == 'put_lms_suffix':
                step_ = 'Induced Sort'
            else:
                step_ = 'Misc'

            # if ((step_ not in items) and (step != algo)):
            #     continue
            spend_time = float(line_list[6])
            if step_ not in m:
                m[step_] = spend_time
            else:
                m[step_] += spend_time
    total_time = 0
    for x in m:
        if x != algo and x != 'Misc':
            total_time += m[x]
    # print(m)
    remain_time = m[algo] - total_time
    m.pop(algo)
    # m.pop('Misc')
    m['Misc'] = remain_time

data = {}
each_step = {}
datat = ["example", 'chm13v2.0', 'mouse', 'zebrafish', 'human_protein_5640', 'mouse_protein_589', 'zebrafish_protein_437']
r_ = ["example", 510, 70, 40, 4, 70, 110, 80, 80, 5, 70]
algos = ['kiss-1', 'kiss-2']
list_threads = [1, 2, 4, 8, 16, 32, 64, 128]

r = {}
for i, d in enumerate(datat):
    r[d] = r_[i]

for d in datat:
    data[d] = {}
    each_step[d] = {}
    for algo in algos:
        data[d][algo] = {}
        each_step[d][algo] = {}
        for num_threads in list_threads:
            data[d][algo][num_threads] = []
            each_step[d][algo][num_threads] = {}
            for item in items:
                each_step[d][algo][num_threads][item] = 0.0

fn = './combined_20240411.csv'
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
        if (k == 256) and (algo == 'kiss-1' or algo == 'kiss-2'):
            times = len(data[test][algo][num_threads]) + 1
            data[test][algo][num_threads].append([t, times])

# calculate the average for each test, algo, k
for d in datat:
    for algo in algos:
        for num_threads in list_threads:
            times = data[d][algo][num_threads]
            if (len(times) == 0):
                continue
            # times = sorted(times)
            if len(times) > 3:
                times = times[-3:]
            # print(d, algo, k, times)
            assert(len(times) == 3)
            for i, p in enumerate(times):
                p[1] = i + 1
            for _, x in times:
                algo_ = ('kiss_1' if algo == 'kiss-1' else 'kiss_2')
                k = 256
                file_name = ('../../build/sa/log_%s_%s_%d_%d_%d.txt' % (algo_, d, k, num_threads, x))
                m1 = {}
                write_data(file_name, m1, algo)
                for step in m1:
                    each_step[d][algo][num_threads][step] += m1[step] / 3

x = np.arange(8)  # the label locations
width = 0.4  # the width of the bars

label = [str(2**i) for i in range(0, 8)]
multiplier = 0

color = seaborn.color_palette("Set3")
f = {}
for i, step in enumerate(items):
    if i == 0:
        i_ = 0
    elif i == 1:
        i_ = 2
    elif i == 2:
        i_ = 5
    else:
        i_ = 4
    f[step] = color[i_]

cnt = 0

plt.style.use(['science', 'nature'])
plt.rcParams.update({"font.size":30})

idx = 0
fig, ax_ = plt.subplots(1, 1, figsize=(17,13))
fig.subplots_adjust(hspace=0.3)
for d in ['chm13v2.0']:
    if d == 'example':
        continue
    ax = ax_
    multiplier = 0
    idx = idx + 1
    for algo in algos:
        bottom = np.zeros(8)
        offset = width * multiplier
        for step in items:
            measurement = np.array([each_step[d][algo][num_threads][step] for num_threads in list_threads])
            ax.bar(x + offset, measurement, width, label=('InducedSort' if step == 'InduceSort' else step), bottom=bottom, color=f[step])
            bottom += measurement
        if idx == 1 and multiplier == 0:
            ax.legend(loc='lower center', ncols=3, bbox_to_anchor=(0.5, -0.35), title = 'Steps', fontsize=30)
        multiplier += 1
    ax.set_ylabel('time (seconds)', fontsize=30)
    ax.set_xlabel('Number of threads', fontsize=30)
    ax.set_title('CHM13v2.0' if d == 'chm13v2.0' else ('frog' if d == 'african_clawed_frog' else d))
    ax.set_xticks(x + 0.5 * width, label, fontsize=28)
    ax.set_yticks(np.arange(0, r[d], (50 if d != 'dna' and d != 'protein' else 1)),
    np.arange(0, r[d], (50 if d != 'dna' and d != 'protein' else 1)), fontsize=28)
plt.savefig('./human_each_step_thread_20240411.png', dpi=300)