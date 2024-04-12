import csv
import matplotlib.pyplot as plt
import scienceplots
import seaborn

plt.style.use(['science', 'nature'])
plt.rcParams.update({"font.size":26})
plt.figure(figsize=(12.75, 9.75))
plt.rc('xtick', labelsize=24)
plt.rc('ytick', labelsize=24)

x = [[[i + 1 for i in range(8)] for _ in range(2)] for _ in range(4)]
y = [[[[] for i in range(8)] for _ in range(2)] for _ in range(4)]

fn = 'fmtree.csv'
testcases = 100000
lengths = [12, 18, 24, 30]
back_search = {}
for i in range(4):
    back_search[lengths[i]] = i

color = seaborn.color_palette("Set1")

with open(fn) as csvFile:
    csvDictReader = csv.DictReader(csvFile)
    for row in csvDictReader:
        length_index = back_search[int(row['pattern_length'])]
        intv_index = int(row['SA_INTV']) - 1
        y[length_index][0][intv_index].append(float(row['fmtree']))
        y[length_index][1][intv_index].append(float(row['naive']))

for length_index in range(0, 4):
    for i in range(0, 8):
        y[length_index][0][i] = sum(y[length_index][0][i][1:4]) / 3 / testcases
        y[length_index][1][i] = sum(y[length_index][1][i][1:4]) / 3 / testcases

line_feed = ['solid', 'dashed', 'dotted',(0, (3, 5, 1, 5))]

for idxx in range(4):
    length = lengths[idxx]
    a = plt.plot(x[idxx][0][1:], y[idxx][0][1:], color=color[idxx], marker='.', linestyle=line_feed[idxx], label='Orginal', markersize=20)
    b = plt.plot(x[idxx][1][1:], y[idxx][1][1:], color=color[idxx], marker='s', linestyle=line_feed[idxx], label='FM-Tree', markersize=10)
    plt.xlabel('Value sampling distance', fontsize=26)
    plt.ylabel('Locating time (seconds)', fontsize=26)
    plt.yscale('log')
plt.legend(["FMtree, Length = 12", "Original, Length = 12", "FMtree, Length = 18", "Original, Length = 18",
"FMtree, Length = 24", "Original, Length = 24", "FMtree, Length = 30", "Original, Length = 30"], loc='lower center', bbox_to_anchor=(0.5, -0.4), ncols=2, fontsize=26)
plt.savefig('fmtree_2.png', dpi=400)
plt.clf()