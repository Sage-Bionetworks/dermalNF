import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import csv
csv.register_dialect('tab', delimiter='\t', quoting=csv.QUOTE_NONE)
avg   = []
dp    = []
sm2id = {}
dat   = None
min   = None
max   = None
with open('allSamples.tab', 'rb') as f:
   reader = csv.reader(f, 'tab')
   i = 0
   for row in reader:
       if row[0]=='SM':
           sm2id[row[4]] = i
           avg.append([i,float(row[1])])
           dp.append([i,float(row[2])])
           i += 1
       elif row[0]=='CN':
           val = 0
           if int(row[2])!=0: val = float(row[1])/int(row[2])
           if not dat:
               dat = [[0]*len(sm2id) for x in xrange(len(sm2id))]
               min = val
               max = val
           id_i = sm2id[row[4]]
           id_j = sm2id[row[5]]
           dat[id_i][id_j] = val
           dat[id_j][id_i] = val
           if min > val: min = val
           if max < val: max = val

if len(sm2id)<=1: exit(1)
if min==max: exit(1)

fig = plt.figure(figsize=(6,7))
gs  = gridspec.GridSpec(2, 1, height_ratios=[1, 1.5])
ax1 = plt.subplot(gs[0])
ax2 = plt.subplot(gs[1])

ax1.plot([x[0] for x in avg],[x[1] for x in avg],'^-', ms=3, color='k')
ax3 = ax1.twinx()
ax3.plot([x[0] for x in dp],[x[1] for x in dp],'^-', ms=3, color='r',mec='r')
for tl in ax3.get_yticklabels():
   tl.set_color('r')
   tl.set_fontsize(9)

im = ax2.imshow(dat,clim=(min),interpolation='nearest',origin='lower')
cb1  = plt.colorbar(im,ax=ax2)
cb1.set_label('Pairwise discordance')
for t in cb1.ax.get_yticklabels(): t.set_fontsize(9)

ax1.tick_params(axis='both', which='major', labelsize=9)
ax1.tick_params(axis='both', which='minor', labelsize=9)
ax2.tick_params(axis='both', which='major', labelsize=9)
ax2.tick_params(axis='both', which='minor', labelsize=9)

ax1.set_title('Sample Discordance Score')
ax2.set_ylabel('Sample ID')
ax2.set_xlabel('Sample ID')
ax3.set_ylabel('Average Depth',color='r')
ax1.set_xlabel('Sample ID')
ax1.set_ylabel('Average discordance')

plt.subplots_adjust(left=0.15,right=0.87,bottom=0.08,top=0.93,hspace=0.25)
plt.savefig('allSamples.png')
plt.close()

