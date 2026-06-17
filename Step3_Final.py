f = open('peak_gene.txt')
pg = f.readlines();f.close()
f = open('dis')
dis = f.readlines();f.close()
f = open('corr')
corr = f.readlines();f.close()

g = open('../Prior/RE_gene_corr_monkey.bed','w')
for i in range(len(pg)):
	pg[i] = pg[i].strip('\n').split('\t')
	g.write(pg[i][0].replace('_','\t')+'\t'+pg[i][1]+'\t'+dis[i].strip('\n')+'\t'+corr[i].strip('\n')+'\n')
g.close()


f = open('SampleNameFile.txt')
a = f.readlines();f.close()
Sample = []
for i in range(len(a)):
        a[i] = a[i].strip('\n')
        Sample.append(a[i])
RPKM = []
for i in range(len(Sample)):
    output_files = './ATAC/'+Sample[i]+'.bed'
    f = open(output_files)
    a = f.readlines();f.close()
    for j in range(len(a)):
        a[j] = a[j].strip('\n').split('\t')
        a[j] = float(a[j][1])
    RPKM.append(a)

import numpy as np

RPKM = np.median(np.array(RPKM),0)
f = open('../../../ATACPeaks/Monkey_dTAC.bed')
a = f.readlines();f.close()
g = open('../Prior/Opn_median_monkey.bed','w')
for i in range(len(a)):
	g.write(a[i].strip('\n')+'\t'+str(RPKM[i])+'\n')
g.close()
