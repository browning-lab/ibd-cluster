
# usage: zcat genotypes.vcf.gz | python3 find_gcs.py myibd.ibdclust.gz minmaf > outfile

# genotypes.vcf.gz is a VCF file of sequence data for one chromosome
# myibd.ibdclust.gz is the output from running ibdcluster.jar on that same VCF file
# minmaf is the minimum minor allele frequency, for example 0.1. It is recommended to use the same minor allele frequency for ibdcluster and this program

# outfile (stdout) has one line per detected allele conversion
# each such line has position (bp), size of subcluster 1, size of subcluster 2, identities of haplotypes in subcluster 1 and then those in subcluster 2.
# each haplotype identity has individual index (zero-based) : haplotype index (0 or 1)


import sys, gzip, time
from collections import Counter

omibdfile = gzip.open(sys.argv[1],mode='rt')
minaf = float(sys.argv[2]) # minimum allele frequency, in order to properly handle multiallelics
maxmaf = 1.0
minclustsize = 4
minpartsize = 2
deletion = 1

ids = omibdfile.readline().split()[3:]

for vline in sys.stdin:
    if vline[1]=="#": continue
    vids = vline.split()[9:]
    break
vindex = [vids.index(thisid)+9 for thisid in ids]

nindiv = len(vindex)
n = 2*nindiv # number of haplotypes

acclustsizes = Counter()

# checks whether the individual that carries hapindex is homozygous
def is_homozygous(valleles,hapindex):
    if hapindex < nindiv: otherindex = hapindex + nindiv
    else: otherindex = hapindex - nindiv
    return valleles[hapindex]==valleles[otherindex]

# checks whether the (heterozygous) individual's other cluster is mixed (so could be fixed with a rephasing)
def is_correctable(oclust,counts,hapindex):
    if hapindex < nindiv: otherindex = hapindex + nindiv
    else: otherindex = hapindex - nindiv
    otherclust = oclust[otherindex]
    thiscounts = counts[otherclust]
    return sum(thiscounts)!=max(thiscounts) # returns true if this is a mixed cluster

# print out the gc events for clusters that are mixed at this locus
def print_gcs(valleles,oclust,counts,pos,clustsize,ndel,acclustsizes):
    nclust = len(oclust)
    for j,x,z in zip(range(nclust),counts,clustsize):
        if z < minclustsize: continue # cluster too small to be mixed
        doubleplus = [y>= minpartsize for y in x]
        if sum(doubleplus) <= 1:
            continue
        if sum(doubleplus) > 2:
            print("more than two alleles in cluster at", pos, "printing first two",file=sys.stderr)
        nhaps = len(oclust)
        nindiv = nhaps/2
        hapindices = [i for i in range(nhaps) if oclust[i]==j]
        homozygous = [is_homozygous(valleles,x) for x in hapindices]
        if deletion and all(homozygous):
            ndel[0] += 1; continue
        whichalleles = [i for i in range(len(x)) if x[i]>1]
        which0 = [i for i in range(len(hapindices)) if valleles[hapindices[i]] == whichalleles[0] ]
        which1 = [i for i in range(len(hapindices)) if valleles[hapindices[i]] == whichalleles[1] ]
        hapindices0 = [hapindices[i] for i in which0]
        hapindices1 = [hapindices[i] for i in which1]
        print(pos, len(hapindices0), len(hapindices1),end=' ')
        for i in hapindices0 + hapindices1:
            print(':'.join([str(int(i % nindiv)), str(int(i/nindiv))]),end=' ')
        print()
        acclustsizes.update([len(hapindices)])
        


oprev = omibdfile.readline().split()
nhaps = (len(oprev) - 3)*2
onew = oprev
ndel = [0]; # counting how many times we invoke deletion hypothesis 
for vline in sys.stdin:
    if vline[0]=='#': continue
    #starttime = time.time()
    vbits = vline.split()
    pos = int(vbits[1])
    # find the closest positions in the multiIBD files
    while pos > int(onew[1]):
        oprev = onew
        onew = omibdfile.readline().split()
        # deal with end of file case
        if len(onew)==0:
            onew = oprev
            break
    othis = onew if abs(int(onew[1])-pos)<abs(int(oprev[1])-pos) else oprev
    # compare vcf variants with clusters
    valleles = [int(x.split('|')[0]) for x in [vbits[i] for i in vindex]] + [int(x.split('|')[1]) for x in [vbits[i] for i in vindex]]
    #n = len(valleles) # should not be needed, already defined above
    uniqalleles = list(set(valleles))
    nalleles = len(uniqalleles)
    if nalleles==1: continue # skip monomorphics
    maxalleleplus = max(uniqalleles)+1
    allelefreqs = [sum([x==i for x in valleles])/n for i in range(maxalleleplus)]
    maf = 1.0 - max(allelefreqs)
    if maf > maxmaf: continue
    if maf < minaf: continue
    oclust = [int(x.split('|')[0]) for x in othis[3:]] + [int(x.split('|')[1]) for x in othis[3:]]
    if len(valleles)!=nhaps:
        print("files don't match in number of individuals",file=sys.stderr)
    maxclustplus = max(oclust)+1
    counts = [[0]*maxalleleplus for j in range(maxclustplus)]
    clustsize = [0]*maxclustplus
    for i in range(nhaps):
        counts[oclust[i]][valleles[i]]+=1
        clustsize[oclust[i]]+=1
    # if truly multiallelic (3 or more alleles with allele count > 1) then remove any alleles with maf < minaf
    if nalleles>2:
        if sum([x>0 for x in allelefreqs])>2:
            for i in range(maxalleleplus):
                if allelefreqs[i] < minaf:
                    for x in counts:
                        x[i] = 0       
    print_gcs(valleles,oclust,counts,pos,clustsize,ndel,acclustsizes)


