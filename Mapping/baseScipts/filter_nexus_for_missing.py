#!/usr/bin/env python3

"""Take a nex or phy or fa file, outputs phylip-relaxed with X missing data

    arguments:
    data  -- extension should be nex/phy/fa
    missing -- the number of species in the alignment allowed to have missing data

    output:
    phylip formatted file ending with _mX.phylip-relaxed where X is the number missing
    """

from __future__ import division
import sys
from os import path
import linecache
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment, AlignInfo
from Bio import AlignIO, SeqIO

def main(alignment_filename, missing_str):

    ######################
    bases = ['A', 'C', 'G', 'T', 'a', 'c', 'g', 't','-']
    formats = {'nex':'nexus', 'phy':'phylip-relaxed', 'fa':'fasta'}
    fformat = formats[alignment_filename.split('.')[-1]]
    missing = int(missing_str)
    data = SeqIO.to_dict(SeqIO.parse(alignment_filename, fformat))

    #Extract loc IDs from nexus
    locline = linecache.getline(alignment_filename, 7)
    locs = locline.split()[1:-1] #Remove brackets from top and bottom

    species = list(data.keys())
    minsp = len(species)-missing
    newlocs_gap = list()
    newlocs_nogap = list()


    for k in data:
        data[k] = list(data[k].seq)

    newdata_gap = {sp: list() for sp in species}
    newdata_nogap = {sp: list() for sp in species}

    for i in range(len(data[species[0]])):
        site = [data[sp][i] for sp in species if data[sp][i] in bases]
        
        if len(set(site)) > 1 and len(site) >= minsp:
            for sp in species:
                newdata_gap[sp].append(data[sp][i])
            newlocs_gap.append(locs[i])
            
            if '-' not in site:
                for sp in species:
                    newdata_nogap[sp].append(data[sp][i])
                newlocs_nogap.append(locs[i])

    # Process gap data
    datalist = []
    for k,v in sorted(newdata_gap.items()):
        seq = SeqRecord(Seq(''.join(v)), id=k)
        datalist.append(seq)

    SeqIO.write(datalist, path.dirname(alignment_filename) + '/'+ path.basename(alignment_filename).split('.')[0] + '_m'+missing_str+'.phylip-relaxed', "phylip-relaxed")
    locfile = open(path.dirname(alignment_filename)+'/'+ path.basename(alignment_filename).replace('.nex','') + '_locs_m'+missing_str+'.txt', 'w')
    locfile.write("\n".join(newlocs_gap))
    locfile.close()
    origLength = len(data[species[0]])
    newLength = len(newlocs_gap)
    print('With '+missing_str+' taxa allowed to be missing, and with gaps allowed, '+str(origLength)+' sites from '+path.basename(alignment_filename)+' ('+str(len(species)-2)+' allowed missing) are reduced to '+str(newLength)+' sites ('+str(origLength-newLength)+' sites or '+str('%.2f' % (((origLength-newLength)/origLength)*100))+'% lost)')

    # Process gapless data
    datalist = []
    for k,v in sorted(newdata_nogap.items()):
        seq = SeqRecord(Seq(''.join(v)), id=k)
        datalist.append(seq)

    SeqIO.write(datalist, path.dirname(alignment_filename) + '/'+ path.basename(alignment_filename).split('.')[0] + '_m'+missing_str+'_nogap.phylip-relaxed', "phylip-relaxed")
    locfile = open(path.dirname(alignment_filename)+'/'+ path.basename(alignment_filename).replace('.nex','') + '_locs_m'+missing_str+'_nogap.txt', 'w')
    locfile.write("\n".join(newlocs_nogap))
    locfile.close()
    origLength = len(data[species[0]])
    newLength = len(newlocs_nogap)
    print('With '+missing_str+' taxa allowed to be missing, and with no gaps allowed, '+str(origLength)+' sites from '+path.basename(alignment_filename)+' ('+str(len(species)-2)+' allowed missing) are reduced to '+str(newLength)+' sites ('+str(origLength-newLength)+' sites or '+str('%.2f' % (((origLength-newLength)/origLength)*100))+'% lost)')


if __name__ == '__main__':
    alignment_filename = sys.argv[1]
    missing_str = sys.argv[2]
    main(alignment_filename, missing_str)

