import gzip
import zipfile
from Modules import Module_2 as mdl
import urllib.request as url
from io import BytesIO


gzfile = url.urlopen('http://uswest.ensembl.org/biomart/martresults/18?file=martquery_0319215120_830.txt.gz')
ens_dict = {}
with gzip.open(gzfile, 'rt') as fi:
    for line in fi.readlines()[1:]:
        items = line.strip('\n').split('\t')
        ens_gen, ens_trn, ens_pro, ens_ent, uni = items
        ens_dict[ens_trn] = ens_dict.get(ens_trn, {'gen': set(), 'pro': set(), 'uni': set(), 'ent': set()})
        if ens_gen:
            ens_dict[ens_trn]['gen'].add(ens_trn)
        if ens_pro:
            ens_dict[ens_trn]['pro'].add(ens_pro)
        if ens_ent:
            ens_dict[ens_trn]['ent'].add(ens_ent)
        if uni:
            ens_dict[ens_trn]['uni'].add(uni)

print(ens_dict)

"""
anizip = url.urlopen('https://www.aniseed.cnrs.fr/aniseed/download/?file=data%2Fci%2FKH-ENS.blast.zip')
kh_to_uni = {}
#kh_to_ens = {}
kh_to_ens_gen = {}
with zipfile.ZipFile(BytesIO(anizip.read())) as zipdir:
    with zipdir.open('KH-ENS.blast', 'r') as fi:
        for line in fi.readlines():
            line = line.decode()
            ens, kh = line.split()[:2]
            kh = mdl.scr_to_gen(kh.split(':')[1])
            try:
                kh_to_ens_gen[kh] = kh_to_ens_gen.get(kh, more_ens[ens])
                kh_to_ens_gen[kh]['trn'] = kh_to_ens_gen[kh].get('trn', set())
                kh_to_ens_gen[kh]['trn'].add(ens)
            except KeyError:
                pass
            try:
                uni = ens_to_uni[ens]
                kh_to_uni[kh] = kh_to_uni.get(kh, set())
                kh_to_uni[kh].add(uni)
            except KeyError:
                pass

with gzip.open('/home/jeff/Desktop/Thesis/Output/Ciona Intestinalis.gene2uniprot.txt.gz', 'wt') as fo:
    fo.write('UniqueID#EntrezGeneID\tUniProtKB_AC\n')
    for kh, uni in kh_to_uni.items():
        uni = '|'.join(uni)
        fo.write('\t'.join([kh, uni]) + '\n')

with gzip.open('/home/jeff/Desktop/Thesis/Output/Ciona Intestinalis.gene2ensembl.txt.gz', 'wt') as fo:
    fo.write('UniqueID#EntrezGeneID\tEnsemblGeneID\tEnsemblTranscriptID\tEnsemblProteinID\n')
    """
"""
    for kh, ens_list in kh_to_ens.items():
        ens_trn, ens_dict = ens_list
        ens_gen, ens_pro = ens_dict['gen'], ens_dict['pro']
        ens_gen = '|'.join(ens_gen)
        ens_pro = '|'.join(ens_pro)
        if ens_gen == '':
            ens_gen = '-'
        if ens_pro == '':
            ens_pro = '-'
        fo.write('\t'.join([kh, ens_gen, ens_trn, ens_pro]) + '\n')
    """
"""    
    for kh, ens_dict in kh_to_ens_gen.items():
        ens_gen, ens_trn, ens_pro = '|'.join(ens_dict['gen']), '|'.join(ens_dict['trn']), '|'.join(ens_dict['pro'])
        if ens_gen == '':
            ens_gen = '-'
        if ens_trn == '':
            ens_trn = '-'
        if ens_pro == '':
            ens_pro = '-'
        fo.write('\t'.join([kh, ens_gen, ens_trn, ens_pro]) + '\n')


gzipf = url.urlopen('https://www.aniseed.cnrs.fr/aniseed/download/?file=data%2Fci%2FCirobu_KHslimTunicate.gaf.gz')
kh_to_name = {}
with gzip.open(gzipf, 'rt') as fi:
    for line in fi.readlines()[5:]:
        items = line.strip('\n').split('\t')
        name, kh = items[2], items[5]
        name = name.split('; ')
        kh = kh.split(':')[1]
        kh_to_name[kh] = kh_to_name.get(kh, set())
        kh_to_name[kh] = kh_to_name[kh].union(name)

gen_set = set(kh_to_name.keys()) | set(kh_to_uni.keys())

with gzip.open('/home/jeff/Desktop/Thesis/Output/Ciona Intestinalis.gene2accession.txt.gz', 'wt') as fo:
    fo.write('UniqueID#EntrezGeneID\tSymbolID\tAccessionID\n')
    for kh in gen_set:
        try:
            name = '|'.join(kh_to_name[kh])
        except KeyError:
            name = '-'
        try:
            uni = '|'.join(kh_to_uni[kh])
        except KeyError:
            uni = '-'
        fo.write('\t'.join([kh, name, uni]) + '\n')
"""
