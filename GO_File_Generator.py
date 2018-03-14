import gzip
import zipfile
import pickle as pkl
from Modules import Module_2 as mdl
import urllib.request as url
from io import BytesIO


gzfile = url.urlopen('http://uswest.ensembl.org/biomart/martresults/18?file=martquery_0313222124_654.txt.gz')
ens_to_uni = {}
with gzip.open(gzfile, 'rt') as fi:
    for line in fi.readlines()[1:]:
        items = line.strip('\n').split(',')
        if items[-1] != '':
            ens_gen, ens_trn, uni = items[0], items[1], items[-1]
            ens_to_uni[ens_trn] = uni

ensgzip = url.urlopen('http://uswest.ensembl.org/biomart/martresults/18?file=martquery_0314202815_922.txt.gz')
more_ens = {}
with gzip.open(ensgzip, 'rt') as fi:
    for line in fi.readlines()[1:]:
        ens_gen, ens_trn, ens_pro = line.strip('\n').split(',')
        more_ens[ens_trn] = more_ens.get(ens_trn, dict())
        more_ens[ens_trn]['gen'] = more_ens[ens_trn].get('gen', set())
        more_ens[ens_trn]['gen'].add(ens_gen)
        more_ens[ens_trn]['pro'] = more_ens[ens_trn].get('pro', set())
        more_ens[ens_trn]['pro'].add(ens_pro)

anizip = url.urlopen('https://www.aniseed.cnrs.fr/aniseed/download/?file=data%2Fci%2FKH-ENS.blast.zip')
kh_to_uni = {}
kh_to_ens = {}
kh_to_ens_gen = {}
with zipfile.ZipFile(BytesIO(anizip.read())) as zipdir:
    with zipdir.open('KH-ENS.blast', 'r') as fi:
        for line in fi.readlines():
            line = line.decode()
            ens, kh = line.split()[:2]
            try:
                kh_to_ens[kh] = [ens, more_ens[ens]]
            except KeyError:
                pass
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

with gzip.open('/home/jeff/Desktop/Thesis/Output/Ciona Robusta.gene2uniprot.txt.gz', 'wt') as fo:
    fo.write('UniqueID#EntrezGeneID\tUniProtKB_AC\n')
    for kh, uni in kh_to_uni.items():
        uni = '|'.join(uni)
        fo.write('\t'.join([kh, uni]) + '\n')

with gzip.open('/home/jeff/Desktop/Thesis/Output/Ciona Robusta.gene2ensembl.txt.gz', 'wt') as fo:
    fo.write('UniqueID#EntrezGeneID\tEnsemblGeneID\tEnsemblTranscriptID\tEnsemblProteinID\n')
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

with gzip.open('/home/jeff/Desktop/Thesis/Output/Ciona Robusta.gene2accession.txt.gz', 'wt') as fo:
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

