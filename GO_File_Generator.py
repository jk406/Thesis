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
            ens, uni = items[1], items[-1]
            ens_to_uni[ens] = uni

anizip = url.urlopen('https://www.aniseed.cnrs.fr/aniseed/download/?file=data%2Fci%2FKH-ENS.blast.zip')
kh_to_uni = {}
with zipfile.ZipFile(BytesIO(anizip.read())) as zipdir:
    with zipdir.open('KH-ENS.blast', 'r') as fi:
        for line in fi.readlines():
            line = line.decode()
            ens, kh = line.split()[:2]
            kh = kh.split(':')[1]
            try:
                uni = ens_to_uni[ens]
                kh = mdl.scr_to_gen(kh)
                kh_to_uni[kh] = kh_to_uni.get(kh, set())
                kh_to_uni[kh].add(uni)
            except KeyError:
                pass

gzipf = url.urlopen('https://www.aniseed.cnrs.fr/aniseed/download/?file=data%2Fci%2FCirobu_KHslimTunicate.gaf.gz')
kh_to_name = {}
with gzip.open(gzipf, 'rt') as fi:
    for line in fi.readlines()[5:]:
        items = line.strip('\n').split('\t')
        name, kh = items[2], items[5]
        name = name.split('; ')
        kh = kh.split(':')[1]
        try:
            uni = kh_to_uni[kh]
            kh_to_name[kh] = kh_to_name.get(kh, set())
            kh_to_name[kh] = kh_to_name[kh].union(name)
        except KeyError:
            pass

out_str = ''

for kh in kh_to_name.keys():
    name = '|'.join(kh_to_name[kh])
    uni = '|'.join(kh_to_uni[kh])
    kh_out = '\t'.join([kh, name, uni]) + '\n'
    out_str += kh_out


with gzip.open('/home/jeff/Desktop/Thesis/Output/Ciona Robusta.gene2accession.txt.gz', 'wt') as fo:
    fo.write('UniqueID#EntrezGeneID\tSymbolID\tAccessionID\n')
    fo.write(out_str)