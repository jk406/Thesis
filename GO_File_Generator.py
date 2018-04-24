import gzip
import zipfile
from Modules import Module_2 as mdl
import urllib.request as url
from io import BytesIO


gzfile = url.urlopen('http://uswest.ensembl.org/biomart/martresults/18?file=martquery_0321234627_448.txt.gz')
ens_dict = {}
with gzip.open(gzfile, 'rt') as fi:
    for line in fi.readlines()[1:]:
        items = line.strip('\n').split('\t')
        ens_gen, ens_trn, ens_pro, ent, uni, name = items
        ens_dict[ens_trn] = ens_dict.get(ens_trn, {'trn': {ens_trn}, 'gen': set(), 'pro': set(), 'uni': set(), 'ent': set(), 'name': set(), 'GO': set()})
        if ens_gen:
            ens_dict[ens_trn]['gen'].add(ens_gen)
        if ens_pro:
            ens_dict[ens_trn]['pro'].add(ens_pro)
        if ent:
            ens_dict[ens_trn]['ent'].add(ent)
        if uni:
            ens_dict[ens_trn]['uni'].add(uni)
        if name:
            ens_dict[ens_trn]['name'].add(name)


anizip = url.urlopen('https://www.aniseed.cnrs.fr/aniseed/download/?file=data%2Fci%2FKH-ENS.blast.zip')
kh_dict = {}
with zipfile.ZipFile(BytesIO(anizip.read())) as zipdir:
    with zipdir.open('KH-ENS.blast', 'r') as fi:
        for line in fi.readlines():
            line = line.decode()
            ens, kh = line.split()[:2]
            kh = mdl.scr_to_gen(kh.split(':')[1])
            try:
                kh_dict[kh] = kh_dict.get(kh, ens_dict[ens])
                for key, val in ens_dict[ens].items():
                    kh_dict[kh][key] = kh_dict[kh][key] | val

            except KeyError:
                pass

gzipf = url.urlopen('https://www.aniseed.cnrs.fr/aniseed/download/?file=data%2Fci%2FCirobu_KHslimTunicate.gaf.gz')
with gzip.open(gzipf, 'rt') as fi:
    for line in fi.readlines()[5:]:
        items = line.strip('\n').split('\t')
        name, GO, kh = items[2], items[4], items[5]
        name = name.split('; ')
        kh = kh.split(':')[1]
        try:
            kh_dict[kh]['name'] = kh_dict[kh]['name'] | set(name)
            kh_dict[kh]['GO'].add(GO)
        except KeyError:
            kh_dict[kh] = {'trn': set(), 'gen': set(), 'pro': set(), 'uni': set(), 'ent': set(), 'name': set(name), 'GO': {GO}}

fo_acc = gzip.open('/home/jeff/Desktop/Thesis/Output/Ciona Intestinalis.gene2accession.txt.gz', 'wt')
fo_acc.write('UniqueID#EntrezGeneID\tKHID\tSymbolID\tAccessionID\tEntrezGeneID\tEnsemblGeneID\tEnsemblTranscriptID\tEnsemblProteinID\tUniProtKB_AC\tGOAnnotation\n')
fo_ens = gzip.open('/home/jeff/Desktop/Thesis/Output/Ciona Intestinalis.gene2ensembl.txt.gz', 'wt')
fo_ens.write('UniqueID#EntrezGeneID\tEnsemblGeneID\tEnsemblTranscriptID\tEnsemblProteinID\n')
fo_uni = gzip.open('/home/jeff/Desktop/Thesis/Output/Ciona Intestinalis.gene2uniprot.txt.gz', 'wt')
fo_uni.write('UniqueID#EntrezGeneID\tUniProtKB_AC\n')
for kh, ens_dict in kh_dict.items():
    out = '{0}\t{0}'.format(kh)
    for key in ['name', 'uni', 'ent', 'gen', 'trn', 'pro', 'uni', 'GO']:
        val = ens_dict[key]
        if val:
            val = '|'.join(val)
        else:
            val = '-'
        out += '\t{0}'.format(val)
    fo_acc.write(out + '\n')

    out = kh
    for key in ['gen', 'trn', 'pro']:
        val = ens_dict[key]
        if val:
            val = '|'.join(val)
        else:
            val = '-'
        out += '\t{0}'.format(val)
    fo_ens.write(out + '\n')

    out = kh
    val = ens_dict['uni']
    if val:
        val = '|'.join(val)
    else:
        val = '-'
    out += '\t{0}'.format(val)
    fo_uni.write(out + '\n')

fo_acc.close()
fo_ens.close()
fo_uni.close()