from Bio import SeqIO
from Bio import motifs
from Bio.Alphabet import IUPAC
import tarfile
import urllib.request as url
import zipfile
from io import BytesIO
import os

#fasta = url.urlopen('http://ghost.zool.kyoto-u.ac.jp/datas/JoinedScaffold.zip')
#with zipfile.ZipFile(BytesIO(fasta.read())) as zip:
with zipfile.ZipFile('/home/jeff/Downloads/JoinedScaffold.zip') as zip:
    fi = zip.extract('JoinedScaffold')
    chrms = SeqIO.parse(fi, 'fasta', IUPAC.unambiguous_dna)

#gene_map = url.urlopen('https://www.aniseed.cnrs.fr/aniseed/download/?file=data%2Fci%2FGeneModel_bbh_mapping.zip')
#with zipfile.ZipFile(BytesIO(gene_map.read())) as zip:
with zipfile.ZipFile('/home/jeff/Downloads/GeneModel_bbh_mapping.zip') as zip:
    with zip.open('GeneModel_Name_mapping.csv', 'r') as fi:
        ani_to_kh = {}
        for line in fi.readlines()[1:]:
            line = line.decode()
            kh, ani = line.split(';')[:2]
            kh = kh.split(':')[1]
            ani = ani.split('.')[1]
            ani_to_kh[ani] = kh

#selex = url.urlopen('https://www.aniseed.cnrs.fr/aniseed/download/?file=data%2Fci%2FSelex_seq_Cirobu_All.zip')
#with zipfile.ZipFile(BytesIO(selex.read())) as zip:
with zipfile.ZipFile('/home/jeff/Downloads/Selex_seq_Cirobu_All.zip') as zip:
    with zip.open('Selex_seq_Cirobu_All.csv', 'r') as list_fi:
        name_dict = {}
        lines = list_fi.readlines()
        for line in lines:
            line = line.decode()
            ani, ac = line.split(';')[:2]
            ani = ani.split('.')[1]
            if ani in ani_to_kh.keys():
                name_dict[ac] = ani_to_kh[ani]

    for fi_name in zip.namelist():
        if 'pwms/' in fi_name:
            ani = fi_name.split('/')[1].split('-')[0]
            kh = name_dict.get(ani, '')
            if kh == 'KH.C8.175':
                break

    name_list, list_list = [], []
    fi = zip.extract(fi_name)
    with tarfile.open(fi, 'r:gz') as fi:
        name_list.append(kh)
        mot_name = fi_name.split('-')[2].split('.')[0]
        thre_list, len_list, mot_list, pssm_list = [], [], [], []
        for tar_info in fi.getmembers():
            mot_name_2 = mot_name + tar_info.name.split('.')[0][len(ani):]
            mot_list.append(mot_name_2)
            pfm_file = fi.extractfile(tar_info)
            motif = motifs.read(pfm_file, "pfm")
            pssm = motif.pssm
            distribution = pssm.distribution()
            threshold = distribution.threshold_fpr(0.0001)
            thre_list.append(threshold)
            len_list.append(len(motif.consensus))
            pssm_list.append(pssm)
        list_list.append([mot_list, thre_list, len_list, pssm_list])

os.remove(fi_name)

print('almost done')

with open('/home/jeff/Desktop/Thesis/Output/TFBS/TFBS_{0}_KhC6.csv'.format(kh), 'w') as fo:
    content_list = ['' for x in range(len(name_list))]
    for x in range(len(name_list)):
        kh = name_list[x]
        for chrm in chrms:
            if chrm.id == 'KhC6':
                for y in range(len(list_list[x][0])):
                    moti = list_list[x][0][y]
                    threshold = list_list[x][1][y]
                    leng = list_list[x][2][y]
                    pssm = list_list[x][3][y]
                    for pos, score in pssm.search(chrm.seq, threshold=threshold):
                        seq = chrm.seq[pos:pos + leng]
                        if pos >= 0:
                            stra = '+'
                        else:
                            stra = '-'
                            pos = len(chrm.seq) + pos + 1
                        fo.write('{0}\t{1}\t{2}\t{3}\t{6}\t{4}\t{5}\n'.format(
                            moti, chrm.id, pos, pos + leng - 1, score, seq, stra))
                break