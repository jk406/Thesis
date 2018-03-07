from Bio import SeqIO
from Bio import motifs
from Bio.Alphabet import IUPAC
import os
import tarfile
from Modules import Module_2 as mdl

chrms = SeqIO.parse('/home/jeff/Desktop/Thesis/Src_Files/JoinedScaffold.fasta', 'fasta', IUPAC.unambiguous_dna)

with open('/home/jeff/Desktop/Thesis/Src_Files/GeneModel_Name_mapping.csv') as fi:
    ani_to_kh = {}
    for line in fi.readlines()[1:]:
        kh, ani = line.split(';')[:2]
        kh = kh.split(':')[1]
        ani = ani.split('.')[1]
        ani_to_kh[ani] = kh

pwm_path = '/home/jeff/Desktop/Thesis/Src_Files/Selex_seq_Cirobu_All/'

with open(pwm_path + '/Selex_seq_Cirobu_All.csv') as list_fi:
    name_dict = {}
    lines = list_fi.readlines()
    for line in lines:
        ani, ac = line.split(';')[:2]
        ani = ani.split('.')[1]
        if ani in ani_to_kh.keys():
            name_dict[ac] = ani_to_kh[ani]

pwm_path += 'pwms/'

tf_set = mdl.get_tf_dict()['MAM']['all']['ox']

name_list, list_list = [], []
for tar_file in os.listdir(pwm_path):
    ani = tar_file.split('-')[0]
    kh = name_dict.get(ani, '')
    if kh == 'KH.C12.675':
        name_list.append(kh)
        mot_name = tar_file.split('-')[2].split('.')[0]
        tar_file = tarfile.open(pwm_path + tar_file)
        thre_list, len_list, mot_list, pssm_list = [], [], [], []
        for tar_info in tar_file.getmembers():
            mot_name_2 = mot_name + tar_info.name.split('.')[0][len(ani):]
            mot_list.append(mot_name_2)
            pfm_file = tar_file.extractfile(tar_info)
            motif = motifs.read(pfm_file, "pfm")
            pssm = motif.pssm
            distribution = pssm.distribution()
            threshold = distribution.threshold_fpr(0.0001)
            thre_list.append(threshold)
            len_list.append(len(motif.consensus))
            pssm_list.append(pssm)
        list_list.append([mot_list, thre_list, len_list, pssm_list])
        tar_file.close()

content_list = ['' for x in range(len(name_list))]
for x in range(len(name_list)):
    kh = name_list[x]
    out = ''
    for chrm in chrms:
        if chrm.id == 'KhC1':
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
                    out += '{0}\t{1}\t{2}\t{3}\t{6}\t{4}\t{5}\n'.format(
                        moti, chrm.id, pos, pos + leng - 1, score, seq, stra)
    with open('/home/jeff/Desktop/Thesis/Output/TFBS/TFBS_{0}_KhC1.csv'.format(kh), 'w') as fo:
        fo.write(out)
