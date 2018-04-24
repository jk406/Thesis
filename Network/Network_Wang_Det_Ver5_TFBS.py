import Modules.Module_2 as mdl
import pickle as pkl
import numpy as np


RNA_Seq = mdl.rna_seq()
th = .88


with open('/home/jeff/Desktop/Thesis/Src_Files/ciona_TF.csv', 'r') as fi:
    tf = set()
    for line in fi.readlines()[4:]:
        items = line.split('\t')
        if len(items) > 1:
            tf.add(items[2])

with open('/home/jeff/Desktop/Thesis/Src_Files/obj/Parser.pkl', 'rb') as fi:
    Parser = pkl.load(fi)

#with open('/home/jeff/Desktop/Thesis/Src_Files/obj/ID_Cnvrt.pkl', 'rb') as fi:
#    ID_Cnvrt = pkl.load(fi)

# set of overexpressed transcripts
gen_expr = mdl.gen_expr()
ox_set = gen_expr['ox']
ux_set = gen_expr['ux']

fo_tf = open('/home/jeff/Desktop/Thesis/Output/TFBS/TF_list_C12_675_{0}.csv'.format(th), 'w')

scr_to_chr = mdl.scr_to_chr()
scr_to_gen = mdl.scr_to_gen()

for con in [con for con in mdl.cons if con != 'raw']:
    with open('/home/jeff/Desktop/Thesis/Src_Files/TFBS/obj/Count_Dict_Con_Det_{0}.pkl'.format(con), 'rb') as fi:
        Count_Dict = pkl.load(fi)

    motif = 'all'
    with open('/home/jeff/Desktop/Thesis/Output/Network/Network_ox_{0}_{1}.txt'.format(con, motif), 'a') as fo:
        #fo.write('source\ttarget\tinteraction\tRPScore\n')
        # Calculation of PC
        tru, pc_list, de_list, pc_dict, de_dict, tf_list = [], [], [], {}, {}, []
        for scr in ox_set:
            pcs, pc, d0 = [], 0, 0
            typ = 'upgap'
            try:
                for reg in Parser[scr_to_chr(scr)][scr_to_gen(scr)][typ][scr]:
                    lo, up = reg
                    d0 += up - lo + 1
                for ord, sites in Count_Dict[motif][scr][typ].items():
                    up_end = []
                    for site in sites:
                        lo, up = site
                        up_end.append(up)
                    di = 1 - max(up_end)
                    pcs.append(np.exp(- (0.5 + 4 * (di / d0))))
                pc += sum(pcs)
            except KeyError:
                pc += 0

            if pc > 0:
                tru.append(scr)
                pc_list.append(pc)
                pc_dict[scr] = pc

        # Calculation of DE
        for scr in tru:
            FDR = RNA_Seq[scr]['FDR']
            PVa = RNA_Seq[scr]['PValue']
            de = FDR/PVa
            de_list.append(de)
            de_dict[scr] = de

        # Calculation of RP
        n = len(tru)
        pc_list, de_list, score_dict = sorted(pc_list), sorted(de_list), {}
        for scr in tru:
            Rpc, Rde = n - pc_list.index(pc_dict[scr]), 1 + de_list.index(de_dict[scr])
            RP = (Rpc / n) * (Rde / n)
            if RP <= th:
                gen = scr_to_gen(scr)
                score_dict[gen] = score_dict.get(gen, [])
                score_dict[gen].append(RP)

        # Calculation of average scor of each gene
        for gen, RP in score_dict.items():
            RP = np.mean(RP)
            score = 1 - RP
            score_dict[gen] = score

        # Writing the output file
        for gen, score in score_dict.items():
            #name = ''.join(ID_Cnvrt['KH'].get(gen, {'name': gen}).get('name', gen)).strip(' ')
            #if name != gen:
            fo.write('KH.C12.675\t{0}\tox\t{1}'.format(gen, score))
            fo.write('\n')
            if gen in tf:
                tf_list.append(gen)

        fo_tf.write('{0}\t{1}\tox\t{2}\t{3}\n'.format(con, motif, str(len(tf_list)), '\t'.join(tf_list)))

    with open('/home/jeff/Desktop/Thesis/Output/Network/Network_ux_{0}_{1}.txt'.format(con, motif), 'a') as fo:
        fo.write('source\ttarget\tinteraction\tRPScore\n')
        tru, pc_list, de_list, pc_dict, de_dict, tf_list = [], [], [], {}, {}, []
        for scr in ux_set:
            d0, pcs, pc = 0, [], 0
            for typ in ['upgap']:
                try:
                    for reg in Parser[scr_to_chr(scr)][scr_to_gen(scr)][typ][scr]:
                        lo, up = reg
                        d0 += up - lo + 1
                    if typ == 'upgap':
                        for ord, sites in Count_Dict[motif][scr][typ].items():
                            up_end = []
                            for site in sites:
                                lo, up = site
                                up_end.append(up)
                            di = 1 - max(up_end)
                            pcs.append(np.exp(- (0.5 + 4 * (di / d0))))
                        pc += sum(pcs)
                    else:
                        for ord, sites in Count_Dict[motif][scr][typ].items():
                            lo_end = []
                            for site in sites:
                                lo, up = site
                                lo_end.append(lo)
                            di = min(lo_end)
                            pcs.append(np.exp(- (0.5 + 4 * (di / d0))))
                        pc += sum(pcs)
                except KeyError:
                    pc += 0

            if pc > 0:
                tru.append(scr)
                pc_list.append(pc)
                pc_dict[scr] = pc

        # Calculation of DE
        for scr in tru:
            FDR = RNA_Seq[scr]['FDR']
            PVa = RNA_Seq[scr]['PValue']
            de = FDR / PVa
            de_list.append(de)
            de_dict[scr] = de

        # Calculation of RP
        n = len(tru)
        pc_list, de_list, score_dict = sorted(pc_list), sorted(de_list), {}
        for scr in tru:
            Rpc, Rde = n - pc_list.index(pc_dict[scr]), 1 + de_list.index(de_dict[scr])
            RP = (Rpc / n) * (Rde / n)
            if RP <= th:
                gen = scr_to_gen(scr)
                score_dict[gen] = score_dict.get(gen, [])
                score_dict[gen].append(RP)

        # Calculation of average score of each gene
        for gen, RP in score_dict.items():
            RP = np.mean(RP)
            score = 1 - RP
            score_dict[gen] = score

        # Writing the output file
        for gen, score in score_dict.items():
            #name = ''.join(ID_Cnvrt['KH'].get(gen, {'name': gen}).get('name', gen)).strip(' ')
            fo.write('KH.C12.675\t{0}\tux\t{1}'.format(gen, score))
            fo.write('\n')
            if gen in tf:
                tf_list.append(gen)

        fo_tf.write('{0}\t{1}\tux\t{2}\t{3}\n'.format(con, motif, str(len(tf_list)), '\t'.join(tf_list)))

    Count_Dict = 0

fo_tf.close()