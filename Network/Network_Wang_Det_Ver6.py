import Modules.Module_2 as mdl
import pickle as pkl
import numpy as np


RNA_Seq = mdl.rna_seq()
# when di = 0
#th = 0.01
# when di = 1
th = 0.009
#th = 1

with open('/home/jeff/Desktop/Thesis/Src_Files/obj/tf.pkl', 'rb') as fi:
    tf = pkl.load(fi)

with open('/home/jeff/Desktop/Thesis/Src_Files/obj/Master.pkl', 'rb') as fi:
    Master = pkl.load(fi)

with open('/home/jeff/Desktop/Thesis/Src_Files/obj/notochord.pkl', 'rb') as fi:
    nc = pkl.load(fi)

#with open('/home/jeff/Desktop/Thesis/Src_Files/obj/ID_Cnvrt.pkl', 'rb') as fi:
#    ID_Cnvrt = pkl.load(fi)

# set of overexpressed transcripts
gen_expr = mdl.gen_expr()
ox_set = gen_expr['ox']
ux_set = gen_expr['ux']

fo_tf = open('/home/jeff/Desktop/Thesis/Output/TF_list_{0}.csv'.format(th), 'w')
fo_nc = open('/home/jeff/Desktop/Thesis/Output/NC_list_{0}.csv'.format(th), 'w')

scr_to_gen = mdl.scr_to_gen
scr_to_chr = mdl.scr_to_chr
kh_to_name = mdl.kh_to_name

for con in [con for con in mdl.cons if con != 'raw']:
    with open('/home/jeff/Desktop/Thesis/Src_Files/obj/Count_Dict_Con_Det_{0}.pkl'.format(con), 'rb') as fi:
        Count_Dict = pkl.load(fi)
    motif = 'all'
    with open('/home/jeff/Desktop/Thesis/Output/Network_ox_{0}_{1}_name.txt'.format(con, motif), 'w') as fo:
        fo.write('source\ttarget\tinteraction\tRPScore\n')
        # Calculation of PC
        pc_list, de_list, pc_dict, de_dict, tf_list = {}, [], {}, {}, []
        for scr in ox_set:
            pcs = {}
            for typ in ['upgap', 'intron', 'downgap']:
                pc, d0 = 0, 0
                try:
                    for reg in Master[scr_to_chr(scr)][scr_to_gen(scr)][scr][typ]:
                        lo, up = reg
                        d0 += up - lo + 1
                    if d0 > 3000:
                        d0 = 3000
                    for ord, sites in Count_Dict[motif][scr][typ].items():
                        up_end = []
                        for site in sites:
                            lo, up = site
                            up_end.append(up)
                        if typ == 'upgap':
                            di = 1 - max(up_end)
                        else:
                            di = 1
                        pc += np.exp(-(0.5 + 4 * (di / d0)))
                except KeyError:
                    pass
                pcs[typ] = pc
            if pcs['upgap'] > 0:
                for typ in pcs:
                    pc_list[typ] = pc_list.get(typ, [])
                    pc_list[typ].append(pcs[typ])
                pc_dict[scr] = pcs

        # Calculation of DE
        for scr in pc_dict.keys():
            FDR = RNA_Seq[scr]['FDR']
            PVa = RNA_Seq[scr]['PValue']
            de = FDR/PVa
            de_list.append(de)
            de_dict[scr] = de

        # Calculation of RP
        n = len(pc_dict.keys())
        for key, value in pc_list.items():
            pc_list[key] = sorted(value)
        de_list, rp_dict = sorted(de_list), {}
        for scr in pc_dict.keys():
            Rpc_up, Rpc_in, Rpc_do, Rde = \
                n - pc_list['upgap'].index(pc_dict[scr]['upgap']), n - pc_list['intron'].index(pc_dict[scr]['intron']), \
                n - pc_list['downgap'].index(pc_dict[scr]['downgap']), 1 + de_list.index(de_dict[scr])
            RP = (Rpc_up / n) * (Rpc_in / n) * (Rpc_do / n) * (Rde / n)
            gen = scr_to_gen(scr)
            rp_dict[gen] = rp_dict.get(gen, [])
            rp_dict[gen].append(RP)
        # Calculation of average score of each gene
        score_dict = {}
        for gen, RP in rp_dict.items():
            RP = np.mean(RP)
            if RP <= th:
                score = 1 - RP
                score_dict[gen] = score

        # Writing the output file
        for gen, score in score_dict.items():
            #gen = ''.join(ID_Cnvrt['KH'].get(gen, {'name': gen}).get('name', gen)).strip(' ')
            fo.write('Pou4\t{0}\tox\t{1}'.format(kh_to_name(gen), score))
            fo.write('\n')
            if gen in tf:
                tf_list.append(gen)
            if gen in nc:
                fo_nc.write('{0}\t{1}\t{2}\t{3}\n'.format(con, gen, kh_to_name(gen), 1 - score_dict[gen]))

        fo_tf.write('{0}\t{1}\tox\t{2}\n'.format(con, motif, str(len(tf_list))))
        for trans in tf_list:
            fo_tf.write('{0}\t{1}\t{2}\n'.format(trans, kh_to_name(trans), score_dict[trans]))
    Count_Dict = 0

fo_tf.close()
fo_nc.close()