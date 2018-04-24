import numpy as np
import pickle as pkl

with open('/home/jeff/Desktop/Thesis/Src_Files/obj/Count_Dict.pkl', 'rb') as fi:
    # Count_Dict[version][rna-seq][motif][transcript_id][type][window] = Count
    Count_Dict = pkl.load(fi)
    ox_set = set(Count_Dict['19']['ox']['all'].keys())
    ux_set = set(Count_Dict['19']['ux']['all'].keys())
    scrs = set(Count_Dict['19']['all']['all'].keys())

with open('/home/jeff/Desktop/Thesis/Src_Files/obj/rna_seq.pkl', 'rb') as fi:
    # rna_seq[transcript][PValue/FDR] = Value
    rna_seq = pkl.load(fi)

with open('/home/jeff/Desktop/Thesis/Src_Files/obj/gen_name.pkl', 'rb') as fi:
    # gen_name[transcript] = transcript name
    gen_name = pkl.load(fi)

for ver in ['20', '19_Con']:
    for reg in ['all', 'ox']:
        gen_pc = {}
        for scr in scrs:
            if scr in rna_seq:
                pc = 0
                d0 = 2000
                for win in range(-2000, 0, 500):
                    try:
                        di, pi = -win - 500, Count_Dict[ver]['all']['all'][scr]['upgap'][win]
                        pc += pi*np.exp(-(0.5 + 4*di/d0))
                    except KeyError:
                        d0 -= 500
                gen_pc[scr] = pc

        tru = []
        pc_lis = []
        de_lis = []
        for scr, pc in gen_pc.items():
            pc_lis.append(gen_pc[scr])
            de_lis.append(rna_seq[scr]['FDR'] / rna_seq[scr]['PValue'])
            if reg == 'all':
                if scr in ox_set or scr in ux_set:
                    if pc > 0:
                        tru.append(scr)
            else:
                if scr in ox_set:
                    if pc > 0:
                        tru.append(scr)

        gen_rp = {}
        for scr in scrs:
            if scr in rna_seq:
                Rpc = len(pc_lis) - sorted(pc_lis).index(gen_pc[scr])
                Rde = sorted(de_lis).index(rna_seq[scr]['FDR']/rna_seq[scr]['PValue']) + 1
                RP = (Rpc/len(tru))*(Rde/len(tru))
                gen_rp[scr] = RP

        with open('/home/jeff/Desktop/Thesis/Output/drctrltn_{1}_{0}.txt'.format(ver, reg), 'w') as fo:
            fo.write('source\ttarget\tinteraction\tscore\n')
            for scr, rp in sorted(gen_rp.items()):
                if rp <= 1:
                    if scr in ox_set:
                        pc = gen_pc[scr]
                        scr = gen_name.get(scr, scr)
                        fo.write('Pou4\t{0}\tpd\t{1}\n'.format(scr, str(pc)))