import pickle as pkl
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.mlab as mlab
from Modules import Module_2 as mdl

with open('/home/jeff/Desktop/Thesis/Src_Files/obj/Master.pkl', 'rb') as fi:
    Master = pkl.load(fi)

scr_to_chr = mdl.scr_to_chr
scr_to_gen = mdl.scr_to_gen

gen_expr = mdl.gen_expr()
ox_set = gen_expr['ox']
ux_set = gen_expr['ux']


label_dict = {'ox': 'Up-regulated', 'ux': 'Down-regulated', 'nc': 'Non-regulated', 'all': 'All', 'upgap': 'Upstream Gap', 'mRNA': 'mRNA', 'five_prime_UTR': "5' UTR", 'intron': 'Intron', 'CDS': 'CDS', 'three_prime_UTR': "3' UTR", 'downgap': 'Downstream Gap'}

motif = 'all'
cons = ['MAM', 'YASS']

leg_list, sam_com, name_com, bin_com = [], [], [], []
typ_sort = {}
for con in ['MAM']:
    typ_sort[con] = {}
    typ_sort[con]['all'] = {}
    with open('/home/jeff/Desktop/Thesis/Src_Files/obj/Count_Dict_Con_Det_{0}.pkl'.format(con), 'rb') as fi:
        Count_Dict = pkl.load(fi)
    di_list = []
    for scr, typ_dict in Count_Dict['all'].items():
        if scr in ox_set:
            reg = 'ox'
        elif scr in ux_set:
            reg = 'ux'
        else:
            reg = 'nc'
        typ_sort[con][reg] = typ_sort[con].get(reg, {})
        for typ in typ_dict.keys():
            typ_sort[con]['all'][typ] = typ_sort[con]['all'].get(typ, list())
            typ_sort[con]['all']['all'] = typ_sort[con]['all'].get('all', list())
            typ_sort[con][reg][typ] = typ_sort[con][reg].get(typ, list())
            typ_sort[con][reg]['all'] = typ_sort[con][reg].get('all', list())
            pc, d0 = 0, 0
            for lot in Master[scr_to_chr(scr)][scr_to_gen(scr)][scr][typ]:
                lo, up = lot
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
                    di_list.append(di)
                else:
                    di = 2970
                    d0 = 2970
                    #di = 0
                pc += np.exp(- (0.5 + 4 * (di / d0)))
            typ_sort[con]['all'][typ].append(pc)
            typ_sort[con]['all']['all'].append(pc)
            typ_sort[con][reg][typ].append(pc)
            typ_sort[con][reg]['all'].append(pc)

    sample, sam_name = [], []
    #reg = 'all'
    #for typ in [typ for typ in typ_sort[con][reg].keys() if typ != 'all' and typ != 'mRNA']:
    #    sam_name.append(label_dict[typ])
    typ = 'intron'
    for reg in [reg for reg in typ_sort[con].keys() if reg != 'all']:
        sam_name.append(label_dict[reg])
        x = typ_sort[con][reg][typ]
        sample.append(x)

    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
    n, bins, patches = plt.hist(sample, alpha=0.5, label=sam_name, density=True, color=colors[:len(sam_name)])
    mu = [np.mean(x_i) for x_i in sample]
    sigma = [np.std(x_i) for x_i in sample]
    for i in range(len(mu)):
        y = mlab.normpdf(bins, mu[i], sigma[i])
        l = plt.plot(bins, y, '--', linewidth=1, color=colors[i])

    plt.legend()
    #plt.axis([0, 2, 0, 5])
    plt.xlabel('Peak Concentration')
    plt.ylabel('Frequency of Genes')
    plt.grid(True)
    plt.savefig('/home/jeff/Desktop/Thesis/Output/intron.png')
    plt.show()
    plt.close()




"""
    plt.figure(2)
    n, bins, patches = plt.hist(sam_com, alpha=0.5, label=name_com, bins=max(bin_com), normed=1)
    mu = [np.mean(x_i) for x_i in sam_com]
    sigma = [np.std(x_i) for x_i in sam_com]
    for i in range(3):
        y = mlab.normpdf(bins, mu[i], sigma[i])
        plt.plot(bins, y, '--', linewidth=1, label='{0} Trendline'.format(cons[i]))
    plt.legend()
    plt.axis([0, 0.01, 0, 800])
    plt.xlabel('Count of Binding Sites / Length of Region')
    plt.ylabel('Frequency of Genes')
    plt.grid(True)
    plt.savefig('/home/jeff/Desktop/Thesis/Output/{0}_histo_trend.png'.format(typ))
    plt.close()
"""