
# coding: utf-8

# # Parse the files and pickle the objects needed for later processes 

# import the necessary module

# In[1]:


import pickle
import urllib.request as url
import zipfile
from io import BytesIO
import tarfile


# store the lower and upper ends of each gene to find the range to look

# In[2]:

save_path = '/home/jeff/Desktop/Thesis/Output/obj/new/'

gffzip = url.urlopen('http://ghost.zool.kyoto-u.ac.jp/datas/KH.KHGene.2012.gff3.zip')
with zipfile.ZipFile(BytesIO(gffzip.read())) as zipdir:
    with zipdir.open('KH.KHGene.2012.gff3', 'r') as fi:
        lower_dict = {}
        upper_dict = {}
        lines = fi.readlines()
        for line in lines:
            words = line.split()
            chrna, typ = words[0], words[2]
            if typ == "gene":
                lo, up = int(words[3]), int(words[4])
                lower_dict[chrna], upper_dict[chrna] = lower_dict.get(chrna, list()), upper_dict.get(chrna, list())
                lower_dict[chrna].append(lo)
                upper_dict[chrna].append(up)
                lower_dict[chrna], upper_dict[chrna] = sorted(lower_dict[chrna]), sorted(upper_dict[chrna])

    # make the master dicitonary with coordinates of genes, mRNA's, CDS's, and UTR's
    with zipdir.open('KH.KHGene.2012.gff3', 'r') as fi:
        Master = {}
        Gen_Info = {}
        lines = fi.readlines()
        for line in lines:
            words = line.split()
            chrna, st_gen, en_gen, typ, stra = words[0], int(words[3]), int(words[4]), words[2], words[6]
            Master[chrna] = Master.get(chrna, dict())
            Gen_Info[chrna] = Gen_Info.get(chrna, dict())
            if typ == "gene":
                for item in words[8].split(";"):
                    if "name=" in item:
                        base = item[5:]
                # lower range end
                x = -1
                while st_gen < upper_dict[chrna][x]:
                    try:
                        lo_gen = upper_dict[chrna][x - 1]
                        x -= 1
                    except IndexError:
                        lo_gen = st_gen - 3000
                        if lo_gen < 0:
                            lo_gen = 0
                        break
                if st_gen - lo_gen > 3000:
                    lo_gen = st_gen - 3000
                # upper range end
                y = 0
                while en_gen > lower_dict[chrna][y]:
                    try:
                        up_gen = lower_dict[chrna][y + 1]
                        y += 1
                    except IndexError:
                        up_gen = en_gen + 3000
                        break
                if up_gen - en_gen > 3000:
                    up_gen = en_gen + 3000
                Gen_Info[chrna][base] = [stra, lo_gen, st_gen, en_gen, up_gen]
            else:
                if typ == "mRNA":
                    base, scrp = words[8].split(";")[0][14:], words[8].split(";")[1][10:]
                    Master[chrna][base] = Master[chrna].get(base, dict())
                    Master[chrna][base][scrp] = Master[chrna][base].get(scrp, dict())
                    Master[chrna][base][scrp][typ] = Master[chrna][base][scrp].get(typ, list())
                    Master[chrna][base][scrp][typ].append([st_gen, en_gen])
                elif typ == "CDS" or typ == "five_prime_UTR" or typ == "three_prime_UTR":
                    scrp = words[8].split(";")[0][14:]
                    base = '.'.join(scrp.split('.')[:3])
                    Master[chrna][base] = Master[chrna].get(base, dict())
                    Master[chrna][base][scrp] = Master[chrna][base].get(scrp, dict())
                    Master[chrna][base][scrp][typ] = Master[chrna][base][scrp].get(typ, list())
                    Master[chrna][base][scrp][typ].append([st_gen, en_gen])


# make coordinates of introns based on given coordinates

# In[4]:


for chrna in Master.keys():
    for base in Master[chrna].keys():
        for scrp in Master[chrna][base].keys():
            ran_loc = {}
            lo_scrp, up_scrp = Master[chrna][base][scrp]['mRNA'][0]
            for typ in Master[chrna][base][scrp].keys():
                if typ != 'mRNA':
                    for ran in Master[chrna][base][scrp][typ]:
                        ran_loc[ran[0]] = ran[1]
            Master[chrna][base][scrp]['intron'] = Master[chrna][base][scrp].get('intron', list())
            keys = sorted(ran_loc.keys())
            if keys[0] > lo_scrp:
                Master[chrna][base][scrp]['intron'].append([lo_scrp, keys[0] - 1])
            for x in range(len(keys)):
                try:
                    if keys[x + 1] - ran_loc[keys[x]] > 1:
                        Master[chrna][base][scrp]['intron'].append([ran_loc[keys[x]] + 1, keys[x + 1] - 1])
                except IndexError:
                    if ran_loc[keys[x]] < up_scrp:
                        Master[chrna][base][scrp]['intron'].append([ran_loc[keys[x]] + 1, up_scrp])


# make gaps based on the data given

# In[5]:


for chrna in Master.keys():
    for base in Master[chrna].keys():
        stra_gene, st_gene, en_gene = Gen_Info[chrna][base][0], Gen_Info[chrna][base][1], Gen_Info[chrna][base][4]
        for scrp in Master[chrna][base].keys():
            lo, up = Master[chrna][base][scrp]['mRNA'][0]
            Master[chrna][base][scrp]['upgap'] = Master[chrna][base][scrp].get('upgap', list())
            Master[chrna][base][scrp]['downgap'] = Master[chrna][base][scrp].get('downgap', list())
            if stra_gene == '+':
                Master[chrna][base][scrp]['upgap'] = [[st_gene, lo - 1]]
                Master[chrna][base][scrp]['downgap'] = [[up + 1, en_gene]]
            elif stra_gene == '-':
                Master[chrna][base][scrp]['downgap'] = [[st_gene, lo - 1]]
                Master[chrna][base][scrp]['upgap'] = [[up + 1, en_gene]]


# pickle the master dictionary

# In[8]:


with open(save_path + 'Master.pkl', 'wb') as fo:
    pickle.dump(Master, fo, pickle.HIGHEST_PROTOCOL)

with open(save_path + 'Gen_Info.pkl', 'wb') as fo:
    pickle.dump(Gen_Info, fo, pickle.HIGHEST_PROTOCOL)


# open the TFBS file and make the dictionary

# In[6]:


Bind_Sites = {}
num = {}
ord = {}
bs_txt = url.urlopen('https://sdsuedu-my.sharepoint.com/personal/jkim6_sdsu_edu/_layouts/15/download.aspx?SourceUrl=%2Fpersonal%2Fjkim6%5Fsdsu%5Fedu%2FDocuments%2FThesis%2FPou4GenomeOutput%2Etxt')
with open(bs_txt, 'r') as fi:
    lines = fi.readlines()
    num['all'], ord['all'] = {}, {}
    for line in lines:
        words = line.split()
        motif = '-'.join(words[0].split('-')[:3])
        chrna, lo, up = words[1], int(words[2]), int(words[3])
        num[motif], ord[motif] = num.get(motif, {}), ord.get(motif, {})
        num['all'][chrna], ord['all'][chrna] = num['all'].get(chrna, 0), ord['all'].get(chrna, 0)
        num[motif][chrna], ord[motif][chrna] = num[motif].get(chrna, 0), ord[motif].get(chrna, 0)

        Bind_Sites['all'] = Bind_Sites.get('all', dict())
        Bind_Sites['all'][chrna] = Bind_Sites['all'].get(chrna, dict())
        Bind_Sites['all'][chrna][ord['all'][chrna]] = Bind_Sites['all'][chrna].get(ord['all'][chrna], [[lo, up]])
        try:
            if Bind_Sites['all'][chrna][ord['all'][chrna]][num['all'][chrna]][0] <= lo <= Bind_Sites['all'][chrna][ord['all'][chrna]][num['all'][chrna]][1] <= up:
                if Bind_Sites['all'][chrna][ord['all'][chrna]][num['all'][chrna]] != [lo, up]:
                    Bind_Sites['all'][chrna][ord['all'][chrna]].append([lo, up])
                    num['all'][chrna] += 1
            elif lo <= Bind_Sites['all'][chrna][ord['all'][chrna]][num['all'][chrna]][0] <= up <= Bind_Sites['all'][chrna][ord['all'][chrna]][num['all'][chrna]][1]:
                if Bind_Sites['all'][chrna][ord['all'][chrna]][num['all'][chrna]] != [lo, up]:
                    Bind_Sites['all'][chrna][ord['all'][chrna]].append([lo, up])
                    num['all'][chrna] += 1
            elif lo == Bind_Sites['all'][chrna][ord['all'][chrna]][num['all'][chrna]][1] + 1:
                Bind_Sites['all'][chrna][ord['all'][chrna]].append([lo, up])
                num['all'][chrna] += 1
            else:
                ord['all'][chrna] += 1
                num['all'][chrna] = 0
                Bind_Sites['all'][chrna][ord['all'][chrna]] = [[lo, up]]
        except IndexError:
            pass

        Bind_Sites[motif] = Bind_Sites.get(motif, dict())
        Bind_Sites[motif][chrna] = Bind_Sites[motif].get(chrna, dict())
        Bind_Sites[motif][chrna][ord[motif][chrna]] = Bind_Sites[motif][chrna].get(ord[motif][chrna], [[lo, up]])
        try:
            if Bind_Sites[motif][chrna][ord[motif][chrna]][num[motif][chrna]][0] <= lo <= Bind_Sites[motif][chrna][ord[motif][chrna]][num[motif][chrna]][1] <= up:
                if Bind_Sites[motif][chrna][ord[motif][chrna]][num[motif][chrna]] != [lo, up]:
                    Bind_Sites[motif][chrna][ord[motif][chrna]].append([lo, up])
                    num[motif][chrna] += 1
            elif lo <= Bind_Sites[motif][chrna][ord[motif][chrna]][num[motif][chrna]][0] <= up <= Bind_Sites[motif][chrna][ord[motif][chrna]][num[motif][chrna]][1]:
                if Bind_Sites[motif][chrna][ord[motif][chrna]][num[motif][chrna]] != [lo, up]:
                    Bind_Sites[motif][chrna][ord[motif][chrna]].append([lo, up])
                    num[motif][chrna] += 1
            elif lo == Bind_Sites[motif][chrna][ord[motif][chrna]][num[motif][chrna]][1] + 1:
                Bind_Sites[motif][chrna][ord[motif][chrna]].append([lo, up])
                num[motif][chrna] += 1
            else:
                ord[motif][chrna] += 1
                num[motif][chrna] = 0
                Bind_Sites[motif][chrna][ord[motif][chrna]] = [[lo, up]]
        except IndexError:
            pass


# pickle the bind site dicitonary

# In[8]:


with open(save_path + 'Bind_Sites.pkl', 'wb') as fo:
    pickle.dump(Bind_Sites, fo, pickle.HIGHEST_PROTOCOL)


# parse the conserved reigon data from the csv file

# In[6]:


con_list = []
# makes a dictionary of each conservative-discovering method
cons_gz = url.urlopen('https://sdsuedu-my.sharepoint.com/personal/jkim6_sdsu_edu/_layouts/15/download.aspx?SourceUrl=%2Fpersonal%2Fjkim6%5Fsdsu%5Fedu%2FDocuments%2FThesis%2FCons%2Etar%2Egz')

with tarfile.open(cons_gz, 'r:gz') as dir_co:
    con_name = dir_co.getnames()
    for con_meth in dir_co:
        con_reg = {}
        chr_con = set()
        with open(con_meth, 'r') as fi_co:
            lines_co = fi_co.readlines()
            num = 1
            tot = len(lines_co[1:])
            for line_co in lines_co[1:]:
                words_co = line_co.split(',')
                chrna_co = words_co[1]
                chr_con.add(chrna_co)
                co_st = int(words_co[2])
                co_en = co_st + int(words_co[3]) - 1
                con_reg[chrna_co] = con_reg.get(chrna_co, dict())
                con_reg[chrna_co][co_st] = co_en
                num += 1
            con_list.append(con_reg)


# organize the conserved regions by each gene and make a dictionary of them

# In[ ]:


for num in range(2):
    con_reg = con_list[num]
    Cons = {}
    for chrna in Master.keys():
        if chrna in con_reg.keys():
            Cons[chrna] = {}
            for gene in Master[chrna].keys():
                Cons[chrna][gene] = {}
                for scr in Master[chrna][gene].keys():
                    for typ in Master[chrna][gene][scr].keys():
                        for raw in Master[chrna][gene][scr][typ]:
                            for reg in con_reg[chrna].keys():
                                end = con_reg[chrna][reg]
                                if reg <= raw[0] and raw[1] <= end:
                                    Cons[chrna][gene][scr] = Cons[chrna][gene].get(scr, dict())
                                    Cons[chrna][gene][scr][typ] = Cons[chrna][gene][scr].get(typ, list())
                                    Cons[chrna][gene][scr][typ].append([raw[0], raw[1]])
                                elif reg <= raw[0]:
                                    if raw[0] <= end <= raw[1]:
                                        Cons[chrna][gene][scr] = Cons[chrna][gene].get(scr, dict())
                                        Cons[chrna][gene][scr][typ] = Cons[chrna][gene][scr].get(typ, list())
                                        Cons[chrna][gene][scr][typ].append([raw[0], end])
                                elif raw[1] <= end:
                                    if raw[0] <= reg <= raw[1]:
                                        Cons[chrna][gene][scr] = Cons[chrna][gene].get(scr, dict())
                                        Cons[chrna][gene][scr][typ] = Cons[chrna][gene][scr].get(typ, list())
                                        Cons[chrna][gene][scr][typ].append([reg, raw[1]])
                                else:
                                    if raw[0] <= reg and end <= raw[1]:
                                        Cons[chrna][gene][scr] = Cons[chrna][gene].get(scr, dict())
                                        Cons[chrna][gene][scr][typ] = Cons[chrna][gene][scr].get(typ, list())
                                        Cons[chrna][gene][scr][typ].append([reg, end])
        name = con_name[num]
        with open(save_path + 'Cons_{}.pkl'.format(name), 'wb') as fo:
            pickle.dump(Cons, fo, pickle.HIGHEST_PROTOCOL)

