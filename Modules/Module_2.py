import pickle as pkl

motifs = ['Pou4F1-SELEX-for', 'Pou4F1-SELEX-rev', 'Pou4F2-Hocomo-For', 'Pou4F2-Hocomo-Rev', 'Pou4F2-Selex-ForA',
          'Pou4F2-Selex-ForB', 'Pou4F2-Selex-RevA', 'Pou4F2-Selex-RevB', 'Pou4F3-Selex-For', 'Pou4F3-Selex-Rev', 'all']
cons = ['raw', 'YASS', 'MAM']

def converter(ip, in_form, out_form):
    """converter function takes transcript ID and desired form as parameters and returns a set of matched ID's.
    Desired forms are 'JG' (JGI), 'EN' (Ensembl), 'Uni' (Uniprot), and 'name' (gene name).
    When there is no match, it returns {'unknown'}.

    example:
        converter('KH.C1.1.v1.A.ND1-1', 'KH', 'Uni') = {'F6Z0S6'}"""
    global op
    with open('/home/jeff/Desktop/Thesis/Src_Files/obj/ID_Cnvrt.pkl', 'rb') as fi:
        ID_Cnvrt = pkl.load(fi)
    op = {}
    op = ID_Cnvrt[in_form].get(ip, {'Undefined'})
    if op != {'Undefinded'}:
        op = ID_Cnvrt[in_form][ip].get(out_form, {'Undefinded'})
    return op

def scr_to_gen(scr):
    """scr_to_gen takes KH transcript ID as a parameter and returns KH gene ID

    example:
        scr_to_gen('KH.C1.1.v1.A.ND1-1') = 'KH.C1.1'"""
    global gen
    gen = '.'.join(scr.split('.')[:3])
    return gen

def scr_to_chr(scr):
    """scr_to_chr function takes KH gene/transcript ID as a parameter and returns chromosome ID

    example:
        scr_to_chr('KH.C1.1.v1.A.ND1-1') = 'KhC1'"""
    global chrna
    chrna = 'Kh' + scr.split('.')[1]
    return chrna

def kh_to_name(ip):
    """kh_to_name function takes KH transcript/gene ID as a parameter and returns a gene name.
    When there is not match, function returns input string.

    example:
        kh_to_name('KH.C1.1') = 'KH.C1.1'"""
    global op
    with open('/home/jeff/Desktop/Thesis/Src_Files/obj/ID_Cnvrt.pkl', 'rb') as fi:
        ID_Cnvrt = pkl.load(fi)
    ip_gen = scr_to_gen(ip)
    op = ID_Cnvrt['KH'].get(ip_gen, ip)
    if op != ip:
        try:
            op = ''.join(ID_Cnvrt['KH'][ip_gen]['name'])
        except KeyError:
            pass
    return op

def gen_expr():
    """gen_expr function returns a dictionary with ox, ux, nc as keys and corresponding sets of KH transcript ID's of that expression as values.

    Example:
    gen_expr() = {'ox': set(.......),...}"""
    global ex_dict
    ox = set()
    ux = set()
    nc = set()
    FDR_th = .05
    with open('/home/jeff/Desktop/Thesis/Src_Files/ACDETopTagsEdgeR_Annotated.csv') as fi_seq:
        lines = fi_seq.readlines()
        for line in lines:
            if 'Gene_Model' not in line:
                words = line.split(',')
                scr = words[0]
                FDR = float(words[5])
                logFC = float(words[2])
                if FDR <= FDR_th:
                    if logFC > 0:
                        ox.add(scr)
                    elif logFC < 0:
                        ux.add(scr)
                    else:
                        nc.add(scr)
                else:
                    nc.add(scr)
    ex_dict = {'ox': ox, 'ux': ux, 'nc': nc}
    return ex_dict

def rna_seq():
    """rna_seq fuction returns a dictionary with transcript as keys and a dictionary with logFC and FDR values of the transcript as values
    Example:
        rna_seq() == dictionary[transcript] = {'logFC': value1, 'FDR': value2}
        """
    with open('/home/jeff/Desktop/Thesis/Src_Files/obj/RNA_Seq.pkl', 'rb') as fi:
        RNA_Seq = pkl.load(fi)
    return RNA_Seq

def count_dict_det(con='raw'):
    global Count_Dict
    """count_dict_det function has conservation as input and returns a corresponding total detailed count dictionary
    Example:
        count_dict('raw') == Count_Dict[motif][transcript][order] = [[lo1, up1], [lo2, up2]]
        """
    if con == 'raw':
        with open('/home/jeff/Desktop/Thesis/Src_Files/obj/Count_Dict_Det.pkl', 'rb') as fi:
            Count_Dict = pkl.load(fi)
            return Count_Dict
    elif con == 'YASS':
        with open('/home/jeff/Desktop/Thesis/Src_Files/obj/Count_Dict_Con_Det_YASS.pkl', 'rb') as fi:
            Count_Dict = pkl.load(fi)
            return Count_Dict
    elif con == 'MAM':
        with open('/home/jeff/Desktop/Thesis/Src_Files/obj/Count_Dict_Con_Det_MAM.pkl', 'rb') as fi:
            Count_Dict = pkl.load(fi)
            return Count_Dict

def get_tf_dict():
    """
    get_tf_dict function returns a dictionary with sets of corresponding transcripiton factor gene names
    :return: tf_dict[Conservation][Motif][Expression] = set()
    """
    global tf_dict
    with open('/home/jeff/Desktop/Thesis/Src_Files/TF_list_0.77.csv') as fi:
        tf_dict = {}
        for line in fi.readlines():
            items = line.split('\t')
            con, motif, expr = items[0:3]
            tf_dict[con] = tf_dict.get(con, {})
            tf_dict[con][motif] = tf_dict[con].get(motif, {})
            for name in items[4:]:
                tf_dict[con][motif][expr] = tf_dict[con][motif].get(expr, set())
                tf_dict[con][motif][expr] |= set(sub_name.strip(' ') for sub_name in name.split('//'))
    return tf_dict

class KH:
    """KH class contains convered ID's of the KH transcript.
        .KH returns a set of KH ID itself
        .JGI returns a set of JGI ID's
        .UNI returns a set of Uniprot ID's
        .ENS returns a set of Ensembl ID's
        .name returns a set of gene names
        If there is no match, returns {'unknown'}"""

    def __init__(self, scr):
        self.KH = {scr}
        self.JGI = converter(scr, 'KH', 'JG')
        self.UNI = converter(scr, 'KH', 'Uni')
        self.ENS = converter(scr, 'KH', 'EN')
        self.name = kh_to_name(scr)

    def __str__(self):
        """returns the matched ID's/names"""
        for kh in self.KH:
            output = '{0}\nJGI:'.format(kh)
        for jg in self.JGI:
            output += '\t{0}'.format(jg)
        output += '\nUniprot:'
        for uni in self.UNI:
            output += '\t{0}'.format(uni)
        output += '\nEnsembl:'
        for en in self.ENS:
            output += '\t{0}'.format(en)
        output += '\nName:'
        for name in self.name:
            output += '\t{0}'.format(name)
        return output