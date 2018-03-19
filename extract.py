import urllib.request as url
import gzip

gzfile = url.urlopen('ftp://ftp.ncbi.nih.gov/gene/DATA/gene_info.gz')
with gzip.open('/home/jeff/Desktop/Thesis/Output/gene.info.txt.gz', 'wt') as fo:
    with open('/home/jeff/Downloads/gene_info', 'r') as fi:
        lines = fi.readlines()
        fo.write(lines[0])
        for line in lines[1:]:
            taxid = line.strip('\n').split(',')[0]
            if taxid == '7719':
                fo.write(line)