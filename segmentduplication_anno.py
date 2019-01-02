import sys
from intervaltree import Interval, IntervalTree

def interTree(filename):   
    t = IntervalTree()
    count = 0
    with open(filename, 'r') as f :       
        for each_line in f:
            line = each_line.strip().split("\t")
            chrname = line[1]   
            posbegin = int(float(line[2]))
            posend = int(float(line[3]))
            strand = line[6]     
            chrname2 = line[7]    
            posbegin2 = int(float(line[8]))     
            posend2 = int(float(line[9]))    
            countname = 'seg' + str(count)
            data = [chrname, chrname2, posbegin2, posend2, strand, countname]   
            t[posbegin:posend] = data
            count += 1
        return(t)
             
class bnd(object):
    def __init__(self,bnd_string):
        self.bnd_string = bnd_string
        self.bnd_string = self.bnd_string.replace(">","")
        self.bnd_string = self.bnd_string.replace("<","")
        if '[' in self.bnd_string:
            first_right_braket = self.bnd_string.index('[')
            next_right_braket = self.bnd_string.index('[', first_right_braket+1)
            if first_right_braket != 0:
                self.stat = "s1"
            else:
                self.stat = "s4"
            self.pos = self.bnd_string[first_right_braket+1:next_right_braket]
        elif ']' in self.bnd_string:
            first_left_braket = self.bnd_string.index(']')
            next_left_braket = self.bnd_string.index(']', first_left_braket+1)
            if first_left_braket != 0:
                self.stat = "s2"
            else:
                self.stat = "s3"
            self.pos = self.bnd_string[first_left_braket+1:next_left_braket]
        self.chrom = self.pos.split(":")[0]
        self.pos_num = self.pos.split(":")[1]
            
class Getvcfinfo():
    def __init__(self, each_line):
        line = each_line.strip().split('\t')
        chrname = line[0]
        posbegin = int(float(line[1]))
        posend = line[7].split(';')[3].split('=')[1]
        chrname2 = line[7].split(';')[2].split('=')[1]
        ID = line[2]
        svtype = line[4]
        if '.' in posend:
            b = bnd(each_line)
            posend = b.pos_num
            chrname2 = b.chrom
        posend = int(float(posend))
        if 'chr' not in chrname:
            chrname = 'chr' + chrname
        if 'chr' not in chrname2:
            chrname2 = 'chr' + chrname2
        self.chrname = chrname
        self.posbegin = posbegin
        self.posend  = posend
        self.chrname2 = chrname2
        self.ID = ID            
        self.svtype = svtype

def calTemp(temp, a_chrname, a_chrname2, a_posbegin, a_posend, pos1, pos2, pos3, pos4):
    if temp:        
        temp = list(temp)
        region1 = ''
        region2 = ''
        #datasetid = ''
        datasetid = []
        #strand_dire = ''
        strand_dire = []
        for i in temp:
            if i[2][0] == a_chrname and (i[0] <= pos1) and (pos2 <= i[1]):
                chr_2 = i[2][1]    
                posbegin_2 = i[2][2]   
                posend_2 = i[2][3] 
                if (chr_2 == a_chrname2) and (posbegin_2 <= pos3) and (pos4 <= posend_2):
                    region1 = '{0}:{1}-{2}'.format(i[2][0], str(i[0]), str(i[1]))
                    region2 = '{0}:{1}-{2}'.format(chr_2, str(posbegin_2), str(posend_2))
                    strand_dire.append(i[2][4])
                    datasetid.append(i[2][5])
                
        return region1,region2, datasetid,strand_dire
    else:
        return '', '', '', ''
        
def calRegion(t, a, each_line2,outfile, bp):
    region1 = ''
    region2 = ''
    datasetid = ''
    strand_dire = ''
    segdup = ''
    if a.chrname == a.chrname2:
        bp = 1
    else:
        bp = bp
    bp = int(bp)
    pos1 = a.posbegin - bp
    pos2 = a.posbegin + bp
    pos3 = a.posend - bp
    pos4 = a.posend + bp
    if t.overlaps(pos1, pos2):
        temp = t.overlap(pos1, pos2)
        region1, region2, datasetid, strand_dire= calTemp(temp, a.chrname, a.chrname2, a.posbegin, a.posend, pos1, pos2, pos3, pos4)
              
    return region1, region2, datasetid, strand_dire

def main():

    filename2 = sys.argv[1]
    outfile = sys.argv[2]
    bp = int(sys.argv[3])
    dataset = sys.argv[4]
    t = interTree(dataset)
     
    flag = 1
    count_vcf = 0 
    count2 = 0
    strand_dire2 = []
    datasetid2 = []

    with open(filename2,'r') as f2:
        for each_line2 in f2:
            if each_line2.strip()[0] == '#':
                continue
            else:
                count2 += 1
                a = Getvcfinfo(each_line2)
                with open(outfile, 'a+') as f3:
                    if flag == 1:
                        print('SVID', 'DatasetID', file = f3)
                        flag = 0
                    region1, region2, datasetid, strand_dire = calRegion(t, a, each_line2, outfile, bp) 
                 
                    if len(a.svtype) == 5  and region1 != '':
                        if 'DUP' in a.svtype or 'DEL' in a.svtype:
                            if len(strand_dire) == 1:
                                if strand_dire[0] == '+': 
                                    print(a.ID, datasetid, file = f3)
                                    count_vcf += 1
                            else:
                                for i in range(len(strand_dire)):
                                    if strand_dire[i] == '+':
                                        datasetid2.append(datasetid[i])
                                        strand_dire2.append(strand_dire[i])
                                if len(strand_dire2):
                                    print(a.ID, datasetid, file = f3)
                                    count_vcf += 1
                                    strand_dire2 = []
                                    datasetid2 = []
                        elif 'INV' in a.svtype:
                            if len(strand_dire) == 1:
                                if strand_dire[0] == "-":
                                    print(a.ID, datasetid, file = f3)
                                    count_vcf += 1
                            else:                                
                                for i in range(len(strand_dire)):                      
                                    if strand_dire[i] == "-":  
                                        datasetid2.append(datasetid[i])
                                        strand_dire2.append(strand_dire[i])
                                if len(strand_dire2):
                                    print(a.ID, datasetid, file = f3)
                                    count_vcf += 1 
                                    strand_dire2 = []
                                    datasetid2 = []
                    else:
                        if (']' in a.svtype or '[' in a.svtype) and region1 != '':
                            print(a.ID, datasetid, file = f3)
                            count_vcf += 1
                        elif (('DUP' in a.svtype and 'INV' in a.svtype) or ('DEL' in a.svtype and 'INV' in a.svtype)) and region1 != '':
                            print(a.ID, dasetid, file = f3)
                            count_vcf += 1
                        else:
                            continue
                        
if __name__ == "__main__":
    main() 

    