

def writeFastqFile(filename,reads):
    fhw=open(filename,"w")
    for read in reads:
        fhw.write("@"+read+"\n")
        fhw.write(reads[read][0]+"\n"+reads[read][1]+"\n"+reads[read][2]+"\n")

def writeFastaFile(filename,seqs):
    fhw=open(filename,"w")
    for id in seqs:
        fhw.write(">"+id+"\n"+seqs[id]+"\n")

def readFastqFile(filename):
    reads={}
    fhr=open(filename,"r")
    while True:
        line=fhr.readline()
        if not line:
            break
        reads[line.split()[0][1:]]=[fhr.readline().strip(),fhr.readline().strip(),fhr.readline().strip()]
        #print(reads[line.split()[0]])
    return reads

def readFastaFile(filename):
    """
    Reads in a fasta file and returns a dictionary
    The keys in the dictionary is same as the fasta header
    for each sequence upto the first space.
    """
    info={}
    fhr=open(filename,"r")
    while(True):
        line=fhr.readline()
        if not line: break
        if(">" in line):
            try:
                info[line.strip()[1:].split()[0]]=fhr.readline().strip()
            except ValueError:
                pass
    return info

