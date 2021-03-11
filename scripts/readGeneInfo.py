


def getGeneInfo(transcriptome_filename,transcript_to_gene_map_filename):
    """
    """
    data={}
    fhr=open(transcript_to_gene_map_filename,"r")
    for line in fhr:
        transcript,gene=line.split()
        data[transcript]={"gene":gene}
    fhr.close()
    
    fhr=open(transcriptome_filename,"r")
    for line in fhr:
        if ">" in line:
            transcript=line[1:].split()[0]
            if "CDS=" in line:
                cds_start,cds_end=line[1:].split("CDS=")[-1].split("-")[0],line[1:].split("CDS=")[-1].split("-")[1]
                data[transcript]["cds_start"]=int(cds_start)
                data[transcript]["cds_end"]=int(cds_end)
            """else:
                data[transcript]["cds_start"]=-1
                #data[transcript]["cds_end"]=len(fhr.readline().strip())
                data[transcript]["cds_end"]=-1"""
    fhr.close()
    
    all_transcripts=list(data.keys())
    for transcript in all_transcripts:
        if "cds_start" not in data[transcript]:
            del data[transcript]
    return data