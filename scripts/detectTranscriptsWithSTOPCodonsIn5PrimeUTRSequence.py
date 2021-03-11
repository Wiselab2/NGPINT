

def allindices(string, sub, offset=0):
    listindex=[]
    i = string.find(sub, offset)
    while i >= 0:
        listindex.append(i)
        i = string.find(sub, i + 1)
    return listindex
 

def findTranscriptsWithSTOPCodonsIn5PrimeUTRSequence(options,gene_info):
    """
    Finds transcripts with STOP codons in 5' UTR sequences
    """
    transcripts_to_in_frame_STOP_codon_locations_in_5_prime_UTR_region={}
    fhr=open(options.transcriptome,"r")
    for line in fhr:
        if ">" in line:
            transcript_id=line.strip().split()[0][1:]
            if transcript_id not in gene_info:continue
            #start,end=int(line.strip().split("CDS=")[-1].split("-")[0]),int(line.strip().split("CDS=")[-1].split("-")[-1])
            start,end=gene_info[transcript_id]["cds_start"],gene_info[transcript_id]["cds_end"]
            seq=fhr.readline().strip()
            five_prime_UTR_region=seq[:start-1]
            if len(five_prime_UTR_region) % 3 == 1:
                five_prime_UTR_region=five_prime_UTR_region[1:]
                dropped=1
            elif len(five_prime_UTR_region) % 3 == 2:
                five_prime_UTR_region=five_prime_UTR_region[2:]
                dropped=2
            else:
                dropped=0
            locations_of_TAA=allindices(five_prime_UTR_region,"TAA")
            locations_of_TAG=allindices(five_prime_UTR_region,"TAG")
            locations_of_TGA=allindices(five_prime_UTR_region,"TGA")
            locations_of_TAA_in_frame_with_CDS_check=[x%3 for x in locations_of_TAA]
            locations_of_TAG_in_frame_with_CDS_check=[x%3 for x in locations_of_TAG]
            locations_of_TGA_in_frame_with_CDS_check=[x%3 for x in locations_of_TGA]
            if dropped==1:
                locations_of_TAA_in_frame_with_CDS=[loc+1+1 for num,loc in enumerate(locations_of_TAA) if locations_of_TAA_in_frame_with_CDS_check[num]==0]
                locations_of_TAG_in_frame_with_CDS=[loc+1+1 for num,loc in enumerate(locations_of_TAG) if locations_of_TAG_in_frame_with_CDS_check[num]==0]
                locations_of_TGA_in_frame_with_CDS=[loc+1+1 for num,loc in enumerate(locations_of_TGA) if locations_of_TGA_in_frame_with_CDS_check[num]==0]
            elif dropped==2:
                locations_of_TAA_in_frame_with_CDS=[loc+2+1 for num,loc in enumerate(locations_of_TAA) if locations_of_TAA_in_frame_with_CDS_check[num]==0]
                locations_of_TAG_in_frame_with_CDS=[loc+2+1 for num,loc in enumerate(locations_of_TAG) if locations_of_TAG_in_frame_with_CDS_check[num]==0]
                locations_of_TGA_in_frame_with_CDS=[loc+2+1 for num,loc in enumerate(locations_of_TGA) if locations_of_TGA_in_frame_with_CDS_check[num]==0]
            else:
                locations_of_TAA_in_frame_with_CDS=[loc+1 for num,loc in enumerate(locations_of_TAA) if locations_of_TAA_in_frame_with_CDS_check[num]==0]
                locations_of_TAG_in_frame_with_CDS=[loc+1 for num,loc in enumerate(locations_of_TAG) if locations_of_TAG_in_frame_with_CDS_check[num]==0]
                locations_of_TGA_in_frame_with_CDS=[loc+1 for num,loc in enumerate(locations_of_TGA) if locations_of_TGA_in_frame_with_CDS_check[num]==0]
            all_locations=[]
            all_locations.extend(locations_of_TAA_in_frame_with_CDS)
            all_locations.extend(locations_of_TAG_in_frame_with_CDS)
            all_locations.extend(locations_of_TGA_in_frame_with_CDS)
            transcripts_to_in_frame_STOP_codon_locations_in_5_prime_UTR_region[transcript_id]=sorted(all_locations)
    fhr.close()
    return transcripts_to_in_frame_STOP_codon_locations_in_5_prime_UTR_region