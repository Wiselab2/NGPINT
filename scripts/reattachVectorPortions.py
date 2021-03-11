

import itertools
import multiprocessing
import os
import sys

from scripts.removeVectors import *


def prepareEachGenomeFileSingleEnded(input):
    """
    """
    genome_alignment_file,options,all_reads_filename,output_genome_alignment_file,genome_alignment_file_round1,logger_proxy,logging_mutex=input
    #print("Input filename",genome_alignment_file)
    #print("Output filename",output_genome_alignment_file)
    #sys.stdout.flush()
    # Add the vector sequences to the trimmed reads
    all_reads=readFastqFile(all_reads_filename)
    """print(all_reads_filename)
    sys.stdout.flush() """ 
    #output_genome_alignment_file=genome_alignment_file[:-3]+"_modified.sam"
    fhr=open(genome_alignment_file,"r")
    fhw=open(output_genome_alignment_file,"w")
    reads_processed=[]
    for line in fhr:
        if line[0]=="@":
            fhw.write(line)
            continue
        line=line.strip().split()
        reads_processed.append(line[0].split("_")[0])
        if "5_prime_end_forward" in line[0]:
            read_id=line[0].split("_")[0]
            start,end=int(line[0].split("_")[-2]),int(line[0].split("_")[-1])
            if start!=1:
                vector_start,vector_end=1,start-1
            else:
                vector_start,vector_end=end+1,len(all_reads[read_id][0])
            vector_seq=all_reads[read_id][0][:vector_end]
            cigar=line[5]
            cigar_values=[int(s) for s in re.findall(r'\d+',cigar)]
            cigar_alphabets=re.findall(r'[A-Z]',cigar)
            #if int(line[1]) in [0,256]:
            if isKthBitSet(int(line[1]), 5)==False:   
                if cigar_alphabets[0]=="S":
                    cigar_values[0]=int(cigar_values[0])+vector_end-vector_start+1
                else:
                    cigar_values.insert(0,str(vector_end-vector_start+1))
                    cigar_alphabets.insert(0,"S")
                cigar_values=list(map(str,cigar_values))
                new_cigar="".join(list(itertools.chain(*zip(cigar_values, cigar_alphabets))))
                line[5]=new_cigar
                line[9]=all_reads[read_id][0]
                line[10]=all_reads[read_id][2]
                line.append("FR:i:1")
                fhw.write("\t".join(line)+"\n")
                #print(cigar,new_cigar)
                sys.stdout.flush()
            else:
                if cigar_alphabets[-1]=="S":
                    cigar_values[-1]=int(cigar_values[-1])+vector_end-vector_start+1
                else:
                    cigar_values.append(str(vector_end-vector_start+1))
                    cigar_alphabets.append("S")
                cigar_values=list(map(str,cigar_values))
                new_cigar="".join(list(itertools.chain(*zip(cigar_values, cigar_alphabets))))
                line[5]=new_cigar
                line[9]=reverseComplement(all_reads[read_id][0])
                line[10]=all_reads[read_id][2][::-1]
                line.append("FR:i:1")
                fhw.write("\t".join(line)+"\n")
                #print(cigar,new_cigar)
                sys.stdout.flush()
        elif "5_prime_end_reverse" in line[0]:
            read_id=line[0].split("_")[0]
            start,end=int(line[0].split("_")[-2]),int(line[0].split("_")[-1])
            if start==1:
                vector_start,vector_end=end+1,len(all_reads[read_id][0])
            else:
                vector_start,vector_end=1,start-1
            #vector_seq=all_reads[read_id][0][vector_start-1:]
            cigar=line[5]
            cigar_values=[int(s) for s in re.findall(r'\d+',cigar)]
            cigar_alphabets=re.findall(r'[A-Z]',cigar)
            if isKthBitSet(int(line[1]), 5)==False:  
                if cigar_alphabets[-1]=="S":
                    cigar_values[-1]=int(cigar_values[-1])+vector_end-vector_start+1
                else:
                    cigar_values.append(str(vector_end-vector_start+1))
                    cigar_alphabets.append("S")
                cigar_values=list(map(str,cigar_values))
                new_cigar="".join(list(itertools.chain(*zip(cigar_values, cigar_alphabets))))
                line[5]=new_cigar
                line[9]=all_reads[read_id][0]
                line[10]=all_reads[read_id][2]
                line.append("FR:i:1")
                fhw.write("\t".join(line)+"\n")
                #print(cigar,new_cigar)
                sys.stdout.flush()
            else:
                if cigar_alphabets[0]=="S":
                    cigar_values[0]=int(cigar_values[0])+vector_end-vector_start+1
                else:
                    cigar_values.insert(0,str(vector_end-vector_start+1))
                    cigar_alphabets.insert(0,"S")
                #temp=[iter(cigar_values),iter(cigar_alphabets)]
                cigar_values=list(map(str,cigar_values))
                new_cigar="".join(list(itertools.chain(*zip(cigar_values, cigar_alphabets))))
                line[5]=new_cigar
                line[9]=reverseComplement(all_reads[read_id][0])
                line[10]=all_reads[read_id][2][::-1]
                line.append("FR:i:1")
                fhw.write("\t".join(line)+"\n")
                """if line[0]=="J00102:13:HTW3WBBXX:1:2205:7253:47436_5_prime_end_reverse_84_149":
                    print(vector_start,vector_end)
                    #print("vector ",vector_seq)
                    print(cigar_values)
                    print(cigar_alphabets)
                    print(new_cigar)
                    print(line)"""
        elif "3_prime_end" in line[0]:
            line.append("FR:i:1")
            fhw.write("\t".join(line)+"\n")
        else:
            line.append("FR:i:0")
            fhw.write("\t".join(line)+"\n")
    fhr.close()
    
    reads_processed=set(reads_processed)
    fhr=open(genome_alignment_file_round1,"r")
    for line in fhr:
        if line[0]=="@":
            continue
        if line.split()[0] not in reads_processed:
            fhw.write(line.strip()+"\tFR:i:0"+"\n")
    fhr.close()
    fhw.close()
    
    # Sort and create indices 
    cmd="samtools view -Sb "+output_genome_alignment_file+"|samtools sort - > "+output_genome_alignment_file+"sorted"
    os.system(cmd)
    cmd="mv "+output_genome_alignment_file+"sorted "+output_genome_alignment_file
    os.system(cmd)
    cmd="samtools index "+output_genome_alignment_file
    os.system(cmd)
    
    with logging_mutex:
        logger_proxy.info("Preparing file for genome browser viewing over for "+output_genome_alignment_file)
    
def prepareEachGenomeFilePairedEnded(input):
    """
    """
    #pprint.pprint(input)
    genome_alignment_file_round2,options,all_reads_filename_mate1,all_reads_filename_mate2,output_genome_alignment_file,genome_alignment_file_round1,fusion_reads_filename_mate1,fusion_reads_filename_mate2,logger_proxy,logging_mutex=input
    # Add the vector sequences to the trimmed reads
    all_reads_mate1=readFastqFile(all_reads_filename_mate1)
    all_reads_mate2=readFastqFile(all_reads_filename_mate2)
    fusion_reads_mate1=readFastqFile(fusion_reads_filename_mate1)
    fusion_reads_mate2=readFastqFile(fusion_reads_filename_mate2)
    #output_genome_alignment_file=genome_alignment_file[:-3]+"_modified.sam"
    fhr=open(genome_alignment_file_round2,"r")
    fhw=open(output_genome_alignment_file,"w")
    reads_processed=[]
    for line in fhr:
        if line[0]=="@":
            fhw.write(line)
            continue
        line=line.strip().split()
        if isKthBitSet(int(line[1]),7)==True:
            if line[0] not in fusion_reads_mate1:
                """print("skipping")
                sys.stdout.flush()"""
                line.append("FR:i:0")
                fhw.write("\t".join(line)+"\n")
                continue
        else:
            if line[0] not in fusion_reads_mate2:
                """print("skipping")
                sys.stdout.flush()"""
                line.append("FR:i:0")
                fhw.write("\t".join(line)+"\n")
                continue
        reads_processed.append(line[0].split("_")[0]+line[1])
        if "5_prime_end_forward" in line[0] and "S" in line[5]:
            read_id=line[0].split("_")[0]
            start,end=int(line[0].split("_")[-2]),int(line[0].split("_")[-1])
            cigar=line[5]
            cigar_values=[int(s) for s in re.findall(r'\d+',cigar)]
            cigar_alphabets=re.findall(r'[A-Z]',cigar)
            """if isKthBitSet(int(line[1]), 7)==True:
                vector_seq=all_reads_mate1[read_id][0][:vector_end]
            else:
                vector_seq=all_reads_mate2[read_id][0][:vector_end]"""
            #if int(line[1]) in [0,256]:
            if isKthBitSet(int(line[1]), 5)==False: 
                if start!=1:
                    vector_start,vector_end=1,start-1
                else:
                    vector_start,vector_end=end+1,len(all_reads_filename_mate1[read_id][0])
                #Forward reads  
                if cigar_alphabets[0]=="S":
                    cigar_values[0]=int(cigar_values[0])+vector_end-vector_start+1
                else:
                    cigar_values.insert(0,str(vector_end-vector_start+1))
                    cigar_alphabets.insert(0,"S")
                old_line=line[:]
                cigar_values=list(map(str,cigar_values))
                new_cigar="".join(list(itertools.chain(*zip(cigar_values, cigar_alphabets))))
                line[5]=new_cigar
                if isKthBitSet(int(line[1]), 7)==True:
                    line[9]=all_reads_mate1[read_id][0]
                    line[10]=all_reads_mate1[read_id][2]
                else:
                    line[9]=all_reads_mate2[read_id][0]
                    line[10]=all_reads_mate2[read_id][2]
                line.append("FR:i:1")
                cigar_values=[int(s) for s in re.findall(r'\d+',new_cigar)]
                cigar_alphabets=re.findall(r'[A-Z]',new_cigar)
                if sum([cigar_values[i] for i in range(len(cigar_alphabets)) if cigar_alphabets[i]!="N"])!=len(line[9]):
                    #print("CIGAR",line[0],line[1],cigar,new_cigar)
                    fhw.write("\t".join(old_line)+"\n")
                    sys.stdout.flush()
                else:
                    fhw.write("\t".join(line)+"\n")
                    """print("Writing new line")
                    sys.stdout.flush()"""
            else:
                # Reverse reads
                if start!=1:
                    vector_start,vector_end=1,start-1
                else:
                    vector_start,vector_end=end+1,len(all_reads_filename_mate2[read_id][0])
                if cigar_alphabets[-1]=="S":
                    cigar_values[-1]=int(cigar_values[-1])+vector_end-vector_start+1
                else:
                    cigar_values.append(str(vector_end-vector_start+1))
                    cigar_alphabets.append("S")
                    
                
                old_line=line[:]
                cigar_values=list(map(str,cigar_values))
                new_cigar="".join(list(itertools.chain(*zip(cigar_values, cigar_alphabets))))
                line[5]=new_cigar
                if isKthBitSet(int(line[1]), 7)==True:
                    line[9]=reverseComplement(all_reads_mate1[read_id][0])
                    line[10]=all_reads_mate1[read_id][2][::-1]
                else:
                    line[9]=reverseComplement(all_reads_mate2[read_id][0])
                    line[10]=all_reads_mate2[read_id][2][::-1]
                
                line.append("FR:i:1")
                fhw.write("\t".join(line)+"\n")
                cigar_values=[int(s) for s in re.findall(r'\d+',new_cigar)]
                cigar_alphabets=re.findall(r'[A-Z]',new_cigar)
                if sum([cigar_values[i] for i in range(len(cigar_alphabets)) if cigar_alphabets[i]!="N"])!=len(line[9]):
                    #print("CIGAR",line[0],line[1],cigar,new_cigar)
                    fhw.write("\t".join(old_line)+"\n")
                    sys.stdout.flush()
                else:
                    fhw.write("\t".join(line)+"\n")
                    """print("Writing new line")
                    sys.stdout.flush()"""
        elif "5_prime_end_reverse" in line[0]  and "S" in line[5]:
            read_id=line[0].split("_")[0]
            start,end=int(line[0].split("_")[-2]),int(line[0].split("_")[-1])
            cigar=line[5]
            cigar_values=[int(s) for s in re.findall(r'\d+',cigar)]
            cigar_alphabets=re.findall(r'[A-Z]',cigar)
            #if int(line[1]) in [0,256]:   
            if isKthBitSet(int(line[1]), 5)==False:  
                if isKthBitSet(int(line[1]), 7)==True:
                    if start==1:
                        vector_start,vector_end=end+1,len(all_reads_mate1[read_id][0])
                    else:
                        vector_start,vector_end=1,start-1
                    #vector_start,vector_end=end+1,len(all_reads_mate1[read_id][0])
                    #vector_seq=all_reads_mate1[read_id][0][vector_start-1:]
                else:
                    if start==1:
                        vector_start,vector_end=end+1,len(all_reads_mate2[read_id][0])
                    else:
                        vector_start,vector_end=1,start-1
                    #vector_start,vector_end=end+1,len(all_reads_mate2[read_id][0])
                    #vector_seq=all_reads_mate2[read_id][0][vector_start-1:]
                
                
                if cigar_alphabets[-1]=="S":
                    cigar_values[-1]=int(cigar_values[-1])+vector_end-vector_start+1
                else:
                    cigar_values.append(str(vector_end-vector_start+1))
                    cigar_alphabets.append("S")
                #temp=[iter(cigar_values),iter(cigar_alphabets)]
                old_line=line[:]
                cigar_values=list(map(str,cigar_values))
                new_cigar="".join(list(itertools.chain(*zip(cigar_values, cigar_alphabets))))
                line[5]=new_cigar
                if isKthBitSet(int(line[1]), 7)==True:
                    line[9]=all_reads_mate1[read_id][0]
                    line[10]=all_reads_mate1[read_id][2]
                else:
                    line[9]=all_reads_mate2[read_id][0]
                    line[10]=all_reads_mate2[read_id][2]
                
                line.append("FR:i:1")
                fhw.write("\t".join(line)+"\n")
                cigar_values=[int(s) for s in re.findall(r'\d+',new_cigar)]
                cigar_alphabets=re.findall(r'[A-Z]',new_cigar)
                if sum([cigar_values[i] for i in range(len(cigar_alphabets)) if cigar_alphabets[i]!="N"])!=len(line[9]):
                    #print("CIGAR",line[0],line[1],cigar,new_cigar)
                    fhw.write("\t".join(old_line)+"\n")
                    sys.stdout.flush()
                else:
                    fhw.write("\t".join(line)+"\n")
                    """print("Writing new line")
                    sys.stdout.flush()"""
            else:
                if isKthBitSet(int(line[1]), 7)==True:
                    if start==1:
                        vector_start,vector_end=end+1,len(all_reads_mate1[read_id][0])
                    else:
                        vector_start,vector_end=1,start-1
                else:
                    if start==1:
                        vector_start,vector_end=end+1,len(all_reads_mate2[read_id][0])
                    else:
                        vector_start,vector_end=1,start-1
                        
                if cigar_alphabets[0]=="S":
                    cigar_values[0]=int(cigar_values[0])+vector_end-vector_start+1
                else:
                    cigar_values.insert(0,str(vector_end-vector_start+1))
                    cigar_alphabets.insert(0,"S")
                #temp=[iter(cigar_values),iter(cigar_alphabets)]
                old_line=line[:]
                cigar_values=list(map(str,cigar_values))
                new_cigar="".join(list(itertools.chain(*zip(cigar_values, cigar_alphabets))))
                line[5]=new_cigar
                if isKthBitSet(int(line[1]), 7)==True:
                    line[9]=reverseComplement(all_reads_mate1[read_id][0])
                    line[10]=all_reads_mate1[read_id][2][::-1]
                else:
                    line[9]=reverseComplement(all_reads_mate2[read_id][0])
                    line[10]=all_reads_mate2[read_id][2][::-1]
                line.append("FR:i:1")
                fhw.write("\t".join(line)+"\n")
                cigar_values=[int(s) for s in re.findall(r'\d+',new_cigar)]
                cigar_alphabets=re.findall(r'[A-Z]',new_cigar)
                if sum([cigar_values[i] for i in range(len(cigar_alphabets)) if cigar_alphabets[i]!="N"])!=len(line[9]):
                    #print("CIGAR",line[0],line[1],cigar,new_cigar)
                    fhw.write("\t".join(old_line)+"\n")
                    sys.stdout.flush()
                else:
                    fhw.write("\t".join(line)+"\n")
                    """print("Writing new line")
                    sys.stdout.flush()"""
        else:
            line.append("FR:i:0")
            fhw.write("\t".join(line)+"\n")
    fhr.close()
    
    reads_processed=set(reads_processed)
    fhr=open(genome_alignment_file_round1,"r")
    for line in fhr:
        if line[0]=="@":
            continue
        if line.split()[0]+line.split()[1] not in reads_processed:
            fhw.write(line.strip()+"\tFR:i:0"+"\n")
    fhr.close()
    fhw.close()
    
    # Sort and create indices 
    cmd="samtools view -Sb "+output_genome_alignment_file+"|samtools sort - > "+output_genome_alignment_file+"sorted"
    os.system(cmd)
    cmd="mv "+output_genome_alignment_file+"sorted "+output_genome_alignment_file
    os.system(cmd)
    cmd="samtools index "+output_genome_alignment_file
    os.system(cmd)
    
    with logging_mutex:
        logger_proxy.info("Preparing file for genome browser viewing over for "+output_genome_alignment_file)

def prepareGenomeFilesForGenomeBrowser(options,logger_proxy,logging_mutex):
    """
    Append reads from round1 which are absent in round2
    """
    pool = multiprocessing.Pool(processes=int(options.CPU))
    allinputs=[]
    for num,eachtype in enumerate([options.selected_sample_N_removed,options.background_sample_N_removed]):
        for file_num,filename in enumerate(eachtype):
            if num==0:
                if options.selected_ended=="SE" and options.background_ended=="SE":
                    pass
                else:
                    if file_num%2==1:
                        pass
                    else:
                        continue
                genome_alignment_file=options.selected_sample_STAR_genome_filename_round2[file_num]
                all_reads_filename=options.selected_sample_N_removed[file_num]
                output_genome_alignment_file=options.selected_sample_genome_browser_per_replicate[file_num]
                genome_alignment_file_round1=options.selected_sample_STAR_genome_filename_round1[file_num]
                #print(genome_alignment_file,all_reads_filename,output_genome_alignment_file,genome_alignment_file_round1)
                #allinputs.append([genome_alignment_file,options,all_reads_filename,output_genome_alignment_file,genome_alignment_file_round1])
                if options.selected_ended=="PE" and options.background_ended=="PE":
                    fusion_reads_filename_mate1=options.selected_sample_fusion_reads[file_num-1]
                    fusion_reads_filename_mate2=options.selected_sample_fusion_reads[file_num]
                    allinputs.append([genome_alignment_file,options,options.selected_sample_N_removed[file_num-1],all_reads_filename,output_genome_alignment_file,genome_alignment_file_round1,fusion_reads_filename_mate1,fusion_reads_filename_mate2,logger_proxy,logging_mutex])
                else:
                    allinputs.append([genome_alignment_file,options,all_reads_filename,output_genome_alignment_file,genome_alignment_file_round1,logger_proxy,logging_mutex])
    if options.selected_ended=="PE" and options.background_ended=="PE":
        pool.map(prepareEachGenomeFilePairedEnded,allinputs)
    else:
        pool.map(prepareEachGenomeFileSingleEnded,allinputs)
    
    # Merge all the replicates together
    if options.selected_ended=="SE" and options.background_ended=="SE":
        cmd="samtools merge -f "+options.selected_sample_genome_browser
        for filename in options.selected_sample_genome_browser_per_replicate:
            cmd+=" "+filename
        os.system(cmd)
        cmd="samtools index "+options.selected_sample_genome_browser
        os.system(cmd)
    else:
        cmd="samtools merge -f "+options.selected_sample_genome_browser
        for num,eachtype in enumerate([options.selected_sample_N_removed,options.background_sample_N_removed]):
            for file_num,filename in enumerate(eachtype):
                if num==0:
                    if file_num%2==1:
                        pass
                    else:
                        continue
                    cmd+=" "+options.selected_sample_genome_browser_per_replicate[file_num]
        os.system(cmd)

