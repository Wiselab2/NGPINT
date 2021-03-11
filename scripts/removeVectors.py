
import multiprocessing
import pickle
import re
from scripts.readWriteOperations import *

def hamming_distance(s1, s2):
    """
    Return the Hamming distance between equal-length sequences
    """
    if len(s1) != len(s2):
        #print(s1)
        #print(s2)
        #raise ValueError("Undefined for sequences of unequal length")
        return len(s1) if len(s1)>len(s2) else len(s2)
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

def isKthBitSet(n, k): 
    if n & (1 << (k - 1)): 
        return True 
    else: 
        return False 

def extractPlasmidSequences(options):
    """
    Returns a dictionary of plasmid sequences
    """
    return readFastaFile(options.plasmid_sequences)

def reverseComplement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return "".join([complement[base] for base in seq])[::-1]

def findVectorInReadPrefix(read_id,read_seq,vector,options,start,end):
    """
    """
    # Check for presence of vector in the 5 prime end of read_seq
    first_start=start
    all_hamming_distances=[]
    """if "read13530682" in read_id:
        print(read_seq)"""
    while start<=end:
        """if "read13530682" in read_id:
            print(start,vector[-start:],read_seq[:start])"""
        all_hamming_distances.append(hamming_distance(vector[-start:],read_seq[:start]))
        start+=1
    """if "read13530682" in read_id:
        print(all_hamming_distances)"""
    if len(all_hamming_distances)==0:
        return 10,10 # Sending a high value of hamming distance so that this read gets eliminated from being characterized as a fusion read
    return min(all_hamming_distances),all_hamming_distances.index(min(all_hamming_distances))+first_start

def findVectorInReadSuffix(read_id,read_seq,vector,options,start,end):
    """
    """
    all_hamming_distances=[]
    first_start=start
    while start<=end:
        #all_hamming_distances.append(hamming_distance(vector[-start:],read_seq[:start]))
        all_hamming_distances.append(hamming_distance(vector[:start],read_seq[-start:]))
        start+=1
    if len(all_hamming_distances)==0:
        return 10,10 # Sending a high value of hamming distance so that this read gets eliminated from being characterized as a fusion read 
    return min(all_hamming_distances),len(read_seq)-all_hamming_distances.index(min(all_hamming_distances))-first_start

def trimVectorFromReadsSingleEnded(inputs):
    """
    Trim vector off reads and stores the results in two different files
    """
    filename,options,reads_to_be_trimmed,trimmed_reads_filename,junction_reads_filename,trimming_stats_filename,reads_discarded=inputs
    fhr=open(filename,"r")
    fhw_trim=open(trimmed_reads_filename,"w")
    fhw_junction=open(junction_reads_filename,"w")
    #fhw_untrimmed_reads=open(filename+"_untrimmed_reads.txt","w")
    reads_to_be_discarded=[]
    for key in reads_discarded:
        reads_to_be_discarded.extend(reads_discarded[key])
    """print(filename,len(reads_to_be_trimmed),len(reads_to_be_discarded))
    sys.stdout.flush()
    return"""
    reads_to_be_discarded=set(reads_to_be_discarded)
    while True:
        line=fhr.readline()
        if not line:break
        if line.strip().split()[0][1:] in reads_to_be_trimmed:
            fhr.readline()
            useless=fhr.readline().strip()
            fhr.readline()
            read_id=line.strip().split()[0][1:]
            retain_from=reads_to_be_trimmed[read_id]["retain_from"]
            retain_to=reads_to_be_trimmed[read_id]["retain_to"]
            
            vector=reads_to_be_trimmed[read_id]["vector"]
            seq=reads_to_be_trimmed[read_id]["seq"]
            qual=reads_to_be_trimmed[read_id]["qual"]
            cigar=reads_to_be_trimmed[read_id]["cigar"]
            
            fhw_junction.write("@"+read_id+"_"+vector+"_"+str(retain_from)+"_"+str(retain_to)+"\n")
            fhw_junction.write(seq+"\n")
            fhw_junction.write(useless+"\n")
            fhw_junction.write(qual+"\n")
            
            fhw_trim.write("@"+read_id+"_"+vector+"_"+str(retain_from)+"_"+str(retain_to)+"\n")
            fhw_trim.write(seq+"\n")
            fhw_trim.write(useless+"\n")
            fhw_trim.write(qual+"\n")
        else:
            if line.strip().split()[0][1:] in reads_to_be_discarded:
                #print("Skipping reads discarded")
                #sys.stdout.flush()
                fhr.readline()
                fhr.readline()
                fhr.readline()
                continue
            fhw_trim.write("@"+line.strip().split()[0][1:]+"\n")
            #fhw_untrimmed_reads.write(line.strip().split()[0][1:]+"\n")
            fhw_trim.write(fhr.readline())
            fhw_trim.write(fhr.readline())
            fhw_trim.write(fhr.readline())
            
    fhr.close()
    fhw_trim.close()
    fhw_junction.close()
    #logger.info("Fusion reads written to "+junction_reads_filename)
    #fhw_untrimmed_reads.close()
    """open(trimming_stats_filename,"a+").write("Reads removed due to small size "+str(reads_removed_small_size)+"\n")
    open(trimming_stats_filename,"a+").write("Reads removed due to polyA "+str(reads_removed_polyA_tail)+"\n")
    open(trimming_stats_filename,"a+").write("num_junction_reads "+str(num_junction_reads)+"\n")
    open(trimming_stats_filename,"a+").write("num_all_reads "+str(num_all_reads)+"\n")"""
    
def trimVectorFromReadsPairedEnded(inputs):
    """
    Trim vector off reads and stores the results in two different files
    """
    filename_mate1,filename_mate2,options,reads_to_be_trimmed,trimmed_reads_filename_mate1,trimmed_reads_filename_mate2,junction_reads_filename_mate1,junction_reads_filename_mate2,trimming_stats_filename,reads_discarded=inputs
    """print("\n".join([filename_mate1,filename_mate2,trimmed_reads_filename_mate1,trimmed_reads_filename_mate2,junction_reads_filename_mate1,junction_reads_filename_mate2,trimming_stats_filename]))
    print("="*150)
    sys.stdout.flush()
    """
    fhr_mate1=open(filename_mate1,"r")
    fhr_mate2=open(filename_mate2,"r")
    fhw_trim_mate1=open(trimmed_reads_filename_mate1,"w")
    fhw_trim_mate2=open(trimmed_reads_filename_mate2,"w")
    fhw_junction_mate1=open(junction_reads_filename_mate1,"w")
    fhw_junction_mate2=open(junction_reads_filename_mate2,"w")
    
    reads_to_be_discarded=[]
    for key in reads_discarded:
        reads_to_be_discarded.extend(reads_discarded[key])
        
    reads_to_be_discarded=set(reads_to_be_discarded)
    """print(filename_mate1,len(reads_to_be_discarded),len(reads_to_be_trimmed))
    sys.stdout.flush()"""
    while True:
        line_mate1=fhr_mate1.readline()
        line_mate2=fhr_mate2.readline()
        if (not line_mate1) or (not line_mate2):break
        if line_mate1.strip().split()[0][1:]+"_1" in reads_to_be_trimmed or line_mate2.strip().split()[0][1:]+"_2" in reads_to_be_trimmed:
            #print("Inside here",line_mate1,line_mate2)
            """fhr_mate1.readline()
            fhr_mate1.readline()
            fhr_mate1.readline()
            
            fhr_mate2.readline()
            fhr_mate2.readline()
            fhr_mate2.readline()
            continue"""
            # progress the file pointers to the next entry in both files
            seq_mate1=fhr_mate1.readline().strip()
            useless_mate1=fhr_mate1.readline().strip()
            qual_mate1=fhr_mate1.readline().strip()
            seq_mate2=fhr_mate2.readline().strip()
            useless_mate2=fhr_mate2.readline().strip()
            qual_mate2=fhr_mate2.readline().strip()
            
            read_id_mate1,retain_from_mate1,retain_to_mate1,read_id_mate2,retain_from_mate2,retain_to_mate2="","","","","",""
            if line_mate1.strip().split()[0][1:]+"_1" in reads_to_be_trimmed:
                read_id_mate1=line_mate1.strip().split()[0][1:]+"_1"
                retain_from_mate1=reads_to_be_trimmed[read_id_mate1]["retain_from"]
                retain_to_mate1=reads_to_be_trimmed[read_id_mate1]["retain_to"]
                
            if line_mate2.strip().split()[0][1:]+"_2" in reads_to_be_trimmed:
                read_id_mate2=line_mate2.strip().split()[0][1:]+"_2"
                retain_from_mate2=reads_to_be_trimmed[read_id_mate2]["retain_from"]
                retain_to_mate2=reads_to_be_trimmed[read_id_mate2]["retain_to"]
            
            if (read_id_mate1=="" and read_id_mate2!="") or (read_id_mate1!="" and read_id_mate2==""):
                # Vector sequence present in either read_id_mate1 or read_id_mate2 
                if read_id_mate1!="" and read_id_mate2=="":
                    read_id,retain_from,retain_to=read_id_mate1,retain_from_mate1,retain_to_mate1
                    fhw_junction=fhw_junction_mate1
                    fhw_trim=fhw_trim_mate1
                    useless=useless_mate1
                elif read_id_mate1=="" and read_id_mate2!="":
                    read_id,retain_from,retain_to=read_id_mate2,retain_from_mate2,retain_to_mate2
                    fhw_junction=fhw_junction_mate2
                    fhw_trim=fhw_trim_mate2
                    useless=useless_mate2
                
                vector=reads_to_be_trimmed[read_id]["vector"]
                seq=reads_to_be_trimmed[read_id]["seq"]
                qual=reads_to_be_trimmed[read_id]["qual"]
                cigar=reads_to_be_trimmed[read_id]["cigar"]
                read_id=read_id[:-2]
                read_header="@"+read_id+"_"+vector+"_"+str(retain_from)+"_"+str(retain_to)
                
                fhw_junction.write(read_header+"\n")
                fhw_trim.write("@"+read_id+"_"+vector+"_"+str(retain_from)+"_"+str(retain_to)+"\n")
                fhw_junction.write(seq+"\n")
                fhw_trim.write(seq+"\n")             
                fhw_junction.write(useless+"\n")
                fhw_trim.write(useless+"\n")
                fhw_junction.write(qual+"\n")
                fhw_trim.write(qual+"\n")
                
                if read_id_mate1!="" and read_id_mate2=="":
                    fhw_junction=fhw_junction_mate2
                    fhw_trim=fhw_trim_mate2
                    useless=useless_mate2
                    line=line_mate2
                    fhr=fhr_mate2
                    seq=seq_mate2
                    qual=qual_mate2
                else:
                    fhw_junction=fhw_junction_mate1
                    fhw_trim=fhw_trim_mate1
                    useless=useless_mate1
                    line=line_mate1
                    fhr=fhr_mate1
                    seq=seq_mate1
                    qual=qual_mate1
                fhw_trim.write(read_header+"\n")
                fhw_trim.write(seq+"\n")
                fhw_trim.write(useless+"\n")
                fhw_trim.write(qual+"\n")
                
                fhw_junction.write(read_header+"\n")
                fhw_junction.write(seq+"\n")
                fhw_junction.write(useless+"\n")
                fhw_junction.write(qual+"\n")
            else:
                # Vector sequence present in both read_id_mate1 and read_id_mate2
                for entry_num,each_entry in enumerate([[read_id_mate1,retain_from_mate1,retain_to_mate1],[read_id_mate2,retain_from_mate2,retain_to_mate2]]):
                    read_id,retain_from,retain_to=each_entry
                    vector=reads_to_be_trimmed[read_id]["vector"]
                    seq=reads_to_be_trimmed[read_id]["seq"]
                    qual=reads_to_be_trimmed[read_id]["qual"]
                    cigar=reads_to_be_trimmed[read_id]["cigar"]
                    
                    read_id=read_id[:-2]
                    if entry_num==0:
                        fhw_junction=fhw_junction_mate1
                        fhw_trim=fhw_trim_mate1
                        useless=useless_mate1
                    else:
                        fhw_junction=fhw_junction_mate2
                        fhw_trim=fhw_trim_mate2
                        useless=useless_mate2
                    
                    fhw_junction.write("@"+read_id+"_"+vector+"_"+str(retain_from)+"_"+str(retain_to)+"\n")
                    fhw_trim.write("@"+read_id+"_"+vector+"_"+str(retain_from)+"_"+str(retain_to)+"\n")
        
                    fhw_junction.write(seq+"\n")
                    fhw_trim.write(seq+"\n")
                    fhw_junction.write(useless+"\n")
                    fhw_trim.write(useless+"\n")
                    fhw_junction.write(qual+"\n")
                    fhw_trim.write(qual+"\n")
        else:
            for entry_num,each_entry in enumerate([[fhr_mate1,fhw_junction_mate1,fhw_trim_mate1,line_mate1],[fhr_mate2,fhw_junction_mate2,fhw_trim_mate2,line_mate2]]):
                fhr,fhw_junction,fhw_trim,line=each_entry
                #print(line)
                fhw_trim.write("@"+line.strip().split()[0][1:]+"\n")
                fhw_trim.write(fhr.readline())
                fhw_trim.write(fhr.readline())
                fhw_trim.write(fhr.readline())
            """fhw_trim_mate1.write("@"+line_mate1.strip().split()[0][1:]+"\n")
            fhw_trim_mate1.write(fhr_mate1.readline())
            fhw_trim_mate1.write(fhr_mate1.readline())
            fhw_trim_mate1.write(fhr_mate1.readline())   
        
            fhw_trim_mate2.write("@"+line_mate2.strip().split()[0][1:]+"\n")
            fhw_trim_mate2.write(fhr_mate2.readline())
            fhw_trim_mate2.write(fhr_mate2.readline())
            fhw_trim_mate2.write(fhr_mate2.readline())  """  
            
    fhr_mate1.close()
    fhr_mate2.close()
    fhw_trim_mate1.close()
    fhw_trim_mate2.close()
    fhw_junction_mate1.close()
    fhw_junction_mate2.close()
    
    #logger.info("Fusion reads written to "+junction_reads_filename_mate1+" and "+junction_reads_filename_mate2)
    """open(trimming_stats_filename,"a+").write("Reads removed due to small size "+str(reads_removed_small_size)+"\n")
    open(trimming_stats_filename,"a+").write("Reads removed due to polyA "+str(reads_removed_polyA_tail)+"\n")
    open(trimming_stats_filename,"a+").write("num_junction_reads "+str(num_junction_reads)+"\n")
    open(trimming_stats_filename,"a+").write("num_all_reads "+str(num_all_reads)+"\n")"""

def checkForVectorsInRead(read_id,read_seq,vectors,cigar,options):
    """
    """
    prefix_min_hamming_dist_5_prime_forward,prefix_min_hamming_dist_5_prime_reverse,prefix_min_hamming_dist_3_prime_forward,prefix_min_hamming_dist_3_prime_reverse=5,5,5,5
    max_limit_end=25
    if re.search(r'^\d*S[A-Z0-9]*\d*S$',cigar):
        start=int(options.min_trimmed_length)
        end=((int(cigar.split("S")[0])+max_limit_end) if (int(cigar.split("S")[0]) > start) else int(options.min_trimmed_length)+max_limit_end)
        if int(cigar.split("S")[0])>=int(options.min_trimmed_length):
            prefix_min_hamming_dist_5_prime_forward,prefix_cut_5_prime_forward=findVectorInReadPrefix(read_id,read_seq,vectors["5_prime_end"],options,start,end)
            prefix_cut_5_prime_forward=prefix_cut_5_prime_forward-int(options.frame_of_TF_fusion)
            start,end=int(options.min_trimmed_length),int(0.80*len(read_seq))
            prefix_min_hamming_dist_3_prime_reverse,prefix_cut_3_prime_reverse=findVectorInReadSuffix(read_id,reverseComplement(read_seq),vectors["3_prime_end"],options,start,end)
        #===================================================================================================================================================================
        
        start=int(options.min_trimmed_length)
        end=([int(s) for s in re.findall(r'-?\d+\.?\d*', cigar)][-1] +max_limit_end) if ([int(s) for s in re.findall(r'-?\d+\.?\d*', cigar)][-1] > start) else (int(options.min_trimmed_length) +max_limit_end)
        if [int(s) for s in re.findall(r'-?\d+\.?\d*', cigar)][-1]>=int(options.min_trimmed_length):
            prefix_min_hamming_dist_5_prime_reverse,prefix_cut_5_prime_reverse=findVectorInReadPrefix(read_id,reverseComplement(read_seq),vectors["5_prime_end"],options,start,end)
            prefix_cut_5_prime_reverse=prefix_cut_5_prime_reverse-int(options.frame_of_TF_fusion)
            start,end=int(options.min_trimmed_length),int(0.80*len(read_seq))
            prefix_min_hamming_dist_3_prime_forward,prefix_cut_3_prime_forward=findVectorInReadSuffix(read_id,read_seq,vectors["3_prime_end"],options,start,end)
    elif re.search(r'\d*S$',cigar):
        start=int(options.min_trimmed_length)
        end=([int(s) for s in re.findall(r'-?\d+\.?\d*', cigar)][-1] +max_limit_end) if ([int(s) for s in re.findall(r'-?\d+\.?\d*', cigar)][-1] > start) else (int(options.min_trimmed_length) +max_limit_end)
        if [int(s) for s in re.findall(r'-?\d+\.?\d*', cigar)][-1]>=int(options.min_trimmed_length):
            prefix_min_hamming_dist_5_prime_reverse,prefix_cut_5_prime_reverse=findVectorInReadPrefix(read_id,reverseComplement(read_seq),vectors["5_prime_end"],options,start,end)
            prefix_cut_5_prime_reverse=prefix_cut_5_prime_reverse-int(options.frame_of_TF_fusion)
            start,end=int(options.min_trimmed_length),int(0.80*len(read_seq))
            prefix_min_hamming_dist_3_prime_forward,prefix_cut_3_prime_forward=findVectorInReadSuffix(read_id,read_seq,vectors["3_prime_end"],options,start,end)    
            prefix_min_hamming_dist_3_prime_reverse,prefix_min_hamming_dist_5_prime_forward=1000,1000 # Assigning dummy values
            
    elif re.search(r'^\d*S',cigar):
        start=int(options.min_trimmed_length)
        end=((int(cigar.split("S")[0])+max_limit_end) if (int(cigar.split("S")[0]) > start) else int(options.min_trimmed_length)+max_limit_end)
        """if "read13530682" in read_id:
            print("startend",start,end)"""
        if int(cigar.split("S")[0])>=int(options.min_trimmed_length):
            prefix_min_hamming_dist_5_prime_forward,prefix_cut_5_prime_forward=findVectorInReadPrefix(read_id,read_seq,vectors["5_prime_end"],options,start,end)
            prefix_cut_5_prime_forward=prefix_cut_5_prime_forward-int(options.frame_of_TF_fusion)
            start,end=int(options.min_trimmed_length),int(0.80*len(read_seq))
            prefix_min_hamming_dist_3_prime_reverse,prefix_cut_3_prime_reverse=findVectorInReadSuffix(read_id,reverseComplement(read_seq),vectors["3_prime_end"],options,start,end)
            prefix_min_hamming_dist_3_prime_forward,prefix_min_hamming_dist_5_prime_reverse=1000,1000 # Assigning dummy values
        
    """if "read13530682" in read_id:
        print(read_id)
        print(cigar)
        print(start,end)
        print(prefix_min_hamming_dist_5_prime_forward,prefix_min_hamming_dist_5_prime_reverse,prefix_min_hamming_dist_3_prime_forward,prefix_min_hamming_dist_3_prime_reverse)
        sys.stdout.flush()"""
    if min([prefix_min_hamming_dist_5_prime_forward,prefix_min_hamming_dist_5_prime_reverse,prefix_min_hamming_dist_3_prime_forward,prefix_min_hamming_dist_3_prime_reverse])>5:
        return -1,-1,-1
    
    
    if prefix_min_hamming_dist_5_prime_forward<min([prefix_min_hamming_dist_5_prime_reverse,prefix_min_hamming_dist_3_prime_forward,prefix_min_hamming_dist_3_prime_reverse]):
        """print("5_prime_end_forward")
        print(cigar)
        print(read_seq)
        print("-"*prefix_cut_5_prime_forward+read_seq[prefix_cut_5_prime_forward:])
        print(prefix_min_hamming_dist_5_prime_forward,prefix_min_hamming_dist_5_prime_reverse,prefix_min_hamming_dist_3_prime_forward,prefix_min_hamming_dist_3_prime_reverse)
        print("="*200)
        """
        return prefix_cut_5_prime_forward+1,len(read_seq),"5_prime_end_forward"
        pass
    elif prefix_min_hamming_dist_5_prime_reverse<min([prefix_min_hamming_dist_5_prime_forward,prefix_min_hamming_dist_3_prime_forward,prefix_min_hamming_dist_3_prime_reverse]):
        """print("5_prime_end_reverse")
        print(cigar)
        print(read_seq)
        print(read_seq[:-prefix_cut_5_prime_reverse]+"-"*prefix_cut_5_prime_reverse)
        print(prefix_min_hamming_dist_5_prime_forward,prefix_min_hamming_dist_5_prime_reverse,prefix_min_hamming_dist_3_prime_forward,prefix_min_hamming_dist_3_prime_reverse)
        print("="*200)
        """
        return 1,len(read_seq)-prefix_cut_5_prime_reverse,"5_prime_end_reverse"
        pass
    elif prefix_min_hamming_dist_3_prime_forward<min([prefix_min_hamming_dist_5_prime_forward,prefix_min_hamming_dist_5_prime_reverse,prefix_min_hamming_dist_3_prime_reverse]):
        """print("3_prime_end_forward")
        print(read_id,cigar)
        print(read_seq)
        print(prefix_cut_3_prime_forward)
        print(read_seq[:prefix_cut_3_prime_forward]+"-"*(len(read_seq)-prefix_cut_3_prime_forward))
        print(prefix_min_hamming_dist_5_prime_forward,prefix_min_hamming_dist_5_prime_reverse,prefix_min_hamming_dist_3_prime_forward,prefix_min_hamming_dist_3_prime_reverse)
        print(1,prefix_cut_3_prime_forward,"3_prime_end_forward")
        print("="*200)"""
        return 1,prefix_cut_3_prime_forward,"3_prime_end_forward"
        pass
    elif prefix_min_hamming_dist_3_prime_reverse<min([prefix_min_hamming_dist_5_prime_forward,prefix_min_hamming_dist_5_prime_reverse,prefix_min_hamming_dist_3_prime_forward]):
        """print("3_prime_end_reverse")
        print(read_id,cigar)
        print(read_seq)
        print(prefix_cut_3_prime_reverse)
        print("-"*(len(read_seq)-prefix_cut_3_prime_reverse)+read_seq[-prefix_cut_3_prime_reverse:])
        print(prefix_min_hamming_dist_5_prime_forward,prefix_min_hamming_dist_5_prime_reverse,prefix_min_hamming_dist_3_prime_forward,prefix_min_hamming_dist_3_prime_reverse)
        
        print(len(read_seq)-prefix_cut_3_prime_reverse+1,len(read_seq),"3_prime_end_reverse")
        print("="*200)"""
        return len(read_seq)-prefix_cut_3_prime_reverse+1,len(read_seq),"3_prime_end_reverse"
        pass
    else:
        """print("MATCHES NOTHING")
        print(prefix_min_hamming_dist_5_prime_forward,prefix_min_hamming_dist_5_prime_reverse,prefix_min_hamming_dist_3_prime_forward,prefix_min_hamming_dist_3_prime_reverse)
        print("="*200)"""
        return -1,-1,-1
    
def findReadsWithVectorSequence(inputs):
    """
    Finds and flags the reads with vector sequence
    """
    vectors,aligned_filename,filename,options,trimming_stats_filename,logfile,logger_proxy,logging_mutex=inputs
    
    reads_discarded={"whole_5_prime_vector":[],
                     "whole_3_prime_vector":[],
                     "incorrect_fragment_amplified_5_prime_vector":[],
                     "incorrect_fragment_amplified_3_prime_vector":[],
                     "polyAtail":[],
                     "too_small_usable_fragment":[],
                     "other":[]}

    five_prime_trim_read_info={}
    three_prime_trim_read_info={}    
    hits_to_yeast_chromomosome={}
    hits_to_yeast_plasmid={}
    
    plasmids=extractPlasmidSequences(options)
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    
    fhw=open(logfile,"w")
    fhw.write("Read_id,explanation,retained_start,retained_end\n")
    
    fhr=open(aligned_filename,"r")
    for line in fhr:
        if line[0]=="@":continue
        line=line.strip()
        #print(line)
        read_id,mapping_orientation_STAR,chromosome,map_pos,score,cigar,useless1,useless2,useless3,read_seq,read_qual=line.split("\t")[:11]
        if options.selected_ended=="PE" and options.background_ended=="PE":
            read_id=read_id+"_"+("1" if isKthBitSet(int(mapping_orientation_STAR),7)==True else "2")
        mapping_orientation="0" if isKthBitSet(int(mapping_orientation_STAR),5)==False else "1"
        # Skipping reads that do not have any soft clips
        if "S" not in cigar:continue
        #if "prime" in chromosome:continue
        if read_id in five_prime_trim_read_info or read_id in three_prime_trim_read_info:continue
        if "prime" in chromosome:
            """if "read14913735" in read_id:
                print(re.search(r'^\d*M[MID0-9]*\d*M$',cigar) )
                print(re.search(r'^\d*M$',cigar))
                sys.stdout.flush()"""
            if re.search(r'^\d*M[MID0-9]*\d*M$',cigar) or re.search(r'^\d*M$',cigar):
                """if "read14913735" in read_id:
                    print("I am inside")
                    sys.stdout.flush()"""
                # Matches cigar of the form 149M --> Matches the whole Vector
                if chromosome=="5_prime_end":
                    reads_discarded["whole_5_prime_vector"].append(read_id)
                    fhw.write(read_id+",whole_5_prime_vector,-1,-1"+"\n")
                else:
                    reads_discarded["whole_3_prime_vector"].append(read_id)
                    fhw.write(read_id+",whole_3_prime_vector,-1,-1"+"\n")
            
            elif re.search(r'^\d*S[a-zA-Z0-9]*\d*S$',cigar) or re.search(r'^(\d*[MID])*\d*S$',cigar):
                # Matches cigar of the form 34S*3S
                if chromosome=="5_prime_end":
                    extract_from=[int(s) for s in re.findall(r'-?\d+\.?\d*', cigar)][-1]
                    if extract_from < options.min_trimmed_length:
                        reads_discarded["too_small_usable_fragment"].append(read_id)
                        fhw.write(read_id+",too_small_fragment_5_prime_vector,-1,-1"+"\n")
                        continue
                    retain_from=len(read_seq)-(extract_from+options.frame_of_TF_fusion-1)
                    retain_to=len(read_seq)
                    five_prime_trim_read_info[read_id]={"retain_from":retain_from,
                                                        "retain_to":retain_to,
                                                        "vector":"5_prime_end_forward" if mapping_orientation=="0" else "5_prime_end_reverse",
                                                        "seq":read_seq[retain_from-1:retain_to] if mapping_orientation=="0" else reverseComplement(read_seq[retain_from-1:retain_to]),
                                                        "qual":read_qual[retain_from-1:retain_to] if mapping_orientation=="0" else read_qual[retain_from-1:retain_to][::-1],
                                                        "cigar":cigar}
                    """if five_prime_trim_read_info[read_id]["seq"][:2]!="GG" and five_prime_trim_read_info[read_id]["seq"][-2:]!="CC":
                        #print(dist_from_5_prime_vector,dist_from_3_prime_vector)
                        print(read_id,five_prime_trim_read_info[read_id]["vector"])
                        print(read_seq)
                        print(retain_from,retain_to)
                        print(five_prime_trim_read_info[read_id]["seq"])
                        print("="*150) 
                        sys.stdout.flush()"""
                    fhw.write(read_id+(",5_prime_end_forward," if mapping_orientation=="0" else ",5_prime_end_reverse,")+str(retain_from)+","+str(retain_to)+"\n")
                else:
                    extract_till=[int(s) for s in re.findall(r'-?\d+\.?\d*', cigar)][0]
                    if extract_till < options.min_trimmed_length:
                        reads_discarded["too_small_usable_fragment"].append(read_id)
                        #fhw.write(read_id+",too_small_fragment_3_prime_vector,-1,-1"+"\n")
                        continue
                    seq_trimmed=read_seq[:extract_till]
                    if seq_trimmed.count('A')/extract_till>0.9 or seq_trimmed.count('T')/extract_till>0.9:
                        reads_discarded["polyAtail"].append(read_id) 
                        #fhw.write(read_id+",polyA_tail,-1,-1"+"\n")
                        continue
                    three_prime_trim_read_info[read_id]={"retain_from":1,
                                                         "retain_to":extract_till,
                                                         "vector":"3_prime_end_forward" if mapping_orientation=="0" else "3_prime_end_reverse",
                                                         "seq":read_seq[:extract_till] if mapping_orientation=="0" else reverseComplement(read_seq[:extract_till]),
                                                         "qual":read_qual[:extract_till] if mapping_orientation=="0" else read_qual[:extract_till][::-1],
                                                         "cigar":cigar}
                    fhw.write(read_id+(",3_prime_end_forward," if mapping_orientation=="0" else ",3_prime_end_reverse,")+str(1)+","+str(extract_till)+"\n")
                    
            elif re.search(r'\d*S$',cigar):
                # Matches cigar of form 130M30S
                if chromosome=="5_prime_end":
                    extract_from=[int(s) for s in re.findall(r'-?\d+\.?\d*', cigar)][-1]
                    if extract_from < options.min_trimmed_length:
                        reads_discarded["too_small_usable_fragment"].append(read_id)
                        fhw.write(read_id+",too_small_fragment_5_prime_vector,-1,-1"+"\n")
                        continue
                    retain_from=len(read_seq)-(extract_from+options.frame_of_TF_fusion-1)
                    retain_to=len(read_seq)
                    five_prime_trim_read_info[read_id]={"retain_from":retain_from,
                                                        "retain_to":retain_to,
                                                        "vector":"5_prime_end_forward" if mapping_orientation=="0" else "5_prime_end_reverse",
                                                        "seq":read_seq[retain_from-1:retain_to] if mapping_orientation=="0" else reverseComplement(read_seq[retain_from-1:retain_to]),
                                                        "qual":read_qual[retain_from-1:retain_to] if mapping_orientation=="0" else read_qual[retain_from-1:retain_to][::-1],
                                                        "cigar":cigar}
                    """if five_prime_trim_read_info[read_id]["seq"][:2]!="GG" and five_prime_trim_read_info[read_id]["seq"][-2:]!="CC":
                        #print(dist_from_5_prime_vector,dist_from_3_prime_vector)
                        print(read_id,five_prime_trim_read_info[read_id]["vector"])
                        print(read_seq)
                        print(retain_from,retain_to)
                        print(five_prime_trim_read_info[read_id]["seq"])
                        print("="*150) 
                        sys.stdout.flush()"""
                    fhw.write(read_id+(",5_prime_end_forward," if mapping_orientation=="0" else ",5_prime_end_reverse,")+str(retain_from)+","+str(retain_to)+"\n")
                elif chromosome=="3_prime_end":
                    reads_discarded["incorrect_fragment_amplified_3_prime_vector"].append(read_id)
                    fhw.write(read_id+",incorrect_fragment_amplified_3_prime_vector,-1,-1"+"\n")
            elif re.search(r'^\d*S',cigar):
                # Matches cigar of the form 30S120M
                if chromosome=="3_prime_end":
                    extract_till=int(cigar.split("S")[0])
                    if extract_till < options.min_trimmed_length:
                        reads_discarded["too_small_usable_fragment"].append(read_id)
                        fhw.write(read_id+",too_small_fragment_3_prime_vector,-1,-1"+"\n")
                        continue
                    seq_trimmed=read_seq[:extract_till]
                    if seq_trimmed.count('A')/extract_till>0.9 or seq_trimmed.count('T')/extract_till>0.9:
                        reads_discarded["polyAtail"].append(read_id) 
                        fhw.write(read_id+",polyA_tail,-1,-1"+"\n")
                        continue
                    three_prime_trim_read_info[read_id]={"retain_from":1,
                                                         "retain_to":extract_till,
                                                         "vector":"3_prime_end_forward" if mapping_orientation=="0" else "3_prime_end_reverse",
                                                         "seq":read_seq[:extract_till] if mapping_orientation=="0" else reverseComplement(read_seq[:extract_till]),
                                                        "qual":read_qual[:extract_till] if mapping_orientation=="0" else read_qual[:extract_till][::-1],
                                                         "cigar":cigar}
                    fhw.write(read_id+(",3_prime_end_forward," if mapping_orientation=="0" else ",3_prime_end_reverse,")+str(1)+","+str(extract_till)+"\n")
                elif chromosome=="5_prime_end":
                    reads_discarded["incorrect_fragment_amplified_5_prime_vector"].append(read_id)
                    fhw.write(read_id+",incorrect_fragment_amplified_5_prime_vector,-1,-1"+"\n")
            else:
                if chromosome=="5_prime_end":
                    reads_discarded["other"].append(read_id)
                    fhw.write(read_id+",other_5_prime_vector,-1,-1"+"\n")
                else:
                    reads_discarded["other"].append(read_id)
                    fhw.write(read_id+",other_3_prime_vector,-1,-1"+"\n")
        elif chromosome in plasmids.keys():
            hits_to_yeast_plasmid[read_id]={"retain_from":0,"retain_to":-1,"seq":read_seq,"qual":read_qual,"cigar":cigar}
            fhw.write(read_id+",plasmid,-1,-1"+"\n")
        else:
            # Most of the read contains transcriptomic sequence
            
            # Inverting the read, the quality values and the CIGAR if the mapping orientation is reversed
            if isKthBitSet(int(mapping_orientation_STAR),5)==True:
                read_seq="".join(complement[base] for base in read_seq)[::-1]
                read_qual=read_qual[::-1]
                cigar_values=[int(s) for s in re.findall(r'\d+',cigar)][::-1]
                cigar_alphabets=re.findall(r'[A-Z]',cigar)
                cigar_alphabets=cigar_alphabets[::-1]
                cigar=""
                for j in range(len(cigar_values)):
                    cigar+=str(cigar_values[j])+cigar_alphabets[j]
            
            retain_from,retain_to,vector=checkForVectorsInRead(read_id,read_seq,vectors,cigar,options)
            if vector==-1:
                continue
            if "5_prime" in vector:
                five_prime_trim_read_info[read_id]={"retain_from":retain_from,
                                                        "retain_to":retain_to,
                                                        "vector":vector,
                                                        "seq":read_seq[retain_from-1:retain_to] ,
                                                        "qual":read_qual[retain_from-1:retain_to],
                                                        "cigar":cigar}
            else:
                three_prime_trim_read_info[read_id]={"retain_from":retain_from,
                                                        "retain_to":retain_to,
                                                        "vector":vector,
                                                        "seq":read_seq[retain_from-1:retain_to] ,
                                                        "qual":read_qual[retain_from-1:retain_to],
                                                        "cigar":cigar}
            fhw.write(read_id+","+vector+","+str(retain_from)+","+str(retain_to)+"\n")
        #print(read_id)
                    
    fhr.close()
    #fhw.close()
    #print("reads_not_trimmed_due_to_small_vector_sequence",len(reads_not_trimmed_due_to_small_vector_sequence))
    length_distros=[[],[],[]]
    fhw=open(filename+"_trimming_info.pkl","wb")
    pickle.dump([filename,length_distros,five_prime_trim_read_info,three_prime_trim_read_info,hits_to_yeast_chromomosome,hits_to_yeast_plasmid,reads_discarded],fhw)
    fhw.close()
    
    fhw=open(trimming_stats_filename,"w")
    five_prime_forward,five_prime_reverse,three_prime_forward,three_prime_reverse=0,0,0,0
    for read_id in five_prime_trim_read_info:
        if "forward" in five_prime_trim_read_info[read_id]["vector"]:
            five_prime_forward+=1
        elif "reverse" in five_prime_trim_read_info[read_id]["vector"]:
            five_prime_reverse+=1
            
    for read_id in three_prime_trim_read_info:
        if "forward" in three_prime_trim_read_info[read_id]["vector"]:
            three_prime_forward+=1
        elif "reverse" in three_prime_trim_read_info[read_id]["vector"]:
            three_prime_reverse+=1
    
    fhw.write("Five Prime Forward "+str(five_prime_forward)+"\n")
    fhw.write("Five Prime Reverse "+str(five_prime_reverse)+"\n")
    fhw.write("Five Prime Vector whole read "+str(len(reads_discarded["whole_5_prime_vector"]))+"\n")
    fhw.write("Three Prime Forward "+str(three_prime_forward)+"\n")
    fhw.write("Three Prime Reverse "+str(three_prime_reverse)+"\n")
    fhw.write("Three Prime Vector whole read "+str(len(reads_discarded["whole_3_prime_vector"]))+"\n")
    fhw.write("Plasmid "+str(len(hits_to_yeast_plasmid))+"\n")
    fhw.write("Yeast Chromosome "+str(len(hits_to_yeast_chromomosome))+"\n")
    fhw.write("Reads removed due to small size "+str(len(reads_discarded["too_small_usable_fragment"]))+"\n")
    fhw.write("Reads removed due to polyA "+str(len(reads_discarded["polyAtail"]))+"\n")
    fhw.write("num_junction_reads "+str(len(five_prime_trim_read_info)+len(three_prime_trim_read_info))+"\n")
    #fhw.write("num_all_reads ",str(len(reads_discarded[""])))
    fhw.close()
    
    with logging_mutex:
        logger_proxy.info("Finding fusion reads for "+filename+" completed")

    """print(trimming_stats_filename)
    for key in reads_discarded:
        print(key,len(reads_discarded[key]))
    print("Five prime vectors",len(five_prime_trim_read_info))
    print("Three prime vectors",len(three_prime_trim_read_info))"""
    #logger.info("Finding fusion reads for "+filename+" completed")   

def findReadsWithVectorSequenceAndTrim(options,logger_proxy,logging_mutex):
    """
    Looks through the STAR alignments to find reads with vector sequence 
    """
    vectors=readFastaFile(options.vector_sequences)
    pool = multiprocessing.Pool(processes=int(options.CPU))
    allinputs=[]
    for num,eachtype in enumerate([options.selected_sample_N_removed,options.background_sample_N_removed]):
        for file_num,filename in enumerate(eachtype):
            if options.selected_ended=="PE" and options.background_ended=="PE":
                if file_num%2==1:
                    pass
                else:
                    continue
            if num==0:
                aligned_filename=options.selected_sample_STAR_genome_filename_round1[file_num]
                trimming_stats_filename=options.selected_sample_trimming_stats[file_num]
                logfile=options.selected_sample_per_read_log[file_num]
            else:
                aligned_filename=options.background_sample_STAR_genome_filename_round1[file_num]
                trimming_stats_filename=options.background_sample_trimming_stats[file_num]
                logfile=options.background_sample_per_read_log[file_num]
            allinputs.append([vectors,aligned_filename,filename,options,trimming_stats_filename,logfile,logger_proxy,logging_mutex])
    #pool.map(findReadsWithVectorSequence,allinputs)
    pool.map(findReadsWithVectorSequence,allinputs)
    
    results=[]
    for num,eachtype in enumerate([options.selected_sample_N_removed,options.background_sample_N_removed]):
        for file_num,filename in enumerate(eachtype):
            if options.selected_ended=="PE" and options.background_ended=="PE":
                if file_num%2==1:
                    pass
                else:
                    continue
            results.append(pickle.load(open(filename+"_trimming_info.pkl","rb")))
    
    #print(results)
    results_dict={}
    for result in results:
        results_dict[result[0]]={"lengths_all_trimmed_reads":result[1][0],
                                 "lengths_all_5_prime_trimmed_reads":result[1][1],
                                 "lengths_all_3_prime_trimmed_reads":result[1][2],
                                 "five_prime_trim_read_info":result[2],
                                 "three_prime_trim_read_info":result[3],
                                 "hits_to_yeast_chromomosome":result[4],
                                 "hits_to_yeast_plasmid":result[5],
                                 "reads_discarded":result[6]
                                 }
        """print("Reading",result[0],"completed")
        sys.stdout.flush()"""
        """print(result[0],len(results_dict[result[0]]["five_prime_trim_read_info"])+len(results_dict[result[0]]["three_prime_trim_read_info"]))
        sys.stdout.flush()"""
    # Clean up the reads and create new files
    pool = multiprocessing.Pool(processes=int(options.CPU)-1)
    allinputs=[]
    for num,eachtype in enumerate([options.selected_sample_N_removed,options.background_sample_N_removed]):
        for file_num,filename in enumerate(eachtype):
            if options.selected_ended=="PE" and options.background_ended=="PE":
                if file_num%2==1:
                    pass
                else:
                    continue
            reads_to_be_trimmed={}
            reads_to_be_trimmed.update(results_dict[filename]["five_prime_trim_read_info"])
            reads_to_be_trimmed.update(results_dict[filename]["three_prime_trim_read_info"])
            #reads_to_be_trimmed.update(results_dict[filename]["hits_to_yeast_chromomosome"])
            #reads_to_be_trimmed.update(results_dict[filename]["hits_to_yeast_plasmid"])
            if num==0:
                trimmed_reads_filename=options.selected_sample_all_reads_vector_trimmed[file_num]
                junction_reads_filename=options.selected_sample_fusion_reads[file_num]
                trimming_stats_filename=options.selected_sample_trimming_stats[file_num]
                if options.selected_ended=="PE" and options.background_ended=="PE":
                    if file_num%2==1:
                        trimmed_reads_filename_mate1=options.selected_sample_all_reads_vector_trimmed[file_num-1]
                        junction_reads_filename_mate1=options.selected_sample_fusion_reads[file_num-1]
                        trimmed_reads_filename_mate2=trimmed_reads_filename
                        junction_reads_filename_mate2=junction_reads_filename
                        filename_mate2=filename
                        filename_mate1=eachtype[file_num-1]
            else:
                trimmed_reads_filename=options.background_sample_all_reads_vector_trimmed[file_num]
                junction_reads_filename=options.background_sample_fusion_reads[file_num]
                trimming_stats_filename=options.background_sample_trimming_stats[file_num]
                if options.selected_ended=="PE" and options.background_ended=="PE":
                    if file_num%2==1:
                        trimmed_reads_filename_mate1=options.background_sample_all_reads_vector_trimmed[file_num-1]
                        junction_reads_filename_mate1=options.background_sample_fusion_reads[file_num-1]
                        trimmed_reads_filename_mate2=trimmed_reads_filename
                        junction_reads_filename_mate2=junction_reads_filename
                        filename_mate2=filename
                        filename_mate1=eachtype[file_num-1]
            if options.selected_ended=="PE" and options.background_ended=="PE":
                allinputs.append([filename_mate1,filename_mate2,options,reads_to_be_trimmed,trimmed_reads_filename_mate1,trimmed_reads_filename_mate2,junction_reads_filename_mate1,junction_reads_filename_mate2,trimming_stats_filename,results_dict[filename]["reads_discarded"]])
            else:
                allinputs.append([filename,options,reads_to_be_trimmed,trimmed_reads_filename,junction_reads_filename,trimming_stats_filename,results_dict[filename]["reads_discarded"]])
    """print(len(allinputs),options.selected_ended, options.background_ended)
    sys.stdout.flush()"""   
    if options.selected_ended=="PE" and options.background_ended=="PE":
        """print("Calling function for trimming")
        sys.stdout.flush()"""
        pool.map(trimVectorFromReadsPairedEnded,allinputs)
    else:
        pool.map(trimVectorFromReadsSingleEnded,allinputs)

