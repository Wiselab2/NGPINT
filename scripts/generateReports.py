
import multiprocessing
import os
import re
from scripts.readWriteOperations import *
from scripts.removeVectors import *

def computeInFrameInfoForEachJunctionRead(inputs):
    """
    Performs the following operations:
        - Takes as input a bamfile for 5' junction reads
        - Computes whether the expression was in-frame
        - If expression was in-frame, then it appends the alignment with portion of the read which is trimmed off
    """
    bamfilename,options,gene_info=inputs
    
    samfilename=bamfilename[:-3]+"sam"
    cmd="samtools view -h "+bamfilename+" > "+samfilename
    os.system(cmd)
    
    number_of_cigars=[[] for i in range(1000)] # ith position stores the number of reads which have i nucleotides soft clipped from the 5' end
    samfilename_junction_info_appended=samfilename[:-4]+"_junction_info_appended.sam"
    fhr=open(samfilename,"r")
    fhw=open(samfilename_junction_info_appended,"w")
    
    for line in fhr:
        bases_soft_clipped=-1
        if line[0]=="@":
            fhw.write(line)
            continue
        #print(len(line.strip().split("\t")),line.strip().split("\t"))
        read_id,orientation,transcript_id,starting_loc,useless1,cigar,useless2,useless3,useless4,read_seq,read_qual=line.strip().split("\t")[:11]
        if transcript_id not in gene_info:continue
        if re.search(r'^\d*M$',cigar) or re.search(r'^\d*M\d*S$',cigar): # No Soft clipping at 5'end
            number_of_cigars[0].append(read_id)
            portion_of_read_unmapped=""
            portion_of_read_mapped=read_seq
        elif re.search(r'^\d*S',cigar):# Some Soft clipping at 5'end
            bases_soft_clipped=[int(s) for s in re.findall(r'-?\d+\.?\d*', cigar)][0]
            number_of_cigars[bases_soft_clipped].append(read_id)
            portion_of_read_unmapped=read_seq[:bases_soft_clipped]
            portion_of_read_mapped=read_seq[bases_soft_clipped:]
        starting_loc=int(starting_loc)
        if bases_soft_clipped==-1 or starting_loc>=gene_info[transcript_id]["cds_end"]:
            line=line.strip()+"\t"+"FN:0"
            fhw.write(line+"\n")
            continue
        if bases_soft_clipped % 3 == 0:
            if (starting_loc - gene_info[transcript_id]["cds_start"]) % 3 == 0:
                line=line.strip()+"\t"+"FI:"+portion_of_read_unmapped+":"+str(starting_loc)
            else:
                line=line.strip()+"\t"+"FN:0"
        elif bases_soft_clipped % 3 == 1:
            if (starting_loc + 2 - gene_info[transcript_id]["cds_start"]) % 3 == 0:
                line=line.strip()+"\t"+"FI:"+portion_of_read_unmapped+":"+str(starting_loc)
            else:
                line=line.strip()+"\t"+"FN:0"
        elif bases_soft_clipped % 3 == 2:
            if (starting_loc + 1 - gene_info[transcript_id]["cds_start"]) % 3 == 0:
                line=line.strip()+"\t"+"FI:"+portion_of_read_unmapped+":"+str(starting_loc)
            else:
                line=line.strip()+"\t"+"FN:0"
        fhw.write(line+"\n")
    fhw.close() 
    fhr.close()

def collectInformationAboutEnrichmentAnalysis(salmon_DGE_filename,salmon_normalized_counts,gene_info):
    """
    """
    info={}
    fhr=open(salmon_normalized_counts,"r")
    for line_num,line in enumerate(fhr):
        if line_num==0:
            line=line.strip().split(",")[1:]
            header=[l.strip("\"") for l in line]
        else:
            gene=line.strip().split(",")[0].strip("\"")
            info[gene]={"norm_counts":{header[i]:val for i,val in enumerate(line.strip().split(",")[1:])}}
            info[gene]["log2FoldChange"]=-1
            info[gene]["padj"]=-1
    fhr.close()
    
    fhr=open(salmon_DGE_filename,"r")
    for line_num,line in enumerate(fhr):
        if line_num==0:continue
        gene=line.strip().split(",")[0].strip("\"")
        info[gene]["log2FoldChange"]=line.strip().split(",")[2]
        info[gene]["padj"]=line.strip().split(",")[-1]
    fhr.close()
    
    genes=list(set([gene_info[transcript]["gene"] for transcript in gene_info]))
    for gene in genes:
        if gene not in info:
            info[gene]={"log2FoldChange":-1,"padj":-1,"norm_counts":{header[i]:-1 for i in range(len(header))}}
    return info

def compileInformationForEachSample(fusion_reads_samfilename,ended,cpu):
    info={}
    
    # Convert samfile to name sorted
    #print(fusion_reads_samfilename)
    #return
    if ended=="PE":
        """cmd=SAMTOOLS+" view -@ "+str(cpu)+" -Sb "
        cmd+=fusion_reads_samfilename+" | "
        cmd+=SAMTOOLS+" sort -@ "+str(cpu)+" -n | "
        cmd+=SAMTOOLS+" view -@ "+str(cpu)+" > "+fusion_reads_samfilename+".namesorted.sam"
        print(cmd)
        return"""
        #os.system(cmd)
        
        fhr=open(fusion_reads_samfilename,"r")
        reads_out_of_frame=[]
        line_num=0
        for line in fhr:
            line_num+=1
            if line[0]=="@":continue
            if "3_prime_end" in line:
                line_num+=1
                continue
            if isKthBitSet(int(line.strip().split()[1]),4)==True:
                #line_num+=1
                continue
            read_id=line.strip().split()[0]
            transcript_id=line.strip().split()[2]
            next_line=fhr.readline().strip()
            next_read=next_line.split()[0]
            next_transcript_id=next_line.split()[2]
            if read_id!=next_read:
                print(line_num,read_id,next_read)
                print("PROBLEM")
            if transcript_id not in info:
                info[transcript_id]={"num_reads_out_of_frame":0,"reads_in_frame":[],"num_reads_in_frame":0}
            if next_transcript_id not in info:
                info[next_transcript_id]={"num_reads_out_of_frame":0,"reads_in_frame":[],"num_reads_in_frame":0}
                
            if "FN:0" in line and "FN:0" in next_line:
                info[transcript_id]["num_reads_out_of_frame"]+=1
            elif ("FI:" in line or "FI:" not in next_line) or ("FI:" not in line or "FI:" in next_line):
                if "FI" in line:
                    info[transcript_id]["reads_in_frame"].append(line.strip().split("FI:")[-1].strip())
                else:
                    info[next_transcript_id]["reads_in_frame"].append(next_line.strip().split("FI:")[-1].strip())
        fhr.close()
        for transcript_id in info:
            info[transcript_id]["num_reads_in_frame"]=len(info[transcript_id]["reads_in_frame"])
            info[transcript_id]["reads_in_frame"]=list(set(info[transcript_id]["reads_in_frame"]))
        
    else:
        fhr=open(fusion_reads_samfilename,"r")
        for line in fhr:
            if line[0]=="@" or "3_prime_end" in line:continue
            read_id=line.strip().split()[0]
            transcript_id=line.strip().split()[2]
            if transcript_id not in info:
                info[transcript_id]={"num_reads_out_of_frame":0,"reads_in_frame":[],"num_reads_in_frame":0}
            if "FN:0" in line:
                info[transcript_id]["num_reads_out_of_frame"]+=1
            else:
                info[transcript_id]["reads_in_frame"].append(line.strip().split("FI:")[-1].strip())
        fhr.close()
        for transcript_id in info:
            info[transcript_id]["num_reads_in_frame"]=len(info[transcript_id]["reads_in_frame"])
            info[transcript_id]["reads_in_frame"]=list(set(info[transcript_id]["reads_in_frame"]))
    
    """for transcript_id in info:
        print(fusion_reads_samfilename.split("/")[-1],transcript_id,info[transcript_id]["num_reads_in_frame"],info[transcript_id]["num_reads_out_of_frame"])"""
    return info
    
def generateReportFile(options,gene_info,transcripts_to_in_frame_STOP_codon_locations_in_5_prime_UTR_region):
    """
    Generates the complete information about all the transcripts and all replicates
    """
    # Check for in-frame fusion reads
    allinputs=[]
    pool = multiprocessing.Pool(processes=int(options.CPU))
    for num,eachtype in enumerate([options.selected_sample,options.background_sample]):
        for file_num,filename in enumerate(eachtype):
            if options.selected_ended=="PE" and options.background_ended=="PE":
                if file_num%2==1:
                    pass
                else:
                    continue
            if num==0:
                junction_reads_filename_bam=options.selected_sample_STAR_transcriptome_bamfilename_round2_fusion_reads[file_num]
            else:
                junction_reads_filename_bam=options.background_sample_STAR_transcriptome_bamfilename_round2_fusion_reads[file_num]
            allinputs.append([junction_reads_filename_bam,options,gene_info])
    pool.map(computeInFrameInfoForEachJunctionRead,allinputs)
    
    # Collect gene-level information about enrichment analysis
    enrichment_info=collectInformationAboutEnrichmentAnalysis(options.salmon_DGE_filename,
                                            options.salmon_normalized_counts,
                                            gene_info)
    
    # Collect all relevant information about each transcript and each replicate
    all_info={}
    for num,eachtype in enumerate([options.selected_sample,options.background_sample]):
        for file_num,filename in enumerate(eachtype):
            if options.selected_ended=="PE" and options.background_ended=="PE":
                if file_num%2==1:
                    pass
                else:
                    continue
            if num==0:
                fusion_reads_samfilename=options.selected_sample_STAR_transcriptome_bamfilename_round2_fusion_reads[file_num][:-4]+"_junction_info_appended.sam"
            else:
                fusion_reads_samfilename=options.background_sample_STAR_transcriptome_bamfilename_round2_fusion_reads[file_num][:-4]+"_junction_info_appended.sam"
            all_info[filename]=compileInformationForEachSample(fusion_reads_samfilename,options.selected_ended,options.CPU)
            for transcript_id in transcripts_to_in_frame_STOP_codon_locations_in_5_prime_UTR_region:
                if transcript_id not in all_info[filename]:
                    all_info[filename][transcript_id]={"num_reads_out_of_frame":-1,"num_reads_in_frame":-1,"reads_in_frame":[]}
    headers=["Transcript_id",
             "Replicate",
             "CDS_start",
             "CDS_end",
             "Length",
             "num_fusion_reads_in_frame_selected_sample",
             "num_fusion_reads_selected_sample",
             "num_fusion_reads_in_frame_background_sample",
             "num_fusion_reads_background_sample",
             "starting_location_of_in_frame_junction_reads",
             "log2FoldChange",
             "padj",
             "norm_counts_selected_sample",
             "norm_counts_background_sample",
             "location_of_STOP_codons_in_5_prime_UTR_region_in_frame_with_CDS"
             ]
    transcriptome=readFastaFile(options.transcriptome)
    
    """# Calculate in-frame scores
    for transcript_id in transcripts_to_in_frame_STOP_codon_locations_in_5_prime_UTR_region:
        gene=gene_info[transcript_id]["gene"]
        for selected_sample_filename_num,selected_sample_filename in enumerate(options.selected_sample):
            if options.selected_ended=="PE" and options.background_ended=="PE":
                if selected_sample_filename_num%2!=0:
                    pass
                else:
                    continue
            background_sample_filename=options.background_sample[selected_sample_filename_num]
            #print(transcript_id,gene,gene in enrichment_info)
            if all_info[selected_sample_filename][transcript_id]["num_reads_in_frame"]<2 or all_info[background_sample_filename][transcript_id]["num_reads_in_frame"]<2:continue
            calculate_in_frame_score(transcript_id,
                                     selected_sample_filename,
                                     all_info[selected_sample_filename][transcript_id]["num_reads_in_frame"],
                                     all_info[selected_sample_filename][transcript_id]["num_reads_out_of_frame"],
                                     all_info[background_sample_filename][transcript_id]["num_reads_in_frame"],
                                     all_info[background_sample_filename][transcript_id]["num_reads_out_of_frame"])
    """
    fhw=open(options.combined_graph_final,"w")
    fhw.write(",".join(headers)+"\n")
    for transcript_id in transcripts_to_in_frame_STOP_codon_locations_in_5_prime_UTR_region:
        gene=gene_info[transcript_id]["gene"]
        for selected_sample_filename_num,selected_sample_filename in enumerate(options.selected_sample):
            if options.selected_ended=="PE" and options.background_ended=="PE":
                if selected_sample_filename_num%2!=0:
                    pass
                else:
                    continue
            background_sample_filename=options.background_sample[selected_sample_filename_num]
            #print(transcript_id,gene,gene in enrichment_info)
            info_about_transcript=[transcript_id,
                  selected_sample_filename,
                  gene_info[transcript_id]["cds_start"],
                  gene_info[transcript_id]["cds_end"],
                  str(len(transcriptome[transcript_id])),
                  all_info[selected_sample_filename][transcript_id]["num_reads_in_frame"],
                  all_info[selected_sample_filename][transcript_id]["num_reads_in_frame"]+all_info[selected_sample_filename][transcript_id]["num_reads_out_of_frame"],
                  all_info[background_sample_filename][transcript_id]["num_reads_in_frame"],
                  all_info[background_sample_filename][transcript_id]["num_reads_in_frame"]+all_info[background_sample_filename][transcript_id]["num_reads_out_of_frame"],
                  ";".join(all_info[selected_sample_filename][transcript_id]["reads_in_frame"]),
                  enrichment_info[gene]["log2FoldChange"],
                  enrichment_info[gene]["padj"],
                  enrichment_info[gene]["norm_counts"][options.selected_sample[selected_sample_filename_num-1]],
                  enrichment_info[gene]["norm_counts"][options.background_sample[selected_sample_filename_num-1]],
                  ";".join(list(map(str,transcripts_to_in_frame_STOP_codon_locations_in_5_prime_UTR_region[transcript_id])))
                  ]
            fhw.write(",".join(list(map(str,info_about_transcript)))+"\n")
    fhw.close()

def generateRunReportFile(options):
    """
    Generates a complete report of the execution of the program
    This will record the duration of execution of each step and 
    also important details about each of the steps
    """
    fhw=open(options.run_details_info_csv,"w")
    all_information={}
    
    all_information[options.output_directory.split("/")[-1]]={}
    all_information[options.output_directory.split("/")[-1]]["genome_browser_viewing_bamfile"]=options.selected_sample_genome_browser
    all_information[options.output_directory.split("/")[-1]]["final_info_file"]=options.combined_graph_final
    all_information[options.output_directory.split("/")[-1]]["salmon_gene_raw_counts"]=options.salmon_gene_counts_matrix
    all_information[options.output_directory.split("/")[-1]]["deseq2_gene_norm_counts"]=options.deseq2_normalized_counts
    all_information[options.output_directory.split("/")[-1]]["run_details_info_csv"]=options.run_details_info_csv
    all_information[options.output_directory.split("/")[-1]]["record_time"]=options.record_time
    fhw.write(",")
    for num_eachtype,eachtype in enumerate([options.selected_sample,options.background_sample]):
        for file_num,filename in enumerate(eachtype):
            if options.selected_ended=="PE" and options.background_ended=="PE":
                if file_num%2==1:
                    pass
                else:
                    continue
            fhw.write(filename+",")
            if filename not in all_information:
                all_information[filename]={}
            
            """
            Read trimming information
            """
            if num_eachtype==0:
                inputfilename=options.selected_sample_adapter_trimmed_error_file[file_num]
            else:
                inputfilename=options.background_sample_adapter_trimmed_error_file[file_num]
            fhr=open(inputfilename,"r")
            for line in fhr:
                if "Input Reads" in line:
                    line=line.strip()
                    total_reads=int(line.split("Input Reads:")[-1].split()[0].strip())
                    reads_dropped=int(line.split("Dropped:")[-1].split()[0].strip())
                    reads_surviving=int(line.split("Surviving:")[-1].split()[0].strip())
                    all_information[filename]["total_reads"]=total_reads
                    all_information[filename]["num_reads_dropped"]=reads_dropped
                    all_information[filename]["num_reads_after_adapter_trimming"]=reads_surviving
            fhr.close()
            
            """
            Number of reads mapped for trimming
            """
            if num_eachtype==0:
                logfilename=options.selected_sample_STAR_prefix_round1[file_num]+"Log.final.out"
            else:
                logfilename=options.background_sample_STAR_prefix_round1[file_num]+"Log.final.out"
            fhr=open(logfilename,"r")
            for line in fhr:
                if "Uniquely mapped reads" in line:
                    umr=line.strip().split("|")[-1].strip()[:-1]
                elif "% of reads mapped to multiple loci" in line:
                    mmr=line.strip().split("|")[-1].strip()[:-1]
                elif "% of reads mapped to too many loci" in line:
                    mmur=line.strip().split("|")[-1].strip()[:-1]
                elif "% of reads unmapped: too many mismatches" in line or "% of reads unmapped: too short" in line or "% of reads unmapped: other" in line:
                    ur=line.strip().split("|")[-1].strip()[:-1]
            fhr.close()
            all_information[filename]["umr"]=umr
            all_information[filename]["mmr"]=mmr
            all_information[filename]["mmur"]=mmur
            all_information[filename]["ur"]=ur
            
            """
            Number of reads having vector sequence
            """
            if num_eachtype==0:
                inputfilename=options.selected_sample_trimming_stats[file_num]
            else:
                inputfilename=options.background_sample_trimming_stats[file_num]
            #print(inputfilename)
            fhr=open(inputfilename,"r")
            five_prime_forward,five_prime_reverse,five_prime_vector_whole_read,three_prime_forward,three_prime_reverse,three_prime_vector_whole_read,yeast_plasmid,yeast_chromosome,Reads_removed_due_to_small_size,Reads_removed_due_to_polyA,num_junction_reads=[int(line.split()[-1]) for line in fhr if " " in line]
            all_information[filename]["five_prime_forward_vector"]=five_prime_forward
            all_information[filename]["five_prime_reverse_vector"]=five_prime_reverse
            all_information[filename]["three_prime_forward_vector"]=three_prime_forward
            all_information[filename]["three_prime_reverse_vector"]=three_prime_reverse
            all_information[filename]["hits_to_yeast_chromosome"]=yeast_chromosome
            all_information[filename]["hits_to_yeast_plasmid"]=yeast_plasmid
            all_information[filename]["five_prime_vector_whole_read"]=five_prime_vector_whole_read
            all_information[filename]["three_prime_vector_whole_read"]=three_prime_vector_whole_read
            all_information[filename]["Reads_removed_due_to_small_size"]=Reads_removed_due_to_small_size
            all_information[filename]["Reads_removed_due_to_polyA"]=Reads_removed_due_to_polyA
            all_information[filename]["num_junction_reads"]=num_junction_reads
            #all_information[filename]["num_all_reads"]=num_all_reads
            fhr.close()
            
            """
            Genome Browser viewing files
            """
            if num_eachtype==0:
                all_information[filename]["genome_browser_viewing_bamfile_per_replicate"]=options.selected_sample_genome_browser_per_replicate[file_num]
            else:
                all_information[filename]["genome_browser_viewing_bamfile_per_replicate"]=""
    #pprint.pprint(all_information)
    """if options.selected_ended=="PE" and options.background_ended=="PE":
                if file_num%2==1:
                    pass
                else:
                    continue"""
    if options.selected_ended=="PE" and options.background_ended=="PE":
        fhw.write("\n")
        fhw.write("Total reads,"+",".join(list(map(str,[all_information[filename]["total_reads"] for num_eachtype,eachtype in enumerate([options.selected_sample,options.background_sample]) for file_num,filename in enumerate(eachtype) if file_num%2==1 ]))))
        fhw.write("\n")
        fhw.write("Reads retained after Adapter Trimming,"+",".join(list(map(str,[all_information[filename]["num_reads_after_adapter_trimming"] for num_eachtype,eachtype in enumerate([options.selected_sample,options.background_sample]) for file_num,filename in enumerate(eachtype) if file_num%2==1 ]))))
        fhw.write("\n")
        fhw.write("% of Uniquely mapped reads,"+",".join(list(map(str,[all_information[filename]["umr"] for num_eachtype,eachtype in enumerate([options.selected_sample,options.background_sample]) for file_num,filename in enumerate(eachtype) if file_num%2==1 ]))))
        fhw.write("\n")
        fhw.write("% of Multi mapped reads,"+",".join(list(map(str,[all_information[filename]["mmr"] for num_eachtype,eachtype in enumerate([options.selected_sample,options.background_sample]) for file_num,filename in enumerate(eachtype) if file_num%2==1]))))
        fhw.write("\n")
        fhw.write("% of reads unmapped due to too many loci,"+",".join(list(map(str,[all_information[filename]["mmur"] for num_eachtype,eachtype in enumerate([options.selected_sample,options.background_sample]) for file_num,filename in enumerate(eachtype) if file_num%2==1]))))
        fhw.write("\n")
        fhw.write("% of reads unmapped,"+",".join(list(map(str,[all_information[filename]["ur"] for num_eachtype,eachtype in enumerate([options.selected_sample,options.background_sample]) for file_num,filename in enumerate(eachtype) if file_num%2==1]))))
        fhw.write("\n")
        fhw.write("Num of fusion reads with 5' vector in forward orientation,"+",".join(list(map(str,[all_information[filename]["five_prime_forward_vector"] for num_eachtype,eachtype in enumerate([options.selected_sample,options.background_sample]) for file_num,filename in enumerate(eachtype) if file_num%2==1]))))
        fhw.write("\n")
        fhw.write("Num of fusion reads with 5' vector in reverse orientation,"+",".join(list(map(str,[all_information[filename]["five_prime_reverse_vector"] for num_eachtype,eachtype in enumerate([options.selected_sample,options.background_sample]) for file_num,filename in enumerate(eachtype) if file_num%2==1]))))
        fhw.write("\n")
        fhw.write("Num of fusion reads with 3' vector in forward orientation,"+",".join(list(map(str,[all_information[filename]["three_prime_forward_vector"] for num_eachtype,eachtype in enumerate([options.selected_sample,options.background_sample]) for file_num,filename in enumerate(eachtype) if file_num%2==1]))))
        fhw.write("\n")
        fhw.write("Num of fusion reads with 3' vector in reverse orientation,"+",".join(list(map(str,[all_information[filename]["three_prime_reverse_vector"] for num_eachtype,eachtype in enumerate([options.selected_sample,options.background_sample]) for file_num,filename in enumerate(eachtype) if file_num%2==1]))))
        fhw.write("\n")
        fhw.write("Num of reads mapped to yeast plasmid,"+",".join(list(map(str,[all_information[filename]["hits_to_yeast_plasmid"] for num_eachtype,eachtype in enumerate([options.selected_sample,options.background_sample]) for file_num,filename in enumerate(eachtype) if file_num%2==1]))))
        fhw.write("\n")
        fhw.write("Num of reads discarded due to entire 5' vector,"+",".join(list(map(str,[all_information[filename]["five_prime_vector_whole_read"] for num_eachtype,eachtype in enumerate([options.selected_sample,options.background_sample]) for file_num,filename in enumerate(eachtype) if file_num%2==1]))))
        fhw.write("\n")
        fhw.write("Num of reads discarded due to entire 3' vector,"+",".join(list(map(str,[all_information[filename]["three_prime_vector_whole_read"] for num_eachtype,eachtype in enumerate([options.selected_sample,options.background_sample]) for file_num,filename in enumerate(eachtype) if file_num%2==1]))))
        fhw.write("\n")
        fhw.write("Num of reads discarded due very small size,"+",".join(list(map(str,[all_information[filename]["Reads_removed_due_to_small_size"] for num_eachtype,eachtype in enumerate([options.selected_sample,options.background_sample]) for file_num,filename in enumerate(eachtype) if file_num%2==1]))))
        fhw.write("\n")
        fhw.write("Num of reads discarded due to polyA,"+",".join(list(map(str,[all_information[filename]["Reads_removed_due_to_polyA"] for num_eachtype,eachtype in enumerate([options.selected_sample,options.background_sample]) for file_num,filename in enumerate(eachtype) if file_num%2==1]))))
        fhw.write("\n")
        fhw.write("Num of fusion reads,"+",".join(list(map(str,[all_information[filename]["num_junction_reads"] for num_eachtype,eachtype in enumerate([options.selected_sample,options.background_sample]) for file_num,filename in enumerate(eachtype) if file_num%2==1]))))
        fhw.write("\n")
        fhw.write("Location of bamfile for viewing on Genome Browser,"+",".join([all_information[filename]["genome_browser_viewing_bamfile_per_replicate"] for num_eachtype,eachtype in enumerate([options.selected_sample,options.background_sample]) for file_num,filename in enumerate(eachtype) if file_num%2==1]))
        fhw.write("\n")
        fhw.write("Location of bamfile for all replicates,"+all_information[options.output_directory.split("/")[-1]]["genome_browser_viewing_bamfile"])
        fhw.write("\n")
        fhw.write("Location of reports file,"+all_information[options.output_directory.split("/")[-1]]["final_info_file"])
        fhw.write("\n")
        fhw.write("Salmon raw counts,"+all_information[options.output_directory.split("/")[-1]]["salmon_gene_raw_counts"])
        fhw.write("\n")
        fhw.write("DESeq2 normalized countsfile,"+all_information[options.output_directory.split("/")[-1]]["deseq2_gene_norm_counts"])
        fhw.write("\n")
    else:
        fhw.write("\n")
        fhw.write("Total reads,"+",".join(list(map(str,[all_information[filename]["total_reads"] for num_eachtype,eachtype in enumerate([options.selected_sample,options.background_sample]) for file_num,filename in enumerate(eachtype) ]))))
        fhw.write("\n")
        fhw.write("Reads retained after Adapter Trimming,"+",".join(list(map(str,[all_information[filename]["num_reads_after_adapter_trimming"] for num_eachtype,eachtype in enumerate([options.selected_sample,options.background_sample]) for file_num,filename in enumerate(eachtype) ]))))
        fhw.write("\n")
        fhw.write("% of Uniquely mapped reads,"+",".join(list(map(str,[all_information[filename]["umr"] for num_eachtype,eachtype in enumerate([options.selected_sample,options.background_sample]) for file_num,filename in enumerate(eachtype) ]))))
        fhw.write("\n")
        fhw.write("% of Multi mapped reads,"+",".join(list(map(str,[all_information[filename]["mmr"] for num_eachtype,eachtype in enumerate([options.selected_sample,options.background_sample]) for file_num,filename in enumerate(eachtype) ]))))
        fhw.write("\n")
        fhw.write("% of reads unmapped due to too many loci,"+",".join(list(map(str,[all_information[filename]["mmur"] for num_eachtype,eachtype in enumerate([options.selected_sample,options.background_sample]) for file_num,filename in enumerate(eachtype) ]))))
        fhw.write("\n")
        fhw.write("% of reads unmapped,"+",".join(list(map(str,[all_information[filename]["ur"] for num_eachtype,eachtype in enumerate([options.selected_sample,options.background_sample]) for file_num,filename in enumerate(eachtype) ]))))
        fhw.write("\n")
        fhw.write("Num of fusion reads with 5' vector in forward orientation,"+",".join(list(map(str,[all_information[filename]["five_prime_forward_vector"] for num_eachtype,eachtype in enumerate([options.selected_sample,options.background_sample]) for file_num,filename in enumerate(eachtype) ]))))
        fhw.write("\n")
        fhw.write("Num of fusion reads with 5' vector in reverse orientation,"+",".join(list(map(str,[all_information[filename]["five_prime_reverse_vector"] for num_eachtype,eachtype in enumerate([options.selected_sample,options.background_sample]) for file_num,filename in enumerate(eachtype) ]))))
        fhw.write("\n")
        fhw.write("Num of fusion reads with 3' vector in forward orientation,"+",".join(list(map(str,[all_information[filename]["three_prime_forward_vector"] for num_eachtype,eachtype in enumerate([options.selected_sample,options.background_sample]) for file_num,filename in enumerate(eachtype) ]))))
        fhw.write("\n")
        fhw.write("Num of fusion reads with 3' vector in reverse orientation,"+",".join(list(map(str,[all_information[filename]["three_prime_reverse_vector"] for num_eachtype,eachtype in enumerate([options.selected_sample,options.background_sample]) for file_num,filename in enumerate(eachtype) ]))))
        fhw.write("\n")
        fhw.write("Num of reads mapped to yeast plasmid,"+",".join(list(map(str,[all_information[filename]["hits_to_yeast_plasmid"] for num_eachtype,eachtype in enumerate([options.selected_sample,options.background_sample]) for file_num,filename in enumerate(eachtype) ]))))
        fhw.write("\n")
        fhw.write("Num of reads discarded due to entire 5' vector,"+",".join(list(map(str,[all_information[filename]["five_prime_vector_whole_read"] for num_eachtype,eachtype in enumerate([options.selected_sample,options.background_sample]) for file_num,filename in enumerate(eachtype) ]))))
        fhw.write("\n")
        fhw.write("Num of reads discarded due to entire 3' vector,"+",".join(list(map(str,[all_information[filename]["three_prime_vector_whole_read"] for num_eachtype,eachtype in enumerate([options.selected_sample,options.background_sample]) for file_num,filename in enumerate(eachtype) ]))))
        fhw.write("\n")
        fhw.write("Num of reads discarded due very small size,"+",".join(list(map(str,[all_information[filename]["Reads_removed_due_to_small_size"] for num_eachtype,eachtype in enumerate([options.selected_sample,options.background_sample]) for file_num,filename in enumerate(eachtype) ]))))
        fhw.write("\n")
        fhw.write("Num of reads discarded due to polyA,"+",".join(list(map(str,[all_information[filename]["Reads_removed_due_to_polyA"] for num_eachtype,eachtype in enumerate([options.selected_sample,options.background_sample]) for file_num,filename in enumerate(eachtype) ]))))
        fhw.write("\n")
        fhw.write("Num of fusion reads,"+",".join(list(map(str,[all_information[filename]["num_junction_reads"] for num_eachtype,eachtype in enumerate([options.selected_sample,options.background_sample]) for file_num,filename in enumerate(eachtype) ]))))
        fhw.write("\n")
        fhw.write("Location of bamfile for viewing on Genome Browser,"+",".join([all_information[filename]["genome_browser_viewing_bamfile_per_replicate"] for num_eachtype,eachtype in enumerate([options.selected_sample,options.background_sample]) for file_num,filename in enumerate(eachtype) ]))
        
        fhw.write("Location of bamfile for all replicates,"+all_information[options.output_directory.split("/")[-1]]["genome_browser_viewing_bamfile"])
        fhw.write("\n")
        fhw.write("Location of reports file,"+all_information[options.output_directory.split("/")[-1]]["final_info_file"])
        fhw.write("\n")
        fhw.write("Salmon raw counts,"+all_information[options.output_directory.split("/")[-1]]["salmon_gene_raw_counts"])
        fhw.write("\n")
        fhw.write("DESeq2 normalized countsfile,"+all_information[options.output_directory.split("/")[-1]]["deseq2_gene_norm_counts"])
        fhw.write("\n")
    
    fhw.write("Running times")
    fhw.write("\n")
    for eachfunction in all_information[options.output_directory.split("/")[-1]]["record_time"]:
        fhw.write(eachfunction+","+all_information[options.output_directory.split("/")[-1]]["record_time"][eachfunction].replace(",",":"))
        fhw.write("\n")
    fhw.close()

