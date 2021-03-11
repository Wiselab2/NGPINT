
import os
import sys

from scripts.readWriteOperations import *

def createTranscriptomeFileForPrimerDesign(options,gene_info):
    """
    """
    transcriptome=readFastaFile(options.transcriptome)
    if os.path.exists(options.functional_annotation)==True:
        functional_annotation={}
        fhr=open(options.functional_annotation,"r")
        for line in fhr:
            functional_annotation[line.strip().split()[0]]=line.strip().split()[1] if len(line.strip().split())>1 else ""
        fhr.close()
        
    #print(options.combined_graph_final)
    sys.stdout.flush()
    fhr=open(options.combined_graph_final,"r")
    fhw=open(options.design_primers_for_transcript,"w")
    for line_num,line in enumerate(fhr):
        if line_num==0:continue
        try:
            Transcript_id,Replicate,CDS_start,CDS_end,Length,num_fusion_reads_in_frame_selected_sample,num_fusion_reads_selected_sample,num_fusion_reads_in_frame_background_sample,num_fusion_reads_background_sample,starting_location_of_in_frame_junction_reads,log2FoldChange,padj,norm_counts_selected_sample,norm_counts_background_sample,location_of_STOP_codons_in_5_prime_UTR_region_in_frame_with_CDS=line.strip().split(",")
        except ValueError:
            print(line)
            return
        header=Transcript_id+"_"+Replicate+"_"+log2FoldChange+"_"+padj+"_"+CDS_start+"_"+CDS_end
        gene=gene_info[Transcript_id]["gene"]
        CDS_start=int(CDS_start)
        CDS_end=int(CDS_end)
        seq=transcriptome[Transcript_id]
        in_frame_fusion_seq=[" " for i in range(len(seq))]
        stop_codons=[" " for i in range(len(seq))]
        if location_of_STOP_codons_in_5_prime_UTR_region_in_frame_with_CDS!="":
            location_of_STOP_codons_in_5_prime_UTR_region_in_frame_with_CDS=list(map(int,location_of_STOP_codons_in_5_prime_UTR_region_in_frame_with_CDS.split(";")))
        else:
            location_of_STOP_codons_in_5_prime_UTR_region_in_frame_with_CDS=""
        if len(starting_location_of_in_frame_junction_reads)>0:
            starting_location_of_in_frame_junction_reads=list(map(int,list(set([ele.split(":")[1] for ele in starting_location_of_in_frame_junction_reads.split(";")]))))
        else:
            starting_location_of_in_frame_junction_reads=-1
        if starting_location_of_in_frame_junction_reads!=-1:
            for eachloc in starting_location_of_in_frame_junction_reads:
                if eachloc >= CDS_end:
                    print(line)
                    print("TROUBLE")
                """print(location_of_STOP_codons_in_5_prime_UTR_region_in_frame_with_CDS)
                sys.stdout.flush()"""
                in_frame_fusion_seq[eachloc-1]="*" if 0 not in list(set([(ele-eachloc)%3  for ele in location_of_STOP_codons_in_5_prime_UTR_region_in_frame_with_CDS if eachloc<ele])) else "X"
        else:
            continue
        in_frame_fusion_seq="".join(in_frame_fusion_seq)
        fhw.write(header+"\n")
        fhw.write(seq[:CDS_start-1].lower()+seq[CDS_start-1:CDS_end]+seq[CDS_end:].lower()+"\n")
        fhw.write(in_frame_fusion_seq[:CDS_start-1].lower()+in_frame_fusion_seq[CDS_start-1:CDS_end]+in_frame_fusion_seq[CDS_end:].lower()+"\n")
    fhr.close()
    fhw.close()
    
