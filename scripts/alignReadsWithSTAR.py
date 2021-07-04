
import os

def alignReadsWithStarForTrimming(options,logger_proxy,logging_mutex):
    """
    Aligns the untrimmed reads to genome and also generates the 
    mappings to transcriptome
    """
    cmd  = f"STAR "
    cmd += f" --genomeLoad Remove "
    os.system(cmd)
    
    for num,eachtype in enumerate([options.selected_sample_N_removed,options.background_sample_N_removed]):
        for file_num,filename in enumerate(eachtype):
            #filename=filename.split("/")[-1].split(".fastq")[0]
            cmd="STAR "
            cmd+=" --runThreadN "+str(options.CPU)+" "
            cmd+=" --genomeDir "+options.star_genome_index
            cmd+=" --genomeLoad LoadAndKeep "
            if options.selected_ended=="SE" and options.background_ended=="SE":
                cmd+=" --readFilesIn "+filename
            else:
                if file_num%2==1:
                    #print(num,file_num,eachtype[file_num-1],filename)
                    cmd+=" --readFilesIn "+eachtype[file_num-1]+" "+filename
                else:
                    continue
            if num==0:
                cmd+=" --outFileNamePrefix "+options.selected_sample_STAR_prefix_round1[file_num]
            else:
                cmd+=" --outFileNamePrefix "+options.background_sample_STAR_prefix_round1[file_num]
            cmd+=" --outSAMtype SAM "
            #cmd+=" --outReadsUnmapped Fastx "
            cmd+=" --outFilterMultimapNmax 500 "
            cmd+=" --limitOutSAMoneReadBytes 1000000 "
            cmd+=" --outFilterScoreMinOverLread 0.30 --outFilterMatchNminOverLread 0.30 "
            cmd+=" --alignIntronMax 10000 "
            cmd+=" --quantMode TranscriptomeSAM "
            cmd+=" --quantTranscriptomeBan Singleend "
            #cmd+=" --seedPerWindowNmax 100 "
            #cmd+=" --seedPerReadNmax 2000 "
            if num==0:
                cmd+=" > "+options.selected_sample_STAR_round1_output[file_num]
                cmd+=" 2> "+options.selected_sample_STAR_round1_error[file_num]
            else:
                cmd+=" > "+options.background_sample_STAR_round1_output[file_num]
                cmd+=" 2> "+options.background_sample_STAR_round1_error[file_num]
            os.system(cmd)
            
            if num==0:
                cmd="rm "+options.selected_sample_STAR_prefix_round1[file_num]+"Log.out "
                cmd+=options.selected_sample_STAR_prefix_round1[file_num]+"Log.progress.out "
                cmd+=options.selected_sample_STAR_prefix_round1[file_num]+"SJ.out.tab "
                cmd+=options.selected_sample_STAR_round1_output[file_num]
            else:
                cmd="rm "+options.background_sample_STAR_prefix_round1[file_num]+"Log.out "
                cmd+=options.background_sample_STAR_prefix_round1[file_num]+"Log.progress.out "
                cmd+=options.background_sample_STAR_prefix_round1[file_num]+"SJ.out.tab "
                cmd+=options.background_sample_STAR_round1_output[file_num]
            os.system(cmd)
            with logging_mutex:
                logger_proxy.info("STAR round1 mapping for "+filename+" completed")
            #logger.info("STAR round1 mapping for "+filename+" completed")
    cmd="STAR "
    cmd+=" --genomeLoad Remove "
    cmd+=" --genomeDir "+options.star_genome_index
    os.system(cmd)


def reAlignReadsMappedToVector(options):
    """
    Selects fusion reads which are not present in the transcriptome file
    and realign those to the genome and update both genome and
    transcriptome files
    """
    
    # Mapping fusion reads to genome
    for num,eachtype in enumerate([options.selected_sample_N_removed,options.background_sample_N_removed]):
        for file_num,filename in enumerate(eachtype):
            #filename=filename.split("/")[-1].split(".fastq")[0]
            if options.selected_ended=="PE" and options.background_ended=="PE":
                if file_num%2==0:continue
            cmd="STAR "
            cmd+=" --runThreadN "+str(options.CPU)+" "
            cmd+=" --genomeDir "+options.star_genome_index
            cmd+=" --genomeLoad LoadAndKeep "
            if num==0:
                if options.selected_ended=="SE" and options.background_ended=="SE":
                    cmd+=" --readFilesIn "+options.selected_sample_fusion_reads[file_num]
                else:
                    cmd+=" --readFilesIn "+options.selected_sample_fusion_reads[file_num-1]+" "+options.selected_sample_fusion_reads[file_num]
                cmd+=" --outFileNamePrefix "+options.selected_sample_STAR_prefix_round2[file_num]
            else:
                if options.selected_ended=="SE" and options.background_ended=="SE":
                    cmd+=" --readFilesIn "+options.background_sample_fusion_reads[file_num]
                else:
                    cmd+=" --readFilesIn "+options.background_sample_fusion_reads[file_num-1]+" "+options.background_sample_fusion_reads[file_num]
                cmd+=" --outFileNamePrefix "+options.background_sample_STAR_prefix_round2[file_num]
            cmd+=" --outSAMtype SAM "
            #cmd+=" --outReadsUnmapped Fastx "
            cmd+=" --outFilterMultimapNmax 500 "
            
            cmd+=" --outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3 "
            cmd+=" --alignIntronMax 10000 "
            #cmd+=" --quantMode TranscriptomeSAM "
            #cmd+=" --quantTranscriptomeBan Singleend "
            if num==0:
                cmd+=" > "+options.selected_sample_STAR_round2_output[file_num]
                cmd+=" 2> "+options.selected_sample_STAR_round2_error[file_num]
            else:
                cmd+=" > "+options.background_sample_STAR_round2_output[file_num]
                cmd+=" 2> "+options.background_sample_STAR_round2_error[file_num]
            os.system(cmd)
            #print(cmd)
            if num==0:
                cmd="rm "+options.selected_sample_STAR_prefix_round2[file_num]+"Log.out "
                cmd+=options.selected_sample_STAR_prefix_round2[file_num]+"Log.progress.out "
                cmd+=options.selected_sample_STAR_prefix_round2[file_num]+"SJ.out.tab "
                cmd+=options.selected_sample_STAR_round2_output[file_num]
            else:
                cmd="rm "+options.background_sample_STAR_prefix_round2[file_num]+"Log.out "
                cmd+=options.background_sample_STAR_prefix_round2[file_num]+"Log.progress.out "
                cmd+=options.background_sample_STAR_prefix_round2[file_num]+"SJ.out.tab "
                cmd+=options.background_sample_STAR_round2_output[file_num]
            os.system(cmd)
            
            if num==0 and os.stat(options.selected_sample_STAR_round2_error[file_num]).st_size == 0:
                cmd="rm "+options.selected_sample_STAR_round2_error[file_num]
            elif num==1 and os.stat(options.background_sample_STAR_round2_error[file_num]).st_size == 0:
                cmd="rm "+options.background_sample_STAR_round2_error[file_num]
            os.system(cmd)
            
            #logger.info("STAR round2 mapping for "+filename+" completed")
            
            """if num==0:
                cmd="cp "+options.selected_sample_STAR_transcriptome_bamfilename_round2[file_num]+" "+options.selected_sample_STAR_transcriptome_bamfilename_round2_fusion_reads[file_num]
                os.system(cmd)
            else:
                cmd="cp "+options.background_sample_STAR_transcriptome_bamfilename_round2[file_num]+" "+options.background_sample_STAR_transcriptome_bamfilename_round2_fusion_reads[file_num]
                os.system(cmd)"""
            
    cmd="STAR "
    cmd+=" --genomeLoad Remove "
    cmd+=" --genomeDir "+options.star_genome_index
    os.system(cmd)
    
    
    # Mapping fusion reads to transcriptome - Need to repeat this step since STAR is buggy 
    for num,eachtype in enumerate([options.selected_sample_N_removed,options.background_sample_N_removed]):
        for file_num,filename in enumerate(eachtype):
            #filename=filename.split("/")[-1].split(".fastq")[0]
            if options.selected_ended=="PE" and options.background_ended=="PE":
                if file_num%2==0:continue
            cmd="STAR "
            cmd+=" --runThreadN "+str(options.CPU)+" "
            #cmd+=" --genomeDir "+options.star_genome_index
            cmd+=" --genomeDir "+options.transcriptome_index
            cmd+=" --genomeLoad LoadAndKeep "
            if num==0:
                if options.selected_ended=="SE" and options.background_ended=="SE":
                    cmd+=" --readFilesIn "+options.selected_sample_fusion_reads[file_num]
                else:
                    cmd+=" --readFilesIn "+options.selected_sample_fusion_reads[file_num-1]+" "+options.selected_sample_fusion_reads[file_num]
                cmd+=" --outFileNamePrefix "+options.selected_sample_STAR_prefix_round2[file_num]+"_transcriptome_"
            else:
                if options.selected_ended=="SE" and options.background_ended=="SE":
                    cmd+=" --readFilesIn "+options.background_sample_fusion_reads[file_num]
                else:
                    cmd+=" --readFilesIn "+options.background_sample_fusion_reads[file_num-1]+" "+options.background_sample_fusion_reads[file_num]
                cmd+=" --outFileNamePrefix "+options.background_sample_STAR_prefix_round2[file_num]+"_transcriptome_"
            cmd+=" --outSAMtype SAM "
            #cmd+=" --outReadsUnmapped Fastx "
            cmd+=" --outFilterMultimapNmax 500 "
            cmd+=" --seedPerReadNmax 5000 "
            cmd+=" --seedPerWindowNmax 100 "
            cmd+=" --outFilterScoreMinOverLread 0.8 --outFilterMatchNminOverLread 0.85 "
            cmd+=" --alignIntronMax 10000 "
            #cmd+=" --quantMode TranscriptomeSAM "
            #cmd+=" --quantTranscriptomeBan Singleend "
            if num==0:
                cmd+=" > "+options.selected_sample_STAR_round2_output[file_num]
                cmd+=" 2> "+options.selected_sample_STAR_round2_error[file_num]
            else:
                cmd+=" > "+options.background_sample_STAR_round2_output[file_num]
                cmd+=" 2> "+options.background_sample_STAR_round2_error[file_num]
            #print(cmd)
            os.system(cmd)
            #print(cmd)
            #continue
            cmd="samtools view -bSh "
            if num==0:
                cmd+=" "+options.selected_sample_STAR_prefix_round2[file_num]+"_transcriptome_Aligned.out.sam "
                cmd+=" > "+options.selected_sample_STAR_prefix_round2[file_num]+"Aligned.toTranscriptome.out.bam"
            else:
                cmd+=" "+options.background_sample_STAR_prefix_round2[file_num]+"_transcriptome_Aligned.out.sam "
                cmd+=" > "+options.background_sample_STAR_prefix_round2[file_num]+"Aligned.toTranscriptome.out.bam"
            os.system(cmd)
            #print(cmd)
            
            if num==0:
                cmd="cp "+options.selected_sample_STAR_transcriptome_bamfilename_round2[file_num]+" "
                cmd+=options.selected_sample_STAR_transcriptome_bamfilename_round2_fusion_reads[file_num]
                #print(cmd)
                os.system(cmd)
            else:
                cmd="cp "+options.background_sample_STAR_transcriptome_bamfilename_round2[file_num]+" "
                cmd+=options.background_sample_STAR_transcriptome_bamfilename_round2_fusion_reads[file_num]
                #print(cmd)
                os.system(cmd)
            
            
            #logger.info("STAR round2 mapping for "+filename+" completed")
    cmd="STAR "
    cmd+=" --genomeLoad Remove "
    cmd+=" --genomeDir "+options.transcriptome_index
    os.system(cmd)
 
