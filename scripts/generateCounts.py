
import os
import subprocess


def runSalmonToGenerateCounts(options,logger_proxy,logging_mutex):
    """
    Use Salmon to generate counts. Then perform differential analysis.
    """
    # Combine Salmon fragments to one
    """cmd="cat lib/salmon-latest_linux_x86_64/bin/salmonaa lib/salmon-latest_linux_x86_64/bin/salmonab lib/salmon-latest_linux_x86_64/bin/salmonac > "+options.output_directory+"/salmon"
    os.system(cmd)
    cmd="chmod 777 "+options.output_directory+"/salmon"
    os.system(cmd)"""
    for num,eachtype in enumerate([options.selected_sample_N_removed,options.background_sample_N_removed]):
        for file_num,filename in enumerate(eachtype):
            if options.selected_ended=="PE" and options.background_ended=="PE":
                if file_num%2==1:
                    pass
                else:
                    continue
            if num==0:  
                inputfilename=options.selected_sample_STAR_transcriptome_bamfilename_round2[file_num]
                outputfilename=options.selected_sample_salmon_counts_outputfile[file_num]
                errorfilename=options.selected_sample_salmon_counts_error[file_num]
            else:
                inputfilename=options.background_sample_STAR_transcriptome_bamfilename_round2[file_num]
                outputfilename=options.background_sample_salmon_counts_outputfile[file_num]
                errorfilename=options.background_sample_salmon_counts_error[file_num]
            cmd="salmon quant "
            #cmd=options.output_directory+"/salmon quant "
            cmd+=" -t "+options.transcriptome
            cmd+=" --libType A "
            cmd+=" -a "+inputfilename
            cmd+=" -p "+str(options.CPU)+" "
            cmd+=" -o "+outputfilename
            cmd+=" -g "+options.transcript_to_gene_map +" "
            cmd+=" -s "
            cmd+=" 2> "+errorfilename
            os.system(cmd)
            with logging_mutex:
                logger_proxy.info("Salmon run completed for "+inputfilename)
            #logger.info("Salmon run completed for "+inputfilename)
            """cmd="lib/samtools/samtools sort -@ "+str(options.CPU)+" "+inputfilename+" > "+inputfilename+"sorted"
            os.system(cmd)
            cmd="mv "+inputfilename+"sorted "+inputfilename
            os.system(cmd)
            cmd="lib/samtools/samtools index "+inputfilename
            os.system(cmd)"""
    
    cmd="paste <(cut -f1 "+outputfilename+"/quant.genes.sf|tail -n +2) "
    for num,eachtype in enumerate([options.selected_sample_N_removed,options.background_sample_N_removed]):
        for file_num,filename in enumerate(eachtype):
            if options.selected_ended=="PE" and options.background_ended=="PE":
                if file_num%2==1:
                    pass
                else:
                    continue
            if num==0:  
                outputfilename=options.selected_sample_salmon_counts_outputfile[file_num]
            else:
                outputfilename=options.background_sample_salmon_counts_outputfile[file_num]
            cmd+=" <(cut -f4 "+outputfilename+"/quant.genes.sf|tail -n +2) "
    cmd+=" > "+options.salmon_gene_counts_matrix
    subprocess.check_call(['bash', '-c', cmd])

