#! /usr/bin/env python

import os

def runTrimmomaticToTrimOffAdapters(options,logger_proxy,logging_mutex):
    """
    No quality trimming will be done. Only removal of 
    adapter sequence is essential
    """
    for num,eachtype in enumerate([options.selected_filename,options.background_sample]):
        for file_num,filename in enumerate(eachtype):
            if num==0:
                if options.selected_path[file_num]=="":
                    inputfilename=filename+".fastq "
                else:
                    inputfilename=options.selected_path[file_num]+"/"+filename+".fastq "
                outputfilename=options.selected_sample_adapter_trimmed[file_num]
                errorfilename=options.selected_sample_adapter_trimmed_error_file[file_num]
            else:
                if options.background_path[file_num]=="":
                    inputfilename=filename+".fastq "
                else:
                    inputfilename=options.background_path[file_num]+"/"+filename+".fastq "
                outputfilename=options.background_sample_adapter_trimmed[file_num]
                errorfilename=options.background_sample_adapter_trimmed_error_file[file_num]
            cmd="trimmomatic SE -phred33 -threads "+str(options.CPU)+" "+inputfilename
            cmd+=" "+outputfilename
            cmd+=" ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 2> "
            cmd+=errorfilename
            os.system(cmd)
            with logging_mutex:
                logger_proxy.info(f"Running cmd - {cmd}")
                logger_proxy.info("Trimmomatic run for "+filename+" completed")
            