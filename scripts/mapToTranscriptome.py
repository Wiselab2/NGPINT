

import os
def prepareTranscriptomeFiles(options,logger_proxy,logging_mutex):
    """
    """
    # Merge the transcriptome files from round1 and round2
    for num,eachtype in enumerate([options.selected_sample_N_removed,options.background_sample_N_removed]):
        for file_num,filename in enumerate(eachtype):
            if num==0:
                if options.selected_ended=="PE" and options.background_ended=="PE":
                    if file_num%2==1:
                        pass
                    else:
                        continue
                cmd="samtools cat -o "+options.selected_sample_STAR_transcriptome_bamfilename_round2[file_num]+"tempfile "
                cmd+=options.selected_sample_STAR_transcriptome_bamfilename_round1[file_num]+" "
                cmd+=options.selected_sample_STAR_transcriptome_bamfilename_round2[file_num]
                os.system(cmd)
                cmd="mv "+options.selected_sample_STAR_transcriptome_bamfilename_round2[file_num]+"tempfile "
                cmd+=options.selected_sample_STAR_transcriptome_bamfilename_round2[file_num]
                os.system(cmd)
                with logging_mutex:
                    logger_proxy.info("Preparing transcriptome files completed for "+options.selected_sample_STAR_transcriptome_bamfilename_round2[file_num])
            else:
                if options.selected_ended=="PE" and options.background_ended=="PE":
                    if file_num%2==1:
                        pass
                    else:
                        continue
                cmd="samtools cat -o "+options.background_sample_STAR_transcriptome_bamfilename_round2[file_num]+"tempfile "
                cmd+=options.background_sample_STAR_transcriptome_bamfilename_round1[file_num]+" "
                cmd+=options.background_sample_STAR_transcriptome_bamfilename_round2[file_num]
                os.system(cmd)
                cmd="mv "+options.background_sample_STAR_transcriptome_bamfilename_round2[file_num]+"tempfile "
                cmd+=options.background_sample_STAR_transcriptome_bamfilename_round2[file_num]
                os.system(cmd)
                with logging_mutex:
                    logger_proxy.info("Preparing transcriptome files completed for "+options.background_sample_STAR_transcriptome_bamfilename_round2[file_num])