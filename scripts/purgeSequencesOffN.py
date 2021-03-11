

import multiprocessing

def runRemoveNFromFastqParallel(inputs):
    """
    Parallel implementation of removeNFromFastq
    """
    filename,outputfilename,logger_proxy,logging_mutex=inputs
    fhr=open(filename,"r")
    fhw=open(outputfilename,"w")
    num_seq_with_N=0
    while True:
        line=fhr.readline().strip()
        if not line:break
        seq=fhr.readline().strip()
        useless=fhr.readline().strip()
        quality=fhr.readline().strip()
        #print(line_num)
        if seq[0]=="N" or seq[-1]=="N":
            rseq=seq.rstrip("N")
            right_trim=len(seq)-len(rseq)
            if right_trim!=0:
                quality=quality[:-right_trim]
            lseq=rseq.lstrip("N")
            left_trim=len(rseq)-len(lseq)
            quality=quality[left_trim:]
            if len(lseq)>0:
                fhw.write(line+"\n")
                fhw.write(lseq+"\n")
                fhw.write(useless+"\n")
                fhw.write(quality+"\n")
                num_seq_with_N+=1
        else:
            if len(seq)>0:
                fhw.write(line+"\n")
                fhw.write(seq+"\n")
                fhw.write(useless+"\n")
                fhw.write(quality+"\n")
    fhw.close()
    fhr.close()
    with logging_mutex:
        logger_proxy.info("Removal of flanked Ns for "+filename+" completed")

def removeNFromFastq(options,logger_proxy,logging_mutex):
    """
    Reads in the raw data files and remove leading and trailing Ns.
    """    
    #install_mp_handler()
    pool = multiprocessing.Pool(processes=int(options.CPU))
    allinputs=[]
    for num,eachtype in enumerate([options.selected_sample,options.background_sample]):
        for file_num,filename in enumerate(eachtype):
            if num==0:
                inputfilename=options.selected_sample_adapter_trimmed[file_num]
                outputfilename=options.selected_sample_N_removed[file_num]
            else:
                inputfilename=options.background_sample_adapter_trimmed[file_num]
                outputfilename=options.background_sample_N_removed[file_num]
            allinputs.append([inputfilename,outputfilename,logger_proxy,logging_mutex])
    pool.map(runRemoveNFromFastqParallel,allinputs)