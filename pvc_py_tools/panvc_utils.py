def splitA2M(A2M_FILENAME):
    f = open(A2M_FILENAME,"r")
    out_num = 0
    for line in f:
        if ";" in line :
            continue
        if ">" in line :
            out_num = out_num + 1
            #close?
            out_file = open('recombinant.n'+str(out_num)+'.gapped', 'w')
            continue
        out_file.write(line)
    return

# what about newlines? shall I remove them?

def sequence_num_to_name(SAMPLES_FILENAME, CHRS_FILENAME, PLOIDITY, CHR_ID, SEQ_ID):
    assert (SEQ_ID >= 1)
    hap_suffix = ""
    if (SEQ_ID != 1) :
        offset = SEQ_ID - 2
        target_line = 1 + (offset//PLOIDITY)
        hap_suffix = "_" + str(offset % PLOIDITY)
    else:
        target_line = 0
        hap_suffix = ""

    with open(SAMPLES_FILENAME, 'r') as f:
        curr_name = f.readline()

    # Open and read file into buffer
    f = open(SAMPLES_FILENAME,"r")
    sample_names = f.readlines()

    total_chroms = sum(1 for line in open(CHRS_FILENAME))
    
    chr_suffix = ""
    if (total_chroms != 1) :
        chr_suffix="_CHR"+str(CHR_ID)

    #we are currently deleting whitespaces from the samples names;
    #if a replacement is needed just change the next line.
    SPACEREPLACEMENT=""
    #answer="PG_REF_"+str(SEQ_ID)+"_CHR"+str(CHR_ID)
    answer = sample_names[target_line].strip().replace(" ", SPACEREPLACEMENT)+hap_suffix+chr_suffix
    return answer
