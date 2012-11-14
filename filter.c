#include "generic.h"

int filter_usage(){
    fprintf(stderr, "\n");
    fprintf(stderr, "Obtain alignment statistics of individual loci of each repeat subfamily, family or class.\n\n");
    fprintf(stderr, "Usage:   iteres filter [options] <chromosome size file> <repeat size file> <rmsk.txt> <bam/sam alignment file>\n\n");
    fprintf(stderr, "Options: -S       input is SAM [off]\n");
    fprintf(stderr, "         -Q       mapping Quality threshold [10]\n");
    fprintf(stderr, "         -N       normalized by number of (1: non-redundant reads, 2: mapped reads, 3: total reads) [1])\n");
    fprintf(stderr, "         -n       use repName (subfamily) as filter [null]\n");
    fprintf(stderr, "         -f       use repFamily as filter [null]\n");
    fprintf(stderr, "         -c       use repClass as filter [null]\n");
    fprintf(stderr, "         -t       only output repeats have more than [1] reads mapped\n");
    fprintf(stderr, "         -r       output the list of reads [off]\n");
    fprintf(stderr, "         -D       remove redundant reads [off]\n");
    fprintf(stderr, "         -C       Add 'chr' string as prefix of reference sequence [off]\n");
    fprintf(stderr, "         -P       Input was paired end data [off]\n");
    fprintf(stderr, "         -I       Insert length [1000] (force -P)\n");
    fprintf(stderr, "         -o       output prefix [basename of input without extension]\n");
    fprintf(stderr, "         -h       help message\n");
    fprintf(stderr, "         -?       help message\n");
    fprintf(stderr, "\n");
    return 1;
}

/* main filter function */
int main_filter(int argc, char *argv[]){
    char *output, *subfam, *out, *outReport;
    unsigned long long int *cnt;
    struct hash *hashRmsk = newHash(0);
    struct hash *hashRep = newHash(0);
    struct hash *hashFam = newHash(0);
    struct hash *hashCla = newHash(0);
    int optSam = 0, optthreshold = 1;
    char *optoutput = NULL, *optname = NULL, *optclass = NULL, *optfamily = NULL;
    unsigned int optreadlist = 0, optQual = 10, optisize = 1000;
    int filterField = 0, c, optDup = 0, optNorm = 1, optaddChr = 0, optPair = 0;
    
    time_t start_time, end_time;
    start_time = time(NULL);
    
    while ((c = getopt(argc, argv, "SQ:N:n:c:f:rDt:CPI:o:h?")) >= 0) {
        switch (c) {
            case 'S': optSam = 1; break;
            case 'Q': optQual = (unsigned int)strtol(optarg, 0, 0); break;
            case 'N': optNorm = (unsigned int)strtol(optarg, 0, 0); break;
            case 't': optthreshold = (unsigned int)strtol(optarg, 0, 0); break;
            case 'r': optreadlist = 1; break;
            case 'D': optDup = 1; break;
            case 'C': optaddChr = 1; break;
            case 'n': optname = strdup(optarg); break;
            case 'c': optclass = strdup(optarg); break;
            case 'f': optfamily = strdup(optarg); break;
            case 'P': optPair = 1; break;
            case 'I': optisize = (unsigned int)strtol(optarg, 0, 0); optPair = 1; break;
            case 'o': optoutput = strdup(optarg); break;
            case 'h':
            case '?': return filter_usage(); break;
            default: return 1;
        }
    }
    if (optind + 4 > argc)
        return filter_usage();

    char *chr_size_file = argv[optind];
    char *rep_size_file = argv[optind+1];
    char *rmsk_file = argv[optind+2];
    char *sam_file = argv[optind+3];
    
    if ( (optname && optclass) || (optname && optfamily) || (optclass && optfamily) || (optname && optclass && optfamily))
        errAbort("Please specify only one filter, either -n, -c or -f.");
    
    int nindex = 0;
    if (optNorm == 0){
        nindex = 4;
    } else if (optNorm == 1){
        nindex = 3;
    } else if (optNorm == 2) {
        nindex = 1;    
    } else if (optNorm == 3) {
        nindex = 0;    
    } else{
        errAbort("Wrong normalization method specified");
    }
    
    subfam = cloneString("ALL");
    if (optname) {
        optclass = NULL;
        optfamily = NULL;
        subfam = cloneString(optname);
        filterField = 10;
    }else if (optclass) {
        optname = NULL;
        optfamily = NULL;
        subfam = cloneString(optclass);
        filterField = 11;
    } else if (optfamily) {
        optname = NULL;
        optclass= NULL;
        subfam = cloneString(optfamily);
        filterField = 12;
    }
    if (sameString(subfam, "ALL")){
        fprintf(stderr, "* You didn't specify any filter, will output all repeats\n");
        filterField = 0;
    }
    
    if(optoutput) {
        output = optoutput;
    } else {
        output = cloneString(get_filename_without_ext(basename(sam_file)));
    }
    
    struct hash *chrHash = hashNameIntFile(chr_size_file);
    struct hash *repHash = hashNameIntFile(rep_size_file);
    
    fprintf(stderr, "* Start to parse the rmsk file\n");
    rmsk2binKeeperHash(rmsk_file, chrHash, repHash, &hashRmsk, &hashRep, &hashFam, &hashCla, filterField, subfam);
    
    //sam file
    fprintf(stderr, "* Start to parse the SAM/BAM file\n");
    if (optPair){
        cnt = PEsamFile2nodupRepbedFile(sam_file, chrHash, hashRmsk, hashRep, hashFam, hashCla, optSam, optQual, 1, optDup, optaddChr, optisize);
    } else {
        cnt = samFile2nodupRepbedFile(sam_file, chrHash, hashRmsk, hashRep, hashFam, hashCla, optSam, optQual, 1, optDup, optaddChr);
    }

    fprintf(stderr, "* Preparing the output file\n");
    if (asprintf(&out, "%s_%s.iteres.loci", output, subfam) < 0)
        errAbort("Preparing output wrong");
    if (asprintf(&outReport, "%s_%s.iteres.reportloci", output, subfam) < 0)
        errAbort("Preparing output wrong");
    
    writeFilterOut(hashRmsk, out, optreadlist, optthreshold, subfam, cnt[nindex]); 
    
    //write report file
    fprintf(stderr, "* Preparing report file\n");
    writeReport(outReport, cnt, optQual, subfam);
    
    hashFree(&chrHash);
    hashFree(&repHash);
    hashFree(&hashRmsk);
    hashFree(&hashRep);
    hashFree(&hashFam);
    hashFree(&hashCla);
    free(out);
    free(outReport);
    
    end_time = time(NULL);
    fprintf(stderr, "* Done, time used %.0f seconds.\n", difftime(end_time, start_time));
    return 0;    
}

