#include "generic.h"

int8_t seq_comp_table[16] = { 0, 8, 4, 12, 2, 10, 9, 14, 1, 6, 5, 13, 3, 11, 7, 15 };

/* definitions of functions */

char *get_filename_without_ext(char *filename) {
    char *s;
    s = malloc(strlen(filename) + 1);
    strcpy(s, filename);
    char *dot = strrchr(s, '.');
    if(!dot || dot == s) return s;
    *dot = '\0';
    return s;
}

char *get_filename_ext(char *filename) {
    char *dot = strrchr(filename, '.');
    if(!dot || dot == filename) return "";
    return dot + 1;
}

double cal_rpkm (unsigned long long int reads_count, unsigned long long int total_length, unsigned long long int mapped_reads_num) {
    return reads_count / (mapped_reads_num * 1e-9 * total_length);
}

double cal_rpm (unsigned long long int reads_count, unsigned long long int mapped_reads_num) {
    return reads_count / (mapped_reads_num * 1e-6 );
}

void writeReport(char *outfile, unsigned long long int *cnt, unsigned int mapQ, char *subfam){
    FILE *f = mustOpen(outfile, "w");
    fprintf(f, "total reads (pair): %llu\n", cnt[0]);
    //fprintf(f, "read ends 1: %llu\n", cnt[0]);
    //fprintf(f, "read ends 2: %llu\n", cnt[1]);
    //fprintf(f, "mapped read ends 1: %llu\n", cnt[2]);
    //fprintf(f, "mapped read ends 2: %llu\n", cnt[3]);
    //fprintf(f, "used read ends 1: %llu\n", cnt[4]);
    //fprintf(f, "used read ends 2: %llu\n", cnt[5]);
    fprintf(f, "mappable reads (pair): %llu\n", cnt[6]);
    //fprintf(f, "non-redundant mappable reads (pair): %llu\n", cnt[8]);
    fprintf(f, "uniquely mapped reads (pair) (mapQ >= %u): %llu\n", mapQ, cnt[7]);
    fprintf(f, "non-redundant uniquely mapped reads (pair): %llu\n", cnt[11]);
    //fprintf(f, "non-redundant reads (pair) overlap with [%s] repeats: %llu\n", subfam, cnt[9]);
    fprintf(f, "non-redundant uniquely mapped reads (pair) overlap with [%s] repeats: %llu\n", subfam, cnt[10]);
    carefulClose(&f);
}

void writeReportDensity(char *outfile, unsigned long long int *cnt, unsigned int mapQ){
    FILE *f = mustOpen(outfile, "w");
    fprintf(f, "total reads (pair): %llu\n", cnt[0]);
    fprintf(f, "  read ends 1: %llu\n", cnt[0]);
    fprintf(f, "  read ends 2: %llu\n", cnt[1]);
    fprintf(f, "  mapped read ends 1: %llu\n", cnt[2]);
    fprintf(f, "  mapped read ends 2: %llu\n", cnt[3]);
    fprintf(f, "  used read ends 1: %llu\n", cnt[4]);
    fprintf(f, "  used read ends 2: %llu\n", cnt[5]);
    fprintf(f, "mappable reads (pair): %llu\n", cnt[6]);
    //fprintf(f, "non-redundant mappable reads (pair): %llu\n", cnt[8]);
    fprintf(f, "uniquely mapped reads (pair) (mapQ >= %u): %llu\n", mapQ, cnt[7]);
    fprintf(f, "non-redundant uniquely mapped reads (pair): %llu\n", cnt[9]);
    carefulClose(&f);
}

void writeWigandStat(struct hash *hash, struct hash *hash1, struct hash *hash2, char *of1, char *of2, char *of3, char *of4, char *of5, unsigned long long int reads_num, unsigned long int reads_num_unique){
    FILE *f1 = mustOpen(of1, "w");
    FILE *f2 = mustOpen(of2, "w");
    FILE *f5 = mustOpen(of5, "w");
    unsigned int m;
    struct hashEl *helr;
    struct hashCookie cookier = hashFirst(hash);
    fprintf(f1, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "#subfamily", "family", "class", "consensus_length", "reads_count", "unique_reads_count", "total_length", "genome_count", "all_reads_RPKM", "all_reads_RPM", "unique_reads_RPKM", "unique_reads_RPM");
    while ( (helr = hashNext(&cookier)) != NULL ) {
        struct rep *or = (struct rep *) (helr->val);
        fprintf(f1, "%s\t%s\t%s\t%u\t%llu\t%llu\t%llu\t%llu\t%.3f\t%.3f\t%.3f\t%.3f\n", or->name, or->fname, or->cname, or->length, or->read_count, or->read_count_unique, or->total_length, or->genome_count, cal_rpkm(or->read_count, or->total_length, reads_num), cal_rpm(or->read_count, reads_num), cal_rpkm(or->read_count_unique, or->total_length, reads_num_unique), cal_rpm(or->read_count_unique, reads_num_unique));
        if (or->length != 0){
            fprintf(f2, "fixedStep chrom=%s start=1 step=1 span=1\n", or->name);
            fprintf(f5, "fixedStep chrom=%s start=1 step=1 span=1\n", or->name);
            for (m = 0; m < or->length; m++){ 
                fprintf(f2, "%u\n", (or->bp_total)[m]);
                fprintf(f5, "%u\n", (or->bp_total_unique)[m]);
            }
        }
    }
    carefulClose(&f2);
    carefulClose(&f1);
    carefulClose(&f5);
    FILE *f3 = mustOpen(of3, "w");
    struct hashEl *hel3;
    struct hashCookie cookier3 = hashFirst(hash1);
    fprintf(f3, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "#family", "class", "reads_count", "unique_reads_count", "total_length", "genome_count", "all_reads_RPKM", "all_reads_RPM", "unique_reads_RPKM", "unique_reads_RPM");
    while ( (hel3 = hashNext(&cookier3)) != NULL) {
        struct repfam *or = (struct repfam *) hel3->val;
        fprintf(f3, "%s\t%s\t%llu\t%llu\t%llu\t%llu\t%.3f\t%.3f\t%.3f\t%.3f\n", or->fname, or->cname, or->read_count, or->read_count_unique, or->total_length, or->genome_count, cal_rpkm(or->read_count, or->total_length, reads_num), cal_rpm(or->read_count, reads_num), cal_rpkm(or->read_count_unique, or->total_length, reads_num_unique), cal_rpm(or->read_count_unique, reads_num_unique));
    }
    carefulClose(&f3);
    FILE *f4 = mustOpen(of4, "w");
    struct hashEl *hel4;
    struct hashCookie cookier4 = hashFirst(hash2);
    fprintf(f4, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",  "#class", "reads_count", "unique_reads_count", "total_length", "genome_count", "all_reads_RPKM", "all_reads_RPM", "unique_reads_RPKM", "unique_reads_RPM");
    while ( (hel4 = hashNext(&cookier4)) != NULL) {
        struct repcla *or = (struct repcla *) hel4->val;
        fprintf(f4, "%s\t%llu\t%llu\t%llu\t%llu\t%.3f\t%.3f\t%.3f\t%.3f\n", or->cname, or->read_count, or->read_count_unique, or->total_length, or->genome_count, cal_rpkm(or->read_count, or->total_length, reads_num), cal_rpm(or->read_count, reads_num), cal_rpkm(or->read_count_unique, or->total_length, reads_num_unique), cal_rpm(or->read_count_unique, reads_num_unique));
    }
    carefulClose(&f4);
}

void MREwriteWigandStat(struct hash *hash, struct hash *hash1, struct hash *hash2, char *of1, char *of2, char *of3, char *of4){
    FILE *f1 = mustOpen(of1, "w");
    FILE *f2 = mustOpen(of2, "w");
    unsigned int m;
    struct hashEl *helr;
    struct hashCookie cookier = hashFirst(hash);
    fprintf(f1, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "#subfamily", "family", "class", "consensus_length", "covered_CpG_sites", "CpG_total_score", "total_length", "genome_count");
    while ( (helr = hashNext(&cookier)) != NULL ) {
        struct rep *or = (struct rep *) (helr->val);
        fprintf(f1, "%s\t%s\t%s\t%u\t%u\t%.4f\t%llu\t%llu\n", or->name, or->fname, or->cname, or->length, or->cpgCount, or->cpgTotalScore, or->total_length, or->genome_count);
        if (or->length != 0){
            fprintf(f2, "fixedStep chrom=%s start=1 step=1 span=1\n", or->name);
            for (m = 0; m < or->length; m++){ 
                fprintf(f2, "%.4f\n", (or->cpgScore)[m]);
            }
        }
    }
    carefulClose(&f2);
    carefulClose(&f1);
    FILE *f3 = mustOpen(of3, "w");
    struct hashEl *hel3;
    struct hashCookie cookier3 = hashFirst(hash1);
    fprintf(f3, "%s\t%s\t%s\t%s\t%s\t%s\n", "#family", "class", "covered_CpG_sites", "CpG_total_score", "total_length", "genome_count");
    while ( (hel3 = hashNext(&cookier3)) != NULL) {
        struct repfam *or = (struct repfam *) hel3->val;
        fprintf(f3, "%s\t%s\t%u\t%.4f\t%llu\t%llu\n", or->fname, or->cname, or->cpgCount, or->cpgTotalScore, or->total_length, or->genome_count);
    }
    carefulClose(&f3);
    FILE *f4 = mustOpen(of4, "w");
    struct hashEl *hel4;
    struct hashCookie cookier4 = hashFirst(hash2);
    fprintf(f4, "%s\t%s\t%s\t%s\t%s\n", "#class", "covered_CpG_sites", "CpG_total_score", "total_length", "genome_count");
    while ( (hel4 = hashNext(&cookier4)) != NULL) {
        struct repcla *or = (struct repcla *) hel4->val;
        fprintf(f4, "%s\t%u\t%.4f\t%llu\t%llu\n", or->cname, or->cpgCount, or->cpgTotalScore, or->total_length, or->genome_count);
    }
    carefulClose(&f4);
}

unsigned long long int *samFile2nodupRepbedFile(char *samfile, struct hash *chrHash, struct hash *hashRmsk, struct hash *hashRep, struct hash *hashFam, struct hash *hashCla, int isSam, unsigned int mapQ, int filter, int rmDup, int addChr) {
    samfile_t *samfp;
    char chr[100], prn[500], key[100], strand;
    unsigned int start, end, cend, rstart, rend;
    unsigned long long int *cnt = malloc(sizeof(unsigned long long int) * 5);
    unsigned long long int mapped_reads_num = 0, reads_num = 0, reads_used = 0, unique_reads = 0, repeat_reads = 0;
    struct hash *nochr = newHash(0), *dup = newHash(0);
    if (isSam) {
        if ( (samfp = samopen(samfile, "r", 0)) == 0) {
            fprintf(stderr, "Fail to open SAM file %s\n", samfile);
            errAbort("Error\n");
        }
    } else {
        if ( (samfp = samopen(samfile, "rb", 0)) == 0) {
            fprintf(stderr, "Fail to open BAM file %s\n", samfile);
            errAbort("Error\n");
        }
    }
    strcpy(prn, "empty");
    bam1_t *b;
    bam_header_t *h;
    h = samfp->header;
    b = bam_init1();
    while ( samread(samfp, b) >= 0) {
        //if ( sameString (bam1_qname(b), prn)) 
        //    continue;
        reads_num++;
        if ((reads_num % 10000) == 0)
            fprintf(stderr, "\r* Processed reads: %llu", reads_num);
        //strcpy(prn, bam1_qname(b));
        //if (b->core.tid < 0)
        if (b->core.flag & BAM_FUNMAP)
            continue;
        if (b->core.qual < mapQ)
            continue;
        mapped_reads_num++;
        //change chr name to chr1, chr2 ...
        strcpy(chr, h->target_name[b->core.tid]);
        if (addChr){
            if (startsWith("GL", h->target_name[b->core.tid])) {
                continue;
            } else if (sameWord(h->target_name[b->core.tid], "MT")) {
                strcpy(chr,"chrM");
            } else if (!startsWith("chr", h->target_name[b->core.tid])) {
                strcpy(chr, "chr");
                strcat(chr, h->target_name[b->core.tid]);
            }
        }
        struct hashEl *he = hashLookup(nochr, chr);
        if (he != NULL)
            continue;
        cend = (unsigned int) (hashIntValDefault(chrHash, chr, 2) - 1);
        if (cend == 1){
            hashAddInt(nochr, chr, 1);
            warn("* Warning: reads mapped to chromosome %s will be discarded as %s not existed in the chromosome size file", chr, chr);
            continue;
        }
        reads_used++;
        start = (unsigned int) b->core.pos;
        int tmpend = b->core.n_cigar? bam_calend(&b->core, bam1_cigar(b)) : b->core.pos + b->core.l_qseq;
        end = min(cend, (unsigned int)tmpend);
        strand = (b->core.flag&BAM_FREVERSE)? '-' : '+';
        //remove dup first
        if (rmDup == 1){
            if (sprintf(key, "%s:%u:%u:%c", chr, start, end, strand) < 0)
                errAbort("Mem ERROR");
            struct hashEl *hel = hashLookup(dup, key);
            if (hel == NULL) {
                hashAddInt(dup, key, 1);
            } else {
                continue;
            }
        }
        unique_reads++;
        //transfer coordinates
        int i, j;
        unsigned int qlen = end - start;
        struct binElement *hitList = NULL, *hit;
        struct hashEl *hel2 = hashLookup(hashRmsk, chr);
        if (hel2 != NULL) {
            struct binKeeper *bs2 = (struct binKeeper *) hel2->val;
            hitList = binKeeperFind(bs2, start, end);
            if(hitList != NULL) {
                for (hit = hitList; hit !=NULL; hit = hit->next) {
                    struct rmsk *ss = (struct rmsk *) hit->val;
                    if (filter == 0){
                        struct hashEl *hel3 = hashLookup(hashRep, ss->name);
                        if (hel3 != NULL){
                            struct rep *rs = (struct rep *) hel3->val;
                            rs->read_count++;
                            if (rs->length != 0){
                                rstart = start - ss->start;
                                rstart = (rstart < 0) ? 0 : rstart;
                                rend = rstart + qlen;
                                rend = (rend < ss->end) ? rend : ss->end;
                                for (i = rstart; i < rend; i++) {
                                    j = i + ss->consensus_start;
                                    if (j >= ss->consensus_end) {
                                        break;
                                    }
                                    if (j >= rs->length) {
                                        break;
                                    }
                                    (rs->bp_total)[j]++;
                                }
                            }
                        }
                        //fill hashFam
                        struct hashEl *hel4 = hashLookup(hashFam, ss->fname);
                        if (hel4 != NULL) {
                            struct repfam *fs = (struct repfam *) hel4->val;
                            fs->read_count++;
                        }
                        //fill hashCla
                        struct hashEl *hel5 = hashLookup(hashCla, ss->cname);
                        if (hel5 != NULL) {
                            struct repcla *cs = (struct repcla *) hel5->val;
                            cs->read_count++;
                        }
                    } else {
                        slNameAddHead(&(ss->sl), bam1_qname(b));
                    }
                    break;
                }
                repeat_reads++;
                slFreeList(hitList);
            }
        }
    }
    fprintf(stderr, "\r* Processed reads: %llu\n", reads_num);
    samclose(samfp);
    bam_destroy1(b);
    freeHash(&nochr);
    freeHash(&dup);
    cnt[0] = reads_num;
    cnt[1] = mapped_reads_num;
    cnt[2] = reads_used;
    cnt[3] = unique_reads;
    cnt[4] = repeat_reads;
    return cnt;
}

float getCov(unsigned int aStart, unsigned int aEnd, unsigned int start, unsigned int end){
    float overlap = positiveRangeIntersection((int)aStart, (int)aEnd, (int)start, (int)end);
    float denominator = (float)(aEnd - aStart);
    float cov = (denominator == 0) ? 0.0 : overlap/denominator;
    return cov;
}

unsigned long long int *samFile2nodupRepbedFileNew(char *samfile, struct hash *chrHash, struct hash *hashRmsk, struct hash *hashRep, struct hash *hashFam, struct hash *hashCla, int isSam, unsigned int mapQ, int filter, int rmDup, int addChr, int discardWrongEnd, unsigned int iSize, unsigned int extension, float minCoverage, int treat, char *outbed, char *outbed_unique) {
    samfile_t *samfp;
    FILE *outbed_f = NULL, *outbed_unique_f = NULL;
    char chr[100], key[100], strand;
    unsigned int start, end, cend, rstart, rend;
    unsigned long long int *cnt = malloc(sizeof(unsigned long long int) * 12);
    //unsigned long long int mapped_reads_num = 0, reads_num = 0, reads_used = 0, unique_reads = 0, repeat_reads = 0;
    unsigned long long int read_end1 = 0, read_end2 = 0;
    unsigned long long int read_end1_mapped = 0, read_end2_mapped = 0;
    unsigned long long int read_end1_used = 0, read_end2_used = 0;
    unsigned long long int reads_nonredundant = 0;
    unsigned long long int reads_nonredundant_unique = 0;
    unsigned long long int reads_mapped = 0;
    unsigned long long int reads_mapped_unique = 0;
    unsigned long long int reads_repeat = 0;
    unsigned long long int reads_repeat_unique = 0;
    struct hash *nochr = newHash(0), *dup = newHash(0);
    if (isSam) {
        if ( (samfp = samopen(samfile, "r", 0)) == 0) {
            fprintf(stderr, "Fail to open SAM file %s\n", samfile);
            errAbort("Error\n");
        }
    } else {
        if ( (samfp = samopen(samfile, "rb", 0)) == 0) {
            fprintf(stderr, "Fail to open BAM file %s\n", samfile);
            errAbort("Error\n");
        }
    }
    //strcpy(prn, "empty");
    if (outbed != NULL)
        outbed_f = mustOpen(outbed, "w");
    if (outbed_unique != NULL)
        outbed_unique_f = mustOpen(outbed_unique, "w");
    bam1_t *b;
    bam_header_t *h;
    h = samfp->header;
    b = bam_init1();
    while ( samread(samfp, b) >= 0) {
        //if ( sameString (bam1_qname(b), prn)) 
        //    continue;
        if (b->core.flag & BAM_FPAIRED) {
            if (b->core.flag & BAM_FREAD1){
                read_end1++;
            }else{
                if(treat)
                    read_end1++;
                else
                    read_end2++;
            }
        }else{
            read_end1++;
        }
        if (((read_end1 + read_end2) % 10000) == 0)
            fprintf(stderr, "\r* Processed read ends: %llu", (read_end1 + read_end2));
        //strcpy(prn, bam1_qname(b));
        //if (b->core.tid < 0)
        if (b->core.flag & BAM_FUNMAP)
            continue;
        //if (b->core.qual < mapQ)
        //    continue;
        if (b->core.flag & BAM_FPAIRED) {
            if (b->core.flag & BAM_FREAD1){
                read_end1_mapped++;
            }else{
                if (treat)
                    read_end1_mapped++;
                else
                    read_end2_mapped++;
            }
        }else{
            read_end1_mapped++;
        }
        //change chr name to chr1, chr2 ...
        strcpy(chr, h->target_name[b->core.tid]);
        if (addChr){
            if (startsWith("GL", h->target_name[b->core.tid])) {
                continue;
            } else if (sameWord(h->target_name[b->core.tid], "MT")) {
                strcpy(chr,"chrM");
            } else if (!startsWith("chr", h->target_name[b->core.tid])) {
                strcpy(chr, "chr");
                strcat(chr, h->target_name[b->core.tid]);
            }
        }
        //check Ref reads mapped to existed in chromosome size file or not
        struct hashEl *he = hashLookup(nochr, chr);
        if (he != NULL)
            continue;
        cend = (unsigned int) (hashIntValDefault(chrHash, chr, 2) - 1);
        if (cend == 1){
            hashAddInt(nochr, chr, 1);
            warn("* Warning: read ends mapped to chromosome %s will be discarded as %s not existed in the chromosome size file", chr, chr);
            continue;
        }
        if (b->core.flag & BAM_FPAIRED) {
            if (b->core.flag & BAM_FREAD1){
                read_end1_used++;
            }else{
                if (treat)
                    read_end1_used++;
                else
                    read_end2_used++;
            }
        }else{
            read_end1_used++;
        }
        //get mapping location for paired-end or single-end
        if (treat){
            reads_mapped++;
            if (b->core.qual >= mapQ)
                reads_mapped_unique++;
            start = (unsigned int) b->core.pos;
            int tmpend = b->core.n_cigar? bam_calend(&b->core, bam1_cigar(b)) : b->core.pos + b->core.l_qseq;
            end = min(cend, (unsigned int)tmpend);
            strand = (b->core.flag&BAM_FREVERSE)? '-' : '+';
            if (extension) {
                if (strand == '+'){
                    end = min(start + extension, cend);
                }else{
                    if (end < extension)
                        start = 0;
                    else
                        start = end - extension;
                    //start = max(end - extension, 0);
                }
            }

        }else{
        if (b->core.flag & BAM_FPAIRED) {
            if (!(b->core.flag & BAM_FMUNMAP)){
                if (b->core.flag & BAM_FREAD1){
                    if (abs(b->core.isize) > iSize || b->core.isize == 0){
                        continue;
                    }else{
                        reads_mapped++;
                        if (b->core.qual >= mapQ)
                            reads_mapped_unique++;
                        if (b->core.isize > 0){
                            start = (unsigned int) b->core.pos;
                            strand = '+';
                            int tmpend = start + b->core.isize;
                            end = min(cend, (unsigned int)tmpend);
                        }else{
                            start = (unsigned int) b->core.pos;
                            strand = '-';
                            int tmpend = start - b->core.isize;
                            end = min(cend, (unsigned int)tmpend);
                        }
                
                    }
                }else{
                    continue;
                }
            }else{
                if (discardWrongEnd){
                    continue;
                }else{
                    reads_mapped++;
                    if (b->core.qual >= mapQ)
                        reads_mapped_unique++;
                    start = (unsigned int) b->core.pos;
                    int tmpend = b->core.n_cigar? bam_calend(&b->core, bam1_cigar(b)) : b->core.pos + b->core.l_qseq;
                    end = min(cend, (unsigned int)tmpend);
                    strand = (b->core.flag&BAM_FREVERSE)? '-' : '+';
                    if (extension) {
                        if (strand == '+'){
                            end = min(start + extension, cend);
                        }else{
                            if (end < extension)
                                start = 0;
                            else
                                start = end - extension;
                            //start = max(end - extension, 0);
                        }
                    }
                }
            }
        }else{
            reads_mapped++;
            if (b->core.qual >= mapQ)
                reads_mapped_unique++;
            start = (unsigned int) b->core.pos;
            int tmpend = b->core.n_cigar? bam_calend(&b->core, bam1_cigar(b)) : b->core.pos + b->core.l_qseq;
            end = min(cend, (unsigned int)tmpend);
            strand = (b->core.flag&BAM_FREVERSE)? '-' : '+';
            if (extension) {
                if (strand == '+'){
                    end = min(start + extension, cend);
                }else{
                    if (end < extension)
                        start = 0;
                    else
                        start = end - extension;
                    //start = max(end - extension, 0);
                }
            }
        }
    }
        //remove dup first
        if (rmDup){
            if (sprintf(key, "%s:%u:%u:%c", chr, start, end, strand) < 0)
                errAbort("Mem ERROR");
            struct hashEl *hel = hashLookup(dup, key);
            if (hel == NULL) {
                hashAddInt(dup, key, 1);
            } else {
                continue;
            }
        }
        reads_nonredundant++;
        if (b->core.qual >= mapQ)
            reads_nonredundant_unique++;
        //output bed
        if (outbed_f)
            fprintf(outbed_f, "%s\t%u\t%u\t%s\t%i\t%c\n", chr, start, end, bam1_qname(b), b->core.qual, strand);
        if (outbed_unique_f){
            if(b->core.qual >= mapQ)
                fprintf(outbed_unique_f, "%s\t%u\t%u\t%s\t%i\t%c\n", chr, start, end, bam1_qname(b), b->core.qual, strand);
        }

        //transfer coordinates
        int i, j;
        int index = 0, tindex = 0;
        float coverage = 0.0, tcoverage = 0.0;
        unsigned int qlen = end - start;
        struct rmsk *ss = NULL;
        struct binElement *hitList = NULL, *hit;
        struct hashEl *hel2 = hashLookup(hashRmsk, chr);
        if (hel2 != NULL) {
            struct binKeeper *bs2 = (struct binKeeper *) hel2->val;
            hitList = binKeeperFind(bs2, start, end);
            if(hitList != NULL) {
                for (hit = hitList; hit !=NULL; hit = hit->next) {
                    index++;
                    struct rmsk *sss = (struct rmsk *) hit->val;
                    float cov = getCov(start, end, sss->start, sss->end);
                    //fprintf(stderr, "coverage: %.2f\n", cov);
                    if (cov > coverage){
                        tindex = index;
                        tcoverage = cov;
                    }
                    coverage = cov;
                }
                if (tcoverage < minCoverage)
                    continue;
                index = 0;
                for (hit = hitList; hit !=NULL; hit = hit->next) {
                    index++;
                    if (index == tindex){
                        ss = (struct rmsk *) hit->val;
                        break;
                    }
                }
                if (filter == 0){
                    struct hashEl *hel3 = hashLookup(hashRep, ss->name);
                    if (hel3 != NULL){
                        struct rep *rs = (struct rep *) hel3->val;
                        rs->read_count++;
                        if (b->core.qual >= mapQ)
                            rs->read_count_unique++;
                        if (rs->length != 0){
                            rstart = start - ss->start;
                            rstart = (rstart < 0) ? 0 : rstart;
                            rend = rstart + qlen;
                            rend = (rend < ss->end) ? rend : ss->end;
                            for (i = rstart; i < rend; i++) {
                                j = i + ss->consensus_start;
                                if (j >= ss->consensus_end) {
                                    break;
                                }
                                if (j >= rs->length) {
                                    break;
                                }
                                (rs->bp_total)[j]++;
                                if (b->core.qual >= mapQ)
                                    (rs->bp_total_unique)[j]++;
                            }
                        }
                    }
                    //fill hashFam
                    struct hashEl *hel4 = hashLookup(hashFam, ss->fname);
                    if (hel4 != NULL) {
                        struct repfam *fs = (struct repfam *) hel4->val;
                        fs->read_count++;
                        if (b->core.qual >= mapQ)
                            fs->read_count_unique++;
                    }
                    //fill hashCla
                    struct hashEl *hel5 = hashLookup(hashCla, ss->cname);
                    if (hel5 != NULL) {
                        struct repcla *cs = (struct repcla *) hel5->val;
                        cs->read_count++;
                        if (b->core.qual >= mapQ)
                            cs->read_count_unique++;
                    }
                } else {
                    slNameAddHead(&(ss->sl), bam1_qname(b));
                    if (b->core.qual >= mapQ)
                        slNameAddHead(&(ss->sl_unique), bam1_qname(b));
                }
                reads_repeat++;
                if (b->core.qual >= mapQ)
                    reads_repeat_unique++;
                slFreeList(hitList);
            }
        }
    }
    fprintf(stderr, "\r* Processed read ends: %llu\n", (read_end1 + read_end2));
    samclose(samfp);
    bam_destroy1(b);
    freeHash(&nochr);
    freeHash(&dup);
    if (outbed_f)
        carefulClose(&outbed_f);
    if (outbed_unique_f)
        carefulClose(&outbed_unique_f);
    cnt[0] = read_end1;
    cnt[1] = read_end2;
    cnt[2] = read_end1_mapped;
    cnt[3] = read_end2_mapped;
    cnt[4] = read_end1_used;
    cnt[5] = read_end2_used;
    cnt[6] = reads_mapped;
    cnt[7] = reads_mapped_unique;
    cnt[8] = reads_nonredundant;
    cnt[9] = reads_repeat;
    cnt[10] = reads_repeat_unique;
    cnt[11] = reads_nonredundant_unique;
    return cnt;
}

void cpgBedGraphOverlapRepeat(char *cpgBedGraphFile, struct hash *hashRmsk, struct hash *hashRep, struct hash *hashFam, struct hash *hashCla, int filter) {
    struct lineFile *infileStream = lineFileOpen(cpgBedGraphFile, TRUE);
    char *row[20], *line;
    unsigned int start, end, rstart, rend, cpgInRepeat = 0, cpglines = 0;
    double score = 0;
    while ( lineFileNextReal(infileStream, &line)) {
        int numFields = chopByWhite(line, row, ArraySize(row));
        if (numFields < 4)
            errAbort("file %s doesn't appear to be in bedGraph format. At least 4 fields required, got %d", cpgBedGraphFile, numFields);
        cpglines++;
        start = (unsigned int)strtol(row[1], NULL, 0);
        end = (unsigned int)strtol(row[2], NULL, 0);
        score = strtod(row[3], NULL);
        //transfer coordinates
        int i, j;
        struct rmsk *ss = NULL;
        struct binElement *hitList = NULL, *hit;
        struct hashEl *hel2 = hashLookup(hashRmsk, row[0]);
        if (hel2 != NULL) {
            struct binKeeper *bs2 = (struct binKeeper *) hel2->val;
            hitList = binKeeperFind(bs2, start, end);
            if(hitList != NULL) {
                for (hit = hitList; hit !=NULL; hit = hit->next) {
                    ss = (struct rmsk *) hit->val;
                    break;
                }
                if (filter){
                    ss->cpgCount++;
                    ss->cpgTotalScore += score;
                }else{
                    struct hashEl *hel3 = hashLookup(hashRep, ss->name);
                    if (hel3 != NULL){
                        struct rep *rs = (struct rep *) hel3->val;
                        rs->cpgCount++;
                        rs->cpgTotalScore += score;
                        if (rs->length != 0){
                            rstart = start - ss->start;
                            rstart = (rstart < 0) ? 0 : rstart;
                            rend = rstart + 2; //CpG site
                            rend = (rend < ss->end) ? rend : ss->end;
                            for (i = rstart; i < rend; i++) {
                                j = i + ss->consensus_start;
                                if (j >= ss->consensus_end) {
                                    break;
                                }
                                if (j >= rs->length) {
                                    break;
                                }
                                (rs->cpgScore)[j] += score;
                            }
                        }
                    }
                    //fill hashFam
                    struct hashEl *hel4 = hashLookup(hashFam, ss->fname);
                    if (hel4 != NULL) {
                        struct repfam *fs = (struct repfam *) hel4->val;
                        fs->cpgCount++;
                        fs->cpgTotalScore += score;
                    }
                    //fill hashCla
                    struct hashEl *hel5 = hashLookup(hashCla, ss->cname);
                    if (hel5 != NULL) {
                        struct repcla *cs = (struct repcla *) hel5->val;
                        cs->cpgCount++;
                        cs->cpgTotalScore += score;
                    }
                }
                slFreeList(hitList);
                cpgInRepeat++;
            }
        }
    }
    fprintf(stderr, "* Processed CpG sites: %u\n", cpglines);
    fprintf(stderr, "* CpG sites in Repeats: %u\n", cpgInRepeat);
    lineFileClose(&infileStream);
}

unsigned long long int *sam2bed(char *samfile, char *outbed, struct hash *chrHash, int isSam, unsigned int mapQ, int rmDup, int addChr, int discardWrongEnd, unsigned int iSize, unsigned int extension, int treat) {
    samfile_t *samfp;
    FILE *outbed_f = mustOpen(outbed, "w");
    char chr[100], key[100], strand;
    unsigned int start, end, cend;
    unsigned long long int *cnt = malloc(sizeof(unsigned long long int) * 10);
    unsigned long long int read_end1 = 0, read_end2 = 0;
    unsigned long long int read_end1_mapped = 0, read_end2_mapped = 0;
    unsigned long long int read_end1_used = 0, read_end2_used = 0;
    unsigned long long int reads_nonredundant = 0;
    unsigned long long int reads_nonredundant_unique = 0;
    unsigned long long int reads_mapped = 0;
    unsigned long long int reads_mapped_unique = 0;
    struct hash *nochr = newHash(0), *dup = newHash(0);
    if (isSam) {
        if ( (samfp = samopen(samfile, "r", 0)) == 0) {
            fprintf(stderr, "Fail to open SAM file %s\n", samfile);
            errAbort("Error\n");
        }
    } else {
        if ( (samfp = samopen(samfile, "rb", 0)) == 0) {
            fprintf(stderr, "Fail to open BAM file %s\n", samfile);
            errAbort("Error\n");
        }
    }
    bam1_t *b;
    bam_header_t *h;
    h = samfp->header;
    b = bam_init1();
    int8_t *buf;
    int max_buf;
    buf = 0;
    max_buf = 0;
    uint8_t *seq;
    while ( samread(samfp, b) >= 0) {
        if (b->core.flag & BAM_FPAIRED) {
            if (b->core.flag & BAM_FREAD1){
                read_end1++;
            }else{
                if(treat)
                    read_end1++;
                else
                    read_end2++;
            }
        }else{
            read_end1++;
        }
        if (((read_end1 + read_end2) % 10000) == 0)
            fprintf(stderr, "\r* Processed read ends: %llu", (read_end1 + read_end2));
        if (b->core.flag & BAM_FUNMAP)
            continue;
        if (b->core.flag & BAM_FPAIRED) {
            if (b->core.flag & BAM_FREAD1){
                read_end1_mapped++;
            }else{
                if (treat)
                    read_end1_mapped++;
                else
                    read_end2_mapped++;
            }
        }else{
            read_end1_mapped++;
        }
        //change chr name to chr1, chr2 ...
        strcpy(chr, h->target_name[b->core.tid]);
        if (addChr){
            if (startsWith("GL", h->target_name[b->core.tid])) {
                continue;
            } else if (sameWord(h->target_name[b->core.tid], "MT")) {
                strcpy(chr,"chrM");
            } else if (!startsWith("chr", h->target_name[b->core.tid])) {
                strcpy(chr, "chr");
                strcat(chr, h->target_name[b->core.tid]);
            }
        }
        //check Ref reads mapped to existed in chromosome size file or not
        struct hashEl *he = hashLookup(nochr, chr);
        if (he != NULL)
            continue;
        cend = (unsigned int) (hashIntValDefault(chrHash, chr, 2) - 1);
        if (cend == 1){
            hashAddInt(nochr, chr, 1);
            warn("* Warning: read ends mapped to chromosome %s will be discarded as %s not existed in the chromosome size file", chr, chr);
            continue;
        }
        if (b->core.flag & BAM_FPAIRED) {
            if (b->core.flag & BAM_FREAD1){
                read_end1_used++;
            }else{
                if (treat)
                    read_end1_used++;
                else
                    read_end2_used++;
            }
        }else{
            read_end1_used++;
        }
        //get mapping location for paired-end or single-end
        if (treat){
            reads_mapped++;
            if (b->core.qual >= mapQ)
                reads_mapped_unique++;
            start = (unsigned int) b->core.pos;
            int tmpend = b->core.n_cigar? bam_calend(&b->core, bam1_cigar(b)) : b->core.pos + b->core.l_qseq;
            end = min(cend, (unsigned int)tmpend);
            strand = (b->core.flag&BAM_FREVERSE)? '-' : '+';
            if (extension) {
                if (strand == '+'){
                    end = min(start + extension, cend);
                }else{
                    if (end < extension)
                        start = 0;
                    else
                        start = end - extension;
                    //start = max(end - extension, 0);
                }
            }

        }else{
        if (b->core.flag & BAM_FPAIRED) {
            if (!(b->core.flag & BAM_FMUNMAP)){
                if (b->core.flag & BAM_FREAD1){
                    if (abs(b->core.isize) > iSize || b->core.isize == 0){
                        continue;
                    }else{
                        reads_mapped++;
                        if (b->core.qual >= mapQ)
                            reads_mapped_unique++;
                        if (b->core.isize > 0){
                            start = (unsigned int) b->core.pos;
                            strand = '+';
                            int tmpend = start + b->core.isize;
                            end = min(cend, (unsigned int)tmpend);
                        }else{
                            start = (unsigned int) b->core.pos;
                            strand = '-';
                            int tmpend = start - b->core.isize;
                            end = min(cend, (unsigned int)tmpend);
                        }
                
                    }
                }else{
                    continue;
                }
            }else{
                if (discardWrongEnd){
                    continue;
                }else{
                    reads_mapped++;
                    if (b->core.qual >= mapQ)
                        reads_mapped_unique++;
                    start = (unsigned int) b->core.pos;
                    int tmpend = b->core.n_cigar? bam_calend(&b->core, bam1_cigar(b)) : b->core.pos + b->core.l_qseq;
                    end = min(cend, (unsigned int)tmpend);
                    strand = (b->core.flag&BAM_FREVERSE)? '-' : '+';
                    if (extension) {
                        if (strand == '+'){
                            end = min(start + extension, cend);
                        }else{
                            if (end < extension)
                                start = 0;
                            else
                                start = end - extension;
                            //start = max(end - extension, 0);
                        }
                    }
                }
            }
        }else{
            reads_mapped++;
            if (b->core.qual >= mapQ)
                reads_mapped_unique++;
            start = (unsigned int) b->core.pos;
            int tmpend = b->core.n_cigar? bam_calend(&b->core, bam1_cigar(b)) : b->core.pos + b->core.l_qseq;
            end = min(cend, (unsigned int)tmpend);
            strand = (b->core.flag&BAM_FREVERSE)? '-' : '+';
            if (extension) {
                if (strand == '+'){
                    end = min(start + extension, cend);
                }else{
                    if (end < extension)
                        start = 0;
                    else
                        start = end - extension;
                }
            }
        }
    }
        //remove dup or not
        if (rmDup){
            if (sprintf(key, "%s:%u:%u:%c", chr, start, end, strand) < 0)
                errAbort("Mem ERROR");
            struct hashEl *hel = hashLookup(dup, key);
            if (hel == NULL) {
                hashAddInt(dup, key, 1);
            } else {
                continue;
            }
        }
        reads_nonredundant++;
        if (b->core.qual >= mapQ)
            reads_nonredundant_unique++;
        //output bed
        int i, qlen = b->core.l_qseq;
        if(b->core.qual >= mapQ){
            fprintf(outbed_f, "%s\t%u\t%u\t", chr, start, end);
            //print read sequence
            if (max_buf < qlen + 1 ) {
                max_buf = qlen + 1;
                kroundup32(max_buf);
                buf = realloc(buf, max_buf);
            }
            buf[qlen] = 0;
            seq = bam1_seq(b);
            for (i = 0; i < qlen; ++i)
                buf[i] = bam1_seqi(seq, i);
            if (b->core.flag & 16) {
                for (i = 0; i < qlen>>1; ++i){
                    int8_t t = seq_comp_table[buf[qlen - 1 - i]];
                    buf[qlen - 1 - i] = seq_comp_table[buf[i]];
                    buf[i] = t;
                }
                if (qlen&1) buf[i] = seq_comp_table[buf[i]];
            }
            for (i = 0; i < qlen; ++i)
                buf[i] = bam_nt16_rev_table[buf[i]];
            fprintf(outbed_f, "%s", (char*)buf);

            fprintf(outbed_f, "\t%i\t%c\n", b->core.qual, strand);
        }
    }
    fprintf(stderr, "\r* Processed read ends: %llu\n", (read_end1 + read_end2));
    samclose(samfp);
    free(buf);
    bam_destroy1(b);
    freeHash(&nochr);
    freeHash(&dup);
    carefulClose(&outbed_f);
    cnt[0] = read_end1;
    cnt[1] = read_end2;
    cnt[2] = read_end1_mapped;
    cnt[3] = read_end2_mapped;
    cnt[4] = read_end1_used;
    cnt[5] = read_end2_used;
    cnt[6] = reads_mapped;
    cnt[7] = reads_mapped_unique;
    cnt[8] = reads_nonredundant;
    cnt[9] = reads_nonredundant_unique;
    return cnt;
}

unsigned long long int *PEsamFile2nodupRepbedFile(char *samfile, struct hash *chrHash, struct hash *hashRmsk, struct hash *hashRep, struct hash *hashFam, struct hash *hashCla, int isSam, unsigned int mapQ, int filter, int rmDup, int addChr, unsigned int iSize) {
    samfile_t *samfp;
    char chr[100], prn[500], key[100], strand;
    unsigned int start, end, cend, rstart, rend;
    unsigned long long int *cnt = malloc(sizeof(unsigned long long int) * 5);
    unsigned long long int mapped_reads_num = 0, reads_num = 0, reads_used = 0, unique_reads = 0, repeat_reads = 0;
    struct hash *nochr = newHash(0), *dup = newHash(0);
    if (isSam) {
        if ( (samfp = samopen(samfile, "r", 0)) == 0) {
            fprintf(stderr, "Fail to open SAM file %s\n", samfile);
            errAbort("Error\n");
        }
    } else {
        if ( (samfp = samopen(samfile, "rb", 0)) == 0) {
            fprintf(stderr, "Fail to open BAM file %s\n", samfile);
            errAbort("Error\n");
        }
    }
    strcpy(prn, "empty");
    bam1_t *b[2];
    int curr, has_prev;
    bam_header_t *h;
    h = samfp->header;
    b[0] = bam_init1();
    b[1] = bam_init1();
    curr = 0;
    has_prev = 0;
    while ( samread(samfp, b[curr]) >= 0) {
        
	bam1_t *cur = b[curr], *pre = b[1-curr];
	if (has_prev) {
	    if (strcmp(bam1_qname(cur), bam1_qname(pre)) == 0) { // identical pair name
		has_prev = 0;
	    }
	} else has_prev = 1;
	curr = 1 - curr;
        
        if (!has_prev) {
            if ( sameString (bam1_qname(cur), prn)) 
                continue;
        }
        reads_num++;
        if ((reads_num % 10000) == 0)
            fprintf(stderr, "\r* Processed reads pair: %llu", reads_num);
        if (!has_prev)
            strcpy(prn, bam1_qname(cur));
        if (has_prev){
            //if (cur->core.tid < 0 || pre->core.tid < 0)
            if (cur->core.flag & BAM_FUNMAP || pre->core.flag & BAM_FUNMAP)
                continue;
        }
        if (has_prev){
            if (cur->core.qual < mapQ || pre->core.qual < mapQ)
                continue;
        }
        mapped_reads_num++;
        if (has_prev) {
            if (cur->core.tid != pre->core.tid)
                continue;
        }
        if (has_prev) {
            if (abs(pre->core.isize) > iSize || pre->core.isize == 0)
            //if (abs(pre->core.isize) > iSize)
                continue;
            strcpy(chr, h->target_name[cur->core.tid]);
            //change chr name to chr1, chr2 ...
            if (addChr){
                if (startsWith("GL", h->target_name[cur->core.tid])) {
                    continue;
                } else if (sameWord(h->target_name[cur->core.tid], "MT")) {
                    strcpy(chr,"chrM");
                } else if (!startsWith("chr", h->target_name[cur->core.tid])) {
                    strcpy(chr, "chr");
                    strcat(chr, h->target_name[cur->core.tid]);
                }
            }
            struct hashEl *he = hashLookup(nochr, chr);
            if (he != NULL)
                continue;
            cend = (unsigned int) (hashIntValDefault(chrHash, chr, 2) - 1);
            if (cend == 1){
                hashAddInt(nochr, chr, 1);
                warn("* Warning: reads mapped to chromosome %s will be discarded as %s not existed in the chromosome size file", chr, chr);
                continue;
            }
            reads_used++;
            int tmpend;
            if (pre->core.isize > 0){
                start = (unsigned int) pre->core.pos;
                strand = '+';
                tmpend = start + pre->core.isize;
            }else{
                start = (unsigned int) cur->core.pos;
                strand = '-';
                tmpend = start - pre->core.isize;
            }
            end = min(cend, (unsigned int)tmpend);
            //remove dup first
            if (rmDup == 1){
                if (sprintf(key, "%s:%u:%u:%c", chr, start, end, strand) < 0)
                    errAbort("Mem ERROR");
                struct hashEl *hel = hashLookup(dup, key);
                if (hel == NULL) {
                    hashAddInt(dup, key, 1);
                } else {
                    continue;
                }
            }
            unique_reads++;
            //transfer coordinates
            int i, j;
            unsigned int qlen = end - start;
            struct binElement *hitList = NULL, *hit;
            struct hashEl *hel2 = hashLookup(hashRmsk, chr);
            if (hel2 != NULL) {
                struct binKeeper *bs2 = (struct binKeeper *) hel2->val;
                hitList = binKeeperFind(bs2, start, end);
                if(hitList != NULL) {
                    for (hit = hitList; hit !=NULL; hit = hit->next) {
                        struct rmsk *ss = (struct rmsk *) hit->val;
                        if (filter == 0){
                            struct hashEl *hel3 = hashLookup(hashRep, ss->name);
                            if (hel3 != NULL){
                                struct rep *rs = (struct rep *) hel3->val;
                                rs->read_count++;
                                if (rs->length != 0){
                                    rstart = start - ss->start;
                                    rstart = (rstart < 0) ? 0 : rstart;
                                    rend = rstart + qlen;
                                    rend = (rend < ss->end) ? rend : ss->end;
                                    for (i = rstart; i < rend; i++) {
                                        j = i + ss->consensus_start;
                                        if (j >= ss->consensus_end) {
                                            break;
                                        }
                                        if (j >= rs->length) {
                                            break;
                                        }
                                        (rs->bp_total)[j]++;
                                    }
                                }
                            }
                            //fill hashFam
                            struct hashEl *hel4 = hashLookup(hashFam, ss->fname);
                            if (hel4 != NULL) {
                                struct repfam *fs = (struct repfam *) hel4->val;
                                fs->read_count++;
                            }
                            //fill hashCla
                            struct hashEl *hel5 = hashLookup(hashCla, ss->cname);
                            if (hel5 != NULL) {
                                struct repcla *cs = (struct repcla *) hel5->val;
                                cs->read_count++;
                            }
                        } else {
                            slNameAddHead(&(ss->sl), bam1_qname(cur));
                        }
                        break;
                    }
                    repeat_reads++;
                    slFreeList(hitList);
                }
            }
        }
    }
    fprintf(stderr, "\r* Processed reads pair: %llu\n", reads_num);
    samclose(samfp);
    bam_destroy1(b[0]);
    bam_destroy1(b[1]);
    freeHash(&nochr);
    freeHash(&dup);
    cnt[0] = reads_num;
    cnt[1] = mapped_reads_num;
    cnt[2] = reads_used;
    cnt[3] = unique_reads;
    cnt[4] = repeat_reads;
    return cnt;
}

void freermsk(struct rmsk *s){
    freeMem(s->name);
    freeMem(s->cname);
    freeMem(s->fname);
    free(s);
}

void rmsk2binKeeperHash(char *rmskfile, struct hash *chrHash, struct hash *repHash, struct hash **hashRmsk, struct hash **hashRep, struct hash **hashFam, struct hash **hashCla, int filterField, char *filterName){
    char strand;
    char *row[17];
    struct lineFile *repeat_stream = lineFileOpen(rmskfile, TRUE);
    int repeat_num = 0;
    struct hash *hash1 = newHash(0); //hashRmsk
    struct hash *hash2 = newHash(0); //hashRep
    struct hash *hash3 = newHash(0); //hashFam
    struct hash *hash4 = newHash(0); //hashCla
    while( lineFileNextRow(repeat_stream, row, ArraySize(row))){
        if (filterField != 0){
            if (strcmp(filterName, row[filterField]) != 0)
                continue;
        }
        struct rmsk *s = malloc(sizeof(struct rmsk));
        repeat_num++;
        s->chr = cloneString(row[5]);
        strand = row[9][0];
        if ( strand == '+' ) {
            s->consensus_start = (unsigned int)strtol(row[13], NULL, 0);    
        } else {
            s->consensus_start = (unsigned int)strtol(row[15], NULL, 0);    
        }
        s->consensus_end = (unsigned int)strtol(row[14], NULL, 0);
        s->start = (unsigned int)strtol(row[6], NULL, 0);
        s->end = (unsigned int)strtol(row[7], NULL, 0);
        s->name = cloneString(row[10]);
        s->cname = cloneString(row[11]);
        s->fname = cloneString(row[12]);
        s->length = s->end - s->start;
        s->sl = NULL;
        s->sl_unique = NULL;
        s->cpgCount = 0;
        s->cpgTotalScore = 0;
        
        struct hashEl *hel = hashLookup(hash1, s->chr);
        if (hel != NULL) {
            struct binKeeper *bk = (struct binKeeper *) hel->val;
            binKeeperAdd(bk, s->start, s->end, s);
        } else {
            int size = hashIntValDefault(chrHash, s->chr, 0);
            if (size == 0) {
                freermsk(s);
                continue;
            }
            struct binKeeper *bk = binKeeperNew(0, size);
            binKeeperAdd(bk, s->start, s->end, s);
            hashAdd(hash1, s->chr, bk);
        }

        if (filterField != 0)
            continue;
        
        struct hashEl *hel2 = hashLookup(hash2, s->name);
        if (hel2 != NULL) {
            struct rep *rr = (struct rep *) hel2->val;
            rr->genome_count++;
            rr->total_length += s->length;
        } else {
            struct rep *ns = malloc(sizeof(struct rep));
            ns->name = cloneString(s->name);
            ns->cname = cloneString(s->cname);
            ns->fname = cloneString(s->fname);
            ns->genome_count = 1;
            ns->total_length = s->length;
            ns->read_count = 0;
            ns->read_count_unique = 0;
            ns->cpgCount = 0;
            ns->cpgTotalScore = 0;
            ns->length = hashIntValDefault(repHash, ns->name, 0);
            ns->bp_total = malloc(sizeof(unsigned int) * ns->length);
            ns->bp_total_unique = malloc(sizeof(unsigned int) * ns->length);
            ns->cpgScore = malloc(sizeof(double) * ns->length);
            int k;
            for (k = 0; k < ns->length; k++){ 
                (ns->bp_total)[k] = 0;
                (ns->bp_total_unique)[k] = 0;
                (ns->cpgScore)[k] = 0;
            }
            hashAdd(hash2, ns->name, ns);
        }
        //fill hash3
        struct hashEl *hel3 = hashLookup(hash3, s->fname);
        if (hel3 != NULL) {
            struct repfam *rr = (struct repfam *) hel3->val;
            rr->genome_count++;
            rr->total_length += s->length;
        } else {
            struct repfam *ns = malloc(sizeof(struct repfam));
            ns->fname = cloneString(s->fname);
            ns->cname = cloneString(s->cname);
            ns->genome_count = 1;
            ns->total_length = s->length;
            ns->read_count = 0;
            ns->read_count_unique = 0;
            ns->cpgCount = 0;
            ns->cpgTotalScore = 0;
            hashAdd(hash3, ns->fname, ns);
        }
        //fill hash4
        struct hashEl *hel4 = hashLookup(hash4, s->cname);
        if (hel4 != NULL) {
            struct repcla *rr = (struct repcla *) hel4->val;
            rr->genome_count++;
            rr->total_length += s->length;
        } else {
            struct repcla *ns = malloc(sizeof(struct repcla));
            ns->cname = cloneString(s->cname);
            ns->genome_count = 1;
            ns->total_length = s->length;
            ns->read_count = 0;
            ns->read_count_unique = 0;
            ns->cpgCount = 0;
            ns->cpgTotalScore = 0;
            hashAdd(hash4, ns->cname, ns);
        }
    }
    lineFileClose(&repeat_stream);
    if (filterField == 0){
        fprintf(stderr, "* Total %d repeats found.\n", repeat_num);
    } else{
        if ( repeat_num <= 0 )
            errAbort("* No repeats found related to [%s], typo? or specify wrong repName/Class/Family filter?", filterName);
        fprintf(stderr, "* Total %d repeats for [%s].\n", repeat_num, filterName);
    }
    *hashRmsk = hash1;
    *hashRep = hash2;
    *hashFam = hash3;
    *hashCla = hash4;
}

struct hash *MREfrag2Hash (char *fragfile, int minlen, int maxlen){
    struct hash *hash = newHash(0);
    char *row[20], *line;
    int start, end, length, rstart, rend;
    char key1[100], key2[100];
    struct lineFile *frag_stream = lineFileOpen(fragfile, TRUE);
    while ( lineFileNextReal(frag_stream, &line)){
        int numFields = chopByWhite(line, row, ArraySize(row));
        if (numFields < 4)
            errAbort("file %s doesn't appear to be in bed format. At least 3 fields required, got %d", fragfile, numFields);
        start = (int) strtol(row[1], NULL, 0);
        end = (int) strtol(row[2], NULL, 0);
        length  = end - start;
        if (length >= minlen && length <= maxlen){
            struct mreFrag *mre1 = malloc(sizeof(struct mreFrag));
            struct mreFrag *mre2 = malloc(sizeof(struct mreFrag));
            if (sameWord(row[3], "CGCG")) { //special case for CGCG
                rstart = start + 2;
                rend = end - 2;
                if (sprintf(key1, "%s:%i:%c", row[0], rstart, '+') < 0)
                    errAbort("Mem ERROR");
                if (sprintf(key2, "%s:%i:%c", row[0], rend, '-') < 0)
                    errAbort("Mem ERROR");
            } else { // case for CCGG, CCGC, GCGC, ACGT
                rstart = start + 1;
                rend = end - 1;
                if (sprintf(key1, "%s:%i:%c", row[0], rstart, '+') < 0)
                    errAbort("Mem ERROR");
                if (sprintf(key2, "%s:%i:%c", row[0], rend, '-') < 0)
                    errAbort("Mem ERROR");
            }
            mre1->head = 1; mre1->start = rstart; mre1->end = rend;
            mre2->head= 0; mre2->start = rstart; mre2->end = rend;
            mre1->reads_count = 0; mre2->reads_count = 0;
            strcpy(mre1->pair, key2); strcpy(mre2->pair, key1);
            strcpy(mre1->site, row[3]);
            strcpy(mre2->site, row[3]);
            strcpy(mre1->chr, row[0]);
            strcpy(mre2->chr, row[0]);
            hashAdd(hash, key1, mre1);
            hashAdd(hash, key2, mre2);
        } else {
            continue;
        }
    }
    lineFileClose(&frag_stream);
    return hash;
}

unsigned long long int *filterReadByMREsite(struct hash *hash, char *inBed, char *outBed, int call){
    FILE *f  = mustOpen(outBed, "w");
    char strand, *row[20], *line;
    struct lineFile *inBedStream = lineFileOpen(inBed, TRUE);
    char key1[100]; //key2[100];
    int start, end, rstart, rend;
    unsigned long long int CCGG = 0, CCGC = 0, GCGC = 0, ACGT = 0, CGCG = 0, unknown = 0;
    unsigned long long int *cnt = malloc(sizeof(unsigned long long int) * 6);
    while( lineFileNextReal(inBedStream, &line)){
        int numFields = chopByWhite(line, row, ArraySize(row));
        if (numFields < 6)
            errAbort("file %s doesn't appear to be in bed format. At least 6 fields required, got %d", inBed, numFields);
        strand = row[5][0];
        start = (int)strtol(row[1], NULL, 0);
        end = (int)strtol(row[2], NULL, 0);
        if (strand == '+'){
            rstart = start - call;
            rend = end;
            if (sprintf(key1, "%s:%i:%c", row[0], rstart, '+') < 0)
                errAbort("Mem ERROR");
            //if (sprintf(key2, "%s:%i:%c", bed->chrom, bed->chromStart - 3, '+') < 0)
            //    errAbort("Mem ERROR");
        } else {
            rstart = start;
            rend = end + call;
            if (sprintf(key1, "%s:%i:%c", row[0], rend, '-') < 0)
                errAbort("Mem ERROR");
            //if (sprintf(key2, "%s:%i:%c", bed->chrom, bed->chromEnd + 3, '-') < 0)
            //    errAbort("Mem ERROR");
        }
        struct hashEl *hel = hashLookup(hash, key1);
        //struct hashEl *hel2 = hashLookup(hash, key2);
        if (hel != NULL) {
            //bedOutputN(bed, 6, f, '\t', '\n');
            fprintf(f, "%s\t%i\t%i\t%s\t%s\t%c\n", row[0], rstart, rend, row[3], row[4], strand);
            struct mreFrag *mre = (struct mreFrag *) hel->val;
            mre->reads_count++;
            if (sameWord(mre->site, "CCGG")) {
                CCGG++;
            } else if (sameWord(mre->site, "CCGC") || sameWord(mre->site, "GCGG")) {
                CCGC++;
            } else if (sameWord(mre->site, "GCGC")){
                GCGC++;
            } else if (sameWord(mre->site, "ACGT")){
                ACGT++;
            } else if (sameWord(mre->site, "CGCG")){
                CGCG++;
            }
        /*} else if(hel2 != NULL){
            bedOutputN(bed, 6, f, '\t', '\n');
            struct mreFrag *mre = (struct mreFrag *) hel2->val;
            mre->reads_count++;
            if (sameWord(mre->site, "CCGG")) {
                CCGG++;
            } else if (sameWord(mre->site, "CCGC") || sameWord(mre->site, "GCGG")) {
                CCGC++;
            } else if (sameWord(mre->site, "GCGC")){
                GCGC++;
            } else if (sameWord(mre->site, "ACGT")){
                ACGT++;
            } else if (sameWord(mre->site, "CGCG")){
                CGCG++;
            }*/
        } else {
            unknown++;
        }
    }
    carefulClose(&f);
    lineFileClose(&inBedStream);
    //bedFreeList(&bedList);
    cnt[0] = CCGG;
    cnt[1] = CCGC;
    cnt[2] = GCGC;
    cnt[3] = ACGT;
    cnt[4] = CGCG;
    cnt[5] = unknown;
    return cnt;
}

double calCpGscore (struct mreFrag *mre, unsigned long long int *cnt){
    unsigned long long int norm = 0;
    if (sameWord(mre->site, "CCGG")) {
        norm = cnt[0];
    //} else if (sameWord(mre->site, "CCGC")) {
    } else if (sameWord(mre->site, "CCGC") || sameWord(mre->site, "GCGG")) {
        norm = cnt[1];
    } else if (sameWord(mre->site, "GCGC")){
        norm = cnt[2];
    } else if (sameWord(mre->site, "ACGT")){
        norm = cnt[3];
    } else if (sameWord(mre->site, "CGCG")){
        norm = cnt[4];
    }
    if (norm > 0)
        return (double) mre->reads_count * 1000000.0 / norm;
    else
        return 0;
}

unsigned long long int CpGscorebedGraph(struct hash *hash, unsigned long long int *cnt, char *outfile){
    /* when same CpG site get two score, add them -- disabled temp -- abled
     * for overlapping site, solve by using single base score - seems no this situation?*/
    FILE *f = mustOpen(outfile, "w");
    struct hashEl *hel, *hel2, *hel3; //*hel4;
    struct hashCookie cookie = hashFirst(hash);
    unsigned long long int c = 0; //count for covered CpG site
    double score = 0;
    char key[100]; //key2[100];
    struct hash *cpgHash = newHash(0);
    while ( (hel = hashNext(&cookie)) != NULL ){
        struct mreFrag *mre = (struct mreFrag *) hel->val;
        if (mre->reads_count > 0){
            score = calCpGscore(mre, cnt);
            if (score > 0){
                struct cpgScore *cgC = malloc(sizeof(struct cpgScore));
                //struct cpgScore *cgG = malloc(sizeof(struct cpgScore));
                cgC->score = score;
                //cgG->score = score;
                strcpy(cgC->chr, mre->chr);
                //strcpy(cgG->chr, mre->chr);
                if (mre->head){
                    //fprintf(f, "%s\t%i\t%i\t%.4f\n", mre->chr, mre->start, mre->start + 2, score);
                    cgC->start = mre->start;
                    //cgG->start = mre->start + 1;
                } else {
                    //fprintf(f, "%s\t%i\t%i\t%.4f\n", mre->chr, mre->end - 2, mre->end, score);
                    cgC->start = mre->end - 2;
                    //cgG->start = mre->end - 1;
                }
                if (sprintf(key, "%s:%i", cgC->chr, cgC->start) < 0)
                    errAbort("Mem ERROR");
                //if (sprintf(key2, "%s:%i", cgG->chr, cgG->start) < 0)
                //    errAbort("Mem ERROR");
                hel2 = hashLookup(cpgHash, key);
                if (hel2 == NULL){
                    hashAdd(cpgHash, key, cgC);
                    c++; //CpG site my coverd by more than one MRE site
                } else {
                    struct cpgScore *cg1 = (struct cpgScore *) hel2->val;
                    cg1->score += score;
                }
                //hel4 = hashLookup(cpgHash, key2);
                //if (hel4 == NULL){
                //    hashAdd(cpgHash, key2, cgG);
                //} else {
                //    struct cpgScore *cg4 = (struct cpgScore *) hel4->val;
                //    cg4->score += score;
                //}
            }
        }
    }
    struct hashCookie cookie2 = hashFirst(cpgHash);
    while ((hel3 = hashNext(&cookie2)) != NULL){
        struct cpgScore *cg2 = (struct cpgScore *) hel3->val;
        fprintf(f, "%s\t%i\t%i\t%.4f\n", cg2->chr, cg2->start, cg2->start + 2, cg2->score);
    }
    carefulClose(&f);
    freeHash(&cpgHash);
    return c;
}

void fragmentStats(struct hash *hash, unsigned long long int *cnt2, unsigned int mapQ, unsigned long long int *cnt, unsigned long long int cnt1, char *outfile, int minlen, int maxlen, int win){
    FILE *f = mustOpen(outfile, "w");
    int highend = 0, lowend = 0, solosite = 0, soloreads = 0, pairsite = 0;
    struct slInt *fraglist = NULL, *frag;
    struct hashEl *hel, *hel2;
    struct hashCookie cookie = hashFirst(hash);
    while( (hel = hashNext(&cookie)) != NULL ){
        struct mreFrag *mre = (struct mreFrag *) hel->val;
        if (mre->reads_count > 0){
            if (mre->head){
                hel2 = hashLookup(hash, mre->pair);
                if (hel2 != NULL){
                    struct mreFrag *mre2 = (struct mreFrag *) hel2->val;
                    if (mre2->reads_count > 0){
                        //a valid fragment
                        pairsite++;
                        frag = slIntNew(mre->end - mre->start);
                        //fprintf(stdout, "%i\n", frag->val);
                        slAddHead(&fraglist, frag);
                        if (mre->reads_count >= mre2->reads_count){
                            highend += mre->reads_count;
                            lowend += mre2->reads_count;
                        }else{
                            highend += mre2->reads_count;
                            lowend += mre->reads_count;
                        }
                    } else{
                        //a head solo
                        soloreads += mre->reads_count;
                        solosite++;
                    }
                } else {
                    errAbort("MRE site Key error.\n");
                }
            } else {
                hel2 = hashLookup(hash, mre->pair);
                if (hel2 != NULL){
                    struct mreFrag *mre2 = (struct mreFrag *) hel2->val;
                    if (mre2->reads_count > 0){
                        //a valid fragment, but already processed
                        continue;
                    } else{
                        //a tail solo
                        soloreads += mre->reads_count;
                        solosite++;
                    }
                } else {
                    errAbort("MRE site Key error.\n");
                }
            }
        }
    }
    slReverse(&fraglist);
    //prepare histgraph
    int bins = (maxlen - minlen) / win;
    int histoMin = 0, histoMax = 0, j;
    struct range *scale = malloc(sizeof(struct range) * bins);
    for (j = 0; j < bins; j++){
        scale[j].s = minlen + j*win;
        scale[j].e = minlen + (j+1)*win;
        scale[j].histo = 0;
    }
    //scale[0]->s = 0; scale[0]->e = minlen; scale[0]->histo = 0;
    //scale[j]->s = maxlen; scale[j]->e = 0; scale[j]->histo = 0;

    for ( frag = fraglist; frag != NULL; frag = frag->next){
        if (frag->val <= minlen){
            histoMin++;
        } else if (frag->val > maxlen){
            histoMax++;
        } else{
            for (j = 0; j < bins; j++){
                if (frag->val > scale[j].s && frag->val <= scale[j].e){
                    (scale[j].histo)++;
                    break;
                }
            }
        }
    }
    
    fprintf(f, "total reads (pair): %llu\n", cnt2[0]);
    fprintf(f, "    read ends 1: %llu\n", cnt2[0]);
    fprintf(f, "    read ends 2: %llu\n", cnt2[1]);
    fprintf(f, "    mapped read ends 1: %llu\n", cnt2[2]);
    fprintf(f, "    mapped read ends 2: %llu\n", cnt2[3]);
    fprintf(f, "    used read ends 1: %llu\n", cnt2[4]);
    fprintf(f, "    used read ends 2: %llu\n", cnt2[5]);
    //fprintf(f, "non-redundant reads (pair): %llu\n\n", cnt2[8]);
    fprintf(f, "mappable reads (pair): %llu\n", cnt2[6]);
    fprintf(f, "unique mapped reads (pair) (mapQ >= %u): %llu\n", mapQ, cnt2[7]);
    fprintf(f, "mre filtered reads:\t%llu\n", cnt[0]+cnt[1]+cnt[2]+cnt[3]+cnt[4]);
    fprintf(f, "    CCGG reads:\t%llu\n", cnt[0]);
    fprintf(f, "    CCGC reads:\t%llu\n", cnt[1]);
    fprintf(f, "    GCGC reads:\t%llu\n", cnt[2]);
    fprintf(f, "    ACGT reads:\t%llu\n", cnt[3]);
    fprintf(f, "    CGCG reads:\t%llu\n", cnt[4]);
    fprintf(f, "    Unkown reads:\t%llu\n", cnt[5]);
    fprintf(f, "Sampled CpG sites:\t%llu\n", cnt1);
    fprintf(f, "solo ends:\t%i\n", solosite);
    fprintf(f, "    reads on solo ends:\t%i\n", soloreads);
    fprintf(f, "fragments:\t%i\n", pairsite);
    fprintf(f, "    reads in higher end of fragments:\t%i\n", highend);
    fprintf(f, "    reads in lower end of fragments:\t%i\n", lowend);
    fprintf(f, "fragment size distribution:\n");
    fprintf(f, "%s\t%10s\t%10s\t%c\n", "Scale", "Count", "Percent", '|');
    fprintf(f, "<=%i\t%10d\t%10.2f\t|%s\n", minlen, histoMin, histoMin*100.0/pairsite, print_bar((int)(histoMin*100.0/pairsite)));
    for (j = 0; j < bins; j++){
        fprintf(f, "%i\t%10d\t%10.2f\t|%s\n", scale[j].e, scale[j].histo, scale[j].histo*100.0/pairsite, print_bar((int)(scale[j].histo*100.0/pairsite)));
    }
    fprintf(f, ">%i\t%10d\t%10.2f\t|%s\n", maxlen, histoMax, histoMax*100.0/pairsite, print_bar((int)(histoMax*100.0/pairsite)));
    
    carefulClose(&f);
}

char *print_bar(int x){
    char *s;
    s = malloc(x);
    int i;
    for (i = 0; i < x; i++){
        s[i] = '*';
    }
    s[i] = '\0';
    return s;
}

void writeFilterOut(struct hash *hash, char *out, int readlist, int threshold, char *subfam, unsigned long long int reads_num){ 
    FILE *out_stream;
    out_stream = mustOpen(out, "w");
    int j = 0;
    struct hashEl *he;
    struct hashCookie cookie = hashFirst(hash);
    if (readlist)
        fprintf(out_stream, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "#chr", "start", "end", "length", "repName", "repClass", "repFamily", "readsCount", "RPKM", "RPM", "readsList");
    else
        fprintf(out_stream, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "#chr", "start", "end", "length", "repName", "repClass", "repFamily", "readsCount", "RPKM", "RPM");
    while ( (he = hashNext(&cookie)) != NULL ) {
        struct binKeeper *bk = (struct binKeeper *) he->val;
        struct binKeeperCookie becookie = binKeeperFirst(bk);
        struct binElement *be;
        while( (be = binKeeperNext(&becookie)) != NULL ){
            struct rmsk *os = (struct rmsk *) (be->val);
            int count = slCount(os->sl);
            if (readlist) {
                if (count >= threshold){
                    j++;
                    slReverse(&(os->sl));
                    char *s = slNameListToString(os->sl, ',');
                    fprintf(out_stream, "%s\t%d\t%d\t%d\t%s\t%s\t%s\t%d\t%.3f\t%.3f\t%s\n", os->chr, os->start, os->end, os->length, os->name, os->cname, os->fname, count, cal_rpkm((unsigned long long int)count,(unsigned long long int)os->length, reads_num), cal_rpm((unsigned long long int)count, reads_num), s);
                    freeMem(s);
                }
            } else {
                if (count >= threshold){
                    j++;
                    fprintf(out_stream, "%s\t%d\t%d\t%d\t%s\t%s\t%s\t%d\t%.3f\t%.3f\n", os->chr, os->start, os->end, os->length, os->name, os->cname, os->fname, count, cal_rpkm((unsigned long long int)count, (unsigned long long int)os->length, reads_num), cal_rpm((unsigned long long int)count, reads_num));
                }
            }
            slFreeList(&(os->sl));
        }
        binKeeperFree(&bk);
    }
    fclose(out_stream);
    fprintf(stderr, "* Total %d [%s] TEs have at least %d reads mapped.\n", j, subfam, threshold);
}

void writeFilterOutMRE(struct hash *hash, char *out, char *subfam, double scoreThreshold){ 
    FILE *out_stream;
    out_stream = mustOpen(out, "w");
    int j = 0;
    struct hashEl *he;
    struct hashCookie cookie = hashFirst(hash);
        fprintf(out_stream, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "#chr", "start", "end", "length", "repName", "repClass", "repFamily", "covered_CpG_site", "total_CpG_score");
    while ( (he = hashNext(&cookie)) != NULL ) {
        struct binKeeper *bk = (struct binKeeper *) he->val;
        struct binKeeperCookie becookie = binKeeperFirst(bk);
        struct binElement *be;
        while( (be = binKeeperNext(&becookie)) != NULL ){
            struct rmsk *os = (struct rmsk *) (be->val);
            double score = os->cpgTotalScore;
            if (score > scoreThreshold){
                j++;
                fprintf(out_stream, "%s\t%d\t%d\t%d\t%s\t%s\t%s\t%d\t%.3f\n", os->chr, os->start, os->end, os->length, os->name, os->cname, os->fname, os->cpgCount, os->cpgTotalScore);
            }
        }
        binKeeperFree(&bk);
    }
    fclose(out_stream);
    fprintf(stderr, "* Total %d [%s] TEs have CpG score larger than %.3f.\n", j, subfam, scoreThreshold);
}

void sortBedfile(char *bedfile) {
    struct lineFile *lf = NULL;
    FILE *f = NULL;
    struct bedLine *blList = NULL, *bl;
    char *line;
    int lineSize;

    lf = lineFileOpen(bedfile, TRUE);
    while (lineFileNext(lf, &line, &lineSize)){
        if (line[0] == '#')
            continue;
        bl = bedLineNew(line);
        slAddHead(&blList, bl);
    }
    lineFileClose(&lf);

    slSort(&blList, bedLineCmp);

    f = mustOpen(bedfile, "w");
    for (bl = blList; bl != NULL; bl = bl->next){
        fprintf(f, "%s\t%s\n", bl->chrom, bl->line);
        if (ferror(f)){
    	    perror("Writing error\n");
	    errAbort("%s is truncated, sorry.", bedfile);
	}
    }
    carefulClose(&f);
    bedLineFreeList(&blList);
}
