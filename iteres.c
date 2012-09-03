#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "binRange.h"
#include "bigWig.h"
#include "obscure.h"
#include "jksql.h"
#include "sam.h"
#define ITERES_VERSION "0.2.6"

//struct hold contens from rmsk line
struct rmsk {
    char *chr;
    unsigned int start, end, consensus_start, consensus_end, length;
    char *name, *fname, *cname;
    struct slName *sl;
};

//struct hold contents for repeat subfamily
struct rep {
    char *name, *fname, *cname;
    unsigned long long int read_count, genome_count, total_length;
    unsigned int length;
    unsigned int *bp_total;
};

//struct hold contents for repeat family
struct repfam {
    char *fname, *cname;
    unsigned long long int read_count, genome_count, total_length;
};

//struct hold contents for repeat class
struct repcla {
    char *cname;
    unsigned long long int read_count, genome_count, total_length;
};

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
    fprintf(f, "Total reads: %llu\n", cnt[0]);
    fprintf(f, "Mapped reads (mapQ >= %u): %llu\n", mapQ, cnt[1]);
    fprintf(f, "Used reads: %llu\n", cnt[2]);
    fprintf(f, "Non-Redundant reads: %llu\n", cnt[3]);
    fprintf(f, "Overlap-[%s]-repeats reads: %llu\n", subfam, cnt[4]);
    carefulClose(&f);
}

void writeWigandStat(struct hash *hash, struct hash *hash1, struct hash *hash2, char *of1, char *of2, char *of3, char *of4, unsigned long long int reads_num){
    FILE *f1 = mustOpen(of1, "w");
    FILE *f2 = mustOpen(of2, "w");
    unsigned int m;
    struct hashEl *helr;
    struct hashCookie cookier = hashFirst(hash);
    fprintf(f1, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "#subfamily", "family", "class", "consensus_length", "reads_count", "total_length", "genome_count", "RPKM", "RPM");
    while ( (helr = hashNext(&cookier)) != NULL ) {
        struct rep *or = (struct rep *) (helr->val);
        fprintf(f1, "%s\t%s\t%s\t%u\t%llu\t%llu\t%llu\t%.3f\t%.3f\n", or->name, or->fname, or->cname, or->length, or->read_count, or->total_length, or->genome_count, cal_rpkm(or->read_count, or->total_length, reads_num), cal_rpm(or->read_count, reads_num));
        if (or->length != 0){
            fprintf(f2, "fixedStep chrom=%s start=1 step=1 span=1\n", or->name);
            for (m = 0; m < or->length; m++) 
                fprintf(f2, "%u\n", (or->bp_total)[m]);
        }
    }
    carefulClose(&f2);
    carefulClose(&f1);
    FILE *f3 = mustOpen(of3, "w");
    struct hashEl *hel3;
    struct hashCookie cookier3 = hashFirst(hash1);
    fprintf(f3, "%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "#family", "class", "reads_count", "total_length", "genome_count", "RPKM", "RPM");
    while ( (hel3 = hashNext(&cookier3)) != NULL) {
        struct repfam *or = (struct repfam *) hel3->val;
        fprintf(f3, "%s\t%s\t%llu\t%llu\t%llu\t%.3f\t%.3f\n", or->fname, or->cname, or->read_count, or->total_length, or->genome_count, cal_rpkm(or->read_count, or->total_length, reads_num), cal_rpm(or->read_count, reads_num));
    }
    carefulClose(&f3);
    FILE *f4 = mustOpen(of4, "w");
    struct hashEl *hel4;
    struct hashCookie cookier4 = hashFirst(hash2);
    fprintf(f4, "%s\t%s\t%s\t%s\t%s\t%s\n",  "#class", "reads_count", "total_length", "genome_count", "RPKM", "RPM");
    while ( (hel4 = hashNext(&cookier4)) != NULL) {
        struct repcla *or = (struct repcla *) hel4->val;
        fprintf(f4, "%s\t%llu\t%llu\t%llu\t%.3f\t%.3f\n", or->cname, or->read_count, or->total_length, or->genome_count, cal_rpkm(or->read_count, or->total_length, reads_num), cal_rpm(or->read_count, reads_num));
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
        if ( sameString (bam1_qname(b), prn)) 
            continue;
        reads_num++;
        if ((reads_num % 10000) == 0)
            fprintf(stderr, "\r* Processed reads: %llu", reads_num);
        strcpy(prn, bam1_qname(b));
        if (b->core.tid < 0)
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
            if (cur->core.tid < 0 || pre->core.tid < 0)
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
            ns->length = hashIntValDefault(repHash, ns->name, 0);
            ns->bp_total = malloc(sizeof(unsigned int) * ns->length);
            int k;
            for (k = 0; k < ns->length; k++) 
                (ns->bp_total)[k] = 0;
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

int stat_usage(){
    fprintf(stderr, "\n");
    fprintf(stderr, "Obtain alignment statistics for each repeat subfamily, family and class.\n\n");
    fprintf(stderr, "Usage:   iteres stat [options] <chromosome size file> <repeat size file> <rmsk.txt> <bam/sam alignment file>\n\n");
    fprintf(stderr, "Options: -S       input is SAM [off]\n");
    fprintf(stderr, "         -Q       mapping Quality threshold [0]\n");
    fprintf(stderr, "         -N       normalized by number of (0: repeat reads, 1: non-redundant reads, 2: mapped reads, 3: total reads) [0])\n");
    fprintf(stderr, "         -D       remove redundant reads [off]\n");
    fprintf(stderr, "         -w       keep the wiggle file [0]\n");
    fprintf(stderr, "         -C       Add 'chr' string as prefix of reference sequence [off]\n");
    fprintf(stderr, "         -P       Input was paired end data [off]\n");
    fprintf(stderr, "         -I       Insert length [1000] (force -P)\n");
    fprintf(stderr, "         -o       output prefix [basename of input without extension]\n");
    fprintf(stderr, "         -h       help message\n");
    fprintf(stderr, "         -?       help message\n");
    fprintf(stderr, "\n");
    return 1;
}

/* main stat funtion */
int main_stat (int argc, char *argv[]) {
    
    char *output, *outReport, *outWig, *outbigWig, *outStat, *outFam, *outCla;
    unsigned long long int *cnt;
    int optSam = 0, optkeepWig = 0, c, optDup = 0, optaddChr = 0, optPair = 0;
    unsigned int optQual = 0, optNorm = 0, optisize = 1000;
    char *optoutput = NULL;
    time_t start_time, end_time;
    struct hash *hashRmsk = newHash(0);
    struct hash *hashRep = newHash(0);
    struct hash *hashFam = newHash(0);
    struct hash *hashCla = newHash(0);
    start_time = time(NULL);
    while ((c = getopt(argc, argv, "SQ:N:DwCo:PI:h?")) >= 0) {
        switch (c) {
            case 'S': optSam = 1; break;
            case 'Q': optQual = (unsigned int)strtol(optarg, 0, 0); break;
            case 'N': optNorm = (unsigned int)strtol(optarg, 0, 0); break;
            case 'D': optDup = 1; break;
            case 'w': optkeepWig = 1; break;
            case 'C': optaddChr = 1; break;
            case 'P': optPair = 1; break;
            case 'I': optisize = (unsigned int)strtol(optarg, 0, 0); optPair = 1; break;
            case 'o': optoutput = strdup(optarg); break;
            case 'h':
            case '?': return stat_usage(); break;
            default: return 1;
        }
    }
    if (optind + 4 > argc)
        return stat_usage();

    char *chr_size_file = argv[optind];
    char *rep_size_file = argv[optind+1];
    char *rmsk_file = argv[optind+2];
    char *sam_file = argv[optind+3];
    
    if(optoutput) {
        output = optoutput;
    } else {
        output = cloneString(get_filename_without_ext(basename(sam_file)));
    }
    
    if(asprintf(&outWig, "%s.iteres.wig", output) < 0)
        errAbort("Mem Error.\n");
    if(asprintf(&outbigWig, "%s.iteres.bigWig", output) < 0)
        errAbort("Mem Error.\n");
    if (asprintf(&outReport, "%s.iteres.report", output) < 0)
        errAbort("Preparing output wrong");
    if (asprintf(&outStat, "%s.iteres.subfamily.stat", output) < 0)
        errAbort("Preparing output wrong");
    if (asprintf(&outFam, "%s.iteres.family.stat", output) < 0)
        errAbort("Preparing output wrong");
    if (asprintf(&outCla, "%s.iteres.class.stat", output) < 0)
        errAbort("Preparing output wrong");
    
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
    
    //if (optPair){
    //    if (optisize == 0)
    //        optisize = 1000;
    //}

    struct hash *chrHash = hashNameIntFile(chr_size_file);
    struct hash *repHash = hashNameIntFile(rep_size_file);
    
    fprintf(stderr, "* Start to parse the rmsk file\n");
    rmsk2binKeeperHash(rmsk_file, chrHash, repHash, &hashRmsk, &hashRep, &hashFam, &hashCla, 0, "ALL");
    
    //sam file
    fprintf(stderr, "* Start to parse the SAM/BAM file\n");
    if (optPair){
        cnt = PEsamFile2nodupRepbedFile(sam_file, chrHash, hashRmsk, hashRep, hashFam, hashCla, optSam, optQual, 0, optDup, optaddChr, optisize);
    } else {
        cnt = samFile2nodupRepbedFile(sam_file, chrHash, hashRmsk, hashRep, hashFam, hashCla, optSam, optQual, 0, optDup, optaddChr);
    }

    fprintf(stderr, "* Writing stats and Wig file\n");
    writeWigandStat(hashRep, hashFam, hashCla, outStat, outWig, outFam, outCla, cnt[nindex]);

    fprintf(stderr, "* Generating bigWig files\n");
    bigWigFileCreate(outWig, rep_size_file, 256, 1024, 0, 1, outbigWig);

    //write report file
    fprintf(stderr, "* Preparing report file\n");
    writeReport(outReport, cnt, optQual, "ALL");
    
    if (!optkeepWig)
        unlink(outWig);
    
    //cleaning
    hashFree(&chrHash);
    hashFree(&repHash);
    hashFree(&hashRmsk);
    hashFree(&hashRep);
    hashFree(&hashFam);
    hashFree(&hashCla);
    free(outWig);
    free(outbigWig);
    free(outReport);
    free(outStat);
    free(outFam);
    free(outCla);
    end_time = time(NULL);
    fprintf(stderr, "* Done, time used %.0f seconds.\n", difftime(end_time, start_time));
    return 0;
}

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

/* main filter funtion */
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

int nearby_usage(){
    fprintf(stderr, "\n");
    fprintf(stderr, "Obtain nearby genes from locations listed in a bed file by querying UCSC database.\n\n");
    fprintf(stderr, "Usage:   iteres nearby [options] <bed file>\n\n");
    fprintf(stderr, "Options: -d       database to query [hg19]\n");
    fprintf(stderr, "         -n       output how many genes each direction [1]\n");
    fprintf(stderr, "         -o       output prefix [basename of input without extension]\n");
    fprintf(stderr, "         -h       help message\n");
    fprintf(stderr, "         -?       help message\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Note: the bed file should contain at least 3 fields which were [chr] [start] [end]\n");
    fprintf(stderr, "      also you need to have an internet connection\n\n");
    return 1;
}

/* main nearby funtion */
int main_nearby(int argc, char *argv[]){
    struct lineFile *infile_stream;
    FILE *out_stream;
    char *output, *outfile, *line, *row[20];
    int c;
    char *optoutput = NULL;
    char *optdb = "hg19";
    unsigned int optthreshold = 1;
    time_t start_time, end_time;
    start_time = time(NULL);
    while ((c = getopt(argc, argv, "d:n:o:h?")) >= 0) {
        switch (c) {
            case 'n': optthreshold = (unsigned int)strtol(optarg, 0, 0); break;
            case 'd': optdb = strdup(optarg); break;
            case 'o': optoutput = strdup(optarg); break;
            case 'h':
            case '?': return nearby_usage(); break;
            default: return 1;
        }
    }
    if (optind + 1 > argc)
        return nearby_usage();
    char *infile = argv[optind];
    if(optoutput) {
        output = optoutput;
    } else {
        output = cloneString(get_filename_without_ext(basename(infile)));
    }
    if(asprintf(&outfile, "%s.iteres.gene", output) < 0)
        errAbort("Mem Error.\n");
    out_stream = mustOpen(outfile, "w");
    infile_stream = lineFileOpen(infile, TRUE);
    fprintf(stderr, "* Making connection to UCSC database %s\n", optdb);
    struct sqlConnection *con = sqlConnectRemote("genome-mysql.cse.ucsc.edu", "genome", "", optdb);
    fprintf(stderr, "* Start to parse the input file\n");
    char *tchr, *query1, *query2;
    int tstart, tend;
    while( lineFileNextReal(infile_stream, &line)){
        int numFields = chopByWhite(line, row, ArraySize(row));
        if (numFields < 3)
            errAbort("file %s doesn't appear to be in bed format. At least 3 fields required, got %d", infile, numFields);
        tchr = cloneString(row[0]);
        tstart = (int)strtol(row[1], NULL, 0);
        tend = (int)strtol(row[2], NULL, 0);
        
        //upstream
        char **row1;
        if (asprintf(&query1, "select e.chrom,e.txStart,e.txEnd,e.alignID,j.geneSymbol FROM \
                   knownGene e, \
                   kgXref j \
                   WHERE e.alignID = j.kgID AND e.chrom='%s' AND e.txEnd < %d \
                   ORDER BY e.txEnd DESC limit %d", tchr, tstart, optthreshold) < 0)
            errAbort("Memory allocating for query string wrong");

        struct sqlResult *sr1 = sqlGetResult(con, query1);
        while ( (row1 = sqlNextRow(sr1)) != NULL ){
            fprintf(out_stream, "%s\t%s\t%s\t%s\t%s\t%s\n", row1[0], row1[1], row1[2], row1[3], row1[4], "upstream");
        }
        sqlFreeResult(&sr1);
        
        //downstream
        char **row2;
        if (asprintf(&query2, "select e.chrom,e.txStart,e.txEnd,e.alignID,j.geneSymbol FROM \
                   knownGene e, \
                   kgXref j \
                   WHERE e.alignID = j.kgID AND e.chrom='%s' AND e.txStart > %d \
                   ORDER BY e.txStart ASC limit %d", tchr, tend, optthreshold) < 0)
            errAbort("Memory allocating for query string wrong");

        struct sqlResult *sr2 = sqlGetResult(con, query2);
        while ( (row2 = sqlNextRow(sr2)) != NULL ){
            fprintf(out_stream, "%s\t%s\t%s\t%s\t%s\t%s\n", row2[0], row2[1], row2[2], row2[3], row2[4], "downstream");
        }
        sqlFreeResult(&sr2);
    }
    lineFileClose(&infile_stream);
    sqlDisconnect(&con);
    fclose(out_stream);
    end_time = time(NULL);
    fprintf(stderr, "* Done, time used %.0f seconds.\n", difftime(end_time, start_time));
    return 0;
}

static int usage() {
    fprintf(stderr, "\n");
    fprintf(stderr, "Program: iteres (Get repeat alignment statistic from BAM/SAM file)\n");
    fprintf(stderr, "Version: %s\n\n", ITERES_VERSION);
    fprintf(stderr, "Usage:   iteres <command> [options]\n\n");
    fprintf(stderr, "Command: stat        get repeat alignment statistics\n");
    fprintf(stderr, "         filter      filter alignment statistic on repName/repFamily/repClass\n");
    fprintf(stderr, "         nearby      fetch nearby genes for locations from bed file\n");
    fprintf(stderr, "\n");
    return 1;
}

int main(int argc, char *argv[]) {
    if (argc < 2) return usage();
    if (strcmp(argv[1], "stat") == 0) return main_stat(argc-1, argv+1);
    else if (strcmp(argv[1], "filter") == 0) return main_filter(argc-1, argv+1);
    else if (strcmp(argv[1], "nearby") == 0) return main_nearby(argc-1, argv+1);
    else {
        fprintf(stderr, "[iteres] unrecognized command '%s'\n", argv[1]);
        return 1;
    }
    return 0;
}
