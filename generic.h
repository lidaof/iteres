#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "binRange.h"
#include "basicBed.h"
#include "bigWig.h"
#include "obscure.h"
#include "jksql.h"
#include "sam.h"

//struct hold contens from rmsk line
struct rmsk {
    char *chr;
    unsigned int start, end, consensus_start, consensus_end, length;
    char *name, *fname, *cname;
    struct slName *sl;
    struct slName *sl_unique;
    unsigned int cpgCount;
    double cpgTotalScore;
};

//struct hold contents for repeat subfamily
struct rep {
    char *name, *fname, *cname;
    unsigned long long int read_count, read_count_unique, genome_count, total_length;
    unsigned int length;
    unsigned int *bp_total, *bp_total_unique;
    unsigned int cpgCount;
    double *cpgScore;
    double cpgTotalScore;
};

//struct hold contents for repeat family
struct repfam {
    char *fname, *cname;
    unsigned long long int read_count, read_count_unique, genome_count, total_length;
    unsigned int cpgCount;
    double cpgTotalScore;
};

//struct hold contents for repeat class
struct repcla {
    char *cname;
    unsigned long long int read_count, read_count_unique, genome_count, total_length;
    unsigned int cpgCount;
    double cpgTotalScore;
};

struct mreFrag {
    char pair[40], chr[20];
    unsigned long long int reads_count;
    int head, start, end;
    char site[5]; //dont forget \0
};

struct range {
    int s, e, histo;
};

struct cpgScore {
    int start;
    double score;
    char chr[20];
};

int stat_usage();
int main_stat(int argc, char *argv[]);
int filter_usage();
int main_filter (int argc, char *argv[]);
int nearby_usage();
int main_nearby(int argc, char *argv[]);
int density_usage();
int main_density(int argc, char *argv[]);
int cpg_usage();
int main_cpg(int argc, char *argv[]);
int cpgstat_usage();
int main_cpgstat(int argc, char *argv[]);
int cpgfilter_usage();
int main_cpgfilter (int argc, char *argv[]);


char *get_filename_without_ext(char *filename);
char *get_filename_ext(char *filename);
double cal_rpkm (unsigned long long int reads_count, unsigned long long int total_length, unsigned long long int mapped_reads_num);
double cal_rpm (unsigned long long int reads_count, unsigned long long int mapped_reads_num);

void writeReport(char *outfile, unsigned long long int *cnt, unsigned int mapQ, char *subfam);
void writeWigandStat(struct hash *hash, struct hash *hash1, struct hash *hash2, char *of1, char *of2, char *of3, char *of4, char *of5, unsigned long long int reads_num, unsigned long int reads_num_unique);
unsigned long long int *samFile2nodupRepbedFile(char *samfile, struct hash *chrHash, struct hash *hashRmsk, struct hash *hashRep, struct hash *hashFam, struct hash *hashCla, int isSam, unsigned int mapQ, int filter, int rmDup, int addChr);
float getCov(unsigned int aStart, unsigned int aEnd, unsigned int start, unsigned int end);
unsigned long long int *samFile2nodupRepbedFileNew(char *samfile, struct hash *chrHash, struct hash *hashRmsk, struct hash *hashRep, struct hash *hashFam, struct hash *hashCla, int isSam, unsigned int mapQ, int filter, int rmDup, int addChr, int discardWrongEnd, unsigned int iSize, unsigned int extension, float minCoverage, int treat, char *outbed, char *outbed_unique);
unsigned long long int *PEsamFile2nodupRepbedFile(char *samfile, struct hash *chrHash, struct hash *hashRmsk, struct hash *hashRep, struct hash *hashFam, struct hash *hashCla, int isSam, unsigned int mapQ, int filter, int rmDup, int addChr, unsigned int iSize);
void freermsk(struct rmsk *s);
void rmsk2binKeeperHash(char *rmskfile, struct hash *chrHash, struct hash *repHash, struct hash **hashRmsk, struct hash **hashRep, struct hash **hashFam, struct hash **hashCla, int filterField, char *filterName);
void writeFilterOut(struct hash *hash, char *out, int readlist, int threshold, char *subfam, unsigned long long int reads_num);
void sortBedfile(char *bedfile);
void writeReportDensity(char *outfile, unsigned long long int *cnt, unsigned int mapQ);
unsigned long long int *sam2bed(char *samfile, char *outbed, struct hash *chrHash, int isSam, unsigned int mapQ, int rmDup, int addChr, int discardWrongEnd, unsigned int iSize, unsigned int extension, int treat);
struct hash *MREfrag2Hash (char *fragfile, int minlen, int maxlen);
unsigned long long int *filterReadByMREsite(struct hash *hash, char *inBed, char *outBed, int call);
double calCpGscore (struct mreFrag *mre, unsigned long long int *cnt);
unsigned long long int  CpGscorebedGraph(struct hash *hash, unsigned long long int *cnt, char *outfile);
void fragmentStats(struct hash *hash, unsigned long long int *cnt2, unsigned int mapQ, unsigned long long int *cnt, unsigned long long int cnt1, char *outfile, int minlen, int maxlen, int win);
char *print_bar(int x);
void cpgBedGraphOverlapRepeat(char *cpgBedGraphFile, struct hash *hashRmsk, struct hash *hashRep, struct hash *hashFam, struct hash *hashCla, int filter);
void MREwriteWigandStat(struct hash *hash, struct hash *hash1, struct hash *hash2, char *of1, char *of2, char *of3, char *of4);
void writeFilterOutMRE(struct hash *hash, char *out, char *subfam, double scoreThreshold); 
