#include "generic.h"

#define ITERES_VERSION "0.2.8"

static int usage() {
    fprintf(stderr, "\n");
    fprintf(stderr, "Program: iteres (Get repeat alignment statistic from BAM/SAM file)\n");
    fprintf(stderr, "Version: %s\n\n", ITERES_VERSION);
    fprintf(stderr, "Usage:   iteres <command> [options]\n\n");
    fprintf(stderr, "Command: stat        get repeat alignment statistics\n");
    fprintf(stderr, "         filter      filter alignment statistic on repName/repFamily/repClass\n");
    fprintf(stderr, "         nearby      fetch nearby genes for locations from bed file\n");
    fprintf(stderr, "         density     generate genome density for ChipSeq data\n");
    fprintf(stderr, "\n");
    return 1;
}

int main(int argc, char *argv[]) {
    if (argc < 2) return usage();
    if (strcmp(argv[1], "stat") == 0) return main_stat(argc-1, argv+1);
    else if (strcmp(argv[1], "filter") == 0) return main_filter(argc-1, argv+1);
    else if (strcmp(argv[1], "nearby") == 0) return main_nearby(argc-1, argv+1);
    else if (strcmp(argv[1], "density") == 0) return main_density(argc-1, argv+1);
    else {
        fprintf(stderr, "[iteres] unrecognized command '%s'\n", argv[1]);
        return 1;
    }
    return 0;
}
