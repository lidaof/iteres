Program: **iteres** (Get repeat alignment statistic from BAM/SAM file)

Usage:   iteres <command> [options]

This program contains 3 modules (commands): 

  * stat:   get repeat alignment statistics for repeat subfamily, family and class
  * filter: filter alignment statistic on individual repeat locus for repeat subfamily, family or class
  * nearby: fetch nearby genes for locations from bed file by querying UCSC database

Version History:

Version 0.2.6
  * Paired-end bam/sam supported

Version 0.2.5
  * -C option added, user selected to add 'chr' as prefix to chromosome name
  * statistics of repeat family and class were available
  * end position of a read mapped to not calculated as start+read_length, using bam_calend function instead
  * some little bug fix

Version: 0.2.4
  * first version pushed to github
