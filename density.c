#include "generic.h"
#include "sqlNum.h"
#include "localmem.h"
#include "dystring.h"
#include "cirTree.h"
#include "sig.h"
#include "zlibFace.h"
#include "bPlusTree.h"
#include "bbiFile.h"
#include "bwgInternal.h"

typedef unsigned int unitSize;

#define MAXCOUNT (unitSize)~0
#define MAXMESSAGE "Overflow of overlap counts. Max is %lu.  Recompile with bigger unitSize or use -max option"
#define INCWOVERFLOW(countArray,x) if(countArray[x] == MAXCOUNT) {if(!doMax) errAbort(MAXMESSAGE,(unsigned long)MAXCOUNT);} else countArray[x]++

boolean doMax = FALSE;   /* if overlap count will overflow, just keep max */
boolean doZero = FALSE;  /* add blocks with 0 counts */
boolean doBed12 = FALSE;  /* expect bed12 and process block by block */
boolean doOutBounds = FALSE;  /* output min/max to stderr */
unitSize overMin = ~1;
unitSize overMax = 0;

static int blockSize = 256;
static int itemsPerSlot = 1024;
static boolean doCompress = TRUE;

struct sectionItem
/* An item in a section of a bedGraph. */
    {
    bits32 start, end;			/* Position in chromosome, half open. */
    float val;				/* Single precision value. */
    };

void writeSections(struct bbiChromUsage *usageList, struct lineFile *lf, 
	int itemsPerSlot, struct bbiBoundsArray *bounds, int sectionCount, FILE *f,
	int resTryCount, int resScales[], int resSizes[], 
	boolean doCompress, bits32 *retMaxSectionSize)
/* Read through lf, chunking it into sections that get written to f.  Save info
 * about sections in bounds. */
{
int maxSectionSize = 0;
struct bbiChromUsage *usage = usageList;
int itemIx = 0, sectionIx = 0;
bits32 reserved32 = 0;
UBYTE reserved8 = 0;
struct sectionItem items[itemsPerSlot];
struct sectionItem *lastB = NULL;
bits32 resEnds[resTryCount];
int resTry;
for (resTry = 0; resTry < resTryCount; ++resTry)
    resEnds[resTry] = 0;
struct dyString *stream = dyStringNew(0);

/* remove initial browser and track lines */
lineFileRemoveInitialCustomTrackLines(lf);

for (;;)
    {
    /* Get next line of input if any. */
    char *row[5];
    int rowSize = lineFileChopNext(lf, row, ArraySize(row));

    /* Figure out whether need to output section. */
    boolean sameChrom = FALSE;
    if (rowSize > 0)
	sameChrom = sameString(row[0], usage->name);
    if (itemIx >= itemsPerSlot || rowSize == 0 || !sameChrom)
        {
	/* Figure out section position. */
	bits32 chromId = usage->id;
	bits32 sectionStart = items[0].start;
	bits32 sectionEnd = items[itemIx-1].end;

	/* Save section info for indexing. */
	assert(sectionIx < sectionCount);
	struct bbiBoundsArray *section = &bounds[sectionIx++];
	section->offset = ftell(f);
	section->range.chromIx = chromId;
	section->range.start = sectionStart;
	section->range.end = sectionEnd;

	/* Output section header to stream. */
	dyStringClear(stream);
	UBYTE type = bwgTypeBedGraph;
	bits16 itemCount = itemIx;
	dyStringWriteOne(stream, chromId);			// chromId
	dyStringWriteOne(stream, sectionStart);		// start
	dyStringWriteOne(stream, sectionEnd);	// end
	dyStringWriteOne(stream, reserved32);		// itemStep
	dyStringWriteOne(stream, reserved32);		// itemSpan
	dyStringWriteOne(stream, type);			// type
	dyStringWriteOne(stream, reserved8);			// reserved
	dyStringWriteOne(stream, itemCount);			// itemCount

	/* Output each item in section to stream. */
	int i;
	for (i=0; i<itemIx; ++i)
	    {
	    struct sectionItem *item = &items[i];
	    dyStringWriteOne(stream, item->start);
	    dyStringWriteOne(stream, item->end);
	    dyStringWriteOne(stream, item->val);
	    }

	/* Save stream to file, compressing if need be. */
	if (stream->stringSize > maxSectionSize)
	    maxSectionSize = stream->stringSize;
	if (doCompress)
	    {
	    size_t maxCompSize = zCompBufSize(stream->stringSize);
	    char compBuf[maxCompSize];
	    int compSize = zCompress(stream->string, stream->stringSize, compBuf, maxCompSize);
	    mustWrite(f, compBuf, compSize);
	    }
	else
	    mustWrite(f, stream->string, stream->stringSize);


	/* If at end of input we are done. */
	if (rowSize == 0)
	    break;

	/* Set up for next section. */
	itemIx = 0;

	if (!sameChrom)
	    {
	    usage = usage->next;
	    assert(usage != NULL);
            if (!sameString(row[0], usage->name))
                errAbort("read %s, expecting %s on line %d in file %s\n", 
                    row[0], usage->name, lf->lineIx, lf->fileName);
	    assert(sameString(row[0], usage->name));
	    lastB = NULL;
	    for (resTry = 0; resTry < resTryCount; ++resTry)
		resEnds[resTry] = 0;
	    }
	}

    /* Parse out input. */
    lineFileExpectWords(lf, 4, rowSize);
    bits32 start = lineFileNeedNum(lf, row, 1);
    bits32 end = lineFileNeedNum(lf, row, 2);
    float val = lineFileNeedDouble(lf, row, 3);

    /* Verify that inputs meets our assumption - that it is a sorted bedGraph file. */
    if (start > end)
        errAbort("Start (%u) after end (%u) line %d of %s", start, end, lf->lineIx, lf->fileName);
    if (lastB != NULL)
        {
	if (lastB->start > start)
	    errAbort("BedGraph not sorted on start line %d of %s", lf->lineIx, lf->fileName);
	if (lastB->end > start)
	    errAbort("Overlapping regions in bedGraph line %d of %s", lf->lineIx, lf->fileName);
	}


    /* Do zoom counting. */
    for (resTry = 0; resTry < resTryCount; ++resTry)
        {
	bits32 resEnd = resEnds[resTry];
	if (start >= resEnd)
	    {
	    resSizes[resTry] += 1;
	    resEnds[resTry] = resEnd = start + resScales[resTry];
	    }
	while (end > resEnd)
	    {
	    resSizes[resTry] += 1;
	    resEnds[resTry] = resEnd = resEnd + resScales[resTry];
	    }
	}

    /* Save values in output array. */
    struct sectionItem *b = &items[itemIx];
    b->start = start;
    b->end = end;
    b->val = val;
    lastB = b;
    itemIx += 1;
    }
assert(sectionIx == sectionCount);

*retMaxSectionSize = maxSectionSize;
}

static struct bbiSummary *writeReducedOnceReturnReducedTwice(struct bbiChromUsage *usageList, 
	struct lineFile *lf, bits32 initialReduction, bits32 initialReductionCount, 
	int zoomIncrement, int blockSize, int itemsPerSlot, boolean doCompress,
	struct lm *lm, FILE *f, bits64 *retDataStart, bits64 *retIndexStart,
	struct bbiSummaryElement *totalSum)
/* Write out data reduced by factor of initialReduction.  Also calculate and keep in memory
 * next reduction level.  This is more work than some ways, but it keeps us from having to
 * keep the first reduction entirely in memory. */
{
struct bbiSummary *twiceReducedList = NULL;
bits32 doubleReductionSize = initialReduction * zoomIncrement;
struct bbiChromUsage *usage = usageList;
struct bbiSummary oneSummary, *sum = NULL;
struct bbiBoundsArray *boundsArray, *boundsPt, *boundsEnd;
boundsPt = AllocArray(boundsArray, initialReductionCount);
boundsEnd = boundsPt + initialReductionCount;

*retDataStart = ftell(f);
writeOne(f, initialReductionCount);
boolean firstRow = TRUE;

struct bbiSumOutStream *stream = bbiSumOutStreamOpen(itemsPerSlot, f, doCompress);

/* remove initial browser and track lines */
lineFileRemoveInitialCustomTrackLines(lf);

for (;;)
    {
    /* Get next line of input if any. */
    char *row[5];
    int rowSize = lineFileChopNext(lf, row, ArraySize(row));

    /* Output last section and break if at end of file. */
    if (rowSize == 0 && sum != NULL)
	{
	bbiOutputOneSummaryFurtherReduce(sum, &twiceReducedList, doubleReductionSize, 
		&boundsPt, boundsEnd, usage->size, lm, stream);
	break;
	}

    /* Parse out row. */
    char *chrom = row[0];
    bits32 start = sqlUnsigned(row[1]);
    bits32 end = sqlUnsigned(row[2]);
    float val = sqlFloat(row[3]);

    /* Update total summary stuff. */
    bits32 size = end-start;
    if (firstRow)
	{
        totalSum->validCount = size;
	totalSum->minVal = totalSum->maxVal = val;
	totalSum->sumData = val*size;
	totalSum->sumSquares = val*val*size;
	firstRow = FALSE;
	}
    else
        {
	totalSum->validCount += size;
	if (val < totalSum->minVal) totalSum->minVal = val;
	if (val > totalSum->maxVal) totalSum->maxVal = val;
	totalSum->sumData += val*size;
	totalSum->sumSquares += val*val*size;
	}

    /* If new chromosome output existing block. */
    if (differentString(chrom, usage->name))
        {
	usage = usage->next;
	bbiOutputOneSummaryFurtherReduce(sum, &twiceReducedList, doubleReductionSize,
		&boundsPt, boundsEnd, usage->size, lm, stream);
	sum = NULL;
	}

    /* If start past existing block then output it. */
    else if (sum != NULL && sum->end <= start)
	{
	bbiOutputOneSummaryFurtherReduce(sum, &twiceReducedList, doubleReductionSize, 
		&boundsPt, boundsEnd, usage->size, lm, stream);
	sum = NULL;
	}

    /* If don't have a summary we're working on now, make one. */
    if (sum == NULL)
        {
	oneSummary.chromId = usage->id;
	oneSummary.start = start;
	oneSummary.end = start + initialReduction;
	if (oneSummary.end > usage->size) oneSummary.end = usage->size;
	oneSummary.minVal = oneSummary.maxVal = val;
	oneSummary.sumData = oneSummary.sumSquares = 0.0;
	oneSummary.validCount = 0;
	sum = &oneSummary;
	}
    
    /* Deal with case where might have to split an item between multiple summaries.  This
     * loop handles all but the final affected summary in that case. */
    while (end > sum->end)
        {
	verbose(3, "Splitting start %d, end %d, sum->start %d, sum->end %d\n", start, end, sum->start, sum->end);
	/* Fold in bits that overlap with existing summary and output. */
	bits32 overlap = rangeIntersection(start, end, sum->start, sum->end);
	sum->validCount += overlap;
	if (sum->minVal > val) sum->minVal = val;
	if (sum->maxVal < val) sum->maxVal = val;
	sum->sumData += val * overlap;
	sum->sumSquares += val*val * overlap;
	bbiOutputOneSummaryFurtherReduce(sum, &twiceReducedList, doubleReductionSize, 
		&boundsPt, boundsEnd, usage->size, lm, stream);
	size -= overlap;

	/* Move summary to next part. */
	sum->start = start = sum->end;
	sum->end = start + initialReduction;
	if (sum->end > usage->size) sum->end = usage->size;
	sum->minVal = sum->maxVal = val;
	sum->sumData = sum->sumSquares = 0.0;
	sum->validCount = 0;
	}

    /* Add to summary. */
    sum->validCount += size;
    if (sum->minVal > val) sum->minVal = val;
    if (sum->maxVal < val) sum->maxVal = val;
    sum->sumData += val * size;
    sum->sumSquares += val*val * size;
    }
bbiSumOutStreamClose(&stream);

/* Write out 1st zoom index. */
int indexOffset = *retIndexStart = ftell(f);
assert(boundsPt == boundsEnd);
cirTreeFileBulkIndexToOpenFile(boundsArray, sizeof(boundsArray[0]), initialReductionCount,
    blockSize, itemsPerSlot, NULL, bbiBoundsArrayFetchKey, bbiBoundsArrayFetchOffset, 
    indexOffset, f);

freez(&boundsArray);
slReverse(&twiceReducedList);
return twiceReducedList;
}

void bedGraphToBigWig(char *inName, char *chromSizes, char *outName)
/* bedGraphToBigWig - Convert a bedGraph program to bigWig.. */
{
verboseTimeInit();
struct lineFile *lf = lineFileOpen(inName, TRUE);
struct hash *chromSizesHash = bbiChromSizesFromFile(chromSizes);
verbose(2, "%d chroms in %s\n", chromSizesHash->elCount, chromSizes);
int minDiff = 0, i;
double aveSize = 0;
bits64 bedCount = 0;
bits32 uncompressBufSize = 0;
struct bbiChromUsage *usageList = bbiChromUsageFromBedFile(lf, chromSizesHash, &minDiff, &aveSize, &bedCount);
verboseTime(2, "pass1");
verbose(2, "%d chroms in %s\n", slCount(usageList), inName);

/* Write out dummy header, zoom offsets. */
FILE *f = mustOpen(outName, "wb");
bbiWriteDummyHeader(f);
bbiWriteDummyZooms(f);

/* Write out dummy total summary. */
struct bbiSummaryElement totalSum;
ZeroVar(&totalSum);
bits64 totalSummaryOffset = ftell(f);
bbiSummaryElementWrite(f, &totalSum);

/* Write out chromosome/size database. */
bits64 chromTreeOffset = ftell(f);
bbiWriteChromInfo(usageList, blockSize, f);

/* Set up to keep track of possible initial reduction levels. */
int resTryCount = 10, resTry;
int resIncrement = 4;
int resScales[resTryCount], resSizes[resTryCount];
int res = minDiff * 2;
if (res > 0)
    {
    for (resTry = 0; resTry < resTryCount; ++resTry)
	{
	resSizes[resTry] = 0;
	resScales[resTry] = res;
	res *= resIncrement;
	}
    }
else
    resTryCount = 0;

/* Write out primary full resolution data in sections, collect stats to use for reductions. */
bits64 dataOffset = ftell(f);
bits64 sectionCount = bbiCountSectionsNeeded(usageList, itemsPerSlot);
writeOne(f, sectionCount);
struct bbiBoundsArray *boundsArray;
AllocArray(boundsArray, sectionCount);
lineFileRewind(lf);
bits32 maxSectionSize = 0;
writeSections(usageList, lf, itemsPerSlot, boundsArray, sectionCount, f,
	resTryCount, resScales, resSizes, doCompress, &maxSectionSize);
verboseTime(2, "pass2");

/* Write out primary data index. */
bits64 indexOffset = ftell(f);
cirTreeFileBulkIndexToOpenFile(boundsArray, sizeof(boundsArray[0]), sectionCount,
    blockSize, 1, NULL, bbiBoundsArrayFetchKey, bbiBoundsArrayFetchOffset, 
    indexOffset, f);
verboseTime(2, "index write");

/* Declare arrays and vars that track the zoom levels we actually output. */
bits32 zoomAmounts[bbiMaxZoomLevels];
bits64 zoomDataOffsets[bbiMaxZoomLevels];
bits64 zoomIndexOffsets[bbiMaxZoomLevels];
int zoomLevels = 0;

/* Write out first zoomed section while storing in memory next zoom level. */
if (minDiff > 0)
    {
    bits64 dataSize = indexOffset - dataOffset;
    int maxReducedSize = dataSize/2;
    int initialReduction = 0, initialReducedCount = 0;

    /* Figure out initialReduction for zoom. */
    for (resTry = 0; resTry < resTryCount; ++resTry)
	{
	bits64 reducedSize = resSizes[resTry] * sizeof(struct bbiSummaryOnDisk);
	if (doCompress)
	    reducedSize /= 2;	// Estimate!
	if (reducedSize <= maxReducedSize)
	    {
	    initialReduction = resScales[resTry];
	    initialReducedCount = resSizes[resTry];
	    break;
	    }
	}
    verbose(2, "initialReduction %d, initialReducedCount = %d\n", 
    	initialReduction, initialReducedCount);

    if (initialReduction > 0)
        {
	struct lm *lm = lmInit(0);
	int zoomIncrement = 4;
	lineFileRewind(lf);
	struct bbiSummary *rezoomedList = writeReducedOnceReturnReducedTwice(usageList, 
		lf, initialReduction, initialReducedCount,
		resIncrement, blockSize, itemsPerSlot, doCompress, lm, 
		f, &zoomDataOffsets[0], &zoomIndexOffsets[0], &totalSum);
	verboseTime(2, "writeReducedOnceReturnReducedTwice");
	zoomAmounts[0] = initialReduction;
	zoomLevels = 1;

	int zoomCount = initialReducedCount;
	int reduction = initialReduction * zoomIncrement;
	while (zoomLevels < bbiMaxZoomLevels)
	    {
	    int rezoomCount = slCount(rezoomedList);
	    if (rezoomCount >= zoomCount)
	        break;
	    zoomCount = rezoomCount;
	    zoomDataOffsets[zoomLevels] = ftell(f);
	    zoomIndexOffsets[zoomLevels] = bbiWriteSummaryAndIndex(rezoomedList, 
	    	blockSize, itemsPerSlot, doCompress, f);
	    zoomAmounts[zoomLevels] = reduction;
	    ++zoomLevels;
	    reduction *= zoomIncrement;
	    rezoomedList = bbiSummarySimpleReduce(rezoomedList, reduction, lm);
	    }
	lmCleanup(&lm);
	verboseTime(2, "further reductions");
	}

    }

/* Figure out buffer size needed for uncompression if need be. */
if (doCompress)
    {
    int maxZoomUncompSize = itemsPerSlot * sizeof(struct bbiSummaryOnDisk);
    uncompressBufSize = max(maxSectionSize, maxZoomUncompSize);
    }

/* Go back and rewrite header. */
rewind(f);
bits32 sig = bigWigSig;
bits16 version = bbiCurrentVersion;
bits16 summaryCount = zoomLevels;
bits16 reserved16 = 0;
bits32 reserved32 = 0;
bits64 reserved64 = 0;

/* Write fixed header */
writeOne(f, sig);
writeOne(f, version);
writeOne(f, summaryCount);
writeOne(f, chromTreeOffset);
writeOne(f, dataOffset);
writeOne(f, indexOffset);
writeOne(f, reserved16);	// fieldCount
writeOne(f, reserved16);	// definedFieldCount
writeOne(f, reserved64);	// autoSqlOffset
writeOne(f, totalSummaryOffset);
writeOne(f, uncompressBufSize);
for (i=0; i<2; ++i)
    writeOne(f, reserved32);
assert(ftell(f) == 64);

/* Write summary headers with data. */
verbose(2, "Writing %d levels of zoom\n", zoomLevels);
for (i=0; i<zoomLevels; ++i)
    {
    verbose(3, "zoomAmounts[%d] = %d\n", i, (int)zoomAmounts[i]);
    writeOne(f, zoomAmounts[i]);
    writeOne(f, reserved32);
    writeOne(f, zoomDataOffsets[i]);
    writeOne(f, zoomIndexOffsets[i]);
    }
/* Write rest of summary headers with no data. */
for (i=zoomLevels; i<bbiMaxZoomLevels; ++i)
    {
    writeOne(f, reserved32);
    writeOne(f, reserved32);
    writeOne(f, reserved64);
    writeOne(f, reserved64);
    }

/* Write total summary. */
fseek(f, totalSummaryOffset, SEEK_SET);
bbiSummaryElementWrite(f, &totalSum);

/* Write end signature. */
fseek(f, 0L, SEEK_END);
writeOne(f, sig);

lineFileClose(&lf);
carefulClose(&f);
}


/* definitions of functions */

static void outputCounts(unitSize *counts, char *chrom, unsigned size, FILE *f){
if (size == 0)
    errAbort("got 0 for size of chrom %s\n", chrom);

if (doOutBounds)
    {
    if (counts[0] < overMin)
	overMin = counts[0];
    if (counts[0] > overMax)
	overMax = counts[0];
    }

int ii;
int prevValue = counts[0];
int startPoint = 0;
for(ii=1; ii < size; ii++)
    {
    if (doOutBounds)
	{
	if (counts[ii] < overMin)
	    overMin = counts[ii];
	if (counts[ii] > overMax)
	    overMax = counts[ii];
	}
    if (counts[ii] != prevValue)
	{
	if (doZero || (prevValue != 0))
	    fprintf(f, "%s\t%u\t%u\t%u\n", chrom, startPoint, ii, prevValue);
	startPoint = ii;
	prevValue = counts[ii];
	}
    }

if (doZero || (prevValue != 0))
    fprintf(f, "%s\t%u\t%u\t%u\n", chrom, startPoint, ii, prevValue);
}

static void bedItemOverlapCount(struct hash *chromHash, char *infile, char *outfile){
unsigned maxChromSize = 0;
unitSize *counts = (unitSize *)NULL;
FILE *f = mustOpen(outfile, "w");
struct hashCookie hc = hashFirst(chromHash);
struct hashEl *hel;
while( (hel = hashNext(&hc)) != NULL) {
    unsigned num = (unsigned) ptToInt(hel->val);
    maxChromSize = max(num, maxChromSize);
}
verbose(2,"#\tmaxChromSize: %u\n", maxChromSize);
if (maxChromSize < 1)
    errAbort("maxChromSize is zero ?");

/*	Allocate just once for the largest chrom and reuse this array */
counts = needHugeMem(sizeof(unitSize) * maxChromSize);

/*	Reset the array to be zero to be reused */
memset((void *)counts, 0, sizeof(unitSize)*(size_t)maxChromSize);

unsigned chromSize = 0;
char *prevChrom = (char *)NULL;
boolean outputToDo = FALSE;
struct hash *seenHash = newHash(5);

    struct lineFile *bf = lineFileOpen(infile , TRUE);
    struct bed *bed = (struct bed *)NULL;
    char *row[12];
    int numFields = doBed12 ? 12 : 3;

    while (lineFileNextRow(bf,row, numFields))
	{
	int i;
	bed = bedLoadN(row, numFields);

	verbose(3,"#\t%s\t%d\t%d\n",bed->chrom,bed->chromStart, bed->chromEnd);

	if (prevChrom && differentWord(bed->chrom,prevChrom)) // End a chr
	    {
	    verbose(2,"#\tchrom %s done, size %d\n", prevChrom, chromSize);
	    if (outputToDo)
		outputCounts(counts, prevChrom, chromSize, f);
	    outputToDo = FALSE;
	    memset((void *)counts, 0,
		sizeof(unitSize)*(size_t)maxChromSize); /* zero counts */
	    freez(&prevChrom); 
	    // prevChrom is now NULL so it will be caught by next if!
	    }
	if ((char *)NULL == prevChrom)  // begin a chr
	    {
	    if (hashLookup(seenHash, bed->chrom))
		errAbort("ERROR:input file not sorted. %s seen before on line %d\n",
		    bed->chrom, bf->lineIx);

	    hashAdd(seenHash, bed->chrom, NULL);
	    prevChrom = cloneString(bed->chrom);
	    chromSize = hashIntVal(chromHash, prevChrom);
	    verbose(2,"#\tchrom %s starting, size %d\n", prevChrom,chromSize);
	    }
	if (bed->chromEnd > chromSize)
	    {
	    // check for circular chrM
	    if (doBed12 || bed->chromStart>=chromSize 
		|| differentWord(bed->chrom,"chrM")) 
		{
		warn("ERROR: %s\t%d\t%d", bed->chrom, bed->chromStart,
		bed->chromEnd);
		errAbort("chromEnd > chromSize ?  %d > %d", 
		    bed->chromEnd,chromSize);
		}

	    for (i = bed->chromStart; i < chromSize; ++i)
		INCWOVERFLOW(counts,i);
	    for (i = 0; i < (bed->chromEnd - chromSize); ++i)
		INCWOVERFLOW(counts,i);
	    }
	else if (doBed12)
	    {
	    int *starts = bed->chromStarts;
	    int *sizes = bed->blockSizes;
	    int *endStarts = &bed->chromStarts[bed->blockCount];

	    for(; starts < endStarts; starts++, sizes++)
		{
		unsigned int end = *starts + *sizes + bed->chromStart;
		for (i = *starts + bed->chromStart; i < end; ++i)
		    INCWOVERFLOW(counts,i);
		}
	    }
	else
	    {
	    for (i = bed->chromStart; i < bed->chromEnd; ++i)
		INCWOVERFLOW(counts, i);
	    }
	outputToDo = TRUE;
	bedFree(&bed); // plug the memory leak
	}

    lineFileClose(&bf);
    // Note, next file could be on same chr!

if (outputToDo)
    outputCounts(counts, prevChrom, chromSize, f);

if (doOutBounds)
    fprintf(stderr, "min %lu max %lu\n", (unsigned long)overMin, (unsigned long)overMax);

verbose(2,"#\tchrom %s done, size %d\n", prevChrom, chromSize);
carefulClose(&f);
freeMem(counts);
freez(&prevChrom);
// hashFreeWithVals(&chromHash, freez);
freeHash(&seenHash);
}

int density_usage(){
    fprintf(stderr, "\n");
    fprintf(stderr, "\nWorks as ChipSeqHandler program.\nPlease noticed that if reads mapped to the chromosomes which didn't existed in size file, this type of reads will be discarded.\n");
    fprintf(stderr, "Usage:   iteres density [options] <chromosome size file> <bam/sam alignment file>\n\n");
    fprintf(stderr, "Options: -S       input is SAM [off]\n");
    fprintf(stderr, "         -Q       unique reads mapping Quality threshold [10]\n");
    fprintf(stderr, "         -r       do NOT remove redundant reads [off]\n");
    fprintf(stderr, "         -T       treat 1 paired-end read as 2 single-end reads [off]\n");
    fprintf(stderr, "         -D       discard if only one end mapped in a paired end reads [off]\n");
    fprintf(stderr, "         -C       Add 'chr' string as prefix of reference sequence [off]\n");
    fprintf(stderr, "         -E       extend reads to represent fragment [150], specify 0 if want no extension\n");
    fprintf(stderr, "         -I       Insert length threshold [500]\n");
    fprintf(stderr, "         -o       output prefix [basename of input without extension]\n");
    fprintf(stderr, "         -h       help message\n");
    fprintf(stderr, "         -?       help message\n");
    fprintf(stderr, "\n");
    return 1;
}

/* main function */
int main_density (int argc, char *argv[]) {
    
    char *output, *outReportfile, *outExtfile, *outbedGraphfile, *outbigWigfile;
    unsigned long long int *cnt;
    int optSam = 0, c, optDup = 1, optaddChr = 0, optDis = 0, optTreat = 0;
    unsigned int optQual = 10, optExt = 150, optisize = 500;
    char *optoutput = NULL;
    time_t start_time, end_time;
    start_time = time(NULL);
    
    while ((c = getopt(argc, argv, "SQ:rTDCo:E:I:h?")) >= 0) {
        switch (c) {
            case 'S': optSam = 1; break;
            case 'Q': optQual = (unsigned int)strtol(optarg, 0, 0); break;
            case 'r': optDup = 0; break;
            case 'T': optTreat = 1; break;
            case 'D': optDis = 1; break;
            case 'C': optaddChr = 1; break;
            case 'E': optExt = (unsigned int)strtol(optarg, 0, 0); break;
            case 'I': optisize = (unsigned int)strtol(optarg, 0, 0); break;
            case 'o': optoutput = strdup(optarg); break;
            case 'h':
            case '?': return density_usage(); break;
            default: return 1;
        }
    }
    if (optind + 2 > argc)
        return density_usage();

    char *chr_size_file = argv[optind];
    char *sam_file = argv[optind+1];
    
    if(optoutput) {
        output = optoutput;
    } else {
        output = cloneString(get_filename_without_ext(basename(sam_file)));
    }
    

    if(asprintf(&outExtfile, "%s.extended.bed", output) < 0)
        errAbort("Mem Error.\n");
    if(asprintf(&outbedGraphfile, "%s.extended.bedGraph", output) < 0)
        errAbort("Mem Error.\n");
    if(asprintf(&outbigWigfile, "%s.bigWig", output) < 0)
        errAbort("Mem Error.\n");
    if (asprintf(&outReportfile, "%s.report", output) < 0)
        errAbort("Preparing output wrong");
    
    struct hash *hash = hashNameIntFile(chr_size_file);
    
    //sam file to bed file
    fprintf(stderr, "* Start to parse the SAM/BAM file ...\n");
    cnt = sam2bed(sam_file, outExtfile, hash, optSam, optQual, optDup, optaddChr, optDis, optisize, optExt, optTreat);
    //sort
    //fprintf(stderr, "\n* Sorting\n");
    //bedSortFile(outBedfile, outBedfile);

    //remove dup
    //fprintf(stderr, "* Removing duplication\n");
    //uniqueBed = removeBedDup(outBedfile, outFilterfile);

    //extend and write extend bed
    //fprintf(stderr, "* Extending to %d and writing extended bed\n", arguments.extend);
    //int extendWarn = extendBed(hash, arguments.extend, outFilterfile, outExtfile);
    //if (extendWarn == 1)
    //    outExtfile = cloneString(outFilterfile);
    //
    //if (extendWarn != 1){
        //sort extend bed
    //    fprintf(stderr, "* Sorting extended bed\n");
    //    bedSortFile(outExtfile, outExtfile);
    //}

    //sort extend bed
    fprintf(stderr, "\n* Sorting extended bed\n");
    sortBedfile(outExtfile);
    
    //bedItemOverlap step
    fprintf(stderr, "* Generating bedGraph\n");
    bedItemOverlapCount(hash, outExtfile, outbedGraphfile);

    //generate bigWig
    fprintf(stderr, "* Generating bigWig\n");
    //bigWigFileCreate(outbedGraphfile, chr_size_file, 256, 1024, 0, 1, outbigWigfile);
    bedGraphToBigWig(outbedGraphfile, chr_size_file, outbigWigfile);

    //write report file
    fprintf(stderr, "* Preparing report file\n");
    writeReportDensity(outReportfile, cnt, optQual);
    
    //cleaning
    hashFree(&hash);
    free(outExtfile);
    free(outbedGraphfile);
    free(outbigWigfile);
    free(outReportfile);
    end_time = time(NULL);
    fprintf(stderr, "* Done, time used %.0f seconds.\n", difftime(end_time, start_time));
    return 0;
}
