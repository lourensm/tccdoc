
/*
 * Current status unclear: I started with a copy of tiffinfo.c:
 * Dollar Id: tiffinfo.c,v 1.22 2013-07-10 00:44:22 fwarmerdam Exp Dollar
 * Maintained by lourensm
 *
 * Copyright (c) 1988-1997 Sam Leffler
 * Copyright (c) 1991-1997 Silicon Graphics, Inc.
 *
 * Permission to use, copy, modify, distribute, and sell this software and 
 * its documentation for any purpose is hereby granted without fee, provided
 * that (i) the above copyright notices and this permission notice appear in
 * all copies of the software and related documentation, and (ii) the names of
 * Sam Leffler and Silicon Graphics may not be used in any advertising or
 * publicity relating to the software without the specific, prior written
 * permission of Sam Leffler and Silicon Graphics.
 * 
 * THE SOFTWARE IS PROVIDED "AS-IS" AND WITHOUT WARRANTY OF ANY KIND, 
 * EXPRESS, IMPLIED OR OTHERWISE, INCLUDING WITHOUT LIMITATION, ANY 
 * WARRANTY OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.  
 * 
 * IN NO EVENT SHALL SAM LEFFLER OR SILICON GRAPHICS BE LIABLE FOR
 * ANY SPECIAL, INCIDENTAL, INDIRECT OR CONSEQUENTIAL DAMAGES OF ANY KIND,
 * OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS,
 * WHETHER OR NOT ADVISED OF THE POSSIBILITY OF DAMAGE, AND ON ANY THEORY OF 
 * LIABILITY, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE 
 * OF THIS SOFTWARE.
 */

/* TODO: consistent structuring of segments and boxes: either x1, x2 or x,width  */
/* TODO: remove copied functionality from tiffinfo.c */
#include "tif_config.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h>

#ifdef HAVE_STRINGS_H
# include <strings.h>
#endif

#ifdef HAVE_UNISTD_H
# include <unistd.h>
#endif

#ifdef NEED_LIBPORT
# include "libport.h"
#endif

#include "tiffiop.h"

#if 0
#define LMIN(a,b) (((a)<(b))?(a):(b))
#define LMAX(a,b) (((a)>(b))?(a):(b))
#endif

#define LMALLOC(x) malloc(x);
#define LFREE(x) free(x)
#define ERROR_EXIT(MSG) error_exit(MSG, __FILE__, __LINE__)
static void error_exit(const char* msg,const char* file, const int line) {
	fprintf(stderr, "***%s:%d  : FATAL ERROR: \"%s\"\n",
		file, line, msg);
	exit(EXIT_FAILURE);
}

static TIFFErrorHandler old_error_handler = NULL;
static int status = 0;                  /* exit status */
static int showdata = 0;		/* show data */
static int doanalyse = 0;		/* start recognizing blobs */
static int test3 = 0;		/* test binary encodings */

static int rawdata = 0;			/* show raw/decoded data */
static int showwords = 0;		/* show data as bytes/words */
static int readdata = 0;		/* read data in file */
static int stoponerr = 1;		/* stop on first read error */
static char* resultdir = "results";
static	void usage(void);
static	void tiffinfo(const char* filename, TIFF*, uint16, long, int);

static void
PrivateErrorHandler(const char* module, const char* fmt, va_list ap)
{
        if (old_error_handler)
                (*old_error_handler)(module,fmt,ap);
	status = 1;
}



int
main(int argc, char* argv[])
{
	int dirnum = -1, multiplefiles, c;
	uint16 order = 0;
	TIFF* tif;
	extern int optind;
	extern char* optarg;
	long flags = 0;
	uint64 diroff = 0;
	int chopstrips = 0;		/* disable strip chopping */

	while ((c = getopt(argc, argv, "f:o:cdDSjilmrstuvwz0123456789")) != -1)
		switch (c) {
		case '0': case '1': case '2': case '3':
		case '4': case '5': case '6': case '7':
		case '8': case '9':
			dirnum = atoi(&argv[optind-1][1]);
			break;
		case 'd':
			showdata++;
			/* fall thru... */
		case 'D':
			readdata++;
			break;
		case 'c':
			flags |= TIFFPRINT_COLORMAP | TIFFPRINT_CURVES;
			break;
		case 'f':		/* fill order */
			if (streq(optarg, "lsb2msb"))
				order = FILLORDER_LSB2MSB;
			else if (streq(optarg, "msb2lsb"))
				order = FILLORDER_MSB2LSB;
			else
				usage();
			break;
		case 'i':
			stoponerr = 0;
			break;
		case 'o':
			diroff = strtoul(optarg, NULL, 0);
			break;
		case 'j':
			flags |= TIFFPRINT_JPEGQTABLES |
				 TIFFPRINT_JPEGACTABLES |
				 TIFFPRINT_JPEGDCTABLES;
			break;
		case 'r':
			rawdata = 1;
			break;
		case 's':
			flags |= TIFFPRINT_STRIPS;
			break;
		case 't':
			doanalyse = 1;
			readdata++;
			break;
		case 'u':
			doanalyse = 1;
			test3 = 1;
			readdata++;
			break;
		case 'w':
			showwords = 1;
			break;
		case 'z':
			chopstrips = 1;
			break;
		case '?':
			usage();
			/*NOTREACHED*/
		}
	if (optind >= argc)
		usage();

	old_error_handler = TIFFSetErrorHandler(PrivateErrorHandler);

	multiplefiles = (argc - optind > 1);
	for (; optind < argc; optind++) {
		const char* filename = argv[optind];
		if (multiplefiles)
			printf("%s:\n", filename);
		tif = TIFFOpen(filename, chopstrips ? "rC" : "rc");
		if (tif != NULL) {
			if (dirnum != -1) {
				if (TIFFSetDirectory(tif, (tdir_t) dirnum))
					tiffinfo(filename, tif, order, flags, 1);
			} else if (diroff != 0) {
				if (TIFFSetSubDirectory(tif, diroff))
					tiffinfo(filename,tif, order, flags, 1);
			} else {
				do {
					toff_t offset=0;

					tiffinfo(filename, tif, order, flags, 1);
					if (TIFFGetField(tif, TIFFTAG_EXIFIFD,
							 &offset)) {
						if (TIFFReadEXIFDirectory(tif, offset)) {
							tiffinfo(filename, tif, order, flags, 0);
						}
					}
				} while (TIFFReadDirectory(tif));
			}
			TIFFClose(tif);
		}
	}
	return (status);
}

char* stuff[] = {
"usage: tiffinfo [options] input...",
"where options are:",
" -D		read data",
" -i		ignore read errors",
" -c		display data for grey/color response curve or colormap",
" -d		display raw/decoded image data",
" -f lsb2msb	force lsb-to-msb FillOrder for input",
" -f msb2lsb	force msb-to-lsb FillOrder for input",
" -j		show JPEG tables",
" -o offset	set initial directory offset",
" -r		read/display raw image data instead of decoded data",
" -s		display strip offsets and byte counts",
" -w		display raw data in words rather than bytes",
" -z		enable strip chopping",
" -#		set initial directory (first directory is # 0)",
NULL
};

static void
usage(void)
{
	char buf[BUFSIZ];
	int i;

	setbuf(stderr, buf);
        fprintf(stderr, "%s\n\n", TIFFGetVersion());
	for (i = 0; stuff[i] != NULL; i++)
		fprintf(stderr, "%s\n", stuff[i]);
	exit(-1);
}

static void
ShowStrip(tstrip_t strip, unsigned char* pp, uint32 nrow, tsize_t scanline)
{
	register tsize_t cc;

	printf("Strip %lu:\n", (unsigned long) strip);
	while (nrow-- > 0) {
		for (cc = 0; cc < scanline; cc++) {
			printf(" %02x", *pp++);
			if (((cc+1) % 24) == 0)
				putchar('\n');
		}
		putchar('\n');
	}
}

void
TIFFReadContigStripData(TIFF* tif)
{
	unsigned char *buf;
	tsize_t scanline = TIFFScanlineSize(tif);

	buf = (unsigned char *)_TIFFmalloc(TIFFStripSize(tif));
	if (buf) {
		uint32 row, h=0;
		uint32 rowsperstrip = (uint32)-1;

		TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h);
		TIFFGetField(tif, TIFFTAG_ROWSPERSTRIP, &rowsperstrip);
		for (row = 0; row < h; row += rowsperstrip) {
			uint32 nrow = (row+rowsperstrip > h ?
			    h-row : rowsperstrip);
			tstrip_t strip = TIFFComputeStrip(tif, row, 0);
			if (TIFFReadEncodedStrip(tif, strip, buf, nrow*scanline) < 0) {
				if (stoponerr)
					break;
			} else if (showdata)
				ShowStrip(strip, buf, nrow, scanline);
		}
		_TIFFfree(buf);
	}
}

void
TIFFReadSeparateStripData(TIFF* tif)
{
	unsigned char *buf;
	tsize_t scanline = TIFFScanlineSize(tif);

	buf = (unsigned char *)_TIFFmalloc(TIFFStripSize(tif));
	if (buf) {
		uint32 row, h=0;
		uint32 rowsperstrip = (uint32)-1;
		tsample_t s, samplesperpixel=0;

		TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h);
		TIFFGetField(tif, TIFFTAG_ROWSPERSTRIP, &rowsperstrip);
		TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL, &samplesperpixel);
		for (row = 0; row < h; row += rowsperstrip) {
			for (s = 0; s < samplesperpixel; s++) {
				uint32 nrow = (row+rowsperstrip > h ?
				    h-row : rowsperstrip);
				tstrip_t strip = TIFFComputeStrip(tif, row, s);
				if (TIFFReadEncodedStrip(tif, strip, buf, nrow*scanline) < 0) {
					if (stoponerr)
						break;
				} else if (showdata)
					ShowStrip(strip, buf, nrow, scanline);
			}
		}
		_TIFFfree(buf);
	}
}

static void
ShowTile(uint32 row, uint32 col, tsample_t sample,
    unsigned char* pp, uint32 nrow, tsize_t rowsize)
{
	uint32 cc;

	printf("Tile (%lu,%lu", (unsigned long) row, (unsigned long) col);
	if (sample != (tsample_t) -1)
		printf(",%u", sample);
	printf("):\n");
	while (nrow-- > 0) {
	  for (cc = 0; cc < (uint32) rowsize; cc++) {
			printf(" %02x", *pp++);
			if (((cc+1) % 24) == 0)
				putchar('\n');
		}
		putchar('\n');
	}
}

void
TIFFReadContigTileData(TIFF* tif)
{
	unsigned char *buf;
	tsize_t rowsize = TIFFTileRowSize(tif);

	buf = (unsigned char *)_TIFFmalloc(TIFFTileSize(tif));
	if (buf) {
		uint32 tw=0, th=0, w=0, h=0;
		uint32 row, col;

		TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &w);
		TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h);
		TIFFGetField(tif, TIFFTAG_TILEWIDTH, &tw);
		TIFFGetField(tif, TIFFTAG_TILELENGTH, &th);
		for (row = 0; row < h; row += th) {
			for (col = 0; col < w; col += tw) {
				if (TIFFReadTile(tif, buf, col, row, 0, 0) < 0) {
					if (stoponerr)
						break;
				} else if (showdata)
					ShowTile(row, col, (tsample_t) -1, buf, th, rowsize);
			}
		}
		_TIFFfree(buf);
	}
}

void
TIFFReadSeparateTileData(TIFF* tif)
{
	unsigned char *buf;
	tsize_t rowsize = TIFFTileRowSize(tif);

	buf = (unsigned char *)_TIFFmalloc(TIFFTileSize(tif));
	if (buf) {
		uint32 tw=0, th=0, w=0, h=0;
		uint32 row, col;
		tsample_t s, samplesperpixel=0;

		TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &w);
		TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h);
		TIFFGetField(tif, TIFFTAG_TILEWIDTH, &tw);
		TIFFGetField(tif, TIFFTAG_TILELENGTH, &th);
		TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL, &samplesperpixel);
		for (row = 0; row < h; row += th) {
			for (col = 0; col < w; col += tw) {
				for (s = 0; s < samplesperpixel; s++) {
					if (TIFFReadTile(tif, buf, col, row, 0, s) < 0) {
						if (stoponerr)
							break;
					} else if (showdata)
						ShowTile(row, col, s, buf, th, rowsize);
				}
			}
		}
		_TIFFfree(buf);
	}
}

void
TIFFReadData(TIFF* tif)
{
	uint16 config = PLANARCONFIG_CONTIG;

	TIFFGetField(tif, TIFFTAG_PLANARCONFIG, &config);
	if (TIFFIsTiled(tif)) {
		if (config == PLANARCONFIG_CONTIG)
			TIFFReadContigTileData(tif);
		else
			TIFFReadSeparateTileData(tif);
	} else {
		if (config == PLANARCONFIG_CONTIG)
			TIFFReadContigStripData(tif);
		else
			TIFFReadSeparateStripData(tif);
	}
}

static void
ShowRawBytes(unsigned char* pp, uint32 n)
{
	uint32 i;

	for (i = 0; i < n; i++) {
		printf(" %02x", *pp++);
		if (((i+1) % 24) == 0)
			printf("\n ");
	}
	putchar('\n');
}

static void
ShowRawWords(uint16* pp, uint32 n)
{
	uint32 i;

	for (i = 0; i < n; i++) {
		printf(" %04x", *pp++);
		if (((i+1) % 15) == 0)
			printf("\n ");
	}
	putchar('\n');
}

void
TIFFReadRawData(TIFF* tif, int bitrev)
{
	tstrip_t nstrips = TIFFNumberOfStrips(tif);
	const char* what = TIFFIsTiled(tif) ? "Tile" : "Strip";
	uint64* stripbc=NULL;

	TIFFGetField(tif, TIFFTAG_STRIPBYTECOUNTS, &stripbc);
	if (nstrips > 0) {
		uint32 bufsize = (uint32) stripbc[0];
		tdata_t buf = _TIFFmalloc(bufsize);
		tstrip_t s;

		for (s = 0; s < nstrips; s++) {
			if (stripbc[s] > bufsize) {
				buf = _TIFFrealloc(buf, (tmsize_t)stripbc[s]);
				bufsize = (uint32) stripbc[s];
			}
			if (buf == NULL) {
				fprintf(stderr,
				   "Cannot allocate buffer to read strip %lu\n",
				    (unsigned long) s);
				break;
			}
			if (TIFFReadRawStrip(tif, s, buf, (tmsize_t) stripbc[s]) < 0) {
				fprintf(stderr, "Error reading strip %lu\n",
				    (unsigned long) s);
				if (stoponerr)
					break;
			} else if (showdata) {
				if (bitrev) {
					TIFFReverseBits(buf, (tmsize_t)stripbc[s]);
					printf("%s %lu: (bit reversed)\n ",
					    what, (unsigned long) s);
				} else
					printf("%s %lu:\n ", what,
					    (unsigned long) s);
				if (showwords)
					ShowRawWords((uint16*) buf, (uint32) stripbc[s]>>1);
				else
					ShowRawBytes((unsigned char*) buf, (uint32) stripbc[s]);
			}
		}
		if (buf != NULL)
			_TIFFfree(buf);
	}
}


typedef struct lsegment {
	int min_x;
	int length;
} LSegment;

typedef struct lcell {
	struct lcell *next;
	/*struct lcell *founder;*/ /* top left touching segment */
	LSegment segment;
	struct lobject* top_object;
} LCell;



typedef struct lbox {
	int x1, y1,x2, y2;
} LBox;

typedef struct lscanline {
	int y;
	struct lcell *first;
} LScanline;

typedef struct lblob {
	struct lobject* objects;
	LBox range;
	int id;
} LBlob;


typedef struct lobject {
	int id;
	struct lblob *origin;
	int x;
	int y;
	struct lobject *blob_next; 
} LObject;


/*
TODO for blob problem detection:
maintain array of object pointers
check when freeing LBlob whether an origin is left with reference to freed.

invariant:
all objects that have a common LBlob should be part of its objects blob_next list.
 */

typedef struct lblobinfo {
	LCell *last_cell_previous_line;
	LCell *last_cell_current_line;
	LObject* last_object;
	LScanline* llines;
} LBlobinfo;


/* areas 
sorted: top first, left first.
Awkward: detect in area seems complex, but one setup per scanline:

Need for scanline: active_areas, sorted wrt startx:
Simplest: look at all entries in active_areas per segment
Next: remove heading entries that are out of range : pixhor >= tophor
*/


typedef struct larea {
	struct larea *next;
	/*double y1f, tophorf, widthf, heightf;*/
	LBox box;
	void* data;
	void (*handlesegment)(LSegment segment, int y, void* data);

	void (*handlepixel)(int x, int y, void* data);
	void (*analyse)(void* data);
} LArea;

typedef struct lxyvalues {
	int z1, z2;
	int lastpospos;
	int* zvalues;
	double z0;
} LXYValues;

typedef struct lpixelrotationdata {
	LBox box, maxbox;
	double max_atan_angle;
	int lastpospos;
	int* xvalues;
	int* yvalues;
	double x0, y0;
	const char * filespec;
} LPixelrotationdata;
/*
Start from the middle of the range.
Add pixel value to array y + offset
n = 0: offset = 0 for all x
n = 1: offset = -1 for 


00000000000000000000000000000000
-1-1-1                1 1 1 1 1 

WRONG: is two dimensional array,
(2*(rdata->lastpospos)+1) * (rdata->box.y2 - rdata->box.y1)

 */
static void handlepixelrotationxy(int z, int z1, int z2, int  oz, double oz0, 
			   int* xyvalues,
			   int lastpospos, double max_atan_angle) {
	int n;
	for (n = -lastpospos;n <= lastpospos;n++) {
		int zz0 = z +
			(int)floor(0.5+(oz - oz0)*n*max_atan_angle/lastpospos);
		int z0d = zz0 - z1;
		if (z0d >=0 && z0d < z2 - z1 + 1) {
			int index = n+lastpospos+ (2*lastpospos+1)*z0d;
			assert(index >= 0);
			assert(index < (z2 -z1 +1)*(2*lastpospos+1));
			xyvalues[n+lastpospos+ (2*lastpospos+1)*z0d]++;
		}
	}
}

static void updatemaxbox(int x, int y, LBox* maxbox) {
	if (x < maxbox->x1) {
		maxbox->x1 = x;
	}
	if (x > maxbox->x2) {
		maxbox->x2 = x;
	}
	if (y < maxbox->y1) {
		maxbox->y1 = y;
	}
	if (y > maxbox->y2) {
		maxbox->y2 = y;
	}
}
static void handlepixelrotation(int x, int y, void* data) {
	LPixelrotationdata* rdata = (LPixelrotationdata*)data;
	int lastpospos = rdata->lastpospos;

	updatemaxbox(x, y, &rdata->maxbox);

	
	
	handlepixelrotationxy(y, (rdata->box).y1, (rdata->box).y2, x, rdata->x0, 
			      rdata->yvalues, lastpospos, rdata->max_atan_angle);
	handlepixelrotationxy(x, (rdata->box).x1, (rdata->box).x2, y, rdata->y0, 
			      rdata->xvalues, lastpospos, rdata->max_atan_angle);
}

static int* setupxyvalues(int z1, int z2, int lastpospos) {
	int *values;
	int yvaluesize = (2*lastpospos+1)*(z2 - z1+1);
	int i;
	assert(z2 >= z1);
	values = (int*)LMALLOC(yvaluesize* sizeof(int));
	for (i = 0; i<yvaluesize;i++) {
		values[i] = 0;
	}
	return values;
}

static LPixelrotationdata* setuppixelrotationdata(const char* filespec, LBox box, double max_atan_angle, int lastpospos) {
	LPixelrotationdata* data = LMALLOC(sizeof (LPixelrotationdata));
	data->yvalues = setupxyvalues(box.y1, box.y2, lastpospos);
	data->xvalues = setupxyvalues(box.x1, box.x2, lastpospos);
	data->box = box;
	data->maxbox.x1= box.x2;
	data->maxbox.x2 = box.x1;
	data->maxbox.y1 = box.y2;
	data->maxbox.y2 = box.y1;
	data->max_atan_angle = max_atan_angle;
	data->lastpospos = lastpospos;
	data->x0 = (box.x2 + box.x1)/2;
	data->y0 = (box.y2 + box.y1)/2;
	data->filespec = filespec;
				      
	return data;
}

struct lcontpt {
	int pos;
	double max;
} LContpt;


static int max_cmp(const void* a, const void* b) {
	const struct lcontpt *a1 = (struct lcontpt *)a;
	const struct lcontpt *b1 = (struct lcontpt *)b;
	return a1->pos - b1->pos;
}

static void analysepixelrotationxy(const char* xory, const char* filespec,
			    int lastpospos,
			    int z1, int z2, int* xyvalues) {
	const int sbufsize = 80;
	char filename[sbufsize];
	int possize;
	double* max2;
	double* max4;
	int z, i1;
	struct lcontpt  maxes[5];
	FILE *f;
	assert(strlen(filespec)+strlen(xory)+strlen("values")+
						 strlen("maxsums")+
						 strlen("__.dat") < sbufsize);
	possize = 2*lastpospos+1;
	max2 = (double*)LMALLOC(possize*sizeof(double));
	max4 = (double*)LMALLOC(possize*sizeof(double));
	{
		int i;
		for (i=0;i < possize;i++) {
			max2[i] = 0;
			max4[i] = 0;
		}
	}
	/* 5 largest yvalues of max2 */

	{
		int i;
		for (i = 0;i < 5;i++) { maxes[i].max = 0;maxes[i].pos = 0; }
	}
	snprintf(filename,sizeof(filename),"%s/%s_%s_values.dat", resultdir,
		 filespec, xory);
	f = fopen(filename, "w");	
	for (z = 0; z < z2 - z1 +1; z++) {
		fprintf(f, "%d", z);
		for (i1 = 0;i1 < possize; i1++) {
			int count = xyvalues[i1 + possize*z];
			double c2;
			fprintf(f, "\t%d", count);
			c2 = count*count;
			max2[i1] += c2;
			max4[i1] += c2*c2;
		}
		fprintf(f, "\n");
	}
	fclose(f);
	printf("Wrote datafile %s\n", filename);
	sprintf(filename, "%s/%s_%s_maxsums.dat", resultdir, filespec, xory);
	f = fopen(filename, "w");
	{
		int i, j;
		for (i = 0;i < possize; i++) {
			fprintf(f, "%d\t%f\t%f\n", i, sqrt(max2[i]), sqrt(sqrt(max4[i])));
			for (j=0;j<5;j++) {
				if (max2[i] > maxes[j].max) {
					int k;
					for (k =4; k > j; k--) {
						maxes[k] = maxes[k-1];
					}
					maxes[j].max = max2[i];
					maxes[j].pos = i;
					break;
				}
			}
		}
	}
	LFREE(max2);
	LFREE(max4);
	fclose(f);
	{
		int mpos = maxes[0].pos;
		int posindex;
		double m1,m2,m3,p1,p2,p3,p0;
		int pl,ph,sl,sh, yl, yh;
		printf("Analysis of %s direction\n", xory);
		qsort(maxes, 5, sizeof(struct lcontpt), max_cmp);
		{
			int i;
			for (i = 0;i<5;i++) {
				if (maxes[i].pos == mpos) {printf("*");posindex = i;} else {printf(" ");}
				printf("i:%d, p:%d, m:%f\n", i, maxes[i].pos, maxes[i].max);
			}
		}
		if (posindex<=0||posindex>=4) {
			fprintf(stderr,"Probably no text in right page\n");
			return;
		}
		assert(posindex>0&&posindex<4);
		m1 = maxes[posindex-1].max;
		m2 = maxes[posindex].max;
		m3 = maxes[posindex+1].max;
		p1 = maxes[posindex-1].pos;
		p2 = maxes[posindex].pos;
		p3 = maxes[posindex+1].pos;
		p0 = (  (m2-m3)*(p1*p1 - p2*p2) - (m1 - m2)*(p2*p2 - p3*p3))
			/
			( (p2 - p3)*(m1 -m2) - (p1 - p2)*(m2-m3))/(-2.0);
		pl = floor(p0);
		ph = pl + 1;
		sl = 0;sh = 0;
		yl = 0;yh = 0;
		/* require at least 12 consecutive nonzeros 
		   if as much consecutive 0s after 
		*/
		{
			int yy;
			for (yy = 0; yy <  z2 - z1 +1; yy++) {
				int countl = xyvalues[pl + possize*yy];
				int counth = xyvalues[ph + possize*yy];
				sl += countl;
				sh += counth;
				if (countl == 0 && sl > 0) sl--;
				if (counth == 0 && sh > 0) sh--;
				if (sl >20) yl = yy;
				if (sh >20) yh = yy;
				printf("yy:%d sl:%d cl:%d\n", yy, sl, countl);
				if (sl>20 && sh >20) break;
			}
		}
		{
			int ylow = yl + (yh - yl)*(p0 - pl)/(ph-pl);
			sl = 0;
			sh = 0;
			{
				int yy;
				for (yy = z2 - z1 ; yy >= 0;yy--) {
					int countl = xyvalues[pl + possize*yy];
					int counth = xyvalues[ph + possize*yy];
					sl += countl;
					sh += counth;
					if (sl >10) yl = yy;
					if (sh >10) yh = yy;
					if (sl>10 && sh >10) break;
				}
			}
			{
				int yhi = yl + (yh - yl)*(p0 - pl)/(ph-pl);
				printf("zl:%d, zh:%d\n", yl, yh);
				printf("maxp:%f, zlow:%d, zhi:%d\n", p0, ylow, yhi);
			}
		}
	}
	/* TODO:
           detect digits of pageno, put font box around them and that defines the
	   x position, and absolute window in case of 13tr, if we know maximum dimensions 
	   of text.
           Requires detection of blobs. and indexing, recognition.
	 */
}

static void   analysepixelrotation(void* data) {
	LPixelrotationdata* rdata = (LPixelrotationdata*)data;
	analysepixelrotationxy("y", rdata->filespec,
			       rdata->lastpospos,
			       (rdata->box).y1, (rdata->box).y2, rdata->yvalues);
	analysepixelrotationxy("x", rdata->filespec,
			       rdata->lastpospos,
			       (rdata->box).x1, (rdata->box).x2, rdata->xvalues);
}

static void setuppixelrotation(LArea* larea, const char* filespec, LBox box,
			double max_atan_angle, int lastpospos) {
	larea->next = NULL;
	larea->box = box;
	larea->data = (void*)setuppixelrotationdata(filespec, box, max_atan_angle, lastpospos);
	larea->handlepixel = &handlepixelrotation;
	larea->analyse = &analysepixelrotation;
	larea->handlesegment = NULL;
}

		
typedef struct lactive_areas {
	struct lactive_areas *next;
	LArea area;
} LActive_areas;
/*
need to partition space:
as rectangular (contiguous? no) disjoint pieces
1)ignore
2)warn if nonwhite occurs
3)analyse
     expect pageno
     determine smallest (rotated < atan(X)) rectangle in which black occurs
     determine rotation:
          based on contour only
	  based on character alignment
	  (rotated) average blackness:
	   should be
     0
     0
       10%
       80%
       10%
     0
1)ignore anything above, below 80%, there should be steep ..

least squares? doesnt seem possible.
a)follow multiple angles:
  atan(.1) as maximum
  atan(.001) as minimum, no, 1% of normal line distance
  and then 100 linearly? 
     analysis result:
          1)are there different fonts detectable?
          2)accuracy of rotation detection

Rotation detection:
1)add to 0 entry
2)divide textarea into 100 pieces,
3)add first 1/100th to f(0,i), to f(-1,i-1), f(-2, i -2) etc..
no, must be fractional...
f(0, i), f(1, i), f(9, i), f(10, i+1),   f(11, i+1),...
second 1/100th = n = 1  f(K, I)
f(0, i),.. f(4,i), f(5, i+1), .. 
I = floor(n*K/L).
split into 1000 and maintain 100, K = -100, 100, L = 1000
8000/2 pixels, 4000 pixels...
K = -100, 100
if (bitset)
for (k in -100, 100) {
     i = floor(n*k/deltan)

zdd deltan = 9909*123/332 dan
lijkt niet te kloppen, leftside, rightside verschillende vullingen
doe n - n0? met n0 in midden van de tekst?
   
}

 */
static void printlbox(const char*txt, LBox box) {
	printf("%s: x1:%d, y1:%d, x2:%d, y2:%d\n",
	       txt, box.x1, box.y1, box.x2, box.y2);
}


static const char* filespec(const char* filename) {
	static const char* called = NULL;
	static char value[80];
	assert(filename != NULL);
	assert(strlen(filename)<80);
	if (called == filename) return value;
	{
		size_t i;
		for (i=0; i<=strlen(filename); i++) {
			value[i] = filename[i];
		}
	}
	assert(called == NULL); /* allow one single call */
	called = filename;
	{
		char* basename1 = basename(value);
		const char* ext = ".tif";
		size_t nl = strlen(basename1), el = strlen(ext);
		if (nl < el || strcmp(basename1 + nl - el, ext)){
			fprintf(stderr, "ERROR:expected .tif extension \"%s\"\n", filename);
			exit(1);
		}
		{
			size_t i;
			for (i = 0;i < nl - el; i++) {
				value[i] = basename1[i];
			}
		}
		value[nl - el] = '\0';
	}
	return value;
}


static void handle_active_areas_per_pixel(const int currentx, const int currenty,
					  LActive_areas** active_areas_left) {
	
	while (*active_areas_left != NULL
	       && (((*active_areas_left)->area).box. y2< currenty)) {
		*active_areas_left = (*active_areas_left)-> next;
					     
	}
	if (*active_areas_left != NULL&&currenty>= ((*active_areas_left)->area).box.y1) {
		LActive_areas* active_areas_now = (*active_areas_left);
		LArea* area_left =  &(*active_areas_left)->area;
		while (active_areas_now != NULL&& currenty<=area_left->box.y2) {
			if (currentx>=area_left->box.x1 && currentx<=area_left->box.x2) {
				area_left->handlepixel(currentx, currenty, area_left->data);
			}
			active_areas_now = active_areas_now->next;
		}	
	}
}

static LActive_areas* setupactionareas(const char* filename, int width, int height) {
	/* two windows for now, detecting text and rotation */
	/* first 1 window.. */
	/*LHS: expected:Hor 23mm, 146mm vert:20mm 202mm
          RHS: 183mm, 306mm
	  (all HOR 332mm, vert:235)
	  allow for rotation of box..first allow around center only
	  latan (x - latan*(235)
	  y - latan
	 */
	const char* filespec1 = filespec(filename);
	double latan = 0.1;
	double hor1 = 23, hor2 = 146, /*hor3 = 183, hor4 = 306,*/ horall = 332;
	double vert1 = 20, vert2 = 202, vertall = 235;
	/*double vertpageno = 213, vertpagenoh = 3;
	  double pagel = 6;*//* 3 digits */
	double scaleh = width/horall;
	double scalev = height/vertall;
	LBox boxl = {(int)round(scaleh*(hor1-latan*( (vert2-vert1)/2))),
		     (int)round(scalev*(vert1 - latan*( (hor2-hor1)/2))),
		     (int)round(scaleh*(hor2+latan*( (vert2-vert1)/2))),
		     (int)round(scalev*(vert2 + latan*( (hor2-hor1)/2)))};
	printlbox("Leftwindow", boxl);
	{
		LActive_areas* r2 = LMALLOC(sizeof (LActive_areas));
		r2->next = NULL;
		setuppixelrotation(&(r2->area),filespec1, boxl, latan, 100);
		return r2;
	}
}


static int lblobcount = 0;
static int nextblobid = 0;
static int nextobject = 0;




static LCell* new_lcell(int startx, int afterx) {
	LCell* res1 = (LCell*)LMALLOC(sizeof(LCell));
	res1->next = NULL;
	(res1->segment).min_x = startx;
	(res1->segment).length = afterx - startx;
	res1->top_object = NULL;

	return res1;
}



/* TODO: link the origin blobs? */
/* TODO: delete the other origin */
static LObject* new_object(LCell* cell, int y) {
	LObject* res;
	LBlob* blob;
	assert(cell->top_object == NULL);
	res = (LObject*)LMALLOC(sizeof(LObject));
	res->id = nextobject;
	nextobject++;
	res->blob_next = NULL;
	cell->top_object = res;
	blob = (LBlob*)LMALLOC(sizeof(LBlob));
	lblobcount++;
	res->y = y;
	res->x = cell->segment.min_x;
	res->origin = blob;
	blob->objects = res;
	
	blob->range.x1 = cell->segment.min_x;
	blob->range.x2 = cell->segment.min_x+cell->segment.length;
	blob->range.y1 = y;
	blob->range.y2 = y;
	blob->id = nextblobid;
	nextblobid++;
	return res;
}


/* Handle new segment! */
/*
  distinguish new object from
  add to object,
  (startnewsegmentx, currentx) (not including)

Maintain current_line last_cell

multiple objects develop, their creation is ordered:
those having smallest y are created earlier, those with lower x after that earlier.

An object could be found to have to be merged with another object:
but that is always an earlier object?
Suppose we wish to have a hierarchy object->joined_next..
we also need object->joined_prev, as we need to find the root.

TOO HARD!!!!
Additional indirection? Have the cell point to a pointer containing the object.
So, all earlier segments point to the same address. And we can change its content in one 
action.
*p = other_object?
real_object_ptr = 

just give cell an ancestor cell pointer,
and then adjust the object of this ancestor cell. 
All descendants of an ancestor cell have NULL or invalid objects.
The descendance does not depend on merging or not.

The ancestor cell points to the current object.
There is a linked list of topcells.

handle_new_segment is incomplete!.
TODO:
make sure there is a link to the founder


How to link merged objects?
1) Leave the objects where they are, but link them together.
   It will be hard then to determine whether two segments really belong to one object.
   We could have objects contain topmost_ancestor objects.
   and there should be an object->next list which traverses the merged objects in some
   preferred way.
2)Is adjusting all cells of to be merged cells an option?
  Remove one of the objects afterwards?

Aspects:
How can we conveniently scan through our object?
top down, left right?
top down better, so top-left should be left.
lowest y value.

How to scan horizontally?
- always possible by doing cell->next until? there is no guarantee.
- find leftmost entry next line?

OK, lets try, LObject reuse:
as root info of each blob-top
as pointer to the real defining object.
That object should have pointers to all constituent objects, ordered in y, x-order.

If we join objects/vertical-blobs that were not linked yet, we should not eat the "embedded"
objects as they are separate, have space on top.
*/
static int cmp_object_yx(LObject* first, LObject* second) {
	assert(first->y>=0);
	assert(second->y>=0);
	assert(first->x>=0);
	assert(second->x>=0);
	if (first->y < second->y) return -1;
	if (second->y < first->y) return 1;
	if (first->x < second->x) return -1;
	if (second->x < first->x) return 1;
	ERROR_EXIT("Shouldnt compare same object\n");
	return -1;
}

static void combine_replace_left_box(LBox *left, LBox* right) {
	if (right->x1 < left->x1) left->x1 = right->x1;
	if (right->x2 > left->x2) left->x2 = right->x2;
	if (right->y1 < left->y1) left->y1 = right->y1;
	if (right->y2 > left->y2) left->y2 = right->y2;
}





static int lobject_length(LObject* first) {
	int res = 0;
	while (first != NULL) {
		res++;
		first = first->blob_next;
	}
	return res;
}
/*
1.No common objects allowed.
2.Ensure that any unneeded allocated ListElements are freed

Left should always have blob origin.

TODO: define primitives of next_left, next_right
 */
static LObject* merge_yx_lists(LObject* left, LObject * right, LBlob* blob) {
	LObject* result = NULL;
	LObject* lastcell = NULL;
	int ln;
	assert(left != NULL);
	assert(right != NULL);
	assert(left->origin == blob);
	ln = lobject_length(left) + lobject_length(right);
	if (left != NULL&& cmp_object_yx(left, right) < 0) {
		result = left;
		left = left->blob_next;
		lastcell = result;
	} else {
		result = right;
		lastcell = result;
		assert(right->origin != blob);
		right->origin = blob;
		if (right != NULL) right = right->blob_next;
	}
	assert(lastcell != NULL);
	assert(lastcell->origin == blob);

	while (1) {
		if (left != NULL&& (right == NULL||cmp_object_yx(left, right) < 0)) {
			lastcell->blob_next = left;
			lastcell = left;
			assert(left->origin == blob);
			left = left->blob_next;
		} else if (right !=NULL && (left == NULL||cmp_object_yx(left, right) > 0)) {
			lastcell->blob_next = right;
			lastcell = right;
			assert(right->origin != blob);
			right->origin = blob;
			right = right->blob_next;
		} else {
			if (left == NULL) {
				lastcell->blob_next = right;
				while (right != NULL) {
					assert(right->origin != blob);
					right->origin = blob;
					right = right->blob_next;	
				}
			} else if (right == NULL) {
				lastcell->blob_next = left;
				while (left != NULL) {
					assert(left->origin == blob);
					left = left->blob_next;	
				}
			} else {
				ERROR_EXIT("merge_yx_lists inconsistency");
			}
			assert(lobject_length(result) == ln);
			return result;
		}
	}
}









static void merge_top_objects_first(LObject* left_top, LObject* right_top) {
	LBlob* toremove = right_top->origin;
	assert(cmp_object_yx(left_top, right_top)<0);
	left_top->origin->objects =
		merge_yx_lists(left_top->origin->objects,
			       right_top->origin->objects, left_top->origin);
	combine_replace_left_box(&(left_top->origin->range), &(right_top->origin->range));
	lblobcount--;
	LFREE((void*)toremove);
}
/*
Which blob to keep? The left-top most.
 */
static void merge_top_objects(LObject* left_top, LObject* right_top) {
	if (left_top->origin == right_top->origin) {
		/* check for embedded/surrounded objects */
	} else {
		if (cmp_object_yx(left_top, right_top)<0) {
			merge_top_objects_first(left_top, right_top);
		} else {
			merge_top_objects_first(right_top, left_top);
		}
	}
}


/* DEBUGGING: kin.tif
we have segment y=14, [24, 30) but should have 21 instead of 24...
Th esegments are read incorrectly:
for line y=14 we should have x:
[4,6][8,11][15,17)[21,30)
earlier: y=10
[4,6)[15,16) instead of [16,17)  hex printout is correct!
 0                                        ",
 1                                        ",
"                                        ",
"                                        ",
"                                        ",
"                                        ",
"                                        ",
07    ..                                  ",
      ..                                  ",
      ..                                  ",
10    ..         ..                       ",
11    ..         ..                       ",
12    ..                                  ",
13    ..   ..    ..    .. ....            ",
14    ..  ...    ..    .........          ",
15    .. ...     ..    ....   ...         ",
16    .....      ..    ...     ..         ",
17    ....       ..    ..      ..         ",
18    ....       ..    ..      ..         ",
19    .....      ..    ..      ..         ",
20    .. ...     ..    ..      ..         ",
      ..  ...    ..    ..      ..         ",
      ..   ...   ..    ..      ..         ",
      ..    ..   ..    ..      ..         ",
                                         ",
"                                        ",
"                                        ",
"                                        ",
"                                        ",
"                                        "
*/
static void printlineblobs(LCell *last_cell_previous_line) {
	LCell* l = last_cell_previous_line;
	while (l !=NULL) {
		printf("s:[%d,l:%d] blobid:%d\n",l->segment.min_x, l->segment.length,
		       l->top_object->origin->id);
		l = l->next;
	}
	printf("ENDprintlineblobs:\n");
	
}

/*
TODO:
if a previous_line_segment does not overlap, check whether there is still
a connection. If not, close the object.

 */

static void handle_new_segment(int startx, int afterx, int y, LBlobinfo* info) {
	LCell* current_cell;
	int first = 1;
	current_cell = new_lcell(startx, afterx);
	if (info->last_cell_current_line == NULL) {
		assert(info->llines[y].first == NULL);
		info->llines[y].first = current_cell;
	} else {
		info->last_cell_current_line->next = current_cell;
	}
	info->last_cell_current_line = current_cell;
	while (info->last_cell_previous_line != NULL &&
	       startx >= info->last_cell_previous_line->segment.min_x +
	       info->last_cell_previous_line->segment.length) {
		info->last_cell_previous_line =
			info->last_cell_previous_line->next;
	}
	while (info->last_cell_previous_line!=NULL&&
	       info->last_cell_previous_line->segment.min_x < afterx -1) {
		if (first) {
			assert(current_cell->top_object == NULL);
			current_cell->top_object  =
				info->last_cell_previous_line->top_object;
			first = 0;
		} else {
			merge_top_objects(current_cell->top_object,
					  info->last_cell_previous_line->top_object);
		}
		if (info->last_cell_previous_line->segment.min_x+
		    info->last_cell_previous_line->segment.length <=
		    afterx -1) {
			info->last_cell_previous_line =
				info->last_cell_previous_line->next;
		} else {
			assert(!first);
			return;
		}
	}
	if (first) {
		/* NEW OBJECT  */
		info->last_object = new_object(current_cell, y);
	}
}


static void analysetiff(const char * filename, TIFF* tif) {

	uint32 imagelength, imagewidth;
        uint32 row;
	uint16 bitspersample;
	tsize_t scanlinesize;
        char* buf;
	char * bufptr;
	LActive_areas* active_areas;
	LActive_areas* active_areas_left;
	int hex = test3 & 1;
	int bin = test3 & 0;
	int currenty = 0;
	LScanline* llines;
	LBlobinfo blobinfo;
	TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &imagewidth);
	TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &bitspersample);
        TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &imagelength);
	llines = (LScanline*)_TIFFmalloc(imagelength*sizeof (LScanline));
	{
		LBlobinfo blobinfo1 = {NULL, NULL, NULL, llines};
		blobinfo = blobinfo1;
	}
	if (bitspersample != 1) {
		fprintf(stderr, "Sorry, only handle 1-bit samples.\n");
		return;
	}
	active_areas = setupactionareas(filename, imagewidth, imagelength);
	active_areas_left= active_areas;
	scanlinesize = TIFFScanlineSize(tif);
        bufptr = buf = (char*)_TIFFmalloc(imagelength*sizeof (LScanline));

	
	{
		unsigned int i;
		for (i = 0; i < imagelength; i++) {
			llines[i].y = i;
			llines[i].first = NULL;
		}
	}
        for (row = 0; row < imagelength; row++) {
		tsize_t charstodo = scanlinesize;
		tsize_t bitstodo = imagewidth;
		int currentx = 0;
		int insegment = 0;
		TIFFReadScanline(tif, buf, row, (tsample_t)0);
		bufptr = buf;

		while (charstodo > 0) {
			unsigned char nextchar = *bufptr;

			tsize_t bitstodochar = 8;
			unsigned char mask = 128;
			int startnewsegmentx;
			bufptr++;
			charstodo --;
			
			if (hex) printf("%02X ", nextchar);
			while (bitstodochar > 0&&bitstodo > 0) {
				if ((nextchar & mask) != 0) {
					if (!insegment) {
						startnewsegmentx = currentx;
						insegment = 1;
					}
					if (bin) printf("x ");
					handle_active_areas_per_pixel(currentx,
								      currenty,
								      &active_areas_left);
				} else {
					if (insegment) {
						handle_new_segment(startnewsegmentx,
								   currentx,
								   currenty,
								   &blobinfo);
						insegment = 0;
					}
					if (bin) printf("  ");
				}
				mask = (mask >> 1);
				bitstodochar --;
				bitstodo --;
				currentx++;
			}
			
		}
		if (test3) printf("|\n");
		blobinfo.last_cell_current_line = NULL;
		blobinfo.last_cell_previous_line =
			blobinfo.llines[currenty].first;
		currenty++;
	}
	active_areas_left= active_areas;
	while (active_areas_left != NULL) {
		if (active_areas_left->area.analyse != NULL) {
			active_areas_left->area.analyse(active_areas_left->area.data);
		}
		active_areas_left = active_areas_left->next;
	}
	printf("Defined %d blobs\n", lblobcount);
	_TIFFfree(buf);
	printf("------------------------------\n");
        /*TIFFClose(tif);*/
}

static void
tiffinfo(const char * filename, TIFF* tif, uint16 order, long flags, int is_image)
{
	
	TIFFPrintDirectory(tif, stdout, flags);
	if (!is_image) return;
	if (!readdata)
		return;
	if (doanalyse) {
		analysetiff(filename, tif);
	}
	if (rawdata) {
		if (order) {
			uint16 o;
			TIFFGetFieldDefaulted(tif,
			    TIFFTAG_FILLORDER, &o);
			TIFFReadRawData(tif, o != order);
		} else
			TIFFReadRawData(tif, 0);
	} else {
		if (order)
			TIFFSetField(tif, TIFFTAG_FILLORDER, order);
		TIFFReadData(tif);
	}
}

/* vim: set ts=8 sts=8 sw=8 noet: */
/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 8
 * fill-column: 78
 * End:
 */
