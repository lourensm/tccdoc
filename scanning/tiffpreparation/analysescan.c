
#include "analysescan.h"
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

#include "tiffcommon.h"
#include "lblobinfo.h"
/*
   Split off analysis functionality is really good!

Next:
-Split off the other part, i.e. the angle determination
-Add generic analysis not per area, but for whole range. Better yet,
 allow generic conditions. And remove from set...based on y ordering?
- Later combine functionality for blobs and angle.
- make segment encoding standardized, x1,x2 vss startx, afterx 


 */


void error_exit(const char* msg,const char* file, const int line) {
	fprintf(stdout, "***%s:%d  : FATAL ERROR: \"%s\"\n",
		file, line, msg);
	exit(EXIT_FAILURE);
}

















typedef struct larea {
	struct larea *next;
	/*double y1f, tophorf, widthf, heightf;*/
	LBox box;
	void* data;
	void (*handlesegment)(int startx, int afterx, int y, void* data);

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
  const char* resultdir;
} LPixelrotationdata;

/*
Start from the middle of the range.
Add pixel value to array y + offset
n = 0: offset = 0 for all x
n = 1: offset = -1 for 
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

	
	
	handlepixelrotationxy(y, (rdata->box).y1, (rdata->box).y2,
			      x, rdata->x0, rdata->yvalues,
			      lastpospos, rdata->max_atan_angle);
	handlepixelrotationxy(x, (rdata->box).x1, (rdata->box).x2, y, rdata->y0, 
			      rdata->xvalues, lastpospos, rdata->max_atan_angle);
}

static int* setupxyvalues(int z1, int z2, int lastpospos) {
	int *values;
	int yvaluesize = (2*lastpospos+1)*(z2 - z1+1);
	int i;
	assert(z2 >= z1);
	values = (int*)malloc(yvaluesize* sizeof(int));
	for (i = 0; i<yvaluesize;i++) {
		values[i] = 0;
	}
	return values;
}


static LPixelrotationdata* setuppixelrotationdata(const char* filespec, LBox box,
						  double max_atan_angle, int lastpospos,
						  const char*resultdir) {
	LPixelrotationdata* data = malloc(sizeof (LPixelrotationdata));
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
	data->resultdir = resultdir;

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
				   int z1, int z2, int* xyvalues,const char* resultdir) {
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
	max2 = (double*)malloc(possize*sizeof(double));
	max4 = (double*)malloc(possize*sizeof(double));
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
	free(max2);
	free(max4);
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
	const char* resultdir = rdata->resultdir;
	analysepixelrotationxy("y", rdata->filespec,
			       rdata->lastpospos,
			       (rdata->box).y1, (rdata->box).y2, rdata->yvalues,
			       resultdir);
	analysepixelrotationxy("x", rdata->filespec,
			       rdata->lastpospos,
			       (rdata->box).x1, (rdata->box).x2, rdata->xvalues,
			       resultdir);
}


/*
  Define 
  - the _box_ in the tif that should be analysed,  
  - the different angles that should be tried: between +- max_atan_angle,
     - using 2*lastpospos + 1 equal spaced increments.
  - the various analysisresults will be stored in files defined by filespec
     - in a directory resultdir.
  This specification is stored in a LArea datastructure. 
*/
static void setuppixelrotation(LArea* larea, const char* filespec, LBox box,
			       double max_atan_angle, int lastpospos,
			       const char*resultdir) {
	larea->next = NULL;
	larea->box = box;
	larea->data = (void*)setuppixelrotationdata(filespec, box, max_atan_angle, lastpospos, resultdir);
	larea->handlepixel = &handlepixelrotation;
	larea->handlesegment = NULL;
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


/* ugly way of constructing strings in C, I am not really good nor interested in this aspect. Use a trick by storing result in one prealloccated static buffer.
 */
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



/*
Do the handlepixel per segment, call it segmentsize times.
 */


static LActive_areas* setupactionareas(const char* filename, int width, int height,
				       const char*resultdir) {
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
	LActive_areas* r2 = malloc(sizeof (LActive_areas));
	r2->next = NULL;
	setuppixelrotation(&(r2->area),filespec1, boxl, latan, 100,
			   resultdir);
	return r2;
}














 


















	      

static void handle_segment_active_areas(int startx, int afterx, int y,
					LActive_areas** active_areas_left) {
	while (*active_areas_left != NULL
	       && (((*active_areas_left)->area).box.y2 < y)) {
		*active_areas_left = (*active_areas_left)-> next;
	}
	if (*active_areas_left != NULL && y>= ((*active_areas_left)->area).box.y1) {
		LActive_areas* active_areas_now = (*active_areas_left);
		LArea* area_left =  &(active_areas_now->area);
		while (active_areas_now != NULL&& y<=area_left->box.y2) {
			if (startx>=area_left->box.x1 && afterx<=area_left->box.x2) {
				if (area_left->handlepixel != NULL) {
					for (int x = startx;x < afterx;x++) {
						area_left->handlepixel(x, y,
								     area_left->data);
					}
				}
				if (area_left->handlesegment != NULL) {
					area_left->handlesegment(startx, afterx,
							  y, area_left->data);
				}
			}
			active_areas_now = active_areas_now->next;
			area_left =  &(active_areas_now->area);
		}	
	}
}

static void handle_segment(int startx, int afterx, int y, LActive_areas** active_areas_left, LBlobinfoPtr info) {
	handle_segment_blobs(startx, afterx, info);
	handle_segment_active_areas(startx, afterx, y, active_areas_left); 
}



LScanOptions* set_scan_options() {
	LScanOptions* res = (LScanOptions*)malloc(sizeof(LScanOptions));
	res->reverse_bits = 1;
	return res;
}






/*
 TShould we handle pixel based analysis and segment based differently?
 */
void analysescan(const char * filename, TIFF* tif, const char* resultdir,
		 LScanOptions* options) {

	uint32 imagelength, imagewidth;
	uint16 bitspersample;
	tsize_t scanlinesize;
	short interpretation;
        unsigned char* buf;
	unsigned char * bufptr;
	LActive_areas* active_areas;
	LActive_areas* active_areas_left;
	int hex = 0;
	int bin = 0;
	int currenty = 0;

	TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &imagewidth);
	TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &bitspersample);
        TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &imagelength);
	TIFFGetField(tif, TIFFTAG_PHOTOMETRIC, &interpretation);
	int invert_bits = 0;
	if (interpretation == 0) invert_bits = ! options->reverse_bits;
	else if (interpretation == 1) invert_bits = options->reverse_bits;
	else ERROR_EXIT("cannot handle some TIFFTAG_PHOTOMETRIC yet");
	LBlobinfoPtr blobinfo = init_lblobinfo(imagelength);
	if (bitspersample != 1) {
		fprintf(stderr, "Sorry, only handle 1-bit samples.\n");
		return;
	}
	active_areas = setupactionareas(filename, imagewidth, imagelength,
					resultdir);
	active_areas_left= active_areas;
	scanlinesize = TIFFScanlineSize(tif);
        bufptr = buf = (unsigned char*)_TIFFmalloc(imagelength*sizeof (unsigned char));

        for (uint32 row = 0; row < imagelength; row++) {
		tsize_t charstodo = scanlinesize;
		tsize_t bitstodo = imagewidth;
		int currentx = 0;
		int insegment = 0;
		TIFFReadScanline(tif, buf, row, (tsample_t)0);
		bufptr = buf;
		while (charstodo > 0) {
			unsigned char nextchar = *bufptr;
			if (invert_bits) nextchar = ~(nextchar);
			unsigned char mask = 128;
			int startnewsegmentx;
			bufptr++;
			charstodo --;			
			if (hex) printf("%02X ", nextchar);
			while (mask != 0&&bitstodo > 0) {
				if ((nextchar & mask) != 0) {
					if (!insegment) {
						startnewsegmentx = currentx;
						insegment = 1;
					}
					if (bin) printf("x ");
				} else {
					if (insegment) {
						handle_segment(startnewsegmentx,
							       currentx, currenty,
							       &active_areas_left,
							       blobinfo);
						insegment = 0;
					}
					if (bin) printf("  ");
				}
				mask = (mask >> 1);
				bitstodo --;
				currentx++;
			}
			
		}
		if (bin||hex) printf("|\n");
		currenty++;
		blobinfo_endline(blobinfo, currenty);

	}
	blobinfo_endimage(blobinfo);
	active_areas_left= active_areas;
	while (active_areas_left != NULL) {
		if (active_areas_left->area.analyse != NULL) {
			active_areas_left->area.analyse(active_areas_left->area.data);
		}
		active_areas_left = active_areas_left->next;
	}
	blobinfo_stats(blobinfo);
	free_lblobinfo(blobinfo);
	_TIFFfree(buf);
	printf("------------------------------\n");
        /*TIFFClose(tif);*/
}






