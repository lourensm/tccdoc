
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

/*
   TODO:
    status?
    on the way to making blobs accessible: e.g. iterator over (x,y) values
    or sequential access to segments, or...

    This requires that all cells are accessible from the blob.
    A blob has a list of objects.
    An object has a pointer to the list of cells.

    After this, I would like code to determine end of blobs.


accessibility probably not available for cell->object_next...
in particular when one previous line segment is touched by multiple current_y


Characterise the blobs on the page.

Possibly: ignore parts of page and group other blobs per page part:
left page, right page, left pageno, right pageno.
 */

#define ERROR_EXIT(MSG) error_exit(MSG, __FILE__, __LINE__)
static void error_exit(const char* msg,const char* file, const int line) {
	fprintf(stdout, "***%s:%d  : FATAL ERROR: \"%s\"\n",
		file, line, msg);
	exit(EXIT_FAILURE);
}

typedef struct lsegment {
	int min_x;
	int length;
} LSegment;


/*
  How to link the LCell's belonging to one LObject? 
At least there is one chain, up and down between LCells.

Give LObject a pointer to first cell and to lastcell? No need to store y values?
Or, only define an object_next pointer...
NO NO NO: joining multiple objects...
 */
typedef struct lcell {
	struct lcell *next;
	struct lcell *object_next;
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
	struct lblob* next;
	struct lblob* prev;
	int open_object_count;
} LBlob;

typedef struct lobject {
	int id;
	struct lblob *origin;
	int x;
	int y;
	LCell* top_cell;
	struct lobject *blob_next; 
} LObject;



/*
last_cell_last_new_object: 
    we have just added a new cell connecting to last_cell_previous_line.
    Could be within handle_segment, or between them.
 */
typedef struct lblobinfo {
	LCell *last_cell_previous_line;
	LCell *last_cell_current_line;
	LObject* last_object;
	LScanline* llines;
	LObject* last_cell_last_new_object;
	int currenty;
} LBlobinfo;




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
	values = (int*)malloc(yvaluesize* sizeof(int));
	for (i = 0; i<yvaluesize;i++) {
		values[i] = 0;
	}
	return values;
}


static LPixelrotationdata* setuppixelrotationdata(const char* filespec, LBox box, double max_atan_angle, int lastpospos, const char*resultdir) {
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
blob determination per area detection? probably, no, probably yes:
one area definition as set of rectangles.
Each complete blob: handle per active area? Store with active area, criteria for 
addition? 
Can we simply detect end of blob? 
*/
/*
Do the handlepixel per segment, call it segmentsize times.

how to organize?
sort per y -> permanently ignore after box passed.
per line? possibly order all active_areas based on endx?
Oh, the ordering of y should also be by lowest last-y 
 */
/*
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
*/
static LActive_areas* setupactionareas(const char* filename, int width, int height, const char*resultdir) {
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
		LActive_areas* r2 = malloc(sizeof (LActive_areas));
		r2->next = NULL;
		setuppixelrotation(&(r2->area),filespec1, boxl, latan, 100,
				   resultdir);
		return r2;
	}
}

static int lblobcount = 0;
static int nextblobid = 0;
static int nextobject = 0;




static LCell* new_lcell(int startx, int afterx) {
	LCell* res1 = (LCell*)malloc(sizeof(LCell));
	res1->next = NULL;
	res1->object_next = NULL;
	(res1->segment).min_x = startx;
	(res1->segment).length = afterx - startx;
	res1->top_object = NULL;

	return res1;
}



static LBlob* last_lblob = NULL;
static LBlob* first_lblob = NULL;
int open_blobs = 0;

static void blob_info(const char* where, int y) {
	return;
	printf("Blob Info at %s, y=%d\n", where, y);
	printf("%d blobs, *%d open blobs\n", lblobcount, open_blobs);
	LBlob* blobptr = first_lblob;
	while(blobptr != NULL) {
		printf("    id:%d x:[%d,%d] y:[%d,%d] openo:%d\n", blobptr->id,
		       blobptr->range.x1, blobptr->range.x2,
		       blobptr->range.y1, blobptr->range.y2,
		       blobptr->open_object_count
		       );
		/*
		printf("    id:%d x:[%d,l=%d), y:[%d, l=%d), open obj:%s\n", blobptr->id,
		       blobptr->range.x1, blobptr->range.x2 - (blobptr->range.x1) + 1,
		       blobptr->range.y1, blobptr->range.y2 - (blobptr->range.y1) + 1,
		       blobptr->open_object_count
		       );
		*/
		blobptr = blobptr->next;
	}
}


static int lobject_length(LObject* first) {
	int res = 0;
	while (first != NULL) {
		res++;
		first = first->blob_next;
	}
	return res;
}

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

/* TODO: standardize x ranges and then join combine_replace_left_box */
static void add_segment_to_blob(int startx, int afterx, int y, LBlob* blob) {
	if (startx < blob->range.x1) {
		blob->range.x1 = startx;
	}
	if (afterx - 1 >  blob->range.x2) {
		blob->range.x2 = afterx - 1;
	}
	if (y < blob->range.y1) {
		blob->range.y1 = y;
	}
	if (y > blob->range.y2) {
		blob->range.y2 = y;
	}
}

/* TODO: link the origin blobs? */
/* TODO: delete the other origin */
static LObject* new_object(LCell* cell, int y, LBlob* blob) {
	LObject* res;
	assert(cell->top_object == NULL);
	assert(cell->object_next == NULL);
	res = (LObject*)malloc(sizeof(LObject));
	res->id = nextobject;
	nextobject++;
	res->blob_next = NULL;
	cell->top_object = res;
	res->top_cell = cell;
	if (blob == NULL) {
		blob = (LBlob*)malloc(sizeof(LBlob));
		lblobcount++;
		blob->objects = res;
		blob->open_object_count = 1;
		blob->range.x1 = cell->segment.min_x;
		blob->range.x2 = cell->segment.min_x+cell->segment.length - 1;
		blob->range.y1 = y;
		blob->range.y2 = y;
		blob->id = nextblobid;
		blob->prev = last_lblob;
		blob->next = NULL;
		open_blobs++;
		if (last_lblob == NULL) {
			first_lblob = blob;
		} else {
			assert(last_lblob->next == NULL);
			last_lblob->next = blob;
		}
		last_lblob = blob;
		nextblobid++;
	} else {
		add_segment_to_blob(cell->segment.min_x,
				    cell->segment.min_x + cell->segment.length, y,
				    blob);
		/* Inefficient.. could define tail of objects */
		blob->objects = merge_yx_lists(blob->objects, res, blob);
		blob->open_object_count++;
	}
	res->y = y;
	res->x = cell->segment.min_x;
	res->origin = blob;

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

handle_segment is incomplete!.
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

static void combine_replace_left_box(LBox *left, LBox* right) {
	if (right->x1 < left->x1) left->x1 = right->x1;
	if (right->x2 > left->x2) left->x2 = right->x2;
	if (right->y1 < left->y1) left->y1 = right->y1;
	if (right->y2 > left->y2) left->y2 = right->y2;
}











static void pre_test(LObject* left_top, const char*where, ...) {
	return;
	if (left_top->origin->open_object_count <= 0) {
		printf("PRETEST:");
		va_list args;
		va_start(args, where);
		vfprintf(stdout, where, args);
		va_end(args);
		printf("\n");
		printf("oid:%d (x:%d,y:%d), lcellhasnext:%p, blobid:%d\n", left_top->id,
		       left_top->x, left_top->y,
		       left_top->top_cell->object_next,
		       left_top->origin->id);
		ERROR_EXIT("merge left object_count 0");
	}
}

static void merge_top_objects_first(LObject* left_top, LObject* right_top) {
	LBlob* toremove = right_top->origin;
	assert(cmp_object_yx(left_top, right_top)<0);
	int new_count = left_top->origin->open_object_count + right_top->origin->open_object_count-1;
	left_top->origin->objects =
		merge_yx_lists(left_top->origin->objects,
			       right_top->origin->objects, left_top->origin);
	combine_replace_left_box(&(left_top->origin->range), &(right_top->origin->range));
	left_top->origin->open_object_count = new_count;
	lblobcount--;
	
	if (toremove->next != NULL) {
		toremove->next->prev = toremove->prev;
	} else {
		assert(toremove == last_lblob);
		last_lblob = toremove->prev;
	}
	if (toremove->prev!= NULL) {
		toremove->prev->next = toremove->next;
	} else {
		assert(toremove == first_lblob);
		first_lblob = toremove->next;
	}
	open_blobs--;
	free((void*)toremove);
}
/*
Which blob to keep? The left-top most.
 */

static void merge_top_objects(LObject* left_top, LObject* right_top) {
	if (left_top->origin == right_top->origin) {
		/* check for embedded/surrounded objects */
		left_top->origin->open_object_count --;
		if(left_top->origin->open_object_count <= 0) {
			printf("ooc:%d, o1:%d o2:%d\n",
			       left_top->origin->open_object_count,
			       left_top->id, right_top->id);
			ERROR_EXIT("merge open_object_count 0");;
		}
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
/*
static void printlineblobs(LCell *last_cell_previous_line) {
	LCell* l = last_cell_previous_line;
	while (l !=NULL) {
		printf("s:[%d,l:%d] blobid:%d\n",l->segment.min_x, l->segment.length,
		       l->top_object->origin->id);
		l = l->next;
	}
	printf("ENDprintlineblobs:\n");
	
}
*/
/*
TODO:
if a previous_line_segment does not overlap, check whether there is still
a connection. If not, close the object.


What do I need to do with blobs?

- Detect the pagenumber
  
- Use that to further improve the alignment of the page

How to detect the pagenumber?
- define the area where to look
- find blobs there
- calculate the contour of the blobs
- match with other result of page alignment

Detect blob in area:
- check at end of blob
- maintain datastructure for blobs which makes access easy: too complex, interesting, but too complex.

So, check at end of blob or check each blob for overlap with all active_areas.

INTERMEZZO:
define rectangle of blob in rotated coordinate system, using the atan_angle.
define rectangle of blobs in rotated coordinate system, using the atan_angle.


Ok, store blobs, make them accessible at end.
Make contour of blob accessible.

Would be better if we store or handle blobs as soon as we know they are done?

Put blobs in linked list... 
 */

static void analyse_blob(LBlob* blob) {
}

static void decrease_blob_object_count(LCell* cell, LBlobinfo* info) {
	assert(cell->top_object->origin->open_object_count > 0);
	cell->top_object->origin->open_object_count --;
	if (cell->top_object->origin->open_object_count == 0) {
		open_blobs--;
	}
}



/*
Invariant should depend on what we have handled of the current line.
So, from current_line first cell to current_line before NOW,
the open_object_counts should reflect the cells that are "active".
Active cells:
before handle_segment, after handle_segment:
 
The y we are dealing with is the y of the current_cells.

invariant open_object_count, before entry:

current line, up to and including last_cell_current_line:
add to open_object_count
previous_line:
starting from last_cell_previous_line, apart from last_cell_previous_line if
last_cell_last_new_object is not NULL (that case is a bit tricky).

 */
typedef struct lblobtest {
	LBlob* blob;
	int open_segments;
	LBox range;
} LBlobtest;

/* REALLY STUPID, should really have used a hashtable but lets ignore as this was only for debugging */
#define BLOBN 80000

static void print_blob_counts(LBlobtest* blobt, int currenty, int start_x, LBlob * suspect, const char* where) {
	for (int i=0;i<BLOBN;i++) {
		LBlob* blob = blobt[i].blob;
		if (blob != NULL) {
			if (blob->open_object_count != blobt[i].open_segments) {
				printf("AT %s\n", where);
				printf("blobid=%d, oc:%d, sc:%d, y(next):%d, startx:%d\n",
				       blob->id, blob->open_object_count,
				       blobt[i].open_segments, currenty, start_x);
				if (suspect!=NULL) {
					if (blob == suspect) {
						ERROR_EXIT("suspect case invariant");
					} else {
						ERROR_EXIT("incomplete suspect case invariant");	
					}
				} else {
					ERROR_EXIT("unexplained open_count mismatch");
				}
			}
			
		}
	}
	printf("TEST INVARIANT OK\n");
}

void print_cell(LCell* cell, LBlobtest* blobt, const char* sp, int first) {
	int blobid = cell->top_object->origin->id;
	char str[15];
	const char* scount;
	if (blobt[blobid].blob == NULL) {
		scount = "N";
	} else {
		sprintf(str, "%d", blobt[blobid].open_segments);
		scount = str;
	}
	if (first) printf("||");
	printf("x:(%d,#%d), b:%d(%d,%s) %s",cell->segment.min_x,
	       cell->segment.length, cell->top_object->origin->id,
	       cell->top_object->origin->open_object_count,
	       scount, sp
	       );
}
/* before some currentx */
void test_invariant_end(LBlobinfo* info) {
	return;
	assert(info->last_cell_last_new_object == NULL);
	assert(info->last_cell_current_line == NULL);
	assert((info->currenty == 0&&info->last_cell_previous_line ==NULL)||
	       info->last_cell_previous_line ==
			info->llines[info->currenty-1].first);
	LBlobtest blobt[BLOBN];
	int deb = 1;
	if (deb) printf("INVARIANT END\n");
	for (int i=0;i<BLOBN;i++) {
		blobt[i].blob = NULL;
		blobt[i].open_segments = 0;
	}
	LCell* cell = info->last_cell_previous_line;
	int first = 1;
	while (cell != NULL) {
		int id = cell->top_object->origin->id;
		if (id < BLOBN) {
			blobt[id].blob = cell->top_object->origin;
			blobt[id].open_segments++;
		} else {
			ERROR_EXIT("cannot test more than BLOBN");
		}
		if (deb) print_cell(cell, blobt, "", first);
		first = 0;
		cell = cell->next;
	}
	print_blob_counts(blobt, info->currenty, -1, NULL, "test_invariant_end");
}


void test_invariant(LBlobinfo* info, int beforestartx, const char* where) {
	return;
	LBlobtest blobt[BLOBN];
	for (int i=0;i<BLOBN;i++) {
		blobt[i].blob = NULL;
		blobt[i].open_segments = 0;
	}
	int deb = 1;
	int first = 1;
	LCell * cell1 = info->currenty == 0?NULL:
		info->llines[((info->currenty) - 1)].first;
	if (deb) printf("CHAIN:previous_line upto last_cell_previous_line:\n");
	while (1) {
		if (cell1 == NULL) break;
		if (deb) print_cell(cell1, blobt, (cell1 == info->last_cell_previous_line&&
						   info->last_cell_last_new_object !=NULL)?"NW":"  ",
				    first);
		first = 0;
		if (cell1 == info->last_cell_previous_line) break;
		cell1 = cell1->next;
	}
	LCell* cell = 
		info->llines[(info->currenty)].first;
	LBlob * suspect = NULL;
	if (deb) printf("\nCHAIN:current_line:\n");
	first = 1;
	while (1) {
		if (cell == NULL) break;
		int id = cell->top_object->origin->id;
		if (id < BLOBN) {
			blobt[id].blob = cell->top_object->origin;
			if (info->last_cell_last_new_object == NULL||
			    cell != info->last_cell_current_line)  {			  
				blobt[id].open_segments++;
			} else {
				assert(suspect == NULL);
				suspect = info->last_cell_last_new_object->origin;
			}
		} else {
			ERROR_EXIT("cannot test more than BLOBN");
		}
		if (deb) print_cell(cell, blobt, suspect!=NULL?"SS":"  ", first);
		first = 0;
		if (cell == info->last_cell_current_line) break;
		cell = cell->next;
	}
	if (deb) {
		if (cell1 == NULL||cell1->next == NULL) {
			printf("\nCHAIN:NO previous_line NEXT");	
		} else {
			first = 1;
			printf("\nCHAIN:previous_line NEXT:\n");
			while (1) {
				cell1 = cell1->next;
				if (cell1 == NULL) break;
				print_cell(cell1, blobt, "  ", first);
				first = 0;
			}
		}
	}
	if (deb) printf("\nCHAIN:end\n");
	print_blob_counts(blobt, info->currenty, beforestartx, suspect, where);
	printf("TEST INVARIANT OK\n");
}

void test_all_open(LCell* first, const char* where, ...) {
	return;
	printf("All cells:\n");
	while (first != NULL) {
		printf("    x:(%d #%d), blobid=%d\n",first->segment.min_x,
		       first->segment.length, first->top_object->origin->id);
		if (first->top_object->origin->open_object_count <= 0) {
			printf("test_all_open:");
			va_list args;
			va_start(args, where);
			vfprintf(stdout, where, args);
			va_end(args);
			printf("\n");
			ERROR_EXIT("some closed cell");
		}
		first = first->next;
	}
}
/*
Focus: long previous_line segment, multiple consecutive new_segments.
Then, only the first new_cell should get to be linked directly, the others should 
introduce new objects.

We need to have info on how the last new_cell got linked to which last_cell_previous_line
We have last_cell_previous_line, we need the new_cell-object to which it has most recently been 
linked, NULL if not by last_cell_previous_line.

If we now have to skip non overlapping, then last_object becomes NULL

startx, afterx : the segment that needs to be incorporated into
info:  LBlobinfo* info

This implies that all segments should be part of a blob. They should be
accessible from the blob.

Segments are linked to a blob such that each segment is linked to parent
segments, until some top_object. A blob then has a list of top_objects that
make up its segments.

With each blob we maintain an oben_object_counter so as to be able to detect when the blob
is finished, i.e. all of its objects and segments no longer have a connection downwards.

handle_segment_blobs links the current segments to segments of the previous line (currenty - 1)

situations:
*current segment is connected to no previous_line segment.
* current segment is connected to exactly one previous_line segment which was not earlier
connected to some current line segment.
* multiple previous_line segments are connected to current_segment 

info->last_cell_last_new_object :
We just linked a current cell to a previousline cell, which is still available for linking to
current cell. We need not decrease its blob->object_count as that has a link already through that
current cell.
 */
static void handle_segment_blobs(int startx, int afterx, 
				 LBlobinfo* info) {
	int deb = 0;
	int y = info->currenty;
	if (deb) printf("HANDLESEGMENTBLOBS: x[%d,%d), y%d, lclno:%p\n", startx,
	       afterx, y,
	       info->last_cell_last_new_object);
	blob_info("handle_segment start", y);
	test_invariant(info, startx, "handle_segment_blobs");

	int first = 1;
	assert(info->last_cell_last_new_object == NULL ||
			       info->last_cell_last_new_object->origin ==
			       info->last_cell_previous_line->top_object->origin);
	LCell* current_cell = new_lcell(startx, afterx);
	if (deb) printf("lcpl:%p\n", info->last_cell_previous_line);
	if (deb) pre_test(info->last_cell_previous_line->top_object, "NU");
	if (info->last_cell_current_line == NULL) {
		assert(info->llines[y].first == NULL);
		assert(y==0||info->last_cell_previous_line == info->llines[y-1].first);
		info->llines[y].first = current_cell;
		info->last_cell_current_line = current_cell;
		if (deb) printf("YY:%p, p:%p\n", info->llines[y].first, info->last_cell_previous_line);
	} else {
		info->last_cell_current_line->next = current_cell;
	}
	info->last_cell_current_line = current_cell;
	/* skip non-overlapping segments previous line */
	while (info->last_cell_previous_line != NULL &&
	       startx >= info->last_cell_previous_line->segment.min_x +
	       info->last_cell_previous_line->segment.length) {
		if (info->last_cell_last_new_object == NULL) {
			if (deb) printf("HSB:DECREASE blob=%d\n",
			       info->last_cell_previous_line->top_object->origin->id);
			decrease_blob_object_count(info->last_cell_previous_line, info);
			if (info->last_cell_previous_line->next != NULL
			    && info->last_cell_previous_line->next->top_object->origin ==
			    info->last_cell_previous_line->top_object->origin)
				{
					assert(info->last_cell_previous_line->top_object->origin->open_object_count >= 1);
				}
		}
		info->last_cell_previous_line =
			info->last_cell_previous_line->next;
		info->last_cell_last_new_object = NULL;
	}

	while (info->last_cell_previous_line!=NULL&&
	       info->last_cell_previous_line->segment.min_x <= afterx - 1) {
		pre_test(info->last_cell_previous_line->top_object, "begin");
		if (first == 1) {  /* Need to assign new_cell to object */
			assert(current_cell->top_object == NULL);
			if (info->last_cell_last_new_object == NULL) {
				/* link to existing object */
				current_cell->top_object  =
					info->last_cell_previous_line->top_object;
				assert(info->last_cell_previous_line->object_next == NULL);
				info->last_cell_previous_line->object_next = current_cell;
				add_segment_to_blob(startx, afterx, y,
						    current_cell->top_object->origin);
				pre_test(current_cell->top_object, "first");
				if (deb) printf("HSB:simple link blob=%d,c=%d, px=[%d,#%d)\n",
				       current_cell->top_object->origin->id,
				       current_cell->top_object->origin->open_object_count,
				       info->last_cell_previous_line->segment.min_x,
				       info->last_cell_previous_line->segment.length);
			}  else {
				assert(info->last_cell_last_new_object->origin ==
				       info->last_cell_previous_line->top_object->origin);
				/* we need to assign new object, but do merging with
                                   info->last_cell_last_new_object
				 */
				info->last_object = new_object(current_cell, y, info->last_cell_last_new_object->origin);
				/* so we create a new blob but it immediately gets merged */
				if (deb) printf("HSB:link blob=%d, c=%d last_new, new_obj blob inc, px=[%d,#%d)\n",
				       current_cell->top_object->origin->id,
				       current_cell->top_object->origin->open_object_count,
				       info->last_cell_previous_line->segment.min_x,
				       info->last_cell_previous_line->segment.length);
			}
			
			first = 0;
		} else {
			int currentblobid = current_cell->top_object->origin->id;
			int otherblobid = info->last_cell_previous_line->top_object->origin->id;
			merge_top_objects(current_cell->top_object,
					  info->last_cell_previous_line->top_object);
			if (deb) printf("HSB:merge to blob=%d,count=%d, decrease blobs=(%d,%d), px=[%d,#%d)\n",
			       current_cell->top_object->origin->id,
			       current_cell->top_object->origin->open_object_count,
			       currentblobid, otherblobid, info->last_cell_previous_line->segment.min_x,
				       info->last_cell_previous_line->segment.length);
		}
		info->last_cell_last_new_object = current_cell->top_object;
		assert(info->last_cell_last_new_object->origin ==
				       info->last_cell_previous_line->top_object->origin);
		if (info->last_cell_previous_line->segment.min_x+
		    info->last_cell_previous_line->segment.length <= afterx -1) {
			info->last_cell_previous_line =
				info->last_cell_previous_line->next;
			info->last_cell_last_new_object = NULL;
		} else {
			assert(!first);
			assert(info->last_cell_last_new_object == NULL ||
			       info->last_cell_last_new_object->origin ==
			       info->last_cell_previous_line->top_object->origin);
			blob_info("handle_segment endlinked", y);
			return;
		}
	}
	if (first) {
		/* NEW OBJECT, no link to earlier segment  */
		info->last_object = new_object(current_cell, y, NULL);
		info->last_cell_last_new_object = NULL;
		if (deb) printf("HSB:blob=%d, NEW\n", current_cell->top_object->origin->id);
	}
	blob_info("handle_segment endnew", y);
}

static void handle_segment_active_areas(int startx, int afterx, int y, LActive_areas** active_areas_left) {
	while (*active_areas_left != NULL
	       && (((*active_areas_left)->area).box. y2< y)) {
		*active_areas_left = (*active_areas_left)-> next;
	}
	if (*active_areas_left != NULL && y>= ((*active_areas_left)->area).box.y1) {
		LActive_areas* active_areas_now = (*active_areas_left);
		LArea* area_left =  &(active_areas_now->area);
		while (active_areas_now != NULL&& y<=area_left->box.y2) {
			if (startx>=area_left->box.x1 && afterx<=area_left->box.x2) {
				if (area_left->handlepixel != NULL) {
					for (int x = startx;x < afterx;x++) {
						area_left->handlepixel(x, y, area_left->data);
					}
				}
				if (area_left->handlesegment != NULL) {
					area_left->handlesegment(startx, afterx, y, area_left->data);
				}
			}
			active_areas_now = active_areas_now->next;
			area_left =  &(active_areas_now->area);
		}	
	}
}

static void handle_segment(int startx, int afterx, int y, LActive_areas** active_areas_left, LBlobinfo* info) {
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
        char* buf;
	char * bufptr;
	LActive_areas* active_areas;
	LActive_areas* active_areas_left;
	int hex = 0;
	int bin = 0;
	int currenty = 0;
	int black_is_zero;
	LScanline* llines;
	LBlobinfo blobinfo;
	TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &imagewidth);
	TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &bitspersample);
        TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &imagelength);
	TIFFGetField(tif, TIFFTAG_PHOTOMETRIC, &interpretation);
	if (interpretation == 0) {
		black_is_zero = 0;
	} else if (interpretation == 1) {
		black_is_zero = 1;
	} else {
		black_is_zero = -1;
		ERROR_EXIT("cannot handle some TIFFTAG_PHOTOMETRIC yet");
	}
	int invert_bits = black_is_zero?options->reverse_bits:! options->reverse_bits;
	llines = (LScanline*)_TIFFmalloc(imagelength*sizeof (LScanline));
	{
		LBlobinfo blobinfo1 = {NULL, NULL, NULL, llines, NULL, currenty};
		blobinfo = blobinfo1;
	}
	if (bitspersample != 1) {
		fprintf(stderr, "Sorry, only handle 1-bit samples.\n");
		return;
	}
	active_areas = setupactionareas(filename, imagewidth, imagelength,
					resultdir);
	active_areas_left= active_areas;
	scanlinesize = TIFFScanlineSize(tif);
        bufptr = buf = (char*)_TIFFmalloc(imagelength*sizeof (LScanline));
	for (unsigned int i = 0; i < imagelength; i++) {
		llines[i].y = i;
		llines[i].first = NULL;
	}
	int deb = 0;
        for (uint32 row = 0; row < imagelength; row++) {
		tsize_t charstodo = scanlinesize;
		tsize_t bitstodo = imagewidth;
		int currentx = 0;
		int insegment = 0;
		TIFFReadScanline(tif, buf, row, (tsample_t)0);
		if (deb) printf("*********%d**********************************\n", currenty);
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
							       &blobinfo);
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
		while (blobinfo.last_cell_previous_line!= NULL) {
			if (blobinfo.last_cell_last_new_object == NULL) {
				decrease_blob_object_count(blobinfo.last_cell_previous_line, &blobinfo);
			} 
			blobinfo.last_cell_previous_line =
				blobinfo.last_cell_previous_line->next;
			blobinfo.last_cell_last_new_object = NULL;
		}
		blobinfo.last_cell_last_new_object = NULL;
		blobinfo.last_cell_current_line = NULL;
		blobinfo.last_cell_previous_line =
			blobinfo.llines[currenty].first;
		test_all_open(blobinfo.last_cell_previous_line, "next_line");
		currenty++;
		blobinfo.currenty = currenty;
		test_invariant_end(&blobinfo);

	}
	/* TODO:Handle open blobs, as if we have an additional 0s line: */
	if (blobinfo.last_cell_previous_line != NULL) {
		while (blobinfo.last_cell_previous_line != NULL) {
			decrease_blob_object_count(blobinfo.last_cell_previous_line, &blobinfo);
			blobinfo.last_cell_previous_line = blobinfo.last_cell_previous_line->next;
		}
	}
	active_areas_left= active_areas;
	while (active_areas_left != NULL) {
		if (active_areas_left->area.analyse != NULL) {
			active_areas_left->area.analyse(active_areas_left->area.data);
		}
		active_areas_left = active_areas_left->next;
	}
	printf("Defined %d blobs (%d blobs still open)\n", lblobcount, open_blobs);
	_TIFFfree(buf);
	printf("------------------------------\n");
        /*TIFFClose(tif);*/
}




/* Lourens...? left it, probably no effect on vim? using emacs anyway */ 
/* vim: set ts=8 sts=8 sw=8 noet: */
/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 8
 * fill-column: 78
 * End:
 */
