/*
 *  pixelrotation
 *
 *
 */

#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include "tiffcommon.h"
#include "activearea.h"
#include "pixelrotation.h"
#include "pixelrotationdev.h"


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

/*
Start from the middle of the range.
Add pixel value to array y + offset
n = 0: offset = 0 for all x
n = 1: offset = -1 for 


Within our rectangle of interest, defined by rdata->box:
Maintain a count in x and y direction of the sum of the pixels int the 
perpendicular direction.

handlepixelrotationxy deals with pixel at coordinate value z (=x or y)
between box extrema z1 and z2.

If z is y direction and oz is x direction, we have offset in y direction.
Then we have rotation relative to horizontal position.

 */
static void handlepixelrotationxy(int z, int z1, int z2, int  oz, double oz0, 
			   int* xyvalues,
			   int lastpospos, double max_atan_angle) {
	int n;
	assert(z >= z1);
	assert(z <= z2);
	for (n = -lastpospos;n <= lastpospos;n++) {
		int zz0 = z - z1 +
			(int)floor(0.5+(oz - oz0)*n*max_atan_angle/lastpospos);
		if (zz0 >=0 && zz0 < z2 - z1 + 1) {
			int index = n+lastpospos+ (2*lastpospos+1)*zz0;
			assert(index >= 0);
			assert(index < (z2 -z1 +1)*(2*lastpospos+1));
			xyvalues[index]++;
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


/*
yvalues:  f(y) = sum x of pixel(x,y)
 */
static void handlepixelrotation(int x, int y, void* data) {
	LPixelrotationdata* rdata = (LPixelrotationdata*)data;
	int lastpospos = rdata->lastpospos;

	updatemaxbox(x, y, &rdata->maxbox);

	
	
	handlepixelrotationxy(y, (rdata->box).y1, (rdata->box).y2,
			      x, rdata->x0, rdata->yvalues,
			      lastpospos, rdata->max_atan_angle);
	handlepixelrotationxy(x, (rdata->box).x1, (rdata->box).x2, y,
			      rdata->y0, rdata->xvalues, lastpospos,
			      rdata->max_atan_angle);
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

static void   freepixelrotationdata(void* data) {
	LPixelrotationdata* rdata = (LPixelrotationdata*)data;
	free(rdata->yvalues);
	free(rdata->xvalues);
}


static LPixelrotationdata* setuppixelrotationdata(const char* filespec, LBox box,
						  double max_atan_angle,
						  int lastpospos,
						  const char*resultdir,
						  const char * description) {
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
	data->description = description;
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


typedef struct lrotresult {
	double zlow, zhi;
	double atanrotation;
} LRotresult;
static void analysepixelrotationxy(const char* xory,
				   int z1, int z2,
				   int* xyvalues, LPixelrotationdata* rdata,
				   LRotresult* result) {
	const char* filespec   = rdata->filespec;
	int         lastpospos = rdata->lastpospos;
	const char* resultdir  = rdata->resultdir;
	const int   sbufsize   = 80;
	char filename[sbufsize];
	int possize;
	double* max2;
	double* max4;
	struct lcontpt  maxes[5];
	FILE *f;
	assert(strlen(filespec)+strlen(xory)+strlen("values")+
						 strlen("maxsums")+
						 strlen("__.dat") < sbufsize);
	possize = 2*lastpospos+1;
	max2 = (double*)malloc(possize*sizeof(double));
	max4 = (double*)malloc(possize*sizeof(double));

	for (int i=0;i < possize;i++) {
		max2[i] = 0;
		max4[i] = 0;
	}

	/* 5 largest yvalues of max2 */
	for (int i = 0;i < 5;i++) { maxes[i].max = 0;maxes[i].pos = 0; }

	snprintf(filename,sizeof(filename),"%s/%s_%s_values.dat", resultdir,
		 filespec, xory);
	f = fopen(filename, "w");	
	for (int z = 0; z < z2 - z1 +1; z++) {
		fprintf(f, "%d", z + z1);
		for (int i1 = 0;i1 < possize; i1++) {
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
	for (int i = 0;i < possize; i++) {
		fprintf(f, "%d\t%f\t%f\n", i, sqrt(max2[i]), sqrt(sqrt(max4[i])));
		for (int j=0;j<5;j++) {
			if (max2[i] > maxes[j].max) {
				for (int k =4; k > j; k--) {
					maxes[k] = maxes[k-1];
				}
				maxes[j].max = max2[i];
				maxes[j].pos = i;
				break;
			}
		}
	}
	free(max2);
	free(max4);
	fclose(f);
	int mpos = maxes[0].pos;
	int posindex;
	double m1,m2,m3,p1,p2,p3,p0;
	int pl,ph,sl,sh, yl, yh;
	printf("Analysis of \"%s\":  %s direction\n", rdata->description, xory);
	qsort(maxes, 5, sizeof(struct lcontpt), max_cmp);
	for (int i = 0;i<5;i++) {
		if (maxes[i].pos == mpos) {printf("*");posindex = i;}
		else {printf(" ");}
		printf("i:%d, p:%d, m:%f\n", i, maxes[i].pos, maxes[i].max);
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
           Detect border of rectangle, if we assume rotation of (pl,ph)
           Awkward ad hoc detection.
	*/
        
	for (int yy = 0; yy <  z2 - z1 +1; yy++) {
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
	int ylow = yl + (yh - yl)*(p0 - pl)/(ph-pl);
	sl = 0;
	sh = 0;
	for (int yy = z2 - z1 ; yy >= 0;yy--) {
		int countl = xyvalues[pl + possize*yy];
		int counth = xyvalues[ph + possize*yy];
		sl += countl;
		sh += counth;
		if (sl >10) yl = yy;
		if (sh >10) yh = yy;
		if (sl>10 && sh >10) break;
	}
	int yhi = yl + (yh - yl)*(p0 - pl)/(ph-pl);
	printf("(zl:%d, zh:%d)\n", yl + z1, yh + z1);
	printf("maxp:%f, %slow:%d, %shi:%d\n", p0, xory, z1 + ylow, xory,z1 + yhi);
	printf("        box:[%d, %d]\n", z1, z2);
	/* TODO:
           detect digits of pageno, put font box around them and that defines the
	   x position, and absolute window in case of 13tr, if we know maximum dimensions 
	   of text.
           Requires detection of blobs. and indexing, recognition.
           ylow, yhi
	 */
	result->zlow = ylow + z1;
	result->zhi = yhi + z1;
	result->atanrotation = (p0- 100.0)*rdata->max_atan_angle/lastpospos;
}

/*
 *    1:y - ylow  = - alpha*(x - xav)
 *    2:x - xlow = alpha*(y - yav)
 *    3:y - yhi  = - alpha*(x - xav)
 *    4:x - xhi  = alpha*(y -yav)
 *
 *    tl=1&2: x = [xlow + alpha*(ylow + alpha*xav - yav)]/(1+alpha*alpha)
 *    br=3&4: x = [xhi  + alpha*(yhi + alpha*xav - yav)]/(1+alpha*alpha)
 *    tr=1&4: x = [xhi  + alpha*(ylow + alpha*xav - yav)]/(1+alpha*alpha)
 *    bl=2&3: x = [xlow + alpha*(yhi + alpha*xav -yav)]/(1 + alpha*alpha)
 * Would be nice if I had better intuition about symmetry.
 * In practice our alpha is -alpha above.
 */
 void   analysepixelrotation(void* data) {
	LPixelrotationdata* rdata = (LPixelrotationdata*)data;
	LBox * box = &(rdata->box);
	LRotresult xr, yr;
	analysepixelrotationxy("y", box->y1, box->y2, rdata->yvalues,rdata,
			       &yr);
	analysepixelrotationxy("x", box->x1, box->x2, rdata->xvalues, rdata,
			       &xr);
	double ylow = yr.zlow;
	double xlow = xr.zlow;
	rdata->xlow = xlow;
	rdata->ylow = ylow;
	double alpha = yr.atanrotation;
	double yav = (box->y2 + box->y1)/2.0;
	double xav = (box->x2 + box->x1)/2.0;

	rdata->xtl = (xlow + alpha*(ylow - yav + alpha*xav))/(1 + alpha*alpha);
	rdata->ytl = ylow  - alpha*(rdata->xtl - xav);
	if (rdata->xtl < box->x1 || rdata->ytl < box->y1) {
		printf("WARNING: top left beyond box limit: %lf,%lf <-> %d,%d\n",
		       rdata->xtl, rdata->ytl, box->x1, box->y1);
	}
	printf("alpha:%lf, xlow:%lf, ylow:%lf\n", alpha, xlow, ylow);
	printf("xtl, ytl:%lf,%lf\n", rdata->xtl, rdata->ytl);

	double yhi = yr.zhi;
	double xhi = xr.zhi;
	rdata->xhi = xhi;
	rdata->yhi = yhi;
	rdata->xbr = (xhi + alpha*(yhi - yav +  alpha*xav))/(1 + alpha*alpha);
	rdata->ybr = yhi  - alpha*(rdata->xbr - xav);
	if (rdata->xbr > box->x2 || rdata->ybr > box->y2) {
		printf("WARNING: bottom right beyond box limit: %lf,%lf <-> %d,%d\n",
		       rdata->xbr, rdata->ybr, box->x2, box->y2);
	}
	printf("xhi:%lf, yhi:%lf\n", xhi, yhi);
	printf("xbr, ybr:%lf,%lf\n", rdata->xbr, rdata->ybr);

	rdata->xtr = (xhi + alpha*(ylow + alpha*xav - yav))/(1+alpha*alpha);
	rdata->ytr = ylow - alpha*(rdata->xtr - xav);
	printf("xtr, ytr:%lf,%lf\n", rdata->xtr, rdata->ytr);

	rdata->xbl = (xlow + alpha*(yhi + alpha*xav - yav))/(1+alpha*alpha);
	rdata->ybl = yhi - alpha*(rdata->xbl - xav);
	printf("xbl, ybl:%lf,%lf\n", rdata->xbl, rdata->ybl);

	rdata->alpha = alpha;
	printf("-draw \"polygon %lf,%lf  %lf,%lf  %lf,%lf %lf,%lf\" -draw \"polygon %lf,%lf  %lf,%lf  %lf,%lf %lf,%lf\"\n", rdata->xtl, rdata->ytl, rdata->xtr, rdata->ytr, rdata->xbr, rdata->ybr, rdata->xbl, rdata->ybl);
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
void setuppixelrotation(LArea* larea, const char* filespec, LBox box,
			       double max_atan_angle, int lastpospos,
			const char*resultdir,
			const char* description) {
	larea->box = box;
	larea->data = (void*)setuppixelrotationdata(filespec, box,
						    max_atan_angle,
						    lastpospos, resultdir,
						    description);
	larea->handlepixel = &handlepixelrotation;
	larea->handlesegment = NULL;
	larea->analyse = &analysepixelrotation;
	larea->handlesegment = NULL;
	larea->handle_endline = NULL;
	larea->free = &freepixelrotationdata;
}


/* wrong file, need some lboxutil file? */
void printlbox(const char*txt, LBox box) {
	printf("%s: x1:%d, y1:%d, x2:%d, y2:%d\n",
	       txt, box.x1, box.y1, box.x2, box.y2);
}




	      


/* vim: set ts=8 sts=8 sw=8 noet: */
/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 8
 * fill-column: 78
 * End:
 */
