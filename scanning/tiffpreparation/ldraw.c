#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <assert.h>
#include <libgen.h>
#include <string.h>

#include "tiffiop.h"
#include "tiffcommon.h"

#include "ldraw.h"

typedef struct lfpoint {
	double x;
	double y;
} LFPoint;

typedef struct segment {
	int x1;
	int x2;
} Segment;

typedef void UNUSED;
/*
y-y1 = (x - x1)*(y2-y1)/(x2-x1)

(y - y1)*(x2 - x1) = (x - x1)*(y2 - y1)
Maintain  DY = (y - y1)*(x2 - x1) and DX = (x - x1)*(y2 - y1)
If DY > DX we should either do y-- or x++

(x,y) such that |DY - DX| minimum
DY - DX:
y++ => DY - DX += x2 - x1
x++ => DY - DX -= y2 - y1

DY - DX > 0 => x++, DY-DX -= y2 - y1
DY - DX < 0 => y++, DY-DX += x2 - x1
DY - DX = 0 => x++;y++; DY-DX -= y2 - y1;DY-DX += x2 - x1

hele lijn moet bedekt zijn met pixels.
Als lijn door top van pixel gaat, dan moet pixel erboven erbij. (Dat gebeurt alleen als coef > 1
Als lijn door rechtervlak gaat, dan pixel rechts.
Ik vermoed dat we de offset van het midden moeten bijhouden.
Voortzetting of horizontaal , of vertikaal.

if (verticaloffset > horizontaloffset) y++
  decrease verticaloffset
else if (verticaloffset , horizontaloffset) x++
  decrease horizontaloffset
else y++, x++

Ok lets give up for now, use:
x >= floor(x1+(y - y1 - 0.5)*(dx/dy))
x <= ceil(x1 + (y - y1 + 0.5)*(dx/dy)-1)

	LFPoint p1 = {10,10};
	LFPoint p2 = {20,11};
 */
#define XLOW1(y) ceil(pylow->x+(y - pylow->y - 0.5)*(deltaxy))
#define XHI1(y) ceil(pylow->x+(y - pylow->y + 0.5)*(deltaxy)-1)
#define XTILDE(y)  pylow->x + (y - pylow->y)*deltaxy
static void test_draw_line(TIFF* out, LFPoint* p1, LFPoint* p2) {
	(UNUSED)out;
	double deltay = p2->y - p1->y;
	double deltax = p2->x - p1->x;
	double deltaxy = deltax/deltay;
	LFPoint *pylow, *pyhi;
	if (p1->y <= p2->y) {
		pylow = p1;
		pyhi = p2;
	} else {
		pylow = p1;
		pyhi = p2;
	}
	int ymax = floor(pyhi->y+0.5);
	int xalow = floor(p1->x<=p2->x?p1->x:p2->x+0.5);
	int xahi = floor(p1->x>=p2->x?p1->x:p2->x+0.5);
	int y = floor(pylow->y+0.5);

	if (deltay == 0) {
		printf("y=%d, x:[%d,%d]\n", y, xalow, xahi);
	} else {
		while (y <= ymax) {
			if (-1 <= deltaxy&& deltaxy <= 1) {
				int x = floor(XTILDE(y)+0.5);
				printf("y=%d, x:%d\n", y, x);
			} else {
				int xlow = XLOW1(y);
				int xhi = XHI1(y);
				if (xlow <= xhi) {
					if (xlow < xalow) xlow = xalow;
					if (xhi > xahi) xhi = xahi;
				} else {
					if (xhi < xalow) xlow = xalow;
					else xlow = xhi;
					if (xlow > xahi) xhi = xahi;
					else xhi = xlow;
				}
				printf("y=%d, x:[%d,%d]\n", y, xlow, xhi);
			} 
			y++;
		}
	}
}
static void test_draw_line1(TIFF* out, LFPoint* p1, LFPoint* p2) {
	(UNUSED)out;
	int y = floor(p1->y+0.5);
	int x1 = p1->x;
	int y1 = p1->y;
	double deltay = p2->y - y1;
	double deltax = p2->x - x1;
	assert(deltay >= 0);
	assert(deltax >= 0);
	while (y <= floor(p2->y+0.5)) {
		int xlow = XLOW(y);
		if (xlow < x1) xlow = x1;
		int xhi = XHI(y);
		if (xhi > p2->x) xhi = p2->x ;
		printf("y=%d, x:[%d,%d]\n", y, xlow, xhi);
		y++;
	}
}
static void test_draw_line2(TIFF* out, LFPoint* p1, LFPoint* p2) {
	(UNUSED)out;
	int x = floor(p1->x+0.5);
	int y = floor(p1->y+0.5);
	double deltay = p2->y - p1->y;
	double deltax = p2->x - p1->x;
	double dyx = (y - p1->y)*deltax - (x - p1->x)*deltay;
	int x1 = x;
	bool insegment = 1;
	while (1) {
		printf("(%d,%d) dyx=%lf\n", x, y, dyx);
		if (dyx > 0) {
			x++;
			dyx -= deltay;
		} else if (dyx < 0) {
			printf("segment([%d,%d],%d)\n",x1,x,y);
			x1 = x;
			y++;
			dyx += deltax;
		} else {
			if (x > floor(p2->x+0.5) &&  y > floor(p2->y+0.5)) {
				break;
			}
			if (deltay == 0) {
				x++;
			} else if (deltax == 0) {
				printf("segment([%d,%d],%d)\n",x1,x,y);
				x1 = x;
				y++;
			} else {
				printf("segment([%d,%d],%d)\n",x1,x,y);
				x++;y++;
				x1 = x;
				dyx -= deltay;
				dyx += deltax;
			}
		}
		printf("dye:%lf\n", dyx);
		printf("x:%d, fx:%d  y:%d, fy:%d\n",x , (int)floor((p2->x)+0.5),
		       y ,(int)floor(p2->y+0.5));
		if  (x >= floor(p2->x+0.5) &&  y >= floor(p2->y+0.5)) {
			break;
		}
	}
	printf("segment([%d,%d],%d)\n",x1,x,y);
}
/*
shapes:
inactive shapes: shapes having firsty > currenty
active shapes: shapes having firsty <= currenty and lasty > currenty
done shapes: rest
for any shape:setup
drawy: list of segments
I hate lists for C!!!
I love interesting results from working around introducing malloc'ed data.
The yshapes are part of one list!

Put them in a list for each y
if handley returns status end of range, remove from active list
 */

typedef enum _shapesegmentstatus {
	SHAPENEXTX,
	SHAPENEXTY,
	SHAPEEND
} ShapeSegmentStatus;

typedef struct yshape {
	void* data;
	struct yshape* next; 
	bool (*handley)(int y, void*data, ShapeSegmentStatus* status, Segment* segment);
} YShape;

typedef struct linedata {
	LFPoint p2;
	double deltax, deltay;
	double dyx;
	int x, y;
} Linedata;

bool handleliney(int y, void* data, ShapeSegmentStatus* status, Segment* segment) {
	Linedata* ldata = (Linedata*)data;
	printf("y=%d,d.y=%d\n", y, ldata->y);
	assert(y == ldata->y);
	segment->x1 = ldata->x;
	segment->x2 = ldata->x;
	if (ldata->dyx > 0) {
		ldata->x++;
		ldata->dyx -= ldata->deltay;
		*status = SHAPENEXTX;
		exit(1);
	} else {
		exit(2);
		*status = SHAPENEXTY;
		ldata->y++;
		ldata->dyx += ldata->deltax;		
		if (ldata->dyx < 0) {
		} else {
			ldata->x++;
			ldata->dyx -= ldata->deltay;
		}
	}
	if  (ldata->x > floor(ldata->p2.x+0.5)
	     &&  ldata->y > floor(ldata->p2.y+0.5)) {
		*status = SHAPEEND;
	}
	return true;
}
void setupline(LFPoint* p1, LFPoint* p2, YShape** yarray) {
	int x = floor(p1->x+0.5);
	int y = floor(p1->y+0.5);
	double deltay = p2->y - p1->y;
	double deltax = p2->x - p1->x;
	assert(deltay >= 0);
	assert(deltax >= 0);
	double dyx = (y - p1->y)*deltax - (x - p1->x)*deltay;
	Linedata ldata = {*p2, deltax, deltay, dyx, x, y};
	Linedata* ldatap = malloc(sizeof(Linedata));
	*ldatap = ldata;
	YShape *res = malloc(sizeof(YShape));
	res->data =(void*)ldatap;
	res->handley = &handleliney;
	YShape* ycont = yarray[y];
	yarray[y] = res;
	res->next = ycont;
}


/*
 * surprisingly hard: x1, x2 bits 8 - x1 till 8 - x2 should be set
 * 0xFF >> x1 0xFF << x2
 *
 */
void outsegment(unsigned char * buf, Segment* segment) {
	int x1 = segment->x1;
	int x2 = segment->x2;
	printf("segment([%d,%d])\n", x1, x2);
	int index1 = x1/8;
	int index2 = x2/8;
	int add = x1%8;
	int add2 = x2%8;
	if (index2 > index1) {
		unsigned char res = 0xFF;
		res >>= add;
		buf[index1] |= res;
		int index2 = x2/8;
		int add2 = x2%8;
		int i;
		for (i = index1+1; i < index2; i++) {
			buf[i] = 0xFF; /* memset */
		}
		if (add2 > 0) {
			res = 0xFF;
			res <<= 7 - add2;
			buf[i] |= res;
		}
	} else {
		unsigned char res = 0xFF;
		res >>= 7 - (x2-x1);  /* remove enough 1's from the left */
		res <<= 7 - x2;  /* add the 0's from the right */
		buf[index1] |= res;
	}
}
/*
 *   active: the chain of active shapes
 *   yarray: the array of chains of shapes that should start at coordinate y
 *
 *   if a shape becomes done, it should be removed from the active chain.
 */
void drawy(unsigned char* buf, int y, YShape** active, YShape** yarray) {
	YShape* a1 = *active;
	YShape* stillactive = NULL;/* this shape is last encountered shape that 
                                      has content left*/
	bool intail = false; /* we have entered the yarray[y] tail */
	bool activeok = false;
	printf("drawy(y=%d)\n", y);
	if (a1 == NULL) {
		a1 = yarray[y];
		intail = true;
		assert(!activeok);
	}
	while (a1 != NULL) {
		
		ShapeSegmentStatus status;
		Segment segment;
		if (a1->handley(y, a1->data, &status, &segment))
			outsegment(buf, &segment);
		while (status == SHAPENEXTX) {
			if (a1->handley(y, a1->data, &status, &segment))
				outsegment(buf, &segment);
		}
		if (status == SHAPENEXTY) {
			stillactive = a1;
			if (!activeok) {
				activeok = true;
				*active = a1;
			}
		} else {
			assert(status == SHAPEEND);
			if (stillactive == NULL) {
				stillactive = a1->next;
				assert(!activeok);
			} else {
				assert(activeok);
				stillactive->next = a1->next;
			}
		}
		if (stillactive == NULL) {
			if (!intail) {
				a1 = yarray[y];
				intail = true;
				assert(!activeok);
			}
		} else if (stillactive->next == NULL) {
			if (!intail) {
				stillactive->next = yarray[y];
				intail = true;
				a1 = stillactive->next;
			} else a1 = NULL;
		} else {
			a1 = stillactive->next;
		}
	}
}

TIFF* setup_tif(const char* file, uint32 imagewidth, uint32 imageheight) {
	TIFF *out = TIFFOpen(file, "w");
	TIFFSetField(out, TIFFTAG_IMAGEWIDTH, imagewidth);
	TIFFSetField(out, TIFFTAG_IMAGELENGTH, imageheight);
	TIFFSetField(out, TIFFTAG_SAMPLESPERPIXEL, 1);
	TIFFSetField(out, TIFFTAG_BITSPERSAMPLE, 1);
        TIFFSetField(out, TIFFTAG_COMPRESSION, 4);
	TIFFSetField(out, TIFFTAG_PHOTOMETRIC, 0);
	printf("DEFRPS:%d\n", TIFFDefaultStripSize(out,0));
	TIFFSetField(out, TIFFTAG_ROWSPERSTRIP, TIFFDefaultStripSize(out,0));
	TIFFSetupStrips(out);
	return out;
}


void drawimage(TIFF* tif, YShape** yarray, uint32 imagewidth,
	       uint32 imageheight) {
	uint32 scanlinesize = TIFFScanlineSize(tif);
	assert(scanlinesize == imagewidth / 8 + (imagewidth % 8 != 0));
	unsigned char* buf = _TIFFmalloc(scanlinesize);
	YShape* active = NULL;
	for (uint32 y = 0; y<imageheight;y++) {
		for (uint32 col = 0; col < scanlinesize; col++) {
			/* USE MEMSET ? */
			buf[col] = 0;
		}
		drawy(buf, y, &active, yarray);
		if (TIFFWriteScanline(tif, buf, y, 0)<0) {
			ERROR_EXIT("problem with writescanline");
		}
		
	}
	_TIFFfree(buf);

}

void input_point(const char * prompt, LFPoint* point) {
	printf("%s",prompt);
	double x, y;
	int res = scanf("%lf,%lf", &x, &y);
	if (res != 2) {
		ERROR_EXIT("problem reading");
	}
	point->x = x;
	point->y = y;
}
void test_drawing(const char * file) {
	/*LFPoint p1 = {10,10}, p2 = {20,20};*/
	LFPoint p1, p2;
		while (true) {
		input_point("point1 x,y:\n",&p1);
		input_point("point2:\n",&p2);
		test_draw_line(NULL, &p1, &p2);
	}
	return;

	uint32 imageheight = 50, imagewidth = 50;
	
	TIFF *out = setup_tif(file, imagewidth, imageheight);

	printf("DEFRPS:%d\n", TIFFDefaultStripSize(out,0));
	
	YShape* yarray[imageheight];  /* Container for shapes, added to their starty */
	for (unsigned int y = 0; y < imageheight;y++) {  /* initialize */
		yarray[y] = NULL;
	}
	/* add shapes : */
	setupline(&p1, &p2, yarray);
	/* handle the pixeldrawing */
	drawimage(out, yarray, imagewidth, imageheight);
	TIFFClose(out);
	printf("Wrote test drawing \"%s\"\n", file);
}





/* vim: set ts=8 sts=8 sw=8 noet: */
/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 8
 * fill-column: 78
 * End:
 */
