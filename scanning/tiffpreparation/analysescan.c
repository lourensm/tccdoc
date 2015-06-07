
#include <stdio.h>
#include <stdlib.h>

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
#include "analysescan.h"
#include "activearea.h"
#include "analysepage.h"
/*
   Split off analysis functionality is really good!

Next:
-Split off the other part, i.e. the angle determination
-Define callbacks.

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



















typedef struct lxyvalues {
	int z1, z2;
	int lastpospos;
	int* zvalues;
	double z0;
} LXYValues;










		







typedef struct lactive_areas {
	struct lactive_areas *next;
	LArea area;
} LActive_areas;

/*
Do the handlepixel per segment, call it segmentsize times.
 */


static void handle_endline(LActive_areasPtr active_areas_left, int currenty) {
	while (active_areas_left != NULL) {
		LArea* area =  &(active_areas_left->area);
		if (currenty>= area->box.y1 && currenty<=area->box.y2) {
			if (area->handle_endline != NULL) {
				area->handle_endline(currenty, area->data);
			}
		}
		active_areas_left = active_areas_left->next;
	}
}


/* touching pixels are connected */
static void handle_segment_active_areas(int startx, int afterx, int y,
					LActive_areas* active_areas_left) {
	while (active_areas_left != NULL) {
		LArea* area =  &active_areas_left->area;
		if (y >= area->box.y1 && y<=area->box.y2 && startx <= area->box.x2
		    && afterx >= area->box.x1) {
			int sx= (startx>= area->box.x1)?startx:area->box.x1;
			int ax= (afterx <= area->box.x2)?afterx:area->box.x2 + 1;
			if (area->handlepixel != NULL) {
				for (int x = sx;x < ax;x++) {
					area->handlepixel(x, y, area->data);
				}
			}
			if (area->handlesegment != NULL) {
				area->handlesegment(sx, ax, y, area->data);
			}
		}
		active_areas_left = active_areas_left->next;
	}
}






static void handle_segment(int startx, int afterx, int y, LActive_areas* active_areas_left) {
	handle_segment_active_areas(startx, afterx, y, active_areas_left); 
}



LScanOptions* set_scan_options() {
	LScanOptions* res = (LScanOptions*)malloc(sizeof(LScanOptions));
	res->reverse_bits = 1;
	res->resultdir = "results";
	return res;
}




LArea* add_active_area(LActive_areasPtr* la) {
	LActive_areas* p = malloc(sizeof (LActive_areas));
	p->next = NULL;
	if (*la == NULL) {
		*la = p;
	} else {
		(*la)->next = p;
	}
	return &(p->area);
}

/*
 TShould we handle pixel based analysis and segment based differently?
 */
void analysescan(const char * filename, TIFF* tif, 
		 LScanOptions* options) {

	uint32 imagelength, imagewidth;
	uint16 bitspersample;
	tsize_t scanlinesize;
	short interpretation;
        unsigned char* buf;
	unsigned char * bufptr;
	LActive_areas* active_areas = NULL;
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
	if (bitspersample != 1) {
		fprintf(stderr, "Sorry, only handle 1-bit samples.\n");
		return;
	}
	setuppageareas(&active_areas, filename, imagewidth, imagelength,
					options);
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
							       active_areas);
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
		handle_endline(active_areas, currenty);
		currenty++;

	}
	LActive_areas* active_areas1 = active_areas;
	while (active_areas1 != NULL) {
		if (active_areas1->area.analyse != NULL) {
			active_areas1->area.analyse(active_areas1->area.data);
		}
		active_areas1 = active_areas1->next;
	}
	while (active_areas != NULL) {
		if (active_areas->area.free != NULL) {
			active_areas->area.free(active_areas->area.data);
		}
		LActive_areas * a1 = active_areas;
		active_areas = active_areas->next;
		free(a1);
	}
	_TIFFfree(buf);
	printf("------------------------------\n");
        /*TIFFClose(tif);*/
}






/* vim: set ts=8 sts=8 sw=8 noet: */
/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 8
 * fill-column: 78
 * End:
 */
