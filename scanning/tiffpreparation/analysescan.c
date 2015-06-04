
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
#include "lblobinfo.h"
#include "activearea.h"
#include "pixelrotation.h"
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










		









/*
Do the handlepixel per segment, call it segmentsize times.
 */





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
	res->resultdir = "results";
	return res;
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
					options);
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






/* vim: set ts=8 sts=8 sw=8 noet: */
/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 8
 * fill-column: 78
 * End:
 */
