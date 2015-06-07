#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <libgen.h>
#include <string.h>

#include "tiffcommon.h"
#include "activearea.h"
#include "pixelrotation.h"
#include "lblobinfo.h"
#include "analysepage.h"

/*
setup:
left and right page for rotation
blob detection.

At analysis time:
iterate over blobs:
in text part do analysis, statistics
detect pagenumbers, left/right.
 */
static void setboxfromranges(int width, double hor1, double hor2, double horall,
			     int height, double vert1, double vert2,
			     double vertall, LBox* box) {
	double latan = 0.1; /* maximal expected misalignment angle */ 
	double scaleh = width/horall;
	double scalev = height/vertall;
	LBox box1 = {(int)round(scaleh*(hor1-latan*( (vert2-vert1)/2))),
		     (int)round(scalev*(vert1 - latan*( (hor2-hor1)/2))),
		     (int)round(scaleh*(hor2+latan*( (vert2-vert1)/2))),
		     (int)round(scalev*(vert2 + latan*( (hor2-hor1)/2)))};
	*box = box1; /* bit wasteful .. */
}





#define ALLOCMAX 20
/* ugly way of constructing strings in C, I am not really good nor interested in 
 * this aspect. 
 * Simply malloc and accept the leak!
 */
static const char* filespec(const char* filename, const char* add) {
	static int called = 0;
	const char* tiffext = ".tif";
	called++;
	if (called > ALLOCMAX) ERROR_EXIT("function filespec urgently needs redesign");
	char* basename1 = basename((char*)filename);
	int sz = strlen(basename1) + (int)strlen(add)+1- (int)strlen(tiffext);
	char * value = malloc(sz* sizeof(char));

	int nl = strlen(basename1), el = strlen(tiffext);
	if (nl < el || strcmp(basename1 + nl - el, tiffext)){
		fprintf(stderr, "ERROR:expected %s extension for \"%s\"\n",
			tiffext, filename);
		exit(1);
	}
	int i;
	for (i = 0;i < nl - el; i++) {
		assert(i<sz);
		value[i] = basename1[i];
	}
	for (size_t j= 0; j < strlen(add); j++) {
		assert(i < sz);
		value[i++] = add[j];
	}
	assert(i<sz);
	value[i] = '\0';
	return value;
}



void setuppageareas(LActive_areasPtr* la, const char* filename,
		    int width, int height, LScanOptions* options) {
	const char* resultdir = options->resultdir;
	/* two windows for now, detecting text and rotation */
	/* first 1 window.. */
	/*LHS: expected:Hor 23mm, 146mm vert:20mm 202mm
          RHS: 183mm, 306mm
	  (all HOR 332mm, vert:235)
	  allow for rotation of box..first allow around center only
	  latan (x - latan*(235)
	  y - latan
	 */
	
	double latan = 0.1;
	double hor1 = 23, hor2 = 146, hor3 = 183, hor4 = 306, horall = 332;
	double vert1 = 20, vert2 = 202, vertall = 235;
	LBox boxl, boxr;
	setboxfromranges(width, hor1, hor2, horall, height, vert1, vert2, vertall, &boxl);
	setboxfromranges(width, hor3, hor4, horall, height, vert1, vert2, vertall, &boxr);
	/*double vertpageno = 213, vertpagenoh = 3;
	  double pagel = 6;*//* 3 digits */
	printlbox("Leftwindow", boxl);
	printlbox("Rightwindow", boxr);
	LArea* p1 = add_active_area(la);
	const char* filespec1 = filespec(filename, "_l");
	setuppixelrotation(p1,filespec1, boxl, latan, 100, resultdir);
	
	LArea* p2 = add_active_area(la);
	filespec1 = filespec(filename, "_r");
	setuppixelrotation(p2,filespec1, boxr, latan, 100, resultdir);
	p1->next = p2;
	LArea* p3 = add_active_area(la);
	LBox blobspace = {0, 0, width, height+10000};
	printlbox("Blobarea", blobspace);
	setupblobarea(p3, blobspace, "wholeimage");
}


/* vim: set ts=8 sts=8 sw=8 noet: */
/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 8
 * fill-column: 78
 * End:
 */
