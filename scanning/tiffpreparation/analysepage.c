#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <assert.h>
#include <libgen.h>
#include <string.h>

#include "tiffiop.h"
#include "tiffcommon.h"
#include "activearea.h"
#include "pixelrotation.h"
#include "pixelrotationdev.h"
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
	double latan = 0.0; /* maximal expected misalignment angle */
	/* TODO: this is wrong, the rotation is not within this box, but 
           over a whole page */
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

typedef struct pagedata {
	const char* filename;
	TIFF* tif;
	unsigned char* buf;
	LScanOptions* options;
	LArea* lefttext;
	LArea* righttext;
	LArea* blobs;
	LArea* lpageno;
	LArea* rpageno;
	uint32 width, height;
} PageData;

PageDataPtr setuppageareas(LActive_areasPtr* la, const char* filename,
		    TIFF* tif, unsigned char* buf, LScanOptions* options) {
	uint32 width, height;

	TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &width);
	TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &height);
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
	LBox boxl, boxr, boxlp, boxrp;

	setboxfromranges(width, hor1, hor2, horall, height, vert1, vert2,
			 vertall, &boxl);
	setboxfromranges(width, hor3, hor4, horall, height, vert1, vert2,
			 vertall, &boxr);
	/*double vertpageno = 213, vertpagenoh = 3;
	  double pagel = 6;
          hor 26, 30, 304, 308, 333.5
          vert 212, 215/235.5
         img012_r.tif big rotation +- 9
         vert: 203, 224  hor 17, 39,   hor 295 313  / 333.5
        *//* 3 digits */
	double pl1 = 17.0, pl2 = 39.0, pr1 = 295.0, pr2 = 313.0, pv1 = 205.0, pv2 = 224.0,
		pla = 333.5, pva = 235.5;
	setboxfromranges(width, pl1, pl2, pla, height, pv1, pv2, pva, &boxlp);
	setboxfromranges(width, pr1, pr2, pla, height, pv1, pv2, pva, &boxrp);	
	printlbox("Leftwindow", boxl);
	printlbox("Rightwindow", boxr);
	printlbox("Leftpageno", boxlp);
	printlbox("Rightpageno", boxrp);
	int leftpage = 1, rightpage = 1, wholeblob = 0, leftpageblob = 1, rightpageblob = 1;
	LArea* lefttext = 0;
	if (leftpage) {
		lefttext = add_active_area(la);
		const char* filespec1 = filespec(filename, "_l");
		setuppixelrotation(lefttext,filespec1, boxl, latan, 100, resultdir,
				   "lefttext");
	}
	LArea* righttext = NULL;
	if (rightpage) {
		righttext = add_active_area(la);
		const char* filespec1 = filespec(filename, "_r");
		setuppixelrotation(righttext,filespec1, boxr, latan, 100, resultdir,
				   "righttext");
	}
	LArea* blobarea = NULL;
	if (wholeblob) {
		blobarea = add_active_area(la);
		LBox blobspace = {0, 0, width, height};
		printlbox("Blobarea", blobspace);
		setupblobarea(blobarea, blobspace, "wholeimage");
	}
	LArea* blobplarea = NULL;
	if (leftpageblob) {
		blobplarea = add_active_area(la);
		setupblobarea(blobplarea, boxlp, "leftpageno");
	}
	LArea* blobprarea = NULL;
	if (rightpageblob) {
		blobprarea = add_active_area(la);
		setupblobarea(blobprarea, boxrp, "rightpageno");
	}
	PageData pd = {filename, tif, buf, options, lefttext, righttext,
		       blobarea, blobplarea, blobprarea, width, height};
	PageData* result = malloc(sizeof(PageData));
	*result = pd;
	return result;
}


/* 
 * The left page: pageno starts at 195/235.5, size 3.5/235.5
 * right page with two digits:ends 124/333 of its leftmargin. two digits are 
 * 4.5/333 wide.
 * Three digits is 7.
 * First calculate the rotated box
 */
void printdrawbox(LBox* box) {
	printf("-draw \"polygon %d,%d  %d,%d  %d,%d %d,%d\" ",
	       box->x1, box->y1, box->x2, box->y1, box->x2, box->y2,
	       box->x1, box->y2);	
}
void printdrawline(double x1,double y1,double x2,double y2) {
	printf("-draw \"line %lf, %lf  %lf, %lf\" ",
	       x1, y1, x2, y2);
}
void printdrawmidlines(LBox* box) {
	double xav = (box->x1 + box->x2)/2;
	printdrawline(xav, box->y1,xav, box->y2);
	double yav = (box->y1 + box->y2)/2;
	printdrawline(box->x1, yav, box->x2, yav);
	       
}
void printdraw(LPixelrotationdata* ldata) {
	printdrawbox(&(ldata->box));
	printdrawmidlines(&(ldata->box));
	return;
	/* y - ylow = -alpha*(x - xav) 
	 * x = x1, x = x2
         * x - xlow = alpha*(y - yav)
	 */
	double xav = (ldata->box.x1 + ldata->box.x2)/2;
	double y1 = ldata->ylow - ldata->alpha*(ldata->box.x1 - xav);
	double y2 = ldata->ylow - ldata->alpha*(ldata->box.x2 - xav);
	printdrawline(ldata->box.x1, y1, ldata->box.x2, y2);
	double yh1 = ldata->yhi - ldata->alpha*(ldata->box.x1 - xav);
	double yh2 = ldata->yhi - ldata->alpha*(ldata->box.x2 - xav);
	printdrawline(ldata->box.x1, yh1, ldata->box.x2, yh2);

	double yav = (ldata->box.y1 + ldata->box.y2)/2;
	double x1 = ldata->xlow + ldata->alpha*(ldata->box.y1 - yav);
	double x2 = ldata->xlow + ldata->alpha*(ldata->box.y2 - yav);
	printdrawline(x1, ldata->box.y1, x2, ldata->box.y2);
	double xh1 = ldata->xhi + ldata->alpha*(ldata->box.y1 - yav);
	double xh2 = ldata->xhi + ldata->alpha*(ldata->box.y2 - yav);
	printdrawline(xh1, ldata->box.y1, xh2, ldata->box.y2);
	return;
	printf("-draw \"polygon %lf,%lf  %lf,%lf  %lf,%lf %lf,%lf\" ",
	       ldata->xtl, ldata->ytl, ldata->xtr, ldata->ytr, ldata->xbr,
	       ldata->ybr, ldata->xbl, ldata->ybl);
	
}


static void outtif(PageData* data) {

	TIFF* tif = data->tif;
	unsigned char* buf = data->buf;
	TIFF *out = TIFFOpen("out.tiff", "w");
	int sampleperpixel;
	uint32 imagelength, imagewidth;
	uint16 bitspersample;
	TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL, &sampleperpixel);
	TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &imagewidth);
	TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &imagelength);
	TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &bitspersample);
	uint32 rowsperstrip;
	TIFFGetField(tif, TIFFTAG_ROWSPERSTRIP, &rowsperstrip);
	printf("ROWSPERSTRIP:%d\n", rowsperstrip);
	short interpretation;
	TIFFGetField(tif, TIFFTAG_PHOTOMETRIC, &interpretation);
	static	uint16 compression;
	TIFFGetField(tif, TIFFTAG_COMPRESSION, &compression);
  
	TIFFSetField(out, TIFFTAG_IMAGEWIDTH, imagewidth);
	TIFFSetField(out, TIFFTAG_IMAGELENGTH, imagelength);
	TIFFSetField(out, TIFFTAG_SAMPLESPERPIXEL, sampleperpixel);
	TIFFSetField(out, TIFFTAG_BITSPERSAMPLE, bitspersample);
	TIFFSetField(out, TIFFTAG_ROWSPERSTRIP, rowsperstrip);
        TIFFSetField(out, TIFFTAG_COMPRESSION, compression);
	TIFFSetField(out, TIFFTAG_PHOTOMETRIC, interpretation);
	TIFFSetupStrips(out);
	uint32 scanlinesize = TIFFScanlineSize(tif);
	for (uint32 col = 0; col < scanlinesize; col++) {
		buf[col] = 0;
	}
	for (uint32 row = 0; row<imagelength;row++) {
		if (TIFFWriteScanline(out, buf, row, 0)<0) {
			ERROR_EXIT("problem with writescanline");
		}
		
	}
	TIFFClose(out);
}
void pageanalysis(PageDataPtr data) {
	uint32 width = data->width, height = data->height;
	/* find left box for pageno area */
	LPixelrotationdata* ldata =(LPixelrotationdata*)(data->lefttext->data);
	LPixelrotationdata* rdata =(LPixelrotationdata*)(data->righttext->data);
	printf("\n\nDRAW:\n");
	printf("convert o.tif -fill none -stroke black -strokewidth 10 ");
	printdraw(ldata);
	printdraw(rdata);
	printf(" o1.tif\n");
	outtif(data);
	free(data);
}

/* vim: set ts=8 sts=8 sw=8 noet: */
/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 8
 * fill-column: 78
 * End:
 */
