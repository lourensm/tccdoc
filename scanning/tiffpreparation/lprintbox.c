#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <stdbool.h>
#include "tiffcommon.h"
#include "activearea.h"
#include "lprintbox.h"


/*
 * TODO
 *  the code really needs debugging, things like correct comparison, nextx, <=8
 *  need to be able to usefully specify part of the whole screen and print out the
 *  selection only, not the 0's before box.x1. BUT SEEMS OK. 
 *  NEED command line option
 */
typedef enum {BIN, HEX} LOutputtype;

typedef struct lprintdata {
	LBox box;
	LOutputtype outtype;
	int nextx;
	short bit;
	unsigned char currentchar;
} LPrintdata;

static void printfree(void* data) {
	free(data);
}

static void singlepixelhex(LPrintdata* pdata, int pixel) {
	assert(pixel ==0||pixel ==1);
	pdata->bit++;
	pdata->currentchar = 2* pdata->currentchar + pixel;
	pdata->nextx++;
	if (pdata->bit >= 8) {
		printf("%02X ", pdata->currentchar);
		pdata->currentchar = 0;
		pdata->	bit = 0;
	}
}
static void singlepixelbin(LPrintdata* pdata, int pixel) {
	assert(pixel ==0||pixel ==1);
	if (pixel == 0) {
		printf(" ");
	} else {
		printf("x");
	}
	pdata->nextx++;
}
static void singlepixel(LPrintdata* pdata, int pixel) {
	switch (pdata->outtype) {
	case BIN:
		singlepixelbin(pdata, pixel);
		break;
	case HEX:
		singlepixelhex(pdata, pixel);
		break;
	default:
	  ERROR_EXIT("Unhandled printtype");
	}
}
static void pixelprint(int x, int y, void* data) {
	(void)y;
	LPrintdata* pdata = (LPrintdata*)data;
	assert(x >= pdata->box.x1);
	assert(x <= pdata->box.x2);
	while (pdata->nextx < x) {
		singlepixel(pdata, 0);
	}
	singlepixel(pdata, 1);
}

static void endlineprint(int y, void* data) {
	(void)y;
	LPrintdata* pdata = (LPrintdata*)data;
	switch (pdata->outtype) {
	case BIN:
		while (pdata->nextx < pdata->box.x2) {
			singlepixelbin(pdata, 0);
		}
		break;
	case HEX:
		while (pdata->nextx < pdata->box.x2) {
			singlepixelhex(pdata, 0);
		}
		while (pdata->bit > 0) {
			singlepixelhex(pdata, 0);
		}
		break;
	default:
		ERROR_EXIT("Unhandled printtype");
	}
}

void setupprintarea(LActive_areasPtr* la, LBox *box, LOutputtype outtype) {
	LArea* printarea = add_active_area(la);
	printarea->box = *box;
	LPrintdata pd = {*box, outtype, box->x1, 0, 0};
	printarea->data = malloc(sizeof(LPrintdata));
	*((LPrintdata*)printarea->data) = pd;
	printarea->handlepixel = &pixelprint;
	printarea->handlesegment = NULL;
	printarea->handle_endline = &endlineprint;
	printarea->free = &printfree;
}

void setupprintareafromoptions(LActive_areasPtr* la, LScanOptions* options,
				int imageheight, int imagewidth) {
	if (options->hexall) {
		LBox box = {0, 0, imagewidth, imageheight};
		setupprintarea(la, &box, HEX);
	}
	if (options->binall) {
		LBox box = {0, 0, imagewidth, imageheight};
		setupprintarea(la, &box, BIN);
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
