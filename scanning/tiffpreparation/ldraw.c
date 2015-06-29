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


void test_drawing(const char * file) {
	TIFF *out = TIFFOpen(file, "w");
	int sampleperpixel;
	uint32 imageheight = 50, imagewidth = 50;
	uint16 bitspersample;
	TIFFSetField(out, TIFFTAG_IMAGEWIDTH, imagewidth);
	TIFFSetField(out, TIFFTAG_IMAGELENGTH, imageheight);
	TIFFSetField(out, TIFFTAG_SAMPLESPERPIXEL, 1);
	TIFFSetField(out, TIFFTAG_BITSPERSAMPLE, 1);
        TIFFSetField(out, TIFFTAG_COMPRESSION, 4);
	TIFFSetField(out, TIFFTAG_PHOTOMETRIC, 0);
	printf("DEFRPS:%d\n", TIFFDefaultStripSize(out,0));
	TIFFSetField(out, TIFFTAG_ROWSPERSTRIP, TIFFDefaultStripSize(out,0));

	TIFFSetupStrips(out);
	uint32 scanlinesize = TIFFScanlineSize(out);
	unsigned char* buf = _TIFFmalloc(scanlinesize);
	for (uint32 col = 0; col < scanlinesize; col++) {
		buf[col] = 0;
	}
	for (uint32 row = 0; row<imageheight;row++) {
		if (TIFFWriteScanline(out, buf, row, 0)<0) {
			ERROR_EXIT("problem with writescanline");
		}
		
	}
	_TIFFfree(buf);
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
