#include "tiffiop.h"

typedef struct lscanoptions {
	int reverse_bits;
} LScanOptions;

extern LScanOptions* set_scan_options();

extern void analysescan(const char * filename, TIFF* tif,
			const char* resultdir, LScanOptions* options);





/* vim: set ts=8 sts=8 sw=8 noet: */
/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 8
 * fill-column: 78
 * End:
 */
