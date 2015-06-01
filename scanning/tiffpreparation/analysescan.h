#include "tiffiop.h"

typedef struct lscanoptions {
	int reverse_bits;
} LScanOptions;

extern LScanOptions* set_scan_options();

extern void analysescan(const char * filename, TIFF* tif,
			const char* resultdir, LScanOptions* options);






