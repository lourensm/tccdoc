/*
 *  buf : so as to be able to do TIFFReadScanline(tif, buf, row, (tsample_t)0);
 */
typedef struct pagedata* PageDataPtr;

PageDataPtr setuppageareas(LActive_areasPtr* la,
		    const char* filename, TIFF* tif, unsigned char* buf,
		    LScanOptions* options);

void pageanalysis(PageDataPtr data);



/* vim: set ts=8 sts=8 sw=8 noet: */
/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 8
 * fill-column: 78
 * End:
 */
