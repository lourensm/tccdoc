typedef struct lpixelrotationdata {
	LBox box, maxbox;
	double max_atan_angle;
	int lastpospos;
	int* xvalues;
	int* yvalues;
	double x0, y0;
	const char * filespec;
	const char* resultdir;
	const char* description;
	double xtl, ytl, xtr, ytr, xbl, ybl, xbr, ybr;
	double xlow, ylow, xhi, yhi;
	double alpha;
} LPixelrotationdata;

/* vim: set ts=8 sts=8 sw=8 noet: */
/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 8
 * fill-column: 78
 * End:
 */
