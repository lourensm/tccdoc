
void error_exit(const char* msg,const char* file, const int line);
#define ERROR_EXIT(MSG) error_exit(MSG, __FILE__, __LINE__)

typedef struct lscanoptions {
	int reverse_bits;
	const char* resultdir;
} LScanOptions;


typedef struct lbox {
	int x1, y1,x2, y2;
} LBox;



/* vim: set ts=8 sts=8 sw=8 noet: */
/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 8
 * fill-column: 78
 * End:
 */
