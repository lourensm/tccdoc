/* REQUIRES #include "tiffcommon.h" */
#ifndef ACTIVEAREA_H
#define ACTIVEAREA_H
typedef struct larea {
	struct larea *next;
	/*double y1f, tophorf, widthf, heightf;*/
	LBox box;
	void* data;
	void (*handlesegment)(int startx, int afterx, int y, void* data);
	void (*handle_endline)(int currenty, void* data);
	void (*handlepixel)(int x, int y, void* data);
	void (*analyse)(void* data);
	void (*free)(void* data);
} LArea;

typedef struct lactive_areas * LActive_areasPtr;
LArea* add_active_area(LActive_areasPtr* la); 






#endif


/* vim: set ts=8 sts=8 sw=8 noet: */
/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 8
 * fill-column: 78
 * End:
 */
