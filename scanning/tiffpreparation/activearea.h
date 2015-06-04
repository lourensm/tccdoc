/* REQUIRES #include "tiffcommon.h" */
#ifndef ACTIVEAREA_H
#define ACTIVEAREA_H
typedef struct larea {
	struct larea *next;
	/*double y1f, tophorf, widthf, heightf;*/
	LBox box;
	void* data;
	void (*handlesegment)(int startx, int afterx, int y, void* data);

	void (*handlepixel)(int x, int y, void* data);
	void (*analyse)(void* data);
} LArea;

typedef struct lactive_areas {
	struct lactive_areas *next;
	LArea area;
} LActive_areas;
#endif
