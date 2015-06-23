
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include "tiffcommon.h"
#include "activearea.h"
#include "lblobinfo.h"
/* TODO:
    on the way to making blobs accessible: e.g. iterator over (x,y) values
    or sequential access to segments, or...

    This requires that all cells are accessible from the blob.
    A blob has a list of objects.
    An object has a pointer to the list of cells.

    After this, I would like code to determine end of blobs.


accessibility probably not available for cell->object_next...
in particular when one previous line segment is touched by multiple current_y


Characterise the blobs on the page.

Possibly: ignore parts of page and group other blobs per page part:
left page, right page, left pageno, right pageno.

TODO?
- blobinfo in separate file? Probably not.
- blobinfo in active_areas code?
- blobinfo currenty??
- remove y from LScanline?
*/
typedef struct lsegment {
	int min_x;
	int length;
} LSegment;

/*
  How to link the LCell's belonging to one LObject? 
At least there is one chain, up and down between LCells.

Give LObject a pointer to first cell and to lastcell? No need to store y values?
Or, only define an object_next pointer...
NO NO NO: joining multiple objects...
 */
typedef struct lcell {
	struct lcell *next;         /* Link to next cell horizontal */
	struct lcell *object_next;  /* Link to cell vertically belonging to one 
                                       LObject */
	LSegment segment;
	struct lobject* top_object;
} LCell;

typedef struct lblob {
	struct lobject* objects;
	LBox range;
	int id;
	struct lblob* next;
	struct lblob* prev;
	int open_object_count;
} LBlob;

typedef struct lobject {
	int id;
	struct lblob *origin;
	int x;
	int y;
	LCell* top_cell;
	struct lobject *blob_next; 
} LObject;
/*
LBlob: connected pixels define a blob. A blob consists of LObjects:
LObject consists of a list of vertically connected LCells (segments).

LBlobinfo: 
Maintaining the information on LBlobs during buildup of Blobs.
 
last_cell_last_new_object: 
    we have just added a new cell connecting to last_cell_previous_line.
    Could be within handle_segment, or between them.

Change name:
It is the object belonging to the last_cell_previous_line that has been linked to
last_cell_current_line.
linked_cell_obj?

last_cell_previous_line     -> lastcell_prevline
last_cell_current_line      -> lastcell_currline
last_cell_last_new_object   -> lastcell_newobj
 */
typedef struct lscanline {
	int y;
	struct lcell *first;
} LScanline;

typedef struct lblobinfo {
	LCell *lastcell_prevline;
	LCell *lastcell_currline;
	LObject* last_object;
	LScanline* llines;
	LObject* lastcell_newobj;
	int lblobcount;
	int nextblobid;
	int nextobjectid;
	LBlob* last_lblob;
	LBlob* first_lblob;
	int open_blobs;
        const char* description;
	LBox box;
} LBlobinfo;

#undef DEBUGBLOB
#ifdef DEBUGBLOB
#include "analysescantests.c"
#endif
static LBlobinfo* init_lblobinfo(const char* description, LBox box) {
	LScanline* llines;
	int ysize = box.y2-box.y1+1;
	llines = (LScanline*)malloc(ysize*sizeof (LScanline));
	LBlobinfo* res = (LBlobinfo*)malloc(sizeof(LBlobinfo));
	LBlobinfo p = {NULL, NULL, NULL, llines, NULL, 0, 0, 0,
		       NULL, NULL, 0, description, box};
	*res = p;
	for (int i = 0; i < ysize; i++) {
		llines[i].y = i+box.y1;
		llines[i].first = NULL;
	}
	return res;
}

void free_lblobinfo(LBlobinfo* info) {
	free(info->llines);
	free(info);
}
static LCell* new_lcell(int startx, int afterx) {
	LCell* res1 = (LCell*)malloc(sizeof(LCell));
	res1->next = NULL;
	res1->object_next = NULL;
	(res1->segment).min_x = startx;
	(res1->segment).length = afterx - startx;
	res1->top_object = NULL;
	return res1;
}

static int lobject_length(LObject* first) {
	int res = 0;
	while (first != NULL) {
		res++;
		first = first->blob_next;
	}
	return res;
}

static int cmp_object_yx(LObject* first, LObject* second) {
	assert(first->y>=0);
	assert(second->y>=0);
	assert(first->x>=0);
	assert(second->x>=0);
	if (first->y < second->y) return -1;
	if (second->y < first->y) return 1;
	if (first->x < second->x) return -1;
	if (second->x < first->x) return 1;
	ERROR_EXIT("Shouldnt compare same object\n");
	return -1;
}
static LObject* merge_yx_lists(LObject* left, LObject * right, LBlob* blob) {
	LObject* result = NULL;
	LObject* lastcell = NULL;
	int ln;
	assert(left != NULL);
	assert(right != NULL);
	assert(left->origin == blob);
	ln = lobject_length(left) + lobject_length(right);
	if (left != NULL&& cmp_object_yx(left, right) < 0) {
		result = left;
		left = left->blob_next;
		lastcell = result;
	} else {
		result = right;
		lastcell = result;
		assert(right->origin != blob);
		right->origin = blob;
		if (right != NULL) right = right->blob_next;
	}
	assert(lastcell != NULL);
	assert(lastcell->origin == blob);

	while (1) {
		if (left != NULL&& (right == NULL||cmp_object_yx(left, right) < 0)) {
			lastcell->blob_next = left;
			lastcell = left;
			assert(left->origin == blob);
			left = left->blob_next;
		} else if (right !=NULL &&
			   (left == NULL||cmp_object_yx(left, right) > 0)) {
			lastcell->blob_next = right;
			lastcell = right;
			assert(right->origin != blob);
			right->origin = blob;
			right = right->blob_next;
		} else {
			if (left == NULL) {
				lastcell->blob_next = right;
				while (right != NULL) {
					assert(right->origin != blob);
					right->origin = blob;
					right = right->blob_next;	
				}
			} else if (right == NULL) {
				lastcell->blob_next = left;
				while (left != NULL) {
					assert(left->origin == blob);
					left = left->blob_next;	
				}
			} else {
				ERROR_EXIT("merge_yx_lists inconsistency");
			}
			assert(lobject_length(result) == ln);
			return result;
		}
	}
}

/* TODO: standardize x ranges and then join combine_replace_left_box */
static void add_segment_to_blob(int startx, int afterx, int y, LBlob* blob) {
	if (startx < blob->range.x1) {
		blob->range.x1 = startx;
	}
	if (afterx - 1 >  blob->range.x2) {
		blob->range.x2 = afterx - 1;
	}
	if (y < blob->range.y1) {
		blob->range.y1 = y;
	}
	if (y > blob->range.y2) {
		blob->range.y2 = y;
	}
}/* TODO: link the origin blobs? */
/* TODO: delete the other origin */
static void new_object(LCell* cell, LBlob* blob, int y, LBlobinfo* info) {
	LObject* res;
	assert(cell->top_object == NULL);
	assert(cell->object_next == NULL);
	res = (LObject*)malloc(sizeof(LObject));
	res->id = info->nextobjectid++;

	res->blob_next = NULL;
	cell->top_object = res;
	res->top_cell = cell;
	if (blob == NULL) {
		blob = (LBlob*)malloc(sizeof(LBlob));
		info->lblobcount++;
		blob->objects = res;
		blob->open_object_count = 1;
		blob->range.x1 = cell->segment.min_x;
		blob->range.x2 = cell->segment.min_x+cell->segment.length - 1;
		blob->range.y1 = y;
		blob->range.y2 = y;
		blob->id = info->nextblobid;
		blob->prev = info->last_lblob;
		blob->next = NULL;
		info->open_blobs++;
		if (info->last_lblob == NULL) {
			info->first_lblob = blob;
		} else {
			assert(info->last_lblob->next == NULL);
			info->last_lblob->next = blob;
		}
		info->last_lblob = blob;
		info->nextblobid++;
	} else {
		add_segment_to_blob(cell->segment.min_x,
				    cell->segment.min_x + cell->segment.length, y,
				    blob);
		/* Inefficient.. could define tail of objects */
		blob->objects = merge_yx_lists(blob->objects, res, blob);
		blob->open_object_count++;
	}
	res->y = y;
	res->x = cell->segment.min_x;
	res->origin = blob;
	info->last_object = res;
}

static void combine_replace_left_box(LBox *left, LBox* right) {
	if (right->x1 < left->x1) left->x1 = right->x1;
	if (right->x2 > left->x2) left->x2 = right->x2;
	if (right->y1 < left->y1) left->y1 = right->y1;
	if (right->y2 > left->y2) left->y2 = right->y2;
}

static void merge_top_objects_first(LObject* left_top, LObject* right_top,
				    LBlobinfo* info) {
	LBlob* toremove = right_top->origin;
	assert(cmp_object_yx(left_top, right_top)<0);
	int new_count = left_top->origin->open_object_count
		      + right_top->origin->open_object_count-1;
	left_top->origin->objects = merge_yx_lists(left_top->origin->objects,
						   right_top->origin->objects,
						   left_top->origin);
	combine_replace_left_box(&(left_top->origin->range),
				 &(right_top->origin->range));
	left_top->origin->open_object_count = new_count;
	info->lblobcount--;
	
	if (toremove->next != NULL) {
		toremove->next->prev = toremove->prev;
	} else {
		assert(toremove == info->last_lblob);
		info->last_lblob = toremove->prev;
	}
	if (toremove->prev!= NULL) {
		toremove->prev->next = toremove->next;
	} else {
		assert(toremove == info->first_lblob);
		info->first_lblob = toremove->next;
	}
	info->open_blobs--;
	free((void*)toremove);
}
/*
Which blob to keep? The left-top most.
 */

static void merge_top_objects(LObject* left_top, LObject* right_top,
			      LBlobinfo* info) {
	if (left_top->origin == right_top->origin) {
		/* check for embedded/surrounded objects */
		left_top->origin->open_object_count --;
		if(left_top->origin->open_object_count <= 0) {
			printf("ooc:%d, o1:%d o2:%d\n",
			       left_top->origin->open_object_count,
			       left_top->id, right_top->id);
			ERROR_EXIT("merge open_object_count 0");;
		}
	} else {
		if (cmp_object_yx(left_top, right_top)<0) {
			merge_top_objects_first(left_top, right_top, info);
		} else {
			merge_top_objects_first(right_top, left_top, info);
		}
	}
}

/*
What do I need to do with blobs?

- Detect the pagenumber
  
- Use that to further improve the alignment of the page

How to detect the pagenumber?
- define the area where to look
- find blobs there
- calculate the contour of the blobs
- match with other result of page alignment

Detect blob in area:
- check at end of blob
- maintain datastructure for blobs which makes access easy: too complex, interesting, but too complex.

So, check at end of blob or check each blob for overlap with all active_areas.

INTERMEZZO:
define rectangle of blob in rotated coordinate system, using the atan_angle.
define rectangle of blobs in rotated coordinate system, using the atan_angle.


Ok, store blobs, make them accessible at end.
Make contour of blob accessible.

Would be better if we store or handle blobs as soon as we know they are done?

Put blobs in linked list... 
 */


static void decrease_blob_object_count(LCell* cell, LBlobinfo* info) {
	assert(cell->top_object->origin->open_object_count > 0);
	cell->top_object->origin->open_object_count --;
	if (cell->top_object->origin->open_object_count == 0) {
		info->open_blobs--;
	}
}



/*
TODO: cleanup of documentation.
Focus: long previous_line segment, multiple consecutive new_segments.
Then, only the first new_cell should get to be linked directly, the others should 
introduce new objects.

We need to have info on how the last new_cell got linked to which 
lastcell_prevline
We have lastcell_prevline, we need the new_cell-object to which it has most 
recently been linked, NULL if not by lastcell_prevline.

If we now have to skip non overlapping, then last_object becomes NULL

startx, afterx : the segment that needs to be incorporated into
info:  LBlobinfo* info

This implies that all segments should be part of a blob. They should be
accessible from the blob.

Segments are linked to a blob such that each segment is linked to parent
segments, until some top_object. A blob then has a list of top_objects that
make up its segments.

With each blob we maintain an open_object_counter so as to be able to detect when the 
blob
is finished, i.e. all of its objects and segments no longer have a connection downwards.

handle_segment_blobs links the current segments to segments of the previous line (currenty - 1)

situations:
*current segment is connected to no previous_line segment.
* current segment is connected to exactly one previous_line segment which was not earlier
connected to some current line segment.
* multiple previous_line segments are connected to current_segment 

info->lastcell_newobj :
We just linked a current cell to a previousline cell, which is still available for linking to
current cell. We need not decrease its blob->object_count as that has a link already through that
current cell.
at entry:
assert(info->lastcell_newobj == NULL ||
			       info->lastcell_newobj->origin ==
			       info->lastcell_prevline->top_object->origin);
 */
void handle_segment_blobs(int startx, int afterx, int y, LBlobinfo* info) {
	int first = 1;
	
	LCell* newcell = new_lcell(startx, afterx);
	if (info->lastcell_currline == NULL) {
		info->llines[y-info->box.y1].first = newcell;
		info->lastcell_currline = newcell;
	} else {
		info->lastcell_currline->next = newcell;
	}
	info->lastcell_currline = newcell;
	/* skip non-overlapping segments previous line */
	while (info->lastcell_prevline != NULL &&
	       startx >= info->lastcell_prevline->segment.min_x +
	       info->lastcell_prevline->segment.length) {
		if (info->lastcell_newobj == NULL) {
			decrease_blob_object_count(info->lastcell_prevline, info);
		}
		info->lastcell_prevline = info->lastcell_prevline->next;
		info->lastcell_newobj = NULL;
	}
	/* The condition of overlap subly different in while condition and return */
	while (info->lastcell_prevline!=NULL&& info->lastcell_prevline->segment.min_x <= afterx - 1) {
		if (first == 1) {  /* Need to assign new_cell to object */
			if (info->lastcell_newobj == NULL) {
				/* link to existing object */
				newcell->top_object  =
					info->lastcell_prevline->top_object;
				assert(info->lastcell_prevline->object_next == NULL);
				info->lastcell_prevline->object_next = newcell;
				add_segment_to_blob(startx, afterx, y,
						    newcell->top_object->origin);

			}  else {
				new_object(newcell,
					   info->lastcell_newobj->origin,
					   y, 
					   info);
			}
			first = 0;
		} else {
			merge_top_objects(newcell->top_object,
					  info->lastcell_prevline->top_object,
					  info);
		}
		info->lastcell_newobj = newcell->top_object;
		if (info->lastcell_prevline->segment.min_x
		    + info->lastcell_prevline->segment.length <= afterx - 1) {
			info->lastcell_prevline = info->lastcell_prevline->next;
			info->lastcell_newobj = NULL;
		} else return;
	}
	if (first) {
		/* NEW OBJECT, no link to earlier segment  */
		new_object(newcell, NULL, y, info);
		info->lastcell_newobj = NULL;
	}
}

void blobinfo_endline(LBlobinfo* info, int currenty) {
	while (info->lastcell_prevline!= NULL) {
		if (info->lastcell_newobj == NULL) {
			decrease_blob_object_count(info->lastcell_prevline,
						   info);
		} 
		info->lastcell_prevline = info->lastcell_prevline->next;
		info->lastcell_newobj = NULL;
	}
	info->lastcell_newobj = NULL;
	info->lastcell_currline = NULL;
	info->lastcell_prevline = info->llines[currenty-info->box.y1].first;
#ifdef DEBUGBLOB
	test_invariant_end(info);
#endif
}
void blobinfo_endimage(LBlobinfo* info) {
	if (info->lastcell_prevline != NULL) {
		while (info->lastcell_prevline != NULL) {
			decrease_blob_object_count(info->lastcell_prevline,info);
			info->lastcell_prevline = info->lastcell_prevline->next;
		}
	}	
}

void blobinfo_stats(LBlobinfo* info) {
	printf("Defined %d blobs (%d blobs still open)\n", info->lblobcount,
	       info->open_blobs);
}




/* LAction_area stuff */






void handleblobsegment(int startx, int afterx, int y, void* data) {
	handle_segment_blobs(startx, afterx, y, (LBlobinfo*)data);
}


static void handleblobendline(int currenty, void* info) {
	blobinfo_endline((LBlobinfo*)info, currenty);
}

static void analyseblob(void* data) {
	blobinfo_endimage((LBlobinfo*)data);
	blobinfo_stats((LBlobinfo*)data);
}

static void freeblob(void*data) {
	free_lblobinfo((LBlobinfo*)data);
}

/* TODO: for now we handle only one blobspace... 
 *     Somehow the segments are linked horizontally, do not really
 *     understand my own code now.
 */ 


void setupblobarea(LArea* larea, LBox box, const char* where) {
	/*	static int called = 0;
	if (called++ > 0) ERROR_EXIT("allow one blobspace for now");
	*/
	larea->next = NULL;
	larea->box = box;
	larea->data = (void*)init_lblobinfo(where, box);
	larea->handlepixel = NULL;
	larea->handlesegment = &handleblobsegment;
	larea->handle_endline = &handleblobendline;
	larea->analyse = &analyseblob;
	larea->free = &freeblob;
	assert(larea->analyse != NULL);
	assert(larea->data != NULL);
}





/* vim: set ts=8 sts=8 sw=8 noet: */
/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 8
 * fill-column: 78
 * End:
 */
