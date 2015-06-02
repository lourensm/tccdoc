
NON COMPILING OLD DEBUGGING FUNCTIONS CUT OUT OF analysescan.c
#define BLOBN 80000


static void pre_test(LObject* left_top, const char*where, ...) {
	return;
	if (left_top->origin->open_object_count <= 0) {
		printf("PRETEST:");
		va_list args;
		va_start(args, where);
		vfprintf(stdout, where, args);
		va_end(args);
		printf("\n");
		printf("oid:%d (x:%d,y:%d), lcellhasnext:%p, blobid:%d\n", left_top->id,
		       left_top->x, left_top->y,
		       left_top->top_cell->object_next,
		       left_top->origin->id);
		ERROR_EXIT("merge left object_count 0");
	}
}

static void print_blob_counts(LBlobtest* blobt, int currenty, int start_x,
			      LBlob * suspect, const char* where, int deb) {
	for (int i=0;i<BLOBN;i++) {
		LBlob* blob = blobt[i].blob;
		if (blob != NULL) {
			if (blob->open_object_count != blobt[i].open_segments) {
				if (deb) printf("AT %s\n", where);
				if (deb) printf("blobid=%d, oc:%d, sc:%d, y(next):%d, startx:%d\n",
				       blob->id, blob->open_object_count,
				       blobt[i].open_segments, currenty, start_x);
				if (suspect!=NULL) {
					if (blob == suspect) {
						ERROR_EXIT("suspect case invariant");
					} else {
						ERROR_EXIT("incomplete suspect case invariant");	
					}
				} else {
					ERROR_EXIT("unexplained open_count mismatch");
				}
			}
			
		}
	}
	if (deb) printf("TEST INVARIANT OK\n");
}

void print_cell(LCell* cell, LBlobtest* blobt, const char* sp, int first) {
	int blobid = cell->top_object->origin->id;
	char str[15];
	const char* scount;
	if (blobt[blobid].blob == NULL) {
		scount = "N";
	} else {
		sprintf(str, "%d", blobt[blobid].open_segments);
		scount = str;
	}
	if (first) printf("||");
	printf("x:(%d,#%d), b:%d(%d,%s) %s",cell->segment.min_x,
	       cell->segment.length, cell->top_object->origin->id,
	       cell->top_object->origin->open_object_count,
	       scount, sp
	       );
}

void test_invariant_end(LBlobinfo* info) {
	
	assert(info->last_cell_last_new_object == NULL);
	assert(info->last_cell_current_line == NULL);
	assert((info->currenty == 0&&info->last_cell_previous_line ==NULL)||
	       info->last_cell_previous_line ==
			info->llines[info->currenty-1].first);
	LBlobtest blobt[BLOBN];
	int deb = 0;
	if (deb) printf("INVARIANT END\n");
	for (int i=0;i<BLOBN;i++) {
		blobt[i].blob = NULL;
		blobt[i].open_segments = 0;
	}
	LCell* cell = info->last_cell_previous_line;
	int first = 1;
	while (cell != NULL) {
		int id = cell->top_object->origin->id;
		if (id < BLOBN) {
			blobt[id].blob = cell->top_object->origin;
			blobt[id].open_segments++;
		} else {
			ERROR_EXIT("cannot test more than BLOBN");
		}
		if (deb) print_cell(cell, blobt, "", first);
		first = 0;
		cell = cell->next;
	}
	print_blob_counts(blobt, info->currenty, -1, NULL, "test_invariant_end",
			  deb);
}



void test_invariant(LBlobinfo* info, int beforestartx, const char* where) {

	LBlobtest blobt[BLOBN];
	for (int i=0;i<BLOBN;i++) {
		blobt[i].blob = NULL;
		blobt[i].open_segments = 0;
	}
	int deb = 1;
	int first = 1;
	LCell * cell1 = info->currenty == 0?NULL:
		info->llines[((info->currenty) - 1)].first;
	if (deb) printf("CHAIN:previous_line upto last_cell_previous_line:\n");
	while (1) {
		if (cell1 == NULL) break;
		if (deb) print_cell(cell1, blobt, (cell1 == info->last_cell_previous_line&&
						   info->last_cell_last_new_object !=NULL)?"NW":"  ",
				    first);
		first = 0;
		if (cell1 == info->last_cell_previous_line) break;
		cell1 = cell1->next;
	}
	LCell* cell = 
		info->llines[(info->currenty)].first;
	LBlob * suspect = NULL;
	if (deb) printf("\nCHAIN:current_line:\n");
	first = 1;
	while (1) {
		if (cell == NULL) break;
		int id = cell->top_object->origin->id;
		if (id < BLOBN) {
			blobt[id].blob = cell->top_object->origin;
			if (info->last_cell_last_new_object == NULL||
			    cell != info->last_cell_current_line)  {			  
				blobt[id].open_segments++;
			} else {
				assert(suspect == NULL);
				suspect = info->last_cell_last_new_object->origin;
			}
		} else {
			ERROR_EXIT("cannot test more than BLOBN");
		}
		if (deb) print_cell(cell, blobt, suspect!=NULL?"SS":"  ", first);
		first = 0;
		if (cell == info->last_cell_current_line) break;
		cell = cell->next;
	}
	if (deb) {
		if (cell1 == NULL||cell1->next == NULL) {
			printf("\nCHAIN:NO previous_line NEXT");	
		} else {
			first = 1;
			printf("\nCHAIN:previous_line NEXT:\n");
			while (1) {
				cell1 = cell1->next;
				if (cell1 == NULL) break;
				print_cell(cell1, blobt, "  ", first);
				first = 0;
			}
		}
	}
	if (deb) printf("\nCHAIN:end\n");
	print_blob_counts(blobt, info->currenty, beforestartx, suspect, where, deb);
	printf("TEST INVARIANT OK\n");
}


void test_all_open(LCell* first, const char* where, ...) {
	return;
	printf("All cells:\n");
	while (first != NULL) {
		printf("    x:(%d #%d), blobid=%d\n",first->segment.min_x,
		       first->segment.length, first->top_object->origin->id);
		if (first->top_object->origin->open_object_count <= 0) {
			printf("test_all_open:");
			va_list args;
			va_start(args, where);
			vfprintf(stdout, where, args);
			va_end(args);
			printf("\n");
			ERROR_EXIT("some closed cell");
		}
		first = first->next;
	}
}
