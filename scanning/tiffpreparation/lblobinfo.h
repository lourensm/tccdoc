struct lblobinfo;
typedef struct lblobinfo * LBlobinfoPtr;
LBlobinfoPtr init_lblobinfo(int imagelength);
void free_lblobinfo(LBlobinfoPtr info);
void handle_segment_blobs(int startx, int afterx, LBlobinfoPtr info);
void blobinfo_endline(LBlobinfoPtr info, int nexty);
void blobinfo_endimage(LBlobinfoPtr info);
void blobinfo_stats(LBlobinfoPtr info);

/* vim: set ts=8 sts=8 sw=8 noet: */
/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 8
 * fill-column: 78
 * End:
 */
