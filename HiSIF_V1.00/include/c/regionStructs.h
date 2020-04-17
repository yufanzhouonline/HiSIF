/**********************************************************
 * Define C structs used during reading process
 *********************************************************/
#if defined (__cplusplus)
extern "C"{
#endif

void *readInteractingRegionsThread(void *args);

// single interaction region
struct interRegion{
	char strand;
	int chr;
	int pos;
	int cuttingSite;
};


// paired interaction regions
struct interRegionPair{
	struct interRegion end1;
	struct interRegion end2;
};


// track next slot, start and end slots
struct regionIndex{
	struct interRegionPair *start;
	struct interRegionPair *end;
	struct interRegionPair *cur;
};

/************thread_args*****************/
struct fragArgs{
		// thread id
	int id;

		// fd to read from, could be pipe
	int readfd;
};

#if defined (__cplusplus)
}
#endif
