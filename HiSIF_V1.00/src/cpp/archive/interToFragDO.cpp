/********************************************************************
 * Perform conversion b/w a interRegionPair struct and a
 * Frag_Info_With_Cutting_Site DO
 *******************************************************************/
#include <stdlib.h>

#include "Frag_Info_With_Cutting_Site.hpp"

// interRegionPair definition
extern "C"{
	#include "lowstructs.h";
}

// convert from struct to DO
struct interRegionPair *interToFragDO(Frag_Info_With_Cutting_Site& frag){
	struct interRegionPair *newInter = 
		(struct interRegionPair)malloc(sizeof(struct interRegionPair));
	
	// copy values over */
	newInter->end1.chr = frag->end1_chr;
	newInter->end1.pos = frag->end1_pos;
	newInter->end1.currintSite = frag->end1_cuttingSite;
	newInter->end1.strand = (char)frag->end1_strand;

	newInter->end2.chr = frag->end2_chr;
	newInter->end2.pos = frag->end2_pos;
	newInter->end2.currintSite = frag->end2_cuttingSite;
	newInter->end2.strand = (char)frag->end2_strand;

	return newInter;
}

// convert from DO to struct
Frag_Info_With_Cutting_Site& fragDOToInter(struct interRegionPair *pair){
	Frag_Info_With_Cutting_Site newFrag = 
		(Frag_Info_With_Cutting_Site)malloc(sizeof(Frag_Info_With_Cutting_Site));

	newFrag->end1_chr = pair->end1.chr;
	newFrag->end1_pos = pair->end1.pos;
	newFrag->end1_cuttingSite = pair->end1.cuttingSite;
	newFrag->end1_strand = (bool)pair->end1.strand;

	newFrag->end2_chr = pair->end2.chr;
	newFrag->end2_pos = pair->end2.pos;
	newFrag->end2_cuttingSite = pair->end2.cuttingSite;
	newFrag->end2_strand = (bool)pair->end2.strand;

	return newFrag;
}
