#include <CPP_Utilities.hpp>

// utility to copy values
//
void fragCopy(cBundle& bundle, const Frag_Info_With_Cutting_Site& frag){
	bundle.frag.end1_chr = frag.end1_chr;
	bundle.frag.end1_pos = frag.end1_pos;
	bundle.frag.end1_cuttingSite = frag.end1_cuttingSite;
	bundle.frag.end1_strand = frag.end1_strand;

	bundle.frag.end2_chr = frag.end2_chr;
	bundle.frag.end2_pos = frag.end2_pos;
	bundle.frag.end2_cuttingSite = frag.end2_cuttingSite;
	bundle.frag.end2_strand = frag.end2_strand;
}

void printcBundle(cBundle& bundle){
	fprintf(stderr, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", bundle.chr, bundle.length,
		bundle.frag.end1_chr, bundle.frag.end1_pos, bundle.frag.end1_cuttingSite,
		bundle.frag.end1_strand, bundle.frag.end2_chr, bundle.frag.end2_pos,
		bundle.frag.end2_cuttingSite, bundle.frag.end2_strand);
}

void cbundleCopy(Frag_Info_With_Cutting_Site& frag, cBundle& bundle){
	// need to make a copy, this is immutable
	frag.end1_chr = bundle.frag.end1_chr;
	frag.end1_pos = bundle.frag.end1_pos;
	frag.end1_cuttingSite = bundle.frag.end1_cuttingSite;
	frag.end1_strand = bundle.frag.end1_strand;

	frag.end2_chr = bundle.frag.end2_chr;
	frag.end2_pos = bundle.frag.end2_pos;
	frag.end2_cuttingSite = bundle.frag.end2_cuttingSite;
	frag.end2_strand = bundle.frag.end2_strand;
}
