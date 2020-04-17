#include <cstdio>
#include <map>
#include <utility>
#include "FragDumper.hpp"

using namespace std;
int main(void){

	map<int, map<int, vector<Frag_Info_With_Cutting_Site> > > b;

	vector<Frag_Info_With_Cutting_Site> vec;
	Frag_Info_With_Cutting_Site frag;

	frag.end1_chr = 9;
	frag.end1_pos = 4;
	frag.end1_cuttingSite = 0;
	frag.end1_strand = 1;

	frag.end2_chr = 1;
	frag.end2_pos = 2;
	frag.end2_cuttingSite = 0;
	frag.end2_strand = 1;

	vec.push_back(frag);

	b[0][0] = vec;
	
	// printf("Value in vec: %d\n", vec.front().end1_strand);
	printf("Value: %d\n", b[0][0].front().end1_strand);
	vec.clear();

	
	FragDumper *f = new FragDumper();

	
	f->dump(b); 
	printf("dumped stuff!\n");
	
	b.clear();

	f->restore(b);
	printf("Got stuff!\n");

	vec = b[0][0];
	
	printf("b[0][0] end1:\n\t%d\t%d\t%d...\n", vec.front().end1_chr, vec.front().end1_pos, vec.front().end1_strand);
	printf("b[0][0] end2:\n\t%d\t%d\t%d...\n", vec.front().end2_chr, vec.front().end2_pos, vec.front().end2_strand);

	/*
	f->saveThing(b);


	// getting things!
	vec.clear();
	f->getVal(vec, 0);
	*/
	

	return 0;
}
