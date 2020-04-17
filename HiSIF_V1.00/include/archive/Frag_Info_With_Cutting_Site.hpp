#ifndef FRAGINFOWITHCUTTINGSITE_H
#define FRAGINFOWITHCUTTINGSITE_H

#include <string>
#include <vector>
#include <map>

using std::pair;
using std::string;
using std::vector;
using std::map;

class Frag_Info_With_Cutting_Site
{
public:
	int 	end1_chr;
	int 	end1_pos;
	int 	end1_cuttingSite;
	bool 	end1_strand;
	
	int 	end2_chr;
	int 	end2_pos;
	int 	end2_cuttingSite;
	bool 	end2_strand;

	Frag_Info_With_Cutting_Site():
	end1_chr(0),
	end1_pos(0),
	end1_cuttingSite(0),
	end1_strand(false),
	end2_chr(0),
	end2_pos(0),
	end2_cuttingSite(0),
	end2_strand(false){}
};

#endif
