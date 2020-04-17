#include <cstdio>
#include <fcntl.h>
#include <map>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <utility>
#include <C_Utilities.hpp>
#include <CPP_Utilities.hpp>

/******************************************************************************
 * Purpose:
 * 	Dump a large nested data structure to a file, and then read from the
 * 	file, to save memory space.
 *
 *		Takes this structure, makes 2 low-c structs from it, and then writes each
 *		to a file, respectively. Since the second file will be filled with
 *		vectors, the first struct will contain the key to find the vector, but
 *		also the size of the vector, to know how many elements to read in
 *		from the second file.
 *
 *	Notes:
 *		Key == [chr#, length, sizeof_vector]
 *****************************************************************************/
// map<int, map<int, vector<Frag_Info_With_Cutting_Site> >

using namespace std;

void funcErr(char*, char*, int);
class FragDumper{
	public:
		// two files to contain the information
		int keyfd, valfd;

		// keySize is 2 for the indexes, 1 for the number of elements
		unsigned int keySize, valSize;

		// how many records do we have? Key:value
		unsigned int totalRecords;

		// open files for rw
		int openRW(char *);
		FragDumper();
		
		// save structure to 2 files
		void dump(const map< int, map<int, vector<Frag_Info_With_Cutting_Site> > >&);

		// restore to the structure
		void restore(map<int, map<int, vector<Frag_Info_With_Cutting_Site> > >&);

		// parse the structure before writing
		// void parseKey(const pair< pair< int, pair<int, int> >, pair< int, pair<int, int> > >&, prox_dump*);
		unsigned int parseVal(const vector<Frag_Info_With_Cutting_Site>&, frag_info**);

		// write each to each respective file
		void saveKeys(vector<unsigned int*>&);
		void saveVal(frag_info[], unsigned int);

		// get keys and values
		// void getKey(prox_dump*, unsigned int);
		void getVal(unsigned int[3], vector<Frag_Info_With_Cutting_Site>&);

		// save things to the file
		// int saveThing(map< pair< pair< int, pair<int, int> >, pair< int, pair<int, int> > >, vector<Frag_Info_With_Cutting_Site> > &b);
};


// open two unused tmp files
FragDumper::FragDumper(){
	keyfd = openRW("keydump");
	valfd = openRW("valdump");
	
	keySize = sizeof(int) * 3;
	valSize = sizeof(frag_info);
	totalRecords = 0;
}


void FragDumper::dump(const map< int, map<int, vector<Frag_Info_With_Cutting_Site> > >& in_map){
	map<int, map<int, vector<Frag_Info_With_Cutting_Site> > >::const_iterator it = in_map.begin();
	map<int, map<int, vector<Frag_Info_With_Cutting_Site> > >::const_iterator it_end = in_map.end();

	unsigned int size;
	frag_info *info;
	vector<unsigned int*> keys;
	unsigned int tmp[3];
	for (; it != it_end; it++){
		map<int, vector<Frag_Info_With_Cutting_Site> >::const_iterator otherit = it->second.begin();
		map<int, vector<Frag_Info_With_Cutting_Site> >::const_iterator otherit_end = it->second.end();

		// first part of the key
		tmp[0] = it->first;
		for (; otherit != otherit_end; otherit++){

			// second part of the key
			tmp[1] = otherit->first;

			// start parsing the values
			size = parseVal(otherit->second, &info);
			tmp[2] = size;

			keys.push_back(tmp);

			// write this large structure
			saveVal(info, size);

			free(info);
			totalRecords++;
		}
	}

	// write all the keys
	saveKeys(keys);
}

void FragDumper::restore(map<int, map<int, vector<Frag_Info_With_Cutting_Site> > > &map){
	int i, j;
	unsigned int key[3];
	vector<Frag_Info_With_Cutting_Site> vec;

	lseek(keyfd, 0, SEEK_SET);
	lseek(valfd, 0, SEEK_SET);
	for (i = 0; i < totalRecords; i++){
		// read key, needed for vector size
		r_read(keyfd, (char*)key, keySize);

		// read in the values
		getVal(key, vec);

		map[key[0]][key[1]] = vec;

		vec.clear();
	}
	
}

void FragDumper::saveKeys(vector<unsigned int*>& vec){
	vector<unsigned int*>::const_iterator it = vec.begin();
	vector<unsigned int*>::const_iterator it_end = vec.end();

	for (; it != it_end; it++){
		unsigned int *ptr = *it;
		// printf("Values: %d %d %d\n", ptr[0], ptr[1], ptr[2]);
		r_write(keyfd, (char*)*it, keySize);
	}
}

void FragDumper::saveVal(frag_info frags[], unsigned int size){
	unsigned int k;
	for (k = 0; k < size; k++){
		r_write(valfd, &(frags[k]), sizeof(frag_info));
	}
}


/*
void FragDumper::getKey(int[3] data, unsigned int index){
	// change to proper offset
	if (index > totalRecords || index < 0){
		funcErr("FragDumper::getKey", "file offset is either negative or past end of file", 1);
	}

	lseek(keyfd, index * keySize, SEEK_SET);
	r_read(keyfd, (char*)data, keySize);
}*/

// read the values in from the vector at the current file offset, used for restore
void FragDumper::getVal(unsigned int key[3], vector<Frag_Info_With_Cutting_Site> &vec){
	Frag_Info_With_Cutting_Site frag;
	prox_dump dump;
	frag_info info;

	int k;
	
	// printf("Getting vector values!\n");
	for (k = 0; k < key[3]; k++){
		r_read(valfd, (char*)&info, sizeof(frag_info));
		
		// copy over to the object
		frag.end1_chr = info.end1[0];
		frag.end1_pos = info.end1[1];
		frag.end1_cuttingSite = info.end1[2];
		frag.end1_strand = info.end1[3];
		
		frag.end2_chr = info.end2[0];
		frag.end2_pos = info.end2[1];
		frag.end2_cuttingSite = info.end2[2];
		frag.end2_strand = info.end2[3];

		vec.push_back(frag);
	}
}


/*
void FragDumper::parseKey(const pair< pair< int, pair<int, int> >, pair< int, pair<int, int> > >& pairs, prox_dump *dump){
	dump->one[0] = pairs.first.first;
	dump->one[1] = pairs.first.second.first;
	dump->one[2] = pairs.first.second.second;

	dump->two[0] = pairs.second.first;
	dump->two[1] = pairs.second.second.first;
	dump->two[2] = pairs.second.second.second;
}*/

int FragDumper::openRW(char *filename){
	int lowfd;
	if (filename == NULL){
		funcErr("FragDumper::openRW", "filename was NULL!", 1);
		return -1;
	}
	
	if ((lowfd = open(filename, O_RDWR | O_CREAT | O_TRUNC, 0666)) < 0){
		funcErr("FragDumper::openRW", "unable to open the file specified", 2);
		return -1;
	}


	// successful opening of the file
	return lowfd;
}

unsigned int FragDumper::parseVal(const vector<Frag_Info_With_Cutting_Site>& vec, frag_info **in_info){
	unsigned int size = vec.size();
	(*in_info) = (frag_info*)malloc(sizeof(frag_info) * size);
	frag_info *info = *in_info;

	frag_info *tmp = info;
	
	Frag_Info_With_Cutting_Site tmpFrag;
	vector<Frag_Info_With_Cutting_Site>::const_iterator it = vec.begin();
	vector<Frag_Info_With_Cutting_Site>::const_iterator end_it = vec.end();

	for (; it != end_it; it++, tmp++){

		tmp->end1[0] = it->end1_chr;
		tmp->end1[1] = it->end1_pos;
		tmp->end1[2] = it->end1_cuttingSite;
		tmp->end1[3] = it->end1_strand;

		tmp->end2[0] = it->end2_chr;
		tmp->end2[1] = it->end2_pos;
		tmp->end2[2] = it->end2_cuttingSite;
		tmp->end2[3] = it->end2_strand;
	}
	
	return size;
}


/*
int FragDumper::saveThing(map< pair< pair< int, pair<int, int> >, pair< int, pair<int, int> > >, vector<Frag_Info_With_Cutting_Site> > &b){
	frag_info *info;
	prox_dump dump;

	// size of each vector
	unsigned int size;

	map< pair< pair< int, pair<int, int> >, pair< int, pair<int, int> > >, vector<Frag_Info_With_Cutting_Site> >::iterator it = b.begin();
	map< pair< pair< int, pair<int, int> >, pair< int, pair<int, int> > >, vector<Frag_Info_With_Cutting_Site> >::iterator it_end = b.end();

	// parse each key/value
	for (; it != it_end; it++){
		size = parseVal(it->second, &info);

		// write this large structure
		saveVal(info, size);

		dump.size = size;
		parseKey(it->first, &dump);

		// write this to a file
		saveKey(&dump);

		// increment how many records
		totalRecords++;
	}

	free(info);
}


*/
