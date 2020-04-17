/********************************************************************
 * Purpose:
 * 	C++ wrapper for C function, uses function getchrlen
 *	
 *
 *	Parameters:
 *		const string&				filepath of .fa file
 *
 *	Note:
 *		Simply converts string to c_str(), and returns value.
 *
 *******************************************************************/
#include <string>
#include <C_Utilities.hpp>

using namespace std;

/*extern "C"{
	int getchrlen(char *filepath);
}*/

int Get_Chr_Len(const string& str){
	return getchrlen((char *)str.c_str());
}
