#include <cstring>
#include <string>
#include <stdlib.h>

using namespace std;
// given a chr20* string, extract 20
// if this contains a capital X or Y, extract that, count as respective
// 23 and 24
int extractChrNum(char *str){

	string stringystr = str;
	// define regular expressions
	// fprintf(stderr, "%d\n%d\n%d\n%d", stringystr.find("X"), stringystr.find("x"), stringystr.find("y"), stringystr.find("Y"));
	if (stringystr.find("X") != string::npos || stringystr.find("x") != string::npos)
		return 23;
	if (stringystr.find("Y") != string::npos || stringystr.find("y") != string::npos)
		return 24;

	string tmpstring;
	
	// otherwise extract the number
	unsigned int i;
	for (i = 0; i < strlen(str); i++){

		// hit the number, gather the characters
		if (isdigit(str[i]+0)){
			// fprintf(stderr, "Found a number @ %c\n", str[i]);
			tmpstring.append(1, str[i]);

			// check for another digit
			if (isdigit(str[i+1] + 0)){
				// fprintf(stderr, "Found a number @ %c\n", str[i+1]);
				tmpstring.append(1, str[i+1]);
				return atoi(tmpstring.c_str());
			}
			else
				return atoi(tmpstring.c_str());
		}
	}

		// no numbers, and no X/x Y/x in this string
	return -1;
}
