// gather the chromosome number
int getChrNum(char *string){
	// fprintf(stderr, "%s == ", string);
	if (string[0] == 'x' || string[0] == 'X')
		return 23;
	else if (string[0] == 'y' || string[0] == 'Y')
		return 24;
	else
		return atoi(string);
}
