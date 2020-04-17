// main program for running the processing of the mapped files
#include <stdio.h>

int main(int argc, char *argv[]){
	if (argc != 4){
		fprintf(stderr, "%s: <indir> <outdir> <-r or -t>\n-t\t\ttraditional 6-column file format\n-r\t\trao format\n\nGenerate chromosome by chromosome files in the specified directory\n\n", argv[0]);
		return -1;
	}

	// input directory, output directory, type of input files
	frag_proc(argv[1], argv[2], argv[3]);
}
