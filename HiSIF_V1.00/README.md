# HiSIF - HiC Significant Interacting Fragments #

----------------------------------


HiSIF was designed to identify characteristic promoter-distal loops from Hi-C data, some of these loops could be the potential promoter enhancer loops which regulate the gene expression. The tool is named from the fact that it is
designed to find the significant interactions from a given sample of Hi-C reads. HiSIF 
supports all of the available three Hi-C protocols (Hi-C, TCC and in situ Hi-C). HiSIF was implemented by using C++/C with parallel processing being written in C. It has been compiled and run exclusively 
on Linux operating systems. This tool only requires the g++ compiler and a reference genome for HG19 or any other kinds of species. Standalone CERN ROOT C++ framework is used 
to extract fit parameters of the CTS interactions. A small C tool is provided to process the initial data from NCBI.


Please refer to the following link for the usuage:

https://github.com/yufanzhouonline/HiSIF


## Citation: ##
Please cite our paper when you use this tool:

Zhou, Y., Cheng, X., Yang, Y., Li, T., Li, J., Huang, T.H., Jin, V.X. (2020) Genome-wide chromatin interactions identify characteristic promoter-distal loops. Genome Medicine. Aug 12;12:69. https://doi.org/10.1186/s13073-020-00769-8

Thank you.
