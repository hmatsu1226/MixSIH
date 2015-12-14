#include <cstdio>
#include <unistd.h>
#include <vector>
#include <string>
#include <map>
#include "MixSIH.h"

using namespace std;

int main(int argc, char* argv[]){
	double alpha = 0.1;
	
	int ch;
	extern char *optarg;
	extern int optind, opterr;
	
	while((ch = getopt(argc, argv, "a:")) != -1){
		switch(ch){
		case 'a':
			alpha = atof(optarg);
			break;
		case ':':
			printf("invalid option\n");
			return 1;
		case '?':
			printf("invalid option\n");
			return 1;
		}
	}
	
	FILE *fin, *fout1, *fout2;
	if((fin=fopen(argv[optind], "r")) == NULL){
		printf("fopenerror\n");
		return 1;
	}
	if((fout1=fopen(argv[optind+1], "w")) == NULL){
		printf("fopenerror\n");
		return 1;
	}
	if((fout2=fopen(argv[optind+2], "w")) == NULL){
		printf("fopenerror\n");
		return 1;
	}

	vector<Haplotype> hap;
	set_fragment(fin, hap);

	for(unsigned int i=0; i<hap.size(); i++){
		hap[i].Assign_fragment_to_SNP();
		hap[i].Init(alpha);
		hap[i].EM();
	}
		
	for(unsigned int i=0; i<hap.size(); i++){
		hap[i].Output_profile(fout1);
		hap[i].Output_haplotype(fout2);
	}
	
	return 0;
}
