//
// Created by parthajit on 23/2/20.
//

#ifndef CPPMET_SYSPARAMS_H
#define CPPMET_SYSPARAMS_H
#include <stdio.h>
#include <string.h>

struct runparams{
    char file1[512];
    char file2[512];
    FILE* ssafp;
    char helix_only;
    char multi_cif;
    int sec_gap;
    int fasta_gap;
    int sec_mat[128][128];
    int fasta_mat[128][128];
    
};

void runparams_init(struct runparams* params);



struct bpfind_params{
	    char cifparam[512];
	    char accnparam[512];
	    char htparam[512];
	    char hdparam[512];
	    char hdvalparam[512];
	    char chainparam[512];
	    char chainvalparam[512];
	    char angparam[512];
	    char angvalparam[512];
	    char chparam[512];
	    char sgparam[512];
	    char evaltypeparam[512];
	    char corparam[50];
	    char nmrparam[50];
        char nmrvalparam[50];
        char accn[512];
        char ext[20];
};
	    void bpfind_params_init(struct bpfind_params* self);
	    void bpfind_params_print(struct bpfind_params* self, FILE* fp);




#endif //CPPMET_SYSPARAMS_H
