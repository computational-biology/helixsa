//
// Created by parthajit on 23/2/20.
//

#include "sysparams.h"

void runparams_init(struct runparams* params){
   params->file1[0] = '\0';
   params->file2[0] = '\0';
   params->ssafp = NULL;
   params->helix_only = 'F';
   params->multi_cif  = 'F';
} 


void bpfind_params_init(struct bpfind_params* self){
		  strcpy(self->cifparam, "-dummyval");
		  strcpy(self->accnparam, "-dummyval");
		  strcpy(self->htparam, "-HT");
		  strcpy(self->hdparam, "-dummyval");
		  strcpy(self->hdvalparam, "-dummyval");
		  strcpy(self->chainparam, "-dummyval");
		  strcpy(self->chainvalparam, "-dummyval");
		  strcpy(self->angparam, "-dummyval");
		  strcpy(self->angvalparam, "-dummyval");
		  strcpy(self->chparam, "-dummyval");
		  strcpy(self->sgparam, "-dummyval");
		  strcpy(self->evaltypeparam, "-dummyval");
		  strcpy(self->corparam, "-NOCOR");
		  strcpy(self->nmrparam, "-dummyval");
          strcpy(self->nmrvalparam, "-dummyval");
          
	    }
	    void bpfind_params_print(struct bpfind_params* self, FILE* fp)
	    {


		  fprintf(fp, "\n---------------------------P A R A M S    O P T E D ---------------------\n");
		    fprintf(fp,"PARAM   ACCN: %s\n", self->accn);
		    if(strcmp(self->ext, ".cif") == 0){
		          fprintf(fp,"PARAM   FILE TYPE: mmCIF\n");
		    }else{
		          fprintf(fp,"PARAM   FILE TYPE: PDB\n");
		    }
		    if(strcmp(self->htparam, "-dummyval") == 0){
			      fprintf(fp,"PARAM   HETATM NOT REQUESTED\n");
			}else{
			      fprintf(fp,"PARAM   HETATM REQUESTED\n");
			}
		    if(strcmp(self->hdparam, "-dummyval") == 0){
			      fprintf(fp,"PARAM   HYDROGEN BOND DIST(FOR BASE PAIR)    3.8A (DEFAULT)\n");
			}else{
			      fprintf(fp,"PARAM   HYDROGEN BOND DIST(FOR BASE PAIR)    %sA (REQUESTED)\n",self->hdvalparam);
			}
/*		    if(strcmp(self->angparam, "-dummyval") == 0){
			      fprintf(fp,"PARAM   ANGLE IN DEGREE(FOR BASE PAIR)    120 (DEFAULT)\n");
			}else{
			      fprintf(fp,"PARAM   ANGLE IN DEGREE(FOR BASE PAIR)    %s (REQUESTED)\n",self->angvalparam);
			}*/
		    if(strcmp(self->chparam, "-dummyval") == 0){
			      fprintf(fp,"PARAM   C-H...O/N MEDIATED BASE PAIR   REQUESTED\n");
			}else{
			      fprintf(fp,"PARAM   C-H...O/N MEDIATED BASE PAIR   NOT REQUESTED\n");
			}
		    if(strcmp(self->sgparam, "-dummyval") == 0){
			      fprintf(fp,"PARAM   SUGAR O2' MEDIATED BASE PAIR   REQUESTED\n");
			}else{
			      fprintf(fp,"PARAM   SUGAR O2' MEDIATED BASE PAIR   NOT REQUESTED\n");
			}
			if(strcmp(self->nmrparam, "-dummyval") != 0){
			      fprintf(fp,"PARAM   NMR MODEL REQUESTED: %s\n", self->nmrvalparam);
			}
			if(strcmp(self->chainparam, "-dummyval") != 0){
			      fprintf(fp,"PARAM   NMR MODEL REQUESTED: %s\n", self->chainvalparam);
			}
		
	    }

