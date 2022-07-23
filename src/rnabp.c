/*
 * =====================================================================================
 *
 *       Filename:  rnabp.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  Thursday 31 March 2022 06:25:34  IST
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  PARTHAJIT ROY (PR), roy.parthajit@gmail.com
 *   Organization:  The University of Burdwan
 *
 * =====================================================================================
 */

#include "rnabp.h"

int is_canonical(char res1, char res2, char* bpname){
       if(strcmp(bpname, "W:WC") != 0) return FALSE;

       if(res1 == 'G' && res2 == 'C') return TRUE;
       if(res1 == 'C' && res2 == 'G') return TRUE;

       if(res1 == 'G' && res2 == 'U') return TRUE;
       if(res1 == 'U' && res2 == 'G') return TRUE;

       if(res1 == 'A' && res2 == 'U') return TRUE;
       if(res1 == 'U' && res2 == 'A') return TRUE;

       return FALSE;
 }

int is_canbp(char* bp, char* res1, char* res2){
      int b1, b2;
      if(strcmp(bp, "W:WC") == 0){
	    if(strcmp(res1, "G") == 0 || strcmp(res1, "DG") == 0){
		  b1 = 2;
	    }else if(strcmp(res1, "A") == 0 || strcmp(res1, "DA") == 0){
		  b1 = 1;
	    }else if(strcmp(res1, "C") == 0 || strcmp(res1, "DC") == 0){
		  b1 = 4;
	    }else if(strcmp(res1, "U") == 0 || strcmp(res1, "DT") == 0){
		  b1 = 3;
	    }else{
		  b1 = 10;
	    } 
	    
	    if(strcmp(res2, "G") == 0 || strcmp(res2, "DG") == 0){
		  b2 = 2;
	    }else if(strcmp(res2, "A") == 0 || strcmp(res2, "DA") == 0){
		  b2 = 1;
	    }else if(strcmp(res2, "C") == 0 || strcmp(res2, "DC") == 0){
		  b2 = 4;
	    }else if(strcmp(res2, "U") == 0 || strcmp(res2, "DT") == 0){
		  b2 = 3;
	    }else{
		  b2 = 10;
	    } 

	    int val = (b1-b2) * (b1 - b2);
	    if( val == 4 ){
		  return TRUE;
	    }else{
		  return FALSE;
	    }
      }else{
	    FALSE;
      }
}

void rnabp_fprint_summary(struct nucbp* self, int numres, char* accn, FILE* fp)
{
      int cano = 0;
      int noncano = 0;
      int multiplet = 0;
      int chain = 0;
      int stdbase = 0;
      int modibase = 0;
      fprintf(fp, "ACCN:%s\n", accn);
      fprintf(fp, "TOTAL RESIDUES:%d\n", numres);
      char bp1[5];
      char bp2[5];
      char bpname[5];
      int othindx;
      for(int i=0; i<numres; ++i){
	    strcpy(bp1, self[i].resname);
	    if(self[i].numbp > 1){
		  multiplet++;
	    }
	    for(int j=0; j<self[i].numbp; ++j){
		  othindx = self[i].oth_base_index[j];
		  strcpy(bp2, self[othindx].resname);
		  strcpy(bpname, self[i].name[j]);
		  if(is_canbp(bpname, bp1, bp2) == TRUE){
			cano++;
		  }else{
			noncano++;
		  }

	    }

	    
      }
      fprintf(fp, "TOTAL CANONICAL BASE-PAIR    : %d\n", cano);
      fprintf(fp, "TOTAL NON CANONICAL BASE-PAIR: %d\n", noncano);
      fprintf(fp, "TOTAL MULTIPLET              : %d\n", multiplet);
}

char get_lw_edge_symb(char e){
      switch(e){
	    case 'W':
	    case 'w':
	    case '+': return 'I';
	    case 'S':
	    case 's':
	    case 'z': return 'S';
	    case 'H':
	    case 'h':
	    case 'g': return 'H';
	    default: return 'X';
      }
}
char get_lw_edge(char e){
      switch(e){
	    case 'W':
	    case 'w':
	    case '+': return 'W';
	    case 'S':
	    case 's':
	    case 'z': return 'S';
	    case 'H':
	    case 'h':
	    case 'g': return 'H';
	    default: return 'X';
      }
}
void to_leontis_westhof_symb(char* lw, char* nuparm)
{
      lw[0] = nuparm[3] == 'C'? '+': '-';
      lw[1] = get_lw_edge_symb(nuparm[0]);
      lw[2] = get_lw_edge_symb(nuparm[2]);
      lw[3] = '\0';
      
}

void to_leontis_westhof(char* lw, char* nuparm)
{
      lw[0] = tolower(nuparm[3]);
      lw[1] = get_lw_edge(nuparm[0]);
      lw[2] = get_lw_edge(nuparm[2]);
      lw[3] = '\0';
      
}

void rnabp_fprint_vertical(struct nucbp* self, int numres, char* accn, FILE* fp)
{



      char* base1 = (char*) malloc ((numres+1) * sizeof(char));
      char* base2 = (char*) malloc ((numres+1) * sizeof(char));
      char* edge1 = (char*) malloc ((numres+1) * sizeof(char));
      char* edge2 = (char*) malloc ((numres+1) * sizeof(char));
      char* ornt1 = (char*) malloc ((numres+1) * sizeof(char));
      fprintf(fp, "ACCN: %s\n", accn);
      fprintf(fp, "TOTAL RESIDUES: %d\n", numres);
      int k = 0;
      char lw[5];
      for(int i=0; i<numres; ++i){
	    base1[i] = self[i].resclass;
            if(self[i].numbp > 0){
                  k = self[i].oth_base_index[0];
		  base2[i] = self[k].resclass;
		  to_leontis_westhof_symb(lw, self[i].name[0]);
		  edge1[i] = lw[1];
		  edge2[i] = lw[2];
		  ornt1[i] = lw[0];

	    }else{

		  base2[i] = ' ';
		  ornt1[i] = ' ';
		  edge1[i] = ' ';
		  edge2[i] = ' ';

	    }
      }
      base1[numres] = '\0';
      base2[numres] = '\0';
      ornt1[numres] = '\0';
      edge1[numres] = '\0';
      edge2[numres] = '\0';

      fprintf(fp, "%s\n", base1);
      fprintf(fp, "%s\n", edge1);
      fprintf(fp, "%s\n", ornt1);
      fprintf(fp, "%s\n", edge2);
      fprintf(fp, "%s\n", base2);
}
void rnabp_fprint(struct nucbp* self, int numres, char* accn, FILE* fp)
{




      fprintf(fp, "ACCN: %s\n", accn);
      fprintf(fp, "TOTAL RESIDUES: %d\n", numres);
      int k = 0;
      char lw[5];
      fprintf(fp, "     RESID INS CHN    RESID INS CHN    BASE-PAIR    ACCN\n");
      fprintf(fp, "     ---------------------------------------------------\n");
      for(int i=0; i<numres; ++i){
            for(int j=0; j<self[i].numbp; ++j){
                  k = self[i].oth_base_index[j];
		  to_leontis_westhof(lw, self[i].name[j]);

                  fprintf(fp, "PAIR %5d %c %3s      %5d %c %3s    %3s-%-3s %4s   %s\n",
                              self[i].cifid, self[i].ins, self[i].chain,
                              self[k].cifid, self[k].ins, self[k].chain,
                              self[i].resname, self[k].resname, lw, accn);

            }
      }
      fprintf(fp, "--------------- END OF RNA STRUCTURE %s --------------\n\n", accn);

}
void rnabp_scan_out(struct nucbp *self, int numres, char *outfile)
{
      FILE* fp = fopen(outfile, "r");
      if(fp == NULL){    /* Exception Handling */
	    fprintf(stderr, "Error in %s at line %d... Unable to open file %s\n", __FILE__, __LINE__, outfile);
	    exit(EXIT_FAILURE);
      }
      char line[256];
      while(fgets(line, 256, fp) != NULL){
	    if(line[0] == '#') continue;
	    break; /* pass the lines */
      }

      int i;
      char sep[] = "\t \n";
      char *token;
      while(line != NULL){
	    if(strlen(line)<=6) break;

	    token = strtok(line, sep);
	    int i = atoi(token);
	    i --;
	    self[i].numbp = 0;

	    token = strtok(NULL, sep);
	    self[i].cifid = atoi(token);

	    token = strtok(NULL, sep);
	    strcpy(self[i].resname, token);
	    self[i].resclass= get_nuc_class(self[i].resname);

	    token = strtok(NULL, sep);
	    if(strlen(token) != 1){
		  fprintf(stderr,"Error... ins-code having more than one characters\n");
		  exit(EXIT_FAILURE);
	    }
	    self[i].ins = token[0];

	    token = strtok(NULL, sep);
	    strcpy(self[i].chain, token);

	    /* For base pairs */
	    while((token = strtok(NULL, sep)) != NULL){

		  int numbp = self[i].numbp;
		  int j = atoi(token);
		  j --;

		  self[i].oth_base_index[numbp] = j;


		  token = strtok(NULL, sep);

		  token = strtok(NULL, sep);

		  token = strtok(NULL, sep);

		  token = strtok(NULL, sep);

		  token = strtok(NULL, sep);
		  strcpy(self[i].name[numbp], token);

		  token = strtok(NULL, sep);
		  strcpy(self[i].type[numbp], token);

		  token = strtok(NULL, sep);
		  self[i].eval[numbp] = atof(token);
		  self[i].numbp++;

	    }
	    fgets(line, 256, fp);

      }
      fclose(fp);
}


int get_numres(char* outfile)
{
      FILE* fp = fopen(outfile, "r");
      if(fp == NULL){    /* Exception Handling */
	    fprintf(stderr, "Error in %s at line %d... Unable to open file %s\n", __FILE__, __LINE__, outfile);
	    exit(EXIT_FAILURE);
      }
      char line[256];
      int flag = 0;
      int nres = -1;
      while(fgets(line, 256, fp) != NULL){
	    if(line[0] != '#') break;
	    if(strncmp(line, "#HEADER   Cleaned number", 24) != 0) continue;
	    nres = atoi(line+38);
	    flag = 1;
	    break;
      }
      if(flag == 0){
	    fprintf(stderr, "Error.... Something wrong in out file. No cleaned residue entry found.\n");
	    fclose(fp);
	    exit(EXIT_FAILURE);
      }
      return nres;
}

