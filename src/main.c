/*
 * =====================================================================================
 *
 *       Filename:  main.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  Thursday 30 June 2022 09:32:55  IST
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  PARTHAJIT ROY (PR), roy.parthajit@gmail.com
 *   Organization:  The University of Burdwan
 *
 * =====================================================================================
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include "helix.h"
#include "rnabp.h"
#include "hlxseq.h"
#include "bioio.h"
#include "editdist.h"
#include "sysparams.h"
extern void callbpfindc(char [],  char [], char [], char [], char [], char [], char [], char [], char [], char [], char [], char[], char[], char[], char[]);

//void needleman_wunsch_sec_seq(int argc, char* argv[]){
void comp_needleman_wunsch(struct multi_fasta* seq1, struct multi_fasta* seq2,  struct multi_fasta* fasta1, struct multi_fasta* fasta2, struct runparams* params, FILE* fp)
{





//     char filepath[512];
//     char basename[64];
//     char ext[20];
//     char infile[512];
//     char seqfile[515];
//     char output_basename[64];
//     int dat_mat[128][128];
//     rna_secseq_scoring(dat_mat);
//     int fasta_mat[128][128];
//     rna_fasta_scoring(dat_mat);
//     int dat_gap = -2;
//     int fasta_gap = -2;
//
//     fname_split(filepath, basename, ext, argv[2]);
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//     fname_join(seqfile, filepath, output_basename, ".ssa");
//
//     FILE	*fp;										/* output-file pointer */
//
//     fp	= fopen( seqfile, "w" );
//     if ( fp == NULL ) {
//         fprintf ( stderr, "couldn't open file '%s'; %s\n",
//     		seqfile, strerror(errno) );
//         exit (EXIT_FAILURE);
//     }
//     printf("%d\n", argc);
//     for(int k=3; k<argc; ++k){  
//         struct hlxinfo rna2;
//         hlxinfo_create(&rna2, rule, argv[k]);
//         fname_split(filepath, basename, ext, argv[k]);
//         fname_join(infile, filepath, basename, ".sec");
//         //strcat(output_basename, basename);
//         char* t_aligned;
//         struct multi_fasta seq2;
//         scan_multifasta(&seq2, infile);
//         fname_join(infile, filepath, basename, ".fasta");
//         struct multi_fasta seq2_fasta;
//         scan_multifasta(&seq2_fasta, infile);   
         for(int i=0; i<seq1->size; ++i){
     	  //manipulate_secseq_coil(seq1.fasta[i].seq, seq1_fasta.fasta[i].seq);
     	  for(int j=0; j<seq2->size; ++j){
     		//manipulate_secseq_coil(seq2.fasta[j].seq, seq2_fasta.fasta[j].seq);
		char* s_aligned = NULL;
		char* t_aligned = NULL;
		
		fprintf(stderr, "Trace..... Executing File %s at line %d.\n", __FILE__, __LINE__);
     		fprintf(stderr, "1: PDB-%s  CHN-%s        2: PDB-%s   CHN-%s\n", seq1->fasta[i].accn, seq1->fasta[i].chain,
     			    seq2->fasta[j].accn, seq2->fasta[j].chain);
		
     		needleman_wunsch(seq1->fasta[i].seq, seq2->fasta[j].seq, params->sec_mat, params->sec_gap, &s_aligned, & t_aligned);
		fprintf(stderr, "Trace..... Executing File %s at line %d.\n", __FILE__, __LINE__);
     		fprintf(stderr, "1: PDB-%s  CHN-%s        2: PDB-%s   CHN-%s\n", seq1->fasta[i].accn, seq1->fasta[i].chain,
     			    seq2->fasta[j].accn, seq2->fasta[j].chain);

     		fprintf(fp, "1: PDB-%s  CHN-%s        2: PDB-%s   CHN-%s\n", seq1->fasta[i].accn, seq1->fasta[i].chain,
     			    seq2->fasta[j].accn, seq2->fasta[j].chain);
     		//fprintf(fp, "SECSEQ\n");
     		fprintf(fp, "%s", s_aligned);
     		fprintf(fp, "%s", t_aligned);
     		//fprintf(fp, "TER\n\n");

		fprintf(stderr, "Trace..... Executing File %s at line %d.\n", __FILE__, __LINE__);
		
     		int score = match_score_simple(s_aligned, t_aligned);
		fprintf(stderr, "Trace..... Executing File %s at line %d.\n", __FILE__, __LINE__);
     		fprintf(fp, "SEC SCORE: %d\n\n", score);
     		free(s_aligned);
     		s_aligned = NULL;
     		free(t_aligned);
     		t_aligned = NULL;



		fprintf(stderr, "Trace..... Executing File %s at line %d.\n", __FILE__, __LINE__);
     		needleman_wunsch(fasta1->fasta[i].seq, fasta2->fasta[j].seq, params->fasta_mat, params->fasta_gap, &s_aligned, & t_aligned);
		fprintf(stderr, "Trace..... Executing File %s at line %d.\n", __FILE__, __LINE__);

     		//fprintf(fp, "FASTA   1: PDB-%s  CHN-%s        2: PDB-%s   CHN-%s\n", seq1_fasta.fasta[i].accn, seq1_fasta.fasta[i].chain,
     		//                                                            seq2_fasta.fasta[j].accn, seq2_fasta.fasta[j].chain);
     		//fprintf(fp, "FASTA\n");
     		fprintf(fp, "%s", s_aligned);
     		fprintf(fp, "%s", t_aligned);
		fprintf(stderr, "Trace..... Executing File %s at line %d.\n", __FILE__, __LINE__);
     		score = match_score_simple(s_aligned, t_aligned);
		fprintf(stderr, "Trace..... Executing File %s at line %d.\n", __FILE__, __LINE__);
     		fprintf(fp, "PRI SCORE: %d\n", score);
     		fprintf(fp, "TER\n\n");
     		free(s_aligned);
     		s_aligned = NULL;
     		free(t_aligned);
     		t_aligned = NULL;


     	  }
         }
//         multifasta_free(&seq2);
//         multifasta_free(&seq2_fasta);
//         hlxinfo_free(&rna2);
//     }
//
//
//     if( fclose(fp) == EOF ) {			/* close output file   */
//         fprintf ( stderr, "couldn't close file '%s'; %s\n",
//     		seqfile, strerror(errno) );
//         exit (EXIT_FAILURE);
//     }
//     multifasta_free(&seq1);
//     multifasta_free(&seq1_fasta);
//     hlxinfo_free(&rna1);
//
}






int main(int argc, char* argv[])
{
      struct bpfind_params bpfindprm;
      bpfind_params_init(&bpfindprm);
      struct runparams params;
      runparams_init(&params);
      params.sec_gap = -2;
      params.fasta_gap = -2;
      
      rna_secseq_scoring(params.sec_mat);
      rna_fasta_scoring(params.fasta_mat);
      int* file_index = (int*) malloc( argc * sizeof(int) );
      int file_count = 0;

      for(int i=1; i<argc; ++i){
	    if(strcmp(argv[i], "-help") == 0){
		  fprintf(stderr, "Under construction\n");
		  exit(EXIT_SUCCESS);
	    }
	    if(strcmp(argv[i], "-pub") == 0){
		  fprintf(stderr, "publication under construction\n");
		  exit(EXIT_SUCCESS);
	    }
	    if(strcmp(argv[i], "-nsg") == 0){
		  strcpy(bpfindprm.sgparam, "-SG");
	    }else if(strcmp(argv[i], "-nch") == 0){
		  strcpy(bpfindprm.chparam, "-CH");
	    }else if(strcmp(argv[i], "-nht") == 0){
		  strcpy(bpfindprm.htparam, "-dummyval");
	    }else if(strcmp(argv[i], "-ho") == 0){
		  params.helix_only = 'T';
	    }else if(strcmp(argv[i], "-m") == 0){
		  params.multi_cif = 'T';
	    }else{
		  file_index[file_count] = i;
		  file_count ++;
	    }
      }
      char file_name[512];
      char file_path[512];
      char file_extn[512];
      char file1_name[512];
      char file1_path[512];
      char file1_extn[512];
      if( file_count == 1 ){    /* Exception Handling */ 
	    fprintf(stderr, "Error in function %s()... only one file supplied for comparison.\n", __func__);
	    exit(EXIT_FAILURE);
      }
      if(params.multi_cif == 'F' && file_count != 2){    /* Exception Handling */ 
	    fprintf(stderr, "Error in function %s()... please supply exactly two structure files.\n", __func__);
	    exit(EXIT_FAILURE);
      }

      fname_split(file1_path, file1_name, file1_extn, argv[file_index[0]]);

      strcpy(bpfindprm.accn, file1_name);
      strcpy(bpfindprm.ext, file1_extn);
      if(strcmp(file1_extn, ".cif") == 0){
	    strcpy(bpfindprm.cifparam, "-CIF");
	    strcpy(bpfindprm.accnparam, argv[file_index[0]]);
      }else if(strcmp(file1_extn, ".pdb") == 0){
	    strcpy(bpfindprm.cifparam, "-dummyval");
	    strcpy(bpfindprm.accnparam, argv[file_index[0]]);
      }else{    /* Exception Handling */ 
	    fprintf(stderr, "Error in function %s()... please supply either mmCIF or PDB file.\n", __func__);
	    exit(EXIT_FAILURE);
      }
      fprintf(stdout, "Please wait... computing base pairs for structure %s", file1_name);
      callbpfindc(bpfindprm.cifparam, bpfindprm.accnparam, bpfindprm.htparam,
		  bpfindprm.hdparam, bpfindprm.hdvalparam, bpfindprm.angparam,
		  bpfindprm.angvalparam, bpfindprm.chparam, bpfindprm.sgparam,
		  bpfindprm.corparam, bpfindprm.evaltypeparam,
		  bpfindprm.chainparam, bpfindprm.chainvalparam,
		  bpfindprm.nmrparam, bpfindprm.nmrvalparam);
      strcpy(params.file1, argv[file_index[0]]);
      struct hlxinfo rna1;


      char rule[6];
      strcpy(rule, "TTTT");

      hlxinfo_create(&rna1, rule, params.file1);
      char infile[512];
      fname_join(infile, file1_path, file1_name, ".sec");
      struct multi_fasta secseq1;
      scan_multifasta(&secseq1, infile);

      fname_join(infile, file1_path, file1_name, ".fasta");
      struct multi_fasta fastaseq1;
      scan_multifasta(&fastaseq1, infile);
      fprintf(stdout, " ...[finished]\n");



      char base_name[512];
      strcpy(base_name, file1_name);

      if(params.multi_cif == 'F'){
	    fname_split(file_path, file_name, file_extn, argv[file_index[1]]);
	    strcat(base_name, "_");
	    strcat(base_name, file_name);
      }

      FILE	*ssafp;										/* output-file pointer */
      char	ssafp_file_name[512];		/* output-file name    */


      fname_join(ssafp_file_name, file1_path, base_name, ".ssa");

      ssafp	= fopen( ssafp_file_name, "w" );
      if ( ssafp == NULL ) {
	    fprintf ( stderr, "couldn't open file '%s'; %s\n",
			ssafp_file_name, strerror(errno) );
	    exit (EXIT_FAILURE);
      }
      params.ssafp = ssafp;

      for(int i=1; i<file_count; ++i){
	    strcpy(params.file2, argv[file_index[i]]);
	    fname_split(file_path, file_name, file_extn, argv[file_index[i]]);

	    strcpy(bpfindprm.accn, file_name);
	    strcpy(bpfindprm.ext, file_extn);
	    if(strcmp(file_extn, ".cif") == 0){
		  strcpy(bpfindprm.cifparam, "-CIF");
		  strcpy(bpfindprm.accnparam, argv[file_index[i]]);
	    }else if(strcmp(file_extn, ".pdb") == 0){
		  strcpy(bpfindprm.cifparam, "-dummyval");
		  strcpy(bpfindprm.accnparam, argv[file_index[i]]);
	    }else{    /* Exception Handling */ 
		  fprintf(stderr, "Error in function %s()... please supply either mmCIF or PDB file.\n", __func__);
		  exit(EXIT_FAILURE);
	    }
	    fprintf(stdout, "Please wait... computing base pairs for structure %s", file_name);
	    callbpfindc(bpfindprm.cifparam, bpfindprm.accnparam, bpfindprm.htparam,
			bpfindprm.hdparam, bpfindprm.hdvalparam, bpfindprm.angparam,
			bpfindprm.angvalparam, bpfindprm.chparam, bpfindprm.sgparam,
			bpfindprm.corparam, bpfindprm.evaltypeparam,
			bpfindprm.chainparam, bpfindprm.chainvalparam,
			bpfindprm.nmrparam, bpfindprm.nmrvalparam);
	    struct hlxinfo rna2;
	    hlxinfo_create(&rna2, rule, params.file2);
	    fname_join(infile, file_path, file_name, ".sec");
	    struct multi_fasta secseq2;
	    scan_multifasta(&secseq2, infile);

	    fname_join(infile, file_path, file_name, ".fasta");
	    struct multi_fasta fastaseq2;
	    scan_multifasta(&fastaseq2, infile);

	    fprintf(stdout, " ...[finished]\n");

//	    needleman_wunsch();

	    comp_needleman_wunsch(&secseq1, &secseq2, &fastaseq1, &fastaseq2, &params, ssafp);
	    hlxinfo_free(&rna2);
	    multifasta_free(&secseq2);
	    multifasta_free(&fastaseq2);
		  
      }

      hlxinfo_free(&rna1);
      multifasta_free(&secseq1);
      multifasta_free(&fastaseq1);


      if( fclose(ssafp) == EOF ) {			/* close output file   */
	    fprintf ( stderr, "couldn't close file '%s'; %s\n",
			ssafp_file_name, strerror(errno) );
	    exit (EXIT_FAILURE);
      }
      params.ssafp = NULL;




      //      needleman_wunsch_sec_seq(argc, argv);
}

