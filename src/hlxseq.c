/*
 * =====================================================================================
 *
 *       Filename:  hlxseq.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  Friday 01 July 2022 02:20:44  IST
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  PARTHAJIT ROY (PR), roy.parthajit@gmail.com
 *   Organization:  The University of Burdwan
 *
 * =====================================================================================
 */

#include "hlxseq.h"

void hlxinfo_free(struct hlxinfo* hlx)
{
      helix_free(hlx->helix);

      free ( hlx->rnabp );
      hlx->rnabp = NULL;
}

void hlxinfo_create_seq(struct hlxinfo* hlx, char* rule)
{
      struct helix* helix = NULL;
      for(int i=0; i<hlx->hlxcount; ++i){
	    helix = hlx->helix + i;
	    helix_create_seq(helix, hlx->rnabp, helix->seq,  rule);
      }
}
void hlxinfo_create(struct hlxinfo* hlx, char* rule, char* infile)
{
      char filepath[512];
      char basename[64];
      char ext[20];
      char outfile[512];
      char hlxfile[515];
      fname_split(filepath, basename, ext, infile);
      fname_join(outfile, filepath, basename, ".out");
      fname_join(hlxfile, filepath, basename, ".hlx");



      hlx->numres = get_numres(outfile);
      if( hlx->numres == 0 ){    /* Exception Handling */ 
	    fprintf(stderr, "Error in function %s()... no rna found for %s\n", __func__, basename);
	    exit(EXIT_FAILURE);
      }
      hlx->rnabp = (struct nucbp*) malloc ( hlx->numres * sizeof(struct nucbp) );
      if ( hlx->rnabp == NULL ) {
	    fprintf ( stderr, "\ndynamic memory allocation failed in function %s()\n" , __func__);
	    exit (EXIT_FAILURE);
      }

      rnabp_scan_out(hlx->rnabp, hlx->numres, outfile);

      helix_init(&hlx->helix, hlx->numres);
      helix_compute(hlx->helix, &hlx->hlxcount, hlx->rnabp, hlx->numres);
      hlxinfo_create_seq(hlx, rule);
      
      FILE* helixfp = fopen(hlxfile, "w");
      if ( helixfp == NULL ) {
	    fprintf ( stderr, "couldn't open file '%s'; %s\n",
			hlxfile, strerror(errno) );
	    exit (EXIT_FAILURE);
      }
      helix_fprint(hlx->helix, hlx->rnabp, hlx->hlxcount, helixfp);
      
      if( fclose(helixfp) == EOF ) {			/* close input file   */
	    fprintf ( stderr, "couldn't close file '%s'; %s\n",
			hlxfile, strerror(errno) );
	    exit (EXIT_FAILURE);
      }
}




void lcsubstr_init(struct lcsubstr* lcs, char* s, char* t)
{
      lcs->s = s;
      lcs->t = t;
      lcs->slen  = strlen(s);
      lcs->tlen  = strlen(t);
      lcs->maxlen = (lcs->slen > lcs->tlen) ? lcs->slen : lcs->tlen;
      lcs->longest = 0;
      

      
//      lcs->sindex = (int*) malloc ((lcs->maxlen + 4) * sizeof(int) );
//      if ( lcs->sindex==NULL ) {
//	    fprintf ( stderr, "\ndynamic memory allocation failed in function %s()\n" , __func__);
//	    exit (EXIT_FAILURE);
//      }
//      
//
//      lcs->tindex = (int*) malloc ( (lcs->maxlen + 4) * sizeof(int) );
//      if ( lcs->tindex==NULL ) {
//	    fprintf ( stderr, "\ndynamic memory allocation failed in function %s()\n" , __func__);
//	    exit (EXIT_FAILURE);
//      }
//

}

void lcsubstr_free(struct lcsubstr* lcs)
{
//      free ( lcs->sindex );
//      lcs->sindex	= NULL;
//
//      free ( lcs->tindex );
//      lcs->tindex	= NULL;
;
}


void lcsubstr_compute(struct lcsubstr* lcs) 
{
	int k=0;
	int longest = 0;
//	int* strs[2];
//	int **strs = (int **)malloc(2 * sizeof(int *));
	
        int strs[2][100];
//	strs[0] = NULL;
//	strs[1] = NULL;
	
	
//	for (int i = 0; i < 2; i++)
//		strs[i] = (int *)calloc(lcs->tlen, sizeof(int));
//	
//	for(int i=0; i<2; ++i){
//	      strs[i] = (int*) malloc ( (lcs->maxlen + 3) * sizeof(int) );
//	      if ( strs[i]==NULL ) {
//		    fprintf ( stderr, "\ndynamic memory allocation failed in function %s()\n" , __func__);
//		    exit (EXIT_FAILURE);
//	      }
//

//	}
	for(int i=0; i<lcs->tlen; ++i){
	      strs[0][i] = 0;
	}
	for(int i=0; i<lcs->tlen; ++i){
	      strs[1][i] = 0;
	}
	
	for (int i = 0; i < lcs->slen; ++i) {
	      for(int j=0; j<lcs->tlen; ++j){
		    strs[0][j] = strs[1][j];
	      }
	      //		memcpy(strs[0], strs[1], lcs->tlen * sizeof(int));
	      for (int j=0; j<lcs->tlen; ++j) {
		    if (lcs->s[i] == lcs->t[j]) {
			  if (i == 0 || j == 0) {
				strs[1][j] = 1;
			  } else {
				strs[1][j] = strs[0][j-1] + 1;
			  }
			  if (strs[1][j] > longest) {
				longest = strs[1][j];
				k = 0;
			  }
			  if (strs[1][j] == longest) {
				lcs->sindex[k] = i;
				lcs->tindex[k] = j;
				k++;
			  }
		    } else {
			  strs[1][j] = 0;
		    }
	      }
	}
	for(int i=0; i<k; ++i){
	      lcs->sindex[i] = lcs->sindex[i] - (longest - 1);
	      lcs->tindex[i] = lcs->tindex[i] - (longest - 1);
	}
	lcs->numstr = k;

	
//	for (int i = 0; i < 2; i++){
////	      free ( strs[i] );
////	      strs[i]	= NULL;
//	}
//	free(strs);

	
	lcs->longest = longest;
}

void lcsubstr_fprint(struct lcsubstr* lcs, struct helix* h1, struct helix* h2, struct nucbp* nbp1,  struct nucbp* nbp2, FILE* fp)
{
      int pos1, pos2;
      int pos3, pos4;
      int pairindx = 0;
      if(lcs->longest <4) return;
      for(int i=0; i<lcs->numstr; ++i){
	    fprintf(fp, "%5d %-4s", nbp1[h1->i[0]].cifid,nbp1[h1->i[0]].chain) ;
	    fprintf(fp, "%5d %-4s", nbp1[h1->j[0]].cifid,nbp1[h1->j[0]].chain) ;

	    fprintf(fp, "%5d %-4s", nbp2[h2->i[0]].cifid,nbp2[h2->i[0]].chain) ;
	    fprintf(fp, "%5d %-4s", nbp2[h2->j[0]].cifid,nbp2[h2->j[0]].chain) ;
	    fprintf(fp, "\n");
	    fprintf(fp, "MATCH1 ");
	    for(int j=0; j<lcs->longest; ++j){
		  pos1 = h1->i[lcs->sindex[i]+ j];
		  pos2 = h1->j[lcs->sindex[i]+ j];
		  pairindx = get_pair_index(nbp1, pos1, pos2);
		  if( pairindx < 0 ){    /* Exception Handling */ 
			fprintf(stderr, "Error in function %s()... invalid base pair.\n", __func__);
			exit(EXIT_FAILURE);
		  }

//		  pos3 = h2->i[lcs->tindex[i]+ j];
//		  pos4 = h2->j[lcs->tindex[i]+ j];
		  fprintf(fp, "%c-%c %s  ", nbp1[pos1].resclass, nbp1[pos2].resclass, nbp1[pos1].name[pairindx] );
//		  pairindx = get_pair_index(nbp2, pos3, pos4);
//		  if( pairindx < 0 ){    /* Exception Handling */ 
//			fprintf(stderr, "Error in function %s()... invalid base pair.\n", __func__);
//			exit(EXIT_FAILURE);
//		  }

//		  fprintf(fp, "%c-%c %s \n", nbp2[pos3].resclass, nbp2[pos4].resclass, nbp2[pos3].name[pairindx] );
	    }
	    fprintf(fp, "\n");
	    fprintf(fp, "MATCH2 ");
	    for(int j=0; j<lcs->longest; ++j){
//		  pos1 = h1->i[lcs->sindex[i]+ j];
//		  pos2 = h1->j[lcs->sindex[i]+ j];
//		  pairindx = get_pair_index(nbp1, pos1, pos2);
//		  if( pairindx < 0 ){    /* Exception Handling */ 
//			fprintf(stderr, "Error in function %s()... invalid base pair.\n", __func__);
//			exit(EXIT_FAILURE);
//		  }

		  pos3 = h2->i[lcs->tindex[i]+ j];
		  pos4 = h2->j[lcs->tindex[i]+ j];
		  pairindx = get_pair_index(nbp2, pos3, pos4);
		  if( pairindx < 0 ){    /* Exception Handling */ 
			fprintf(stderr, "Error in function %s()... invalid base pair.\n", __func__);
			exit(EXIT_FAILURE);
		  }

		  fprintf(fp, "%c-%c %s  ", nbp2[pos3].resclass, nbp2[pos4].resclass, nbp2[pos3].name[pairindx] );
	    }
	    fprintf(fp, "\n");
	    fprintf(fp, "\n%5d %-4s", nbp1[h1->i[h1->size - 1]].cifid,nbp1[h1->i[h1->size - 1]].chain) ;
	    fprintf(fp, "%5d %-4s", nbp1[h1->j[h1->size - 1]].cifid,nbp1[h1->j[h1->size - 1]].chain) ;

	    fprintf(fp, "%5d %-4s", nbp2[h2->i[h2->size - 1]].cifid,nbp2[h2->i[h2->size - 1]].chain) ;
	    fprintf(fp, "%5d %-4s", nbp2[h2->j[h2->size - 1]].cifid,nbp2[h2->j[h2->size - 1]].chain) ;
	    fprintf(fp, "\n");
	    fprintf(fp, "----------\n");
//	        printf("%d-%d\n", lcs.sindex[i], lcs.tindex[i]);
//		printf("%.*s\n", lcs.longest, &argv[1][lcs.sindex[i]]);
//		printf("%.*s\n", lcs.longest, &argv[2][lcs.tindex[i]]);
	}
	    fprintf(fp, "==========\n");
}
