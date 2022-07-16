//
// Created by parthajit on 15/8/20.
//

#ifndef CONSHELIX_EDITDIST_H
#define CONSHELIX_EDITDIST_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "exception.h"
#include "ndarray.h"


void manipulate_secseq_coil(char* s, char* t);
void rna_secseq_scoring(int SMAT[][128]);
void rna_fasta_scoring(int SMAT[][128]);
double match_score_simple(char* a, char* b);


void needleman_wunsch(char* s, char* t, int NWSim[][128], int gap_penalty, char** s_aligned, char** t_aligned);




void smith_waterman_seq_align(char* s, char* t, char** s_aligned, char** t_aligned);






void hirschberg_seq_align(char* s, char* t, char** s_aligned, char** t_aligned);











int edit_dist(const char* source, const char* target, char** s_aligned, char** t_aligned);



int Wagner_Fischer_weighted_edit_dist(const char* source, const char* target, char* action_seq, char* action_char);





/* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
 *
 *  KMP_match_first(...)
 *
 *  This algorithm matches the first occurance of the pattern in the text.
 *  The return value is 0 or +ve value to indicate the match start index in the text array.
 *  The return value -1 indicates no match found.
 *
 *  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */

int KMP_match_first(char* pattern, int p_len, char* text, int t_len);


/* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
 *
 *    lc_subseq(...)
 *
 *    This function computes the longest common subsequence of two strings. The Implementation  *
 *    is done based on Coremen et al. book on Algorithm. The Algorithm uses dynamic programming *
 *    and DP table for its execution.
 *
 *    INPUT: A is a string of size na. (array size is na+1 for accomodationof \0)
 *    INPUT: B is a string of size nb. (array size is nb+1 for accomodationof \0)
 *    OUTPUT: lcs_out is a pointer to char for storing the LCS. (size is min(na, nb) + 1
 *    OUTPUT: n is a pointer to long for storing the LCS size.
 *    OUTPUT: A1 is the modified A array with all non LCS positions as '.'. Some thing like
 *    if LCS is BAB, the A is like ..B.A..B.. (size is na+1)
 *    OUTPUT: B1 is same as A1 but for B.     (size is nb+1)
 *
 *    RETURN: It returns the length of LCS.
 *
 *    Note: In all the cases of lcs_out, A1 and B1 the memory SHOULD NOT be allocated by the user.
 *    The algorithm DOSE allocate all the memory for them. The user needs to call
 *    lc_subseq_free(...) when the task is done.
 *
 *    $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$*/

void lc_subseq_free(char* lcs, char* a, char* b);

int lc_subseq(const char* A, const char* B, char** lcs, char** a, char** b);



#endif //CONSHELIX_EDITDIST_H
