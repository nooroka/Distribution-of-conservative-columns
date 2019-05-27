#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h> 
#include <math.h>



int cmp(const void * x1, const void * x2);

struct alignment {
    int numseq;
    struct seq* sequence;
};

struct seq {
    char* name;
    char* description;
    char* sequence;
    int  length;
};

struct intervals {
      float* dist;
      int number;
      int n;
      float* E;
      float* F; 
      double scale;
      float* M;
      float D;
};

struct distribution {
      int* znach;
      int* x;
      int masslen;
      float* fun;
};


struct seq seqread(FILE *infastafile);
void seqwrite(FILE *outfastafile, struct seq sequence, int k);
struct intervals distancepro (struct seq chain, char residue);
struct distribution empirical (struct seq chain, char residue, struct intervals distancepro);
void free_intervals(struct intervals *x);
void free_seq(struct seq *x);
double probks(double alam);
void ksone(float *data, int n, float scale, double *d, double *prob);
/* *mean1(struct intervals irr1);
double *variance1(struct intervals irr1);*/

/*float probks(float alam);*/
/*void ksone(struct intervals irr, float (*func)(float), float *d,
	float *prob);*/


