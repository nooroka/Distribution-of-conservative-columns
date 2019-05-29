#include "proteins.h"
/*write destructors for releasing memory*/
#define NMAX 2000 
#define EPS1 0.001
#define EPS2 1.0e-8 

static float maxarg1, maxarg2;

int cmp(const void * x1, const void * x2)
{
  return ( *(int*)x1 - *(int*)x2 );              
}

struct seq seqread(FILE *infastafile) {
    char sym;
    int maxlen = 1000, m = 800;
    int maxnamelen = 200;
    int maxdesclen = 1000;
    int maxseqlen = 9000;
    int t1 = 0,  t2 = 0,  t3 = 0;
    struct seq result;
    long int tmp;

    result.name = (char *)malloc(maxnamelen*sizeof(char ));
    result.description = (char *)malloc(maxdesclen*sizeof(char )); 
    result.sequence = (char *) malloc(maxseqlen*sizeof(char ));  
    
    result.length = 0;
    sym = fgetc(infastafile);
    if (sym == '>'){
        t1 = 0;
        while ( (sym != '\n') && (sym != ' ') && (sym != '\r') ) {
            sym = fgetc(infastafile);
            result.name[t1] = sym;
            t1++;
            if (t1 > maxnamelen){
                maxnamelen += 200;
                result.name = (char *)realloc(result.name, maxnamelen*sizeof(char)); 
            }
        }
        result.name[t1 - 1] = '\0';

        if (sym == ' ') { /* there is a description */
            t2 = 0;
            while( (sym != '\n') && (sym != '\r') ){
                sym = fgetc(infastafile);
                result.description[t2] = sym;
                t2++;  
                if (t2 > maxdesclen){
                    maxdesclen += 200;
                    result.description = (char *)realloc(result.description, maxdesclen*sizeof(char)); 
                }          
           }
           result.description[t2 - 1] = '\0';
       }
       else { /* empty description */
           result.description[0] = '\0';
       }
    } 
    else {
      /*fprintf(stderr, "First symbol is %c\n", sym);
      perror("\nNot fasta!\n");*/
      exit(1);
    }

    t3 = 0;
    sym = fgetc(infastafile);
   /* result.sequence[t3] = *(int *) malloc(mama*sizeof(int));*/
    while (!feof(infastafile) && sym != '>' ) { 
       /* printf("%c",sym);*/
        if ( isalpha(sym) || sym == '-' || sym == '~' || sym == '.' ) {
           /* result.sequence[t3] = (char )malloc(mama*sizeof(char ));*/
            result.sequence[t3] = sym;
            t3++;
            if (t3 >= maxseqlen){
                maxseqlen += 1200;
                result.sequence = (char *)realloc(result.sequence, maxseqlen*(sizeof(char ))); 
           }  
        }
        sym = fgetc(infastafile);
        if (sym == '>') {
            tmp = ftell(infastafile);
         /* printf("%ld\n", tmp);*/
            fseek(infastafile, tmp - 1, SEEK_SET);
        /*  printf("%ld\n", tmp);*/
        }
   }
   result.sequence[t3] = '\0';
   result.length = strlen(result.sequence);
  
   return result;
}

void seqwrite(FILE *outfastafile, struct seq sequence, int k){
    fprintf(outfastafile, ">");
    fprintf(outfastafile, "%s", sequence.name);
    if ( strlen(sequence.description) > 0 ) {
        fprintf(outfastafile, " %s\n", sequence.description);
    }
    else {
        fprintf(outfastafile, "\n");
    }

    for (int l = 0; l < strlen(sequence.sequence); l++) {
        fprintf(outfastafile, "%c", sequence.sequence[l]);
        if ( (l + 1) % k == 0 ){
            fprintf(outfastafile, "\n");
        }
   }
   fprintf(outfastafile, "\n");
   return;
}

struct intervals distancepro (struct seq chain, char residue) {
    
    struct intervals result;
    int temp = 0;
    int maxnum = 9000;
    int w;
    int i;
    int maxl = 0;
    int j = 0;
    int maxnumber = 1000;
    
    result.dist = (float *) malloc(sizeof(float) * maxnum);
    i = 0;
    for (w = 0; w < chain.length; w++) {
        if ( toupper(chain.sequence[w]) == toupper(residue) ) {
            result.dist[i] = (float)temp;
            printf(".%4f",result.dist[i]);
            temp = 0;
            i++;
            if ( i >= maxnum ) {
                maxnum += 200;
                result.dist = (float *) realloc(result.dist, sizeof(float) * maxnum);
            }
        }
        else {
            temp++;
        }
    }
    result.dist[i] = temp;
    result.number = i + 1;
 
    result.n = i;

/*    printf("Non sorted %d intervals between given residues  ", result.number); 
    while (j < result.number) {
        printf("%.1f ", result.dist[j]);
        j++;
    } 
    printf("\n"); */

    qsort(result.dist, i + 1, sizeof(int), cmp );
    
    result.E = (float *) malloc(sizeof(float) * result.number);
    result.F = (float *) malloc(sizeof(float) * result.number);
    result.M = (float *) malloc(sizeof(float) * result.number);
    result.scale = ((double)(chain.length - i))/(i + 1); 
    for (int d = 0; d < result.number; d++) {
        result.E[d] = ((float)(d+1))/result.number;       
        result.F[d] = (float)(1 - exp(-result.dist[d]/result.scale));
        
        if (d > 0) {
/*           result.M[d] = fmax(fabs(result.E[d] - result.F[d]), fabs(result.E[d - 1] - result.F[d])); */
             if (fabs(result.E[d] - result.F[d]) > fabs(result.E[d - 1] - result.F[d])) result.M[d] = fabs(result.E[d] - result.F[d]);
             else result.M[d] =  fabs(result.E[d - 1] - result.F[d]);
        }
        else  {
           if (fabs(result.E[d] - result.F[d]) > result.F[d]) result.M[d] = fabs(result.E[d] - result.F[d]); 
           else result.M[d] = result.F[d];
        }
    }
    for (int b = 0; b < result.number; b++) {
        if (result.M[b] > result.M[maxl]) {
            maxl = b;
        }
    }
    result.D = result.M[maxl];    
    return result;
}

double probks(double alam) {
  int j;
  float a2 = 0.0;
  float fac;
  double sum = 0.0;
  double term, termbf = 0.0;

  fac = 2.0;
  a2 = -fac * alam * alam;
  for (j = 1; j <= 1000; j++) {
      term = fac * exp(a2 * j * j);
      sum += term;
      if ( fabs(term) < EPS1 * termbf || fabs(term) < EPS2 * sum ) return sum;
      fac = -fac;
      termbf = fabs(term);
  }
  return 1.0;
}

void ksone(float *data, int n, float scale, double *d, double *prob)
/* Given an array data[1..n], and given a user-supplied function of a single variable func which
is a cumulative distribution function ranging from 0 (for smallest values of its argument) to 1
(for largest values of its argument), this routine returns the K�S statistic d, and the significance
level prob. Small values of prob showtha t the cumulative distribution function of data is
significantly different from func. The array data is modified by being sorted into ascending
order. */
{
  unsigned long j;
  float dt, en, ff, fn, fo=0.0;
  double sqren;
  en = (float)(n);
  sqren = sqrt(en);
  *d=0.0;
  for (j=1; j<=n; j++) { // Loop over the sorted data points.
    fn = j/en; // Data�s c.d.f. after this step.
    ff= 1.0 - exp(-data[j-1]/scale); // Compare to the user-supplied function.
    if (fabs(fo - ff) > fabs(fn - ff)) dt = fabs(fo - ff); else dt = fabs(fn - ff); // Maximum distance.
    if (dt > *d) *d = dt;
    fo = fn;
  }
  *prob = probks((sqren + 0.12 + 0.11/sqren) * (*d)); // Compute significance.
  return;
}

void free_intervals(struct intervals *x) {
  if (x -> dist != NULL) {
    free(x -> dist);
    x -> dist = NULL;  
  }

  if (x -> E != NULL) { // doesn't do outfile 
    free(x -> E);
    x -> E = NULL;  
  }
  
  if (x -> F != NULL) { 
    free(x -> F);
    x -> F = NULL;  
  }
  if (x -> M != NULL) {
    free(x -> M);
    x -> M = NULL;  
  }
   
  x -> number = 0;
  x -> n = 0;
  x -> scale = 0;
  return;
}

void free_seq(struct seq *x) {
  if (x -> name != NULL) {
    free(x -> name);
    x -> name = NULL;  
  }
  if (x -> description != NULL) {
  free(x -> description);
  x -> description = NULL;  
  }
  if (x -> sequence != NULL) {
  free(x -> sequence);
  x -> sequence = NULL;  
  }
  x -> length = 0; 
  return;
}

