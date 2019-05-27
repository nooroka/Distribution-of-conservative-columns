#include "sequence3.h"
#define SIZE	800
int main(int argc, char **argv) {
    /*arg1 is infile, arg2 is amino acid, arg3 is outfile*/
    struct seq mysequence1;
    
    FILE *myfile;
//    FILE *outfile;
//    FILE *outfile2;
    FILE *outresult;
    FILE *Dd;
    struct intervals res1;
    struct distribution resk;
    int i, j, k;
    int p, b, r;
    int violent;
    int sqc = 0;
    int sqcall;
    int nap = 0;
    double alam;
    double d;
    double prob;
    int maxlen = 1000, m = 800;
    int maxname = 20;
    int maxdesc = 100;
    int maxseq = 900;
    int all = 900;
    struct seq seqall1[SIZE];
   
    if (argc < 2){
        printf("Not enough arguments!");
        exit(0);
    }
    myfile = fopen (argv[1], "r");
    if (myfile == NULL){
        fprintf(stderr, "Cannot open file %s!\n", argv[1]);
        exit(1);
    }
    outresult = fopen (argv[2], "a");
    if (outresult == NULL){
        fprintf(stderr, "Cannot open file %s!\n", argv[3]);
        exit(1);
    }
 
   
/*    outfile2 = fopen ("outR.txt", "w"); */
    
    while (!feof(myfile)  ) {
       mysequence1 = seqread(myfile, sqc, seqall);    
       violent = mysequence1.length;
       sqc+=1;   
       sqcall = sqc;
       printf("\nSqc is cycle: %d\n",sqc);
       /*free_seq(&mysequence1);*/
    } /* while (!feof(myfile)) */

    
    printf("\nSqc is: %d\n",sqc);
    res1 = distancepro(violent, seqall, sqcall); 
    free_seq(&mysequence1);
    printf("\nNumber of conservative columns: %d\n", res1.n);
    if (res1.n == 0){
        printf("\nNo  conservative columns:(");
        exit (1);
    }

  /*  while (i < res1.number) {
           printf("%.2f ", res1.dist[i]);
           i++;          
    } */
    printf("\n");
    if (res1.n > 4 ){
    /*       printf ("%s\n", mysequence1.name);
           printf ("%s\n", mysequence1.description);
           printf ("%s\n", mysequence1.sequence);
           printf("Sorted %d intervals between given residues  ", res1.number); */
           printf("D manual ");
           printf("%.6f",res1.D);
           printf("\n");
           ksone(res1.dist, res1.number, res1.scale, &d, &prob);
           printf("%.2e\n", d);
           Dd = fopen ("p-value1.txt", "a");
           fprintf(Dd, "%.2e\n",prob);       
           fclose(Dd);                    
           printf("%.2e\n", prob);
           fprintf(outresult,"%s  ", argv[1]); 
           fprintf(outresult,"%.2e\n", prob); 
           free_intervals(&res1);
    }    
    fclose(myfile);
    fclose(outresult);
    return 0;
} 
