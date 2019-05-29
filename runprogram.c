#include "proteins.h"

int main(int argc, char **argv) {
    /*arg1 is infile, arg2 is amino acid, arg3 is outfile*/
    struct seq mysequence1;
    FILE *myfile;
//    FILE *outfile;
//    FILE *outfile2;
    FILE *outresult;
    char residue;
    struct intervals res1;
    struct distribution resk;
    int i, j, k;
    int p, b, r;
    int count = 0;
    int f1 = 0;
    double alam;
    double d;
    double prob;
    double threshold = 0;
    double divide;
    if (argc < 4){
        printf("Not enough arguments!");
        exit(0);
    }
    myfile = fopen (argv[1], "r");
    if (myfile == NULL){
        fprintf(stderr, "Cannot open file %s!\n", argv[1]);
        exit(1);
    }
    outresult = fopen (argv[3], "w");
    if (outresult == NULL){
        fprintf(stderr, "Cannot open file %s!\n", argv[3]);
        exit(1);
    }
 
    residue = argv[2][0];
    threshold = atof(argv[4]);
    printf("%.2f\n",threshold);
/*    outfile2 = fopen ("outR.txt", "w"); */
    while (!feof(myfile)) {
       mysequence1 = seqread(myfile);

/*
        outfile = fopen ("out2.txt", "w"); 
        if (outfile == NULL){
            fprintf(stderr, "Cannot open file out2.txt for one sequence!");
            exit(1);
       }
       mysequence1 = seqread(myfile);
       seqwrite(outfile, mysequence1, 72);
       fclose(outfile);
       outfile = fopen ("out2.txt", "r");
       mysequence2 = seqread(outfile);
       printf("\nSequence length is %d\n", mysequence2.length);
       res1 = distancepro(mysequence2, residue,outfile2); */

/*       printf("\nSequence length is %d\n", mysequence1.length); */
       res1 = distancepro(mysequence1, residue);
       if (res1.n > 5 ){
          /* i = 0;*/
           f1+=1;
           printf ("%s\n", mysequence1.name);
           printf ("%s\n", mysequence1.description);
           printf ("%s\n", mysequence1.sequence);
   /*      printf("Sorted %d intervals between given residues  ", res1.number); */
           printf("D manual ");
           printf("%.6f",res1.D);
           printf("\n");
           ksone(res1.dist, res1.number, res1.scale, &d, &prob);
           printf("%.2e\n", d);
           printf("%.2e\n", prob); 
           if (prob < threshold){
              fprintf(outresult,"%.2e\n", prob); 
              count +=1;
           }
       } 
       free_intervals(&res1);
       free_seq(&mysequence1);
    } /* while (!feof(myfile)) */
    divide = (double) count/f1;
    fprintf(outresult,"p-value < %.2f: %d\n", threshold, count); 
    fprintf(outresult,"All is %d\n", f1); 
    fprintf(outresult,"Divide is %4f\n", divide); 
    fclose(myfile);
    fclose(outresult);
   
    return 0;
} /* main */
