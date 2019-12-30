These programs can be compiled with "compile.sh". <br>
"proteins.h", "alignments_abs.h", "alignments_not_abs.h" are header files,
"razbor_proteins.c", "razbor_alignments_abs.c", "razbor_alignments_not_abs.c" contain functions for reading  
a set of sequences and a set of alignments. These functions calculate the intervals between given amino acids
or  conservative column and compare empirical and theoretical distributions according to the Kolmogorov-Smirnov criterion.
For alignments, an additional function was written to check the column for conservatism.
"runprogram.c", "runprogram_abs.c", "runprogram_not_abs.c" - run these functions for a given set of sequences or alignments.


The conservativeness threshold is not set via the command line; it is set manually as a function of checking for conservatism.
Program launch: <br>
1. for amino acid sequences, a command  like ./yourprogram of yours_array_specified_specified_amino acid_file_with_a result of your_por_p-value is launched
2. for alignments - a bash-program "runsh.sh" was created.  It passes through the folder with files with alignments and for each file with alignments
computes D-statistic and p-value. Then, using the "lookatyourres.py", you can calculate the percentage of alignments, the p-value for which
less than the specified threshold. 
