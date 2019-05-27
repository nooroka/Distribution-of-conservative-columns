for a in `ls -1 ess_aln/* | sort`
do
echo $a
./myprog $a resultess_aln.txt > outputess_aln.txt
done
