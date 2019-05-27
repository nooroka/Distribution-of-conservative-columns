for a in `ls -1 filesacids/* | sort`
do
echo $a
python look.py $a
done
