import sys
with open(sys.argv[1], "r") as f1:
    last_line = f1.readlines()[-1]
#print (sys.argv[1])
print(last_line)