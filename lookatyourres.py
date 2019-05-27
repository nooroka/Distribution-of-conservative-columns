import sys
w = open("resultall.txt","a")
f = open(sys.argv[1],"r")
c_all = 0
c = 0
for line in f:
    line = line.split()
  #  print(line[1])
    if float(line[1]) < sys.argv[2]:
       c+=1
    c_all+=1
print (c)
print(c_all)
k = float(c)/float(c_all)
print("%.3f"%k)
w.write("\n")
w.write(str(sys.argv[1][6:]))
w.write(" ")
w.write(str(k))
w.close()