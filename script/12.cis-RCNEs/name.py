import sys

dict = {}
with open(sys.argv[1],"r") as f:
	for i in f:
		i = i.strip().split("\t")
		x =("\t").join(i[0:3])
		dict[x]=i[3]
		#print(x)
with open(sys.argv[2],"r") as f:
	for i in f:
		i = i.strip().split("\t")
		x =("\t").join(i[0:3])
		#print(x)
		for x in dict:
			print(x+"\t"+dict[x]+"\t"+"0"+"\t"+"+")

