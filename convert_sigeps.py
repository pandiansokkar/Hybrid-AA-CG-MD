import os
import re

out = open("output.txt","w")
nbfile = open("nonbond.tmp","r")
lines = nbfile.readlines()
for i in lines:
	line=i.split()
	atom1=line[0]
	atom2=line[1]
	func=line[2]
	c6=line[3]
	c12=line[4]
	cmd="gmx sigeps -c6 " + c6 + " -cn " + c12 + "> temp"
	os.system(cmd)
	os.system("rm ./#*")
	parfile = open("temp","r")
	parlines=parfile.readlines()
	for j in parlines:
		if j.startswith("sigma"):
			par = re.split(",",j)
               		par1 = par[0].split()
                	sig = par1[2]
                	par2 = par[1].split()
                	eps = par2[2]
			
	parfile.close()
	
	if (sig == "-nan"):
		sig = 0.47
	if (eps == "-nan"):
		eps = 0.0

	print('{0}\t{1}\t{2}\t{3}\t{4}\n'.format(atom1, atom2, func, sig, eps))

	out.write('{0}\t{1}\t{2}\t{3}\t{4}\n'.format("v"+atom1, "v"+atom2, func, sig, eps))

out.close()
					
