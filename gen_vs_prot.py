# This code is used to generate virtual sites to perform multiscale AA/CG simulations with GROMACS
# to run this script,
# python gen_vs_prot.py -h
# The above command will show you the help menu
# This script works only for single chain of a protein
# Compatibility with martini2.2 and martini2.2p forcefields
# For multichain proteins and other molecules, stay tuned


# This code was developed as part of DST-SERB, SRG/2019/002156 project.
# If you need any help with the setup of multiscale MD, please write to pandian.sokkar@gmail.com

import re
import numpy as np
import argparse
import sys
import os
from array import *
from shutil import copyfile

# global variables
anum = []
resid = []
resname = {}
aname = {}
resatoms = {}
vsfunc = 2
lastatom = 0
coords = {}
allpdbinfo = {}
Deb = "yes"
chains=[]
NTER_res=[]
CTER_res=[]
# main program to call other functions
def main():
	args = ParseOptions()
	readPDB(args.inputfile,args.debug)
	#print(allpdbinfo)
	Deb=args.debug
	if args.debug=="full":
		print('resid:	',resid)
		print('resname:	',resname)
		print('atoms:	',aname)
		print('resatoms:',resatoms)
		print('atom num:',anum)

	nres=len(resid)
	nresnam=len(resname)
	natms=len(aname)
	nresatms=len(resatoms)
	nchains = len(chains)
	nNTER = len(NTER_res)
	#print("summary \n")
	#print('no res: {0}, no resnames: {1}, no atoms: {2}, no resatoms {3}'.format(nres,nresnam,natms,nresatms))
	#print('test printing...')

	#print('test {0}'.format(resatoms[1]))
	#copyfile(args.inputfile,args.vsitepdb)
	getVS(args.outputfile,args.forcefield,args.vsitepdb,args.inputfile)
	get_VSindex(args.vsitepdb, args.indexfile)

# block to parse arguments	
def ParseOptions():
	parser = argparse.ArgumentParser()
	parser.add_argument("-if", "--inputfile", help='coordinates file (PDB)', action="store", default='input.pdb', metavar='input.pdb')
	parser.add_argument("-vspdb", "--vsitepdb", help='coordinates file (PDB)', action="store", default='viste.pdb', metavar='vsite.pdb')
	parser.add_argument("-of", "--outputfile", help='name of the output ilfe (TXT)', action="store", default='vs.out', metavar='vs.out')
	parser.add_argument("-cid", "--chainid", help='chain id to get vsites (A, B, C or all', action="store", default='all', metavar='all')
	parser.add_argument("-ff", "--forcefield", help='forcefield to use (martini22, martini22p)', action="store", default='all', metavar='all')
	parser.add_argument("-deb", "--debug", help='Perform debug (yes/no)', action="store", default='no', metavar='no')
	parser.add_argument("-ndx", "--indexfile", help='index file (ndx)', action="store", default='index.ndx', metavar='index.ndx')
	
	args = parser.parse_args()
	if not os.path.exists(args.inputfile):
		print('\n file {0} not found ..\n' .format(args.inputfile))
		parser.print_help()
		sys.exit(1)
	
	return args

# read pdb and get necessary info 
def readPDB(FileName,Deb):
	try:
		pdbfile = open(FileName,'r')
	except:
		raise IOError ('\n cannot open {0}...\n'.format(FileName))
	
	
	pdblines=pdbfile.readlines()
	
	
	sno=0
	prev_nter_res=0
	prev_cter_res=0
	for x in pdblines:
		if x.startswith('CRYST'):
			allpdbinfo['Title']=x
		if x.startswith('ATOM'):
			prev_sno=sno
			sno=int(x[6:11].strip())
			anm=x[11:17].strip()
			pdb_resname=x[17:21].strip()
			chn=x[21].strip()
			resno=int(x[22:26].strip())
			ax=x[32:38].strip()
			ay=x[40:46].strip()
			az=x[48:54].strip()
			coords[sno] = [ax,ay,az]

			alphacol=x[56:60]
			betacol=x[62:67]
			elemid=x[78].strip()
			
			
			#find out if the residue is NTER or CTER or MID
			
			if chn not in chains:
				nter="yes"
				NTER_res.append(resno)
				prev_nter_res = resno
			else:
				if prev_nter_res == resno:
					nter = "yes"
				else:
					nter = "no"
			
			if anm == "OT1"  or anm == "OT2" or anm == "OXT":
				cter = "yes"
				prev_cter_res = resno
			else:
				cter= "no"
		#	break
			
			allpdbinfo[sno]={'AtomName':anm, 'ResName':pdb_resname, 'Chain':chn, 'ResNo': resno, 'X':ax, 'Y':ay, 'Z':az, 'Alpha':alphacol, 'Beta':betacol, 'ElemId':elemid, 'NTER':nter, 'CTER':cter}

			chains.append(chn) if chn not in chains else chains			
			# append to the global variables
			if resno not in resid:
				resid.append(resno)
				resname[resno] = pdb_resname

			aname[sno]=anm
			resatoms[sno]=resno
			anum.append(sno)
			
			
			if Deb == "full":
				print('SNo. {0}: Atm {1}: Resn {2}: Chn {3}: resno {4}'.format(int(sno), anm, pdb_resname, chn, int(resno)))
	
	# set all 
	pdbfile.close()

	
# collect virtual site info from other functions	
def getVS(outname,ff,outpdb,inputfile):
	vsout=open(outname,"w")
	vsout.write(";; the following should be copied into a atomistic topology file in suitable places\n")
	vsout.write('[ virtual_sitesn ]'+"\n")
	vspdb=open(outpdb,"w+")
	#vspdb.write(allpdbinfo['Title']
	
	with open(inputfile,'r') as aapdb:
		for x in aapdb:
			if x.startswith('CRYST') or x.startswith('ATOM'):
				vspdb.write(x)
		
	aapdb.close()

	lastatom = anum[-1]
	n = 1
	vsnum = lastatom + n
	#print(resatoms)
	vstoplines = ""

	for i in resid:
		aa=resname[i]
		serlist=[k for k,v in resatoms.items() if v == i]
		atomlist=[aname[k] for k in serlist]
		#print('last atom number: {0}'.format(vsnum))
		sitedef=VSdef(aa,atomlist,i,atomlist,serlist)
		#print('{0} - {1}'.format(aa,sitedef))
		#print('residues: {0} {1}'.format(i,aa))
		#print('serial: {0} and name {1}'.format(serlist,atomlist))
		n=sitedef['vsnum']
		VSwriteOut=VSwrite(i,sitedef,vsnum)
		vsout.write(VSwriteOut['vsitesnm'])
		vspdb.write(VSwriteOut['pdblines'])
		vstoplines += VSwriteOut['vstoplines']
		#print(lines)
		vsnum += n
	print("virtual sites generated!!!! \n")
	print('the files {0} and {1} are written' . format(outname,outpdb))
	vsout.write(";; the following section is for vsite topologies\n")
	vsout.write(vstoplines)	
	vsout.close()
	vspdb.close()

	#


	


# virtual site definitions are done here
def VSdef(AAres,pdblist,AAid,AAlist,SERlist):
	
	AAsidechain = {
		"ALA" : ["CB"],
		"ARG" : ["CB","CG","CD","NE","CZ","NH1","NH2"],
		"ASN" : ["CB","CG","OD1","ND2"],
		"ASP" : ["CB","CG","OD1","OD2"],
		"CYS" : ["CB","SG"],
		"GLU" : ["CB","CG","CD","OE1","OE2"],
		"GLN" : ["CB","CG","CD","OE1","NE2"],
		"GLY" : [],
		"HIS" : ["CB","CG","ND1","CE1","NE2","CD2"],
		"HSD" : ["CB","CG","ND1","CE1","NE2","CD2"],
		"HSE" : ["CB","CG","ND1","CE1","NE2","CD2"],
		"ILE" : ["CB","CG2","CG1","CD"],
		"LEU" : ["CB","CG","CD1","CD2"],
		"LYS" : ["CB","CG","CD","CE","NZ"],
		"LSN" : ["CB","CG","CD","CE","NZ"],
		"MET" : ["CB","CG","SD","CE"],
		"PHE" : ["CB","CG","CD1","CE1","CZ","CD2","CE2"],
		"PRO" : ["CB","CG"],
		"SER" : ["CB","OG"],
		"THR" : ["CB","OG1","CG2"],
		"TRP" : ["CB","CG","CD2","CD1","NE1","CE2","CE3","CZ3","CZ2","CH2"],
		"TYR" : ["CB","CG","CD1","CE1","CZ","OH","CD2","CE2"],
		"VAL" : ["CB","CG1","CG2"]
	}
	AAbb = ["N","CA","C","O","OT1","OT2"]
	
	#print('atoms from def: {0}'.format(AAfull))
	#print('atoms from pdb: {0}'.format(pdblist))
	#print('warning {0} are excluded'.format(excluded))


	vbb=[x for x in AAbb if x in pdblist]
	side=AAsidechain.get(AAres)
	vsc1=side
	vsc2=[]
	vsc3=[]
	vsc4=[]
	if AAres == "ALA":
		vbb = vbb + side
		vsc1 = []
	if AAres == "ARG":
		vsc1 = side[0:4]
		vsc2 = side[4:]
	if AAres == "HIS" or AAres=="HSD" or AAres=="HSE" or AAres =="HSP":
		vsc1 = side[0:2]
		vsc2 = side[2:4]
		vsc3 = side[4:]
	if AAres == "LYS" or AAres=="LSN":
		vsc1 = side[0:3]
		vsc2 = side[3:]
	if AAres == "PHE":
		vsc1 = side[0:3]
		vsc2 = side[3:5]
		vsc3 = side[5:]
	if AAres == "TRP":
		vsc1 = side[0:3]
		vsc2 = side[3:6]
		vsc3 = side[6:8]
		vsc4 = side[8:]
	if AAres == "TYR":
		vsc1 = side[0:3]
		vsc2 = side[6:]
		vsc3 = side[3:6]


	vsdef = {}
	vsdef = { "vBB" : vbb, "vSC1" : vsc1, "vSC2" : vsc2, "vSC3" : vsc3, "vSC4" : vsc4 }

	vsdef2 = {}
	servbb=[]
	for x in vbb:
		servbb.append(SERlist[AAlist.index(x)])

	servsc1=[]
	for x in vsc1:
		servsc1.append(SERlist[AAlist.index(x)])
	servsc2=[]
	for x in vsc2:
		servsc2.append(SERlist[AAlist.index(x)])
	servsc3=[]
	for x in vsc3:
		servsc3.append(SERlist[AAlist.index(x)])
	servsc4=[]
	for x in vsc4:
		servsc4.append(SERlist[AAlist.index(x)])
	
	count = 0
	if len(servbb) != 0:
		count += 1
	
	#print(AAres)
	#print('index of vBB: {0} {1} {2}'.format(AAlast,vsfunc,servbb))

	if len(servsc1) != 0:
		count += 1
	#print('index of vSC1: {0} {1} {2}'.format(AAlast,vsfunc,servsc1))

	if len(servsc2) != 0:
		count += 1
	#print('index of vSC2: {0} {1} {2}'.format(AAlast,vsfunc,servsc2))

	if len(servsc3) != 0:
		count += 1
	#print('index of vSC3: {0} {1} {2}'.format(AAlast,vsfunc,servsc3))

	if len(servsc4) != 0:
		count += 1
	#print('index of vSC4: {0} {1} {2}'.format(AAlast,vsfunc,servsc4))
	
	VSdefOut = { "vBB" : servbb, "vSC1" : servsc1, "vSC2" : servsc2, "vSC3" : servsc3, "vSC4" : servsc4, "vsnum" : count}
	#print(side)
	return VSdefOut

def VSwrite(AAres,VSites,num):
	count=VSites['vsnum']
	VSwriteOut={}
	lines=""
	pdblines=""
	vstoplines=""
	
	for i in range(0,count):
		VSser=num+i
		if i == 0:
			key = "vBB"
		elif i == 1:
			key = "vSC1"
		elif i == 2:
			key = "vSC2"
		elif i == 3:
			key = "vSC3"
		elif i == 4:
			key = "vSC4"
		vsatoms = ""
		VSAtoms = []
		for x in VSites[key]:
			vsatoms += " " + str(x)
			VSAtoms.append(x)

		lines += str(VSser) + "\t" + str(vsfunc) +"\t"+ vsatoms + "\n"
		pdblines += PDBwrite(VSser, VSAtoms, key)
		#print ("sending vstoplines += str(VStop(VSser, VSAtoms, key)) for {0} {1} {2}".format(VSser, VSAtoms, key))
		vstoplines += VStop(VSser, VSAtoms, key)
		
		
	VSwriteOut = {'vsitesnm':lines, 'pdblines':pdblines, 'vstoplines':vstoplines}	
	return VSwriteOut	


def PDBwrite(VSser,VSatoms, VSname):
	#print(allpdbinfo)
	#print("so far good")
	#print("VSatoms: ",VSatoms)
	xx,yy,zz=[],[],[]
	for i in VSatoms:
		#print ("i=",i)
		ii=int(i)
		#print (allpdbinfo[ii]['ResName'])
		xx.append(allpdbinfo[ii]['X'])
		yy.append(allpdbinfo[ii]['Y'])
		zz.append(allpdbinfo[ii]['Z'])
	
	#print ("x of vsite: ",VSser,VSname,":",xx) 	
	vscoord=com_calc(xx,yy,zz)
	vsx=vscoord[0]
	vsy=vscoord[1]
	vsz=vscoord[2]
	
	PDBwriteOut=format('%-6s' % 'ATOM')+ format('%5s' % str(VSser))+VSname.center(6,' ')+format('%-4s' % allpdbinfo[i]['ResName'])+allpdbinfo[i]['Chain']+format('%4s' % allpdbinfo[i]['ResNo'])+format('%12.3f' % vsx)+format('%8.3f' % vsy)+format('%8.3f' % vsz)+format('%6s' % '1.00')+format('%6s' % '00.00')+format('%12s' % 'X')+"\n"
	
	#print line
	
	return PDBwriteOut
		
	

def com_calc(xx,yy,zz):
	
	#x,y,z=0.0,0.0,0.0
	if Deb == "full":
		print ("com_calc module\n")
		print ("xx is: ",xx)
		print ("yy is: ",yy)
		print ("zz is: ",zz)

	com=[]
	#for i,j,k in xx,yy,zz:
	#	x+=float(i[0])
	#	y+=float(j[1])
	#	z+=float(k[2])

	count = len(xx)
	comx = sum(float(x) for x in xx if x)/count
	comy = sum(float(x) for x in yy if x)/count
	comz = sum(float(x) for x in zz if x)/count
	com.append(comx)
	com.append(comy)
	com.append(comz)
	#print com
	return com

def VStop(VSser,VSatoms, VSname):
	
	i=VSatoms[len(VSatoms)-1]
	
	getinfo=allpdbinfo[i]
	vs_resname=getinfo['ResName']
	vs_resno=getinfo['ResNo']
	vs_mass = 0
	vs_charge = 0
	chn = getinfo['Chain']
	
	#print ("NTER residues: ", NTER_res)
	for x in VSatoms:
		#print ("VStop x {}".format(x))
		resloc=""
		vs_atomname = allpdbinfo[x]['AtomName']
		if allpdbinfo[x]['NTER'] == "yes":
			resloc = "NTER"
			break
		if allpdbinfo[x]['CTER'] == "yes":
			resloc = "CTER"
			break

		#if vs_atomname == "OT1"  or vs_atomname == "OT2" or vs_atomname == "OXT":
		#	resloc = "CTER"
		#	break
		
	#print ("VS residue {0} is {1}" . format(vs_resno,resloc))
	# resloc tells if the residue is N-terminal or C-terminal
	
	vs_type = get_VStype(VSname, vs_resname, resloc)
	# virtual sites topology
	# S.NO	Type	ResID	Resname	AtomName	ChrgGrpNo	Mass	Charge
	vstop_arr = [VSser, vs_type, vs_resno, vs_resname, VSname, VSser, vs_mass, vs_charge]
	vstop_out = ""
	for x in vstop_arr:
		vstop_out += str(x) + "\t"

	vstop_out += "\n"
		
	return vstop_out
	

def get_VStype(VSname, vs_resname, resloc):
	
	vstype_def = {
	"ALA" : {'vBB': 'vP4'},
	"ARG" : {'vBB': 'vP5', 'vSC1':'vN0', 'vSC2':'vQd'},
	"ASN" : {'vBB': 'vP5', 'vSC1':'vNda'},
	"ASP" : {'vBB': 'vP5', 'vSC1':'vQa'},
	"CYS" : {'vBB': 'vP5', 'vSC1':'vC5'},
	"GLU" : {'vBB': 'vP5', 'vSC1':'vQa'},
	"GLN" : {'vBB': 'vP5', 'vSC1':'vNda'},
	"GLY" : {'vBB': 'vP5'},
	"HIS" : {'vBB': 'vP5', 'vSC1':'vSC4', 'vSC2':'vSP1', 'vSC3':'vSP1'},
	"HSD" : {'vBB': 'vP5', 'vSC1':'vSC4', 'vSC2':'vSP1', 'vSC3':'vSP1'},
	"HSE" : {'vBB': 'vP5', 'vSC1':'vSC4', 'vSC2':'vSP1', 'vSC3':'vSP1'},
	"ILE" : {'vBB': 'vP5', 'vSC1':'vC1'},
	"LEU" : {'vBB': 'vP5', 'vSC1':'vC1'},
	"LYS" : {'vBB': 'vP5', 'vSC1':'vC3', 'vSC2':'vQd'},
	"LSN" : {'vBB': 'vP5', 'vSC1':'vC3', 'vSC2':'vP1'},
	"MET" : {'vBB': 'vP5', 'vSC1':'vC5'},
	"PHE" : {'vBB': 'vP5', 'vSC1':'vSC5', 'vSC2':'vSC5', 'vSC3':'vSC5'},
	"PRO" : {'vBB': 'vP5', 'vSC1':'vC5'},
	"SER" : {'vBB': 'vP5', 'vSC1':'vN0'},
	"THR" : {'vBB': 'vP5', 'vSC1':'vNda'},
	"TRP" : {'vBB': 'vP5', 'vSC1':'vSC4', 'vSC2':'vSNd', 'vSC3':'vSC5', 'vSC4':'vSC5'},
	"TYR" : {'vBB': 'vP5', 'vSC1':'vSC4', 'vSC2':'vSC4', 'vSC3':'vSP1'},
	"VAL" : {'vBB': 'vP5', 'vSC1':'vC2'},
	}
	
	if vs_resname not in vstype_def:
		print ("Error: Resname  {0} not found in the CG def" . format(vs_resname))
		exit()

	if VSname not in vstype_def[vs_resname]:
		print ("Error: VS name  {0} not found in the CG def" . format(VSname))
		exit()
	
	vstype=vstype_def[vs_resname][VSname]
		
	if resloc == "NTER" and VSname == 'vBB':
		vstype = 'vQd'
		
	elif resloc == "CTER" and VSname == 'vBB':
		vstype = 'vQa'
		
	#print ("{0} {1} - {2} ({3})".format(vs_resname, VSname, vstype,resloc)) 
	return vstype

def get_VSindex(pdb,idx):
	groups={}
	with open(pdb,"r") as input:
    		for line in input:
			if line.startswith('ATOM'):
				serial=line[6:11].strip()
				atomname=line[11:17].strip()
				if 'System' in groups:
					groups['System'].append(serial)
				else:
					groups ['System'] = [serial]


				if atomname.startswith('v'):
					if 'VS' in groups:
						groups['VS'].append(serial)
					else:
						groups['VS'] = [serial]
				else:
					if 'AA' in groups:
						groups['AA'].append(serial)
					else:
						groups['AA'] = [serial]	 
					
    
	with open(idx,'w') as output:
		for x in groups:
			output.write("[ {} ]\n".format(x))
			ind=""
        		count = 0
        		for i in groups[x]:
				ind += i + " "
				count += 1
				if count==10:
					ind += "\n"
					count=0
			output.write(ind)
			output.write("\n")

	input.close()
	output.close()


def list_diff(list1, list2):
	out1 = [item for item in list1 if not item in list2]
	out2 = [item for item in list2 if not item in list1]
	return out1+out2

if __name__=="__main__":
	main()

	 
