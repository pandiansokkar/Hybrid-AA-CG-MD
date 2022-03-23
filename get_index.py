# this script can be used to generate indices for gromacs MD simulation.
# especially for virtual sites generated for AA/CG simulation

pdb = 'vsite.pdb'
idx = 'index.ndx'
with open(pdb,"r") as input:
    data={}
    seglist = ['System']

    segments = {}
    segments['System']=[]
    for line in input:
        all=line.split()
        if all[0] == "ATOM":
            serial=all[1]
            segments['System'].append(serial)
            segname=all[len(all)-2]
            if segname in segments:
                segments[segname].append(serial)
            else:
                segments[segname] = [serial]
            
            if segname not in seglist:
                seglist.append(segname)

    #print seglist
    #print segments['HETA']
    
with open(idx,'w') as output:

    for x in segments:
        
        output.write("[ {} ]\n".format(x))
        ind=""
        count = 0
        for i in segments[x]:
            ind += i + " "
            count += 1
            if count==10:
                ind += "\n"
                count=0
        output.write(ind)
        output.write("\n")


input.close()
output.close()

