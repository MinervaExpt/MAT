import sys,os

auddir = sys.argv[1]
print "find %s -type f -name '*.log'"%(auddir)
files = os.popen("find %s -type f -name '*.log'|sort"%(auddir)).readlines()
summ = open("Summary_%s.txt"%(os.path.basename(auddir)),"w")

for f in files:
    filename = f.rstrip("\n")
    ff = open(filename,"r").readlines()
    if(len(ff)==0):
        summ.write("%s BAD EMPTY LOG\n"%(os.path.basename(filename)))
        print "%s BAD EMPTY LOG"%(os.path.basename(filename))
        continue
    if(ff[-1].find("GOOD")==-1):
        summ.write("%s BAD"%(os.path.basename(filename)))
        print os.path.basename(filename)

    
