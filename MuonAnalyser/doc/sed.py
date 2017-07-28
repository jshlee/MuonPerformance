import os, sys

if sys.argv[1:] == []: filename = [sys.argv[1]]
else: filename = sys.argv[1:]
hostname = os.environ['HOSTNAME']

if 'sdfarm' in hostname: outstr = "root:\/\/cms-xrdr.sdfarm.kr:1094\/\/\/xrd\/"
if 'lxplus' in hostname: outstr = "root:\/\/cms-xrd-global.cern.ch\/\/"

for name in filename:
    f = open(name, "r")
    line = f.readlines()[0]
    instr = line[:line.find("store")]
    instr = instr.replace("/","\/")

    os.system("sed -i 's/%s/%s/g' %s"%(instr,outstr,name))
    #print "sed -i 's/%s/%s/g' %s"%(instr,outstr,name)

