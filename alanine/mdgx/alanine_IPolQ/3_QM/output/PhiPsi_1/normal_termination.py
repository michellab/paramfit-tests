import os,sys,glob


allvacu = glob.glob("*/*/*out.vacu")
allsolv = glob.glob("*/*/*out.solv")


problems = open("problem.dat","w")
for vacu in allvacu:
    with open(vacu) as f:
        last_line = f.readlines()[-1]
        if "Normal termination" in last_line:
            pass
        else:
            problems.write("%s" % vacu)

for solv in allsolv:
    with open(solv) as f:
        last_line = f.readlines()[-1]
        if "Normal termination" in last_line:
            pass
        else:
            problems.write("%s" % vacu)
