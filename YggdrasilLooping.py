from yggdrasil import runner
import argparse

parser = argparse.ArgumentParser(description='Run an integration.')
parser.add_argument('yamlfile', nargs='+',help='One or more yaml specification files.')
args1 = parser.parse_args(["./ePhotosynthesis/EphotosynthesisOnly.yml"])
args2 = parser.parse_args(["yggrasil_ODE_FBA_testing.yaml"])

#ensure additional chloroplastic ATP consumption rate (J_ATPase) starts at 0
f1 = open("./ePhotosynthesis/InputATPCost.txt","w")
f1.write("ATPCost 0")
f1.close()

ModelConverged = False

i=0

while not ModelConverged:
    print("Loop "+str(i))

    runner.run(args1.yamlfile)

    #store J_ATPase value
    f2 = open("./ePhotosynthesis/InputATPCost.txt","r")
    line = f2.readline()
    J_ATPase1 = line.replace("\t"," ").split(" ")[1]
    f2.close()

    runner.run(args2.yamlfile)

    #store J_ATPase value
    f2 = open("./ePhotosynthesis/InputATPCost.txt","r")
    line = f2.readline()
    J_ATPase2 = line.replace("\t"," ").split(" ")[1]
    f2.close()

    i=i+1
    print("ODE ATPase "+str(J_ATPase1))
    print("FBA ATPase "+str(J_ATPase2))

    if round(J_ATPase1,3)==round(J_ATPase2,3):
        ModelConverged=True

print("Models converged at "+str(J_ATPase1))
