from yggdrasil import runner
import argparse

parser = argparse.ArgumentParser(description='Run an integration.')
parser.add_argument('yamlfile', nargs='+',help='One or more yaml specification files.')
args1 = parser.parse_args(["../ePhotosynthesis/EphotosynthesisOnly.yml"])
args2 = parser.parse_args(["../FBA/yggrasil_ODE_FBA_testing.yaml"])
args3 = parser.parse_args(["../FBA/yggrasil_ODE_FBA_night.yaml"])

Data_PPFD = {0:9.278350515,2.2:732.4345757,4.5:1486.597938,7:1694.528152,
             9.1:1366.058684,11.2:1181.919112,13.4:18.47739889}

def ProcessPPFD(data):
    PPFD_list = list()
    values = list(Data_PPFD.values())
    keys = list(Data_PPFD.keys())
    for i in range(1, len(Data_PPFD)):
        m = (values[i]-values[i-1])/(keys[i]-keys[i-1])
        c = values[i] - (m*(keys[i]))
        for j in range(0,24):
            y = (m*j)+c
            if j < keys[i-1] or j >= keys[i]:
                continue
            if y < 100:
                #print(str(y)+"<100")
                continue
            PPFD_list.append(y)
    return PPFD_list

PPFD = ProcessPPFD(Data_PPFD)
Vc = list()
Vo = list()
Vpga = list()
Vt3p = list()
Vstarch = list()
Vglycerate = list()
Vglycolate = list()

for p in PPFD:

    #ensure additional chloroplastic ATP consumption rate (J_ATPase) starts at 0
    f1 = open("../ePhotosynthesis/InputATPCost.txt","w")
    f1.write("ATPCost 0")
    f1.close()

    F_weather = open("../ePhotosynthesis/InputEvn.txt","w")
    F_weather.write("CO2 372\nPPFD "+str(p)+"\nSucPath 0"+"\ndaylength "+str(len(PPFD)))
    F_weather.close()

    ModelConverged = False

    i=0

    while not ModelConverged:
        print("Loop "+str(i))

        runner.run(args1.yamlfile)

        #store J_ATPase value
        f2 = open("../ePhotosynthesis/InputATPCost.txt","r")
        line = f2.readline()
        J_ATPase1 = float(line.replace("\t"," ").split(" ")[1])
        f2.close()

        runner.run(args2.yamlfile)

        #store J_ATPase value
        f2 = open("../ePhotosynthesis/InputATPCost.txt","r")
        line = f2.readline()
        J_ATPase2 = float(line.replace("\t"," ").split(" ")[1])
        f2.close()

        i=i+1
        print("ODE ATPase "+str(J_ATPase1))
        print("FBA ATPase "+str(J_ATPase2))

        if round(J_ATPase1,2)==round(J_ATPase2,2):
            ModelConverged=True

    print("Models converged at "+str(J_ATPase1))

    import os
    os.rename("./Daytime_flux.csv","./Daytime_flux_FACE254_372.csv")

    F_fluxes = open("../ePhotosynthesis/OutputRate.txt")
    lines = F_fluxes.readlines()
    Vc.append(float(lines[1].split(",")[1]))
    Vo.append(float(lines[1].split(",")[2]))
    Vpga.append(float(lines[1].split(",")[3]))
    Vt3p.append(float(lines[1].split(",")[4]))
    Vstarch.append(float(lines[1].split(",")[5]))
    Vglycerate.append(float(lines[1].split(",")[6]))
    Vglycolate.append(float(lines[1].split(",")[7]))
    F_fluxes.close()

F_fluxes = open("../ePhotosynthesis/OutputRate.txt","w")
PPFD_avg = sum(PPFD)/len(PPFD)
Vc_avg = sum(Vc)/len(Vc)
Vo_avg = sum(Vo)/len(Vo)
Vpga_avg = sum(Vpga)/len(Vpga)
Vt3p_avg = sum(Vt3p)/len(Vt3p)
Vstarch_avg = sum(Vstarch)/len(Vstarch)
Vglycerate_avg = sum(Vglycerate)/len(Vglycerate)
Vglycolate_avg = sum(Vglycolate)/len(Vglycolate)
F_fluxes.write("Light intensity,Vc,Vo,VPGA,VT3P,Vstarch,Vt_glycerate,Vt_glycolate\n")
F_fluxes.write(str(PPFD_avg)+","+str(Vc_avg)+","+str(Vo_avg)+","+str(Vpga_avg)+","+str(Vt3p_avg)+","+str(Vstarch_avg)+","+str(Vglycerate_avg)+","+str(Vglycolate_avg))
F_fluxes.close()
runner.run(args3.yamlfile)

import os
os.rename("./../FBA/Nighttime_flux.csv","./../FBA/Daytime_flux_FACE254_372.csv")