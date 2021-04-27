from yggdrasil import runner
import argparse

parser = argparse.ArgumentParser(description='Run an integration.')
parser.add_argument('yamlfile', nargs='+',help='One or more yaml specification files.')
args1 = parser.parse_args(["../ePhotosynthesis/EphotosynthesisOnly.yml"])


PPFD = [50,100,200,300,400,500,600,800,1000,1200,1500,2000]
Vc = list()
Vo = list()
Vpga = list()
Vt3p = list()
Vstarch = list()
Vglycerate = list()
Vglycolate = list()

j=0
for ppfd in PPFD:
    j=j+1

    #ensure additional chloroplastic ATP consumption rate (J_ATPase) starts at 0
    f1 = open("../ePhotosynthesis/InputATPCost.txt","w")
    f1.write("ATPCost 0")
    f1.close()

    F_weather = open("../ePhotosynthesis/InputEvn.txt","w")
    F_weather.write("CO2 400\nPPFD "+str(ppfd)+"\nSucPath 1"+"\ndaylength 12")
    F_weather.close()

    f3 = open("../ePhotosynthesis/InputNADPHCost.txt","w")
    f3.write("NADPHCost 0")
    f3.close()

    runner.run(args1.yamlfile)

    import os
    os.rename("./../ePhotosynthesis/OutputFluxT.txt","./../Validations/ODEonly_fluxes_Fig2B_"+str(ppfd)+".csv")
