from yggdrasil import runner
import argparse

parser = argparse.ArgumentParser(description='Run an integration.')
parser.add_argument('yamlfile', nargs='+',help='One or more yaml specification files.')
args1 = parser.parse_args(["../FBA/yggrasil_ODE_FBA_dielFBA_only.yaml"])


PPFD = [50,100,200,300,400,500,600,800,1000,1200,1500,2000]
Vc = list()
Vo = list()
Vpga = list()
Vt3p = list()
Vstarch = list()
Vglycerate = list()
Vglycolate = list()

j=0
for p in PPFD:
    j=j+1

    F_weather = open("../ePhotosynthesis/InputEvn.txt","w")
    F_weather.write("CO2 400\nPPFD "+str(p)+"\nSucPath 1"+"\ndaylength 12")
    F_weather.close()

    runner.run(args1.yamlfile)

    import os
    os.rename("./../FBA/Diel_flux.csv","./../Validations/FBAfluxes_Fig2B_"+str(p)+".csv")
