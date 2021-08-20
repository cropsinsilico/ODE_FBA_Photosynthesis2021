from yggdrasil import runner
import argparse

parser = argparse.ArgumentParser(description='Run an integration.')
parser.add_argument('yamlfile', nargs='+',help='One or more yaml specification files.')
args1 = parser.parse_args(["../../ePhotosynthesis/EphotosynthesisOnly.yml"])
args2 = parser.parse_args(["./yggrasil_ODE_FBA_testing.yaml"])
args3 = parser.parse_args(["./yggrasil_ODE_FBA_night.yaml"])


Vc = list()
Vo = list()
Vpga = list()
Vt3p = list()
Vstarch = list()
Vglycerate = list()
Vglycolate = list()

for MALmult in [0.50,1,1.5]:
	for CITmult in [0.5,1,1.5]:
		f0 = open("./MALmultiplier.txt","w")
		f0.write(str(MALmult))
		f0.close()

		f0 = open("./CITmultiplier.txt","w")
		f0.write(str(CITmult))
		f0.close()

		#ensure additional chloroplastic ATP consumption rate (J_ATPase) starts at 0
		f1 = open("../../ePhotosynthesis/InputATPCost.txt","w")
		f1.write("ATPCost 0")
		f1.close()

		F_weather = open("../../ePhotosynthesis/InputEvn.txt","w")
		F_weather.write("CO2 400\nPPFD 500\nSucPath 1\ndaylength 12")
		F_weather.close()

		f3 = open("../../ePhotosynthesis/InputNADPHCost.txt","w")
		f3.write("NADPHCost 0")
		f3.close()

		ModelConverged = False

		i=0

		while not ModelConverged:
		    print("Loop "+str(i))

		    runner.run(args1.yamlfile)

		    #store J_ATPase value
		    f2 = open("../../ePhotosynthesis/InputATPCost.txt","r")
		    line = f2.readline()
		    J_ATPase1 = float(line.replace("\t"," ").split(" ")[1])
		    f2.close()
		    f4 = open("../../ePhotosynthesis/InputNADPHCost.txt","r")
		    line = f4.readline()
		    J_NADPHox1 = float(line.replace("\t"," ").split(" ")[1])
		    f4.close()

		    runner.run(args2.yamlfile)

		    #store J_ATPase value
		    f2 = open("../../ePhotosynthesis/InputATPCost.txt","r")
		    line = f2.readline()
		    J_ATPase2 = float(line.replace("\t"," ").split(" ")[1])
		    f2.close()
		    f4 = open("../../ePhotosynthesis/InputNADPHCost.txt","r")
		    line = f4.readline()
		    J_NADPHox2 = float(line.replace("\t"," ").split(" ")[1])
		    f4.close()

		    i=i+1
		    print("ODE ATPase "+str(J_ATPase1))
		    print("FBA ATPase "+str(J_ATPase2))
		    print("ODE NADPHox "+str(J_NADPHox1))
		    print("FBA NADPHox "+str(J_NADPHox2))

		    if round(J_ATPase1,2)==round(J_ATPase2,2) and round(J_NADPHox1,2)==round(J_NADPHox2,2):
		        ModelConverged=True
		        print("Models converged at "+str(J_ATPase1))
		        print("Models converged at "+str(J_NADPHox1))
		    if J_ATPase2 == 0 and J_NADPHox1 == 0:
		        print("Breaking Loop to avoid infinite loop")
		        break

		import os

		col = str(int(MALmult*100))+"_"+str(int(CITmult*100))
		os.rename("./Daytime_flux.csv","./TC_FBAfluxes_"+str(col)+"_lowLight.csv")
		os.rename("./../../ePhotosynthesis/OutputFluxT.txt","./TC_ODEfluxes_Fig2A_"+str(col)+"_lowLight.csv")

		F_fluxes = open("../../ePhotosynthesis/OutputRate.txt")
		lines = F_fluxes.readlines()
		Vc.append(float(lines[1].split(",")[1]))
		Vo.append(float(lines[1].split(",")[2]))
		Vpga.append(float(lines[1].split(",")[3]))
		Vt3p.append(float(lines[1].split(",")[4]))
		Vstarch.append(float(lines[1].split(",")[5]))
		Vglycerate.append(float(lines[1].split(",")[6]))
		Vglycolate.append(float(lines[1].split(",")[7]))
		F_fluxes.close()

		runner.run(args3.yamlfile)
		os.rename("./Nighttime_flux.csv","./TC_FBAfluxes_"+str(col)+"_night_lowLight.csv")
