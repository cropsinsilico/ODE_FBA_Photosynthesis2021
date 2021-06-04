
def custom_pFBA(model):
    sol = model.optimize()
    for rxn in model.reactions:
        if rxn.objective_coefficient != 0:
            rxn.lower_buond = round(sol.fluxes[rxn.id],3)
            rxn.upper_bound = round(sol.fluxes[rxn.id],3)
    from sweetlovegroup import FBA
    Irrev_model = FBA.rev2irrev(model)
    for rxn in Irrev_model.reactions:
        if rxn.upper_bound > 0:
            rxn.objective_coefficient = -1
        else:
            rxn.objective_coefficient = 1
    sol2 = Irrev_model.optimize()
    rxnSet=set()
    fluxDict = dict()
    for rxn in Irrev_model.reactions.query("_reverse"):
        rxnSet.add(rxn.id)
        rxnSet.add(rxn.id.replace("_reverse",""))
        fluxDict[rxn.id.replace("_reverse","")]=sol2.fluxes[rxn.id]+sol2.fluxes[rxn.id.replace("_reverse","")]
    for rxn in Irrev_model.reactions:
        if rxn.id in rxnSet:
            continue
        else:
            fluxDict[rxn.id]=sol2.fluxes[rxn.id]
    sol3 = sol2
    sol3.fluxes = fluxDict
    return sol3



def remove_metabolite_from_reaction(rxn,mets):
    '''
    This functions removes a list of metabolites from a reaction
    '''
    for met in mets:
        if met in rxn.metabolites.keys():
            coeff = rxn.metabolites.get(met)
            rxn.add_metabolites({met:-1*coeff})
        else:
            print("Metabolite "+met.id+" not present in reaction "+rxn.id)
    return rxn

import pandas as pd
from io import BytesIO
import sys
import logging
logging.basicConfig()
logger = logging.getLogger('logger')

# Import classes for input/output channels
from yggdrasil.interface.YggInterface import YggInput, YggOutput
# Initialize input/output channels
in_channel = YggInput('input1')

# Loop until there is no longer input or the queues are closed
while True:

    # Receive input from input channel
    # If there is an error, the flag will be False
    flag, msg = in_channel.recv()
    if not flag:
        print("No more input.")
        break
    else:
        df2 = pd.read_csv(BytesIO(msg))
        print(df2.to_string())


# Read day-length from Environemnt file
f1 = open("../ePhotosynthesis/InputEvn.txt")
weather = dict()
for line in f1:
    parts = line.split(" ")
    weather[parts[0]]=float(parts[1])
    print(parts[0])
    print(parts[1])
f1.close()

from cobra import io,flux_analysis
from cobra.core import Reaction, Metabolite

#import model. Update file name and location in the next line
cobra_model = io.sbml.read_sbml_model("./../Data/PlantCoreMetabolism_v2_0_0.xml")

#Remove all metabolites except sucrose from Phloem
rxn = cobra_model.reactions.get_by_id("Phloem_output_tx")
mets2remove = list()

#for met in rxn.metabolites.keys():
    #if "SUCROSE" in met.id:# or "GLC" in met.id or "FRU" in met.id:
    #    continue
    #else:
    #    mets2remove.append(met)
#    mets2remove.append(met)

#remove_metabolite_from_reaction(rxn,mets2remove)
#rxn.add_metabolites({cobra_model.metabolites.get_by_id("sSUCROSE_b"):-1})
#rxn.add_metabolites({cobra_model.metabolites.get_by_id("GAP_c"):-1})
#coeff = sum(rxn.metabolites.values())
#rxn.add_metabolites({cobra_model.metabolites.get_by_id("PROTON_c"):-1*coeff,cobra_model.metabolites.get_by_id("PROTON_e"):coeff})

#no external sucrose or glucose
cobra_model.reactions.get_by_id("Sucrose_tx").lower_bound = 0
cobra_model.reactions.get_by_id("Sucrose_tx").upper_bound = 0
cobra_model.reactions.get_by_id("GLC_tx").lower_bound = 0
cobra_model.reactions.get_by_id("GLC_tx").upper_bound = 0
cobra_model.reactions.get_by_id("NH4_tx").lower_bound=0
cobra_model.reactions.get_by_id("NH4_tx").upper_bound=0

#no external light energy
cobra_model.reactions.get_by_id("Photon_tx").lower_bound = 0
cobra_model.reactions.get_by_id("Photon_tx").upper_bound = 0

#set export of sugars as objective
cobra_model.reactions.get_by_id("Phloem_output_tx").objective_coefficient=0

#add source reaction for TP
rxn = Reaction("STARCH_p_accumulation",name = "Starch store")
rxn.add_metabolites({cobra_model.metabolites.get_by_id("STARCH_p"):1})
rxn.upper_bound = 1000
rxn.lower_bound = 0
cobra_model.add_reaction(rxn)


#remove mGS and cGS
cobra_model.reactions.get_by_id("GLUTAMINESYN_RXN_m").lower_bound = 0
cobra_model.reactions.get_by_id("GLUTAMINESYN_RXN_m").upper_bound = 0
cobra_model.reactions.get_by_id("GLUTAMINESYN_RXN_c").lower_bound = 0
cobra_model.reactions.get_by_id("GLUTAMINESYN_RXN_c").upper_bound = 0

#add malate and citrate accumulation reactions
rxn = Reaction("MAL_v_accumulation")
rxn.add_metabolites({cobra_model.metabolites.MAL_v:-0.7,cobra_model.metabolites.aMAL_v:-0.3})
rxn.lower_bound = -1000
rxn.upper_bound = 1000
cobra_model.add_reaction(rxn)
rxn = Reaction("CIT_v_accumulation")
rxn.add_metabolites({cobra_model.metabolites.CIT_v:-0.5,cobra_model.metabolites.aCIT_v:-0.5})
rxn.lower_bound = -1000
rxn.upper_bound = 1000
cobra_model.add_reaction(rxn)


#settingup model
model = cobra_model.copy()
model.metabolites.get_by_id("aL_1_PHOSPHATIDYL_GLYCEROL_P_m").remove_from_model()
model.metabolites.get_by_id("aL_1_PHOSPHATIDYL_GLYCEROL_P_p").remove_from_model()

model.reactions.get_by_id("PGPPHOSPHA_RXN_m").add_metabolites({model.metabolites.get_by_id("L_1_PHOSPHATIDYL_GLYCEROL_P_m"):-0.03,
                                                               model.metabolites.get_by_id("PROTON_m"):-0.03})
model.reactions.get_by_id("PHOSPHAGLYPSYN_RXN_m").add_metabolites({model.metabolites.get_by_id("L_1_PHOSPHATIDYL_GLYCEROL_P_m"):0.03,
                                                                   model.metabolites.get_by_id("PROTON_m"):0.03})
model.reactions.get_by_id("PGPPHOSPHA_RXN_p").add_metabolites({model.metabolites.get_by_id("L_1_PHOSPHATIDYL_GLYCEROL_P_p"):-0.03,
                                                               model.metabolites.get_by_id("PROTON_p"):-0.03})
model.reactions.get_by_id("PHOSPHAGLYPSYN_RXN_p").add_metabolites({model.metabolites.get_by_id("L_1_PHOSPHATIDYL_GLYCEROL_P_p"):0.03,
                                                                   model.metabolites.get_by_id("PROTON_p"):0.03})
model.reactions.get_by_id("LPG_biosynthesis_c").add_metabolites({model.metabolites.get_by_id("L_1_PHOSPHATIDYL_GLYCEROL_P_p"):-0.03})

import pandas as pd

df = pd.read_csv("./../Data/biomass_soy.csv")

FA=["PALMITATE_p","CPD_9245_p","CPD_17412_p","CPD_17291_p","STEARIC_ACID_p","OLEATE_CPD_p",
    "Octadecadienoate_p","LINOLENIC_ACID_p","ARACHIDIC_ACID_p","CPD_16709_p","DOCOSANOATE_p"]
FACP = {"PALMITATE_p":"Palmitoyl_ACPs_p",
        "CPD_9245_p":"Palmitoleoyl_ACP_p",
        "CPD_17412_p":"hexadecadienoate_ACP_p",
        "CPD_17291_p":"hexadecatrienoate_ACP_p",
        "STEARIC_ACID_p":"Stearoyl_ACPs_p",
        "OLEATE_CPD_p":"Oleoyl_ACPs_p",
        "Octadecadienoate_p":"Octadecadienoyl_ACP_p",
        "LINOLENIC_ACID_p":"Octadecatrienoyl_ACP_p",
        "ARACHIDIC_ACID_p":"Arachidoyl_ACPs_p",
        "CPD_16709_p":"Eicosenoyl_ACP_p",
        "DOCOSANOATE_p":"Behenoyl_ACPs_p"}


PLs = ["ACYL_SN_GLYCEROL_3P_p",
       "L_PHOSPHATIDATE_p","L_PHOSPHATIDATE_m","DIACYLGLYCEROL_p",
       "DIACYLGLYCEROL_r","Triacylglycerols_p","PHOSPHATIDYL_CHOLINE_r",
       "L_1_PHOSPHATIDYL_ETHANOLAMINE_r","L_1_PHOSPHATIDYL_GLYCEROL_p",
       "L_1_PHOSPHATIDYL_GLYCEROL_P_p","L_1_PHOSPHATIDYL_GLYCEROL_P_m",
       "L_1_PHOSPHATIDYL_GLYCEROL_m","2_Lysophosphatidylcholines_r",
       "Lysophosphatidylglycerols_r","CDPDIACYLGLYCEROL_p","CDPDIACYLGLYCEROL_m",
       "D_Galactosyl_12_diacyl_glycerols_p","Galactosyl_galactosyl_diacyl_glycerols_p"]


for met in PLs:
    met=model.metabolites.get_by_id(met)
    met.formula=""

def generateMissingFormula(model,debug=False):
    loop_counter = 0
    former = 0
    for met in model.metabolites:
        if met.formula == "" or met.formula == "NA":
            former = former +1
    latter = 1
    while True:
        loop_counter = loop_counter+1
        if debug:
            print("Loop = "+str(loop_counter))
        former = latter
        for rxn in model.reactions:
            count = 0
            for met in rxn.metabolites:
                if met.formula=="" or met.formula=="NA" or met.formula == None:
                    if met.formula == "NA" or met.formula == None:
                        met.formula = ""
                    count = count + 1
                    coeff = rxn.metabolites[met]
            if count == 1:
                unb = rxn.check_mass_balance()
                eqn = rxn.reaction
                eqn = " "+eqn+" "
                for met in rxn.metabolites.keys():
                    formula = met.formula
                    if formula == None:
                        formula = "0"
                        NF_list.add(rxn.id)
                    eqn=eqn.replace(" "+met.id+" ","("+formula+")")
                if debug:
                    print(eqn)
                    print(unb)
                for met in rxn.metabolites:
                    if met.formula == "":
                        tempForm = ""
                        for a in sorted(unb.keys()):
                            if a=="charge" or round(unb[a],2)==0:
                                continue
                            num = float(abs(unb[a]))/abs(coeff)
                            if str(round(num))==str(num):
                                tempForm = tempForm+a+str(int(round(num)))
                            else:
                                tempForm = tempForm+a+str(num)
                                #print(a)
                                #print(round(num)==num)
                                #print(round(num))
                                #print(num)
                                #print(type(round(num)))
                                #print(type(num))
                        met.formula = tempForm
                        if debug:
                            print(met.id)
                            print(tempForm)
        latter = 0
        for met in model.metabolites:
            if met.formula == "" or met.formula == "NA":
                latter = latter +1
        if former == latter:
            break


leaf_model = model.copy()

k = "leaf"
RXN1 = Reaction("Fatty_acid_mix_"+k)
RXN2 = Reaction("Fatty_acid_ACP_"+k)
tot = 0
for met in df["Unnamed: 0"]:
    #print met
    if met in FA:
        RXN1.add_metabolites({leaf_model.metabolites.get_by_id(met):-1*float(df[df["Unnamed: 0"]==met][k])})
        RXN2.add_metabolites({leaf_model.metabolites.get_by_id(FACP[met]):-1*float(df[df["Unnamed: 0"]==met][k])})
        tot = tot+(float(df[df["Unnamed: 0"]==met][k]))
if tot==0:
    RXN1.add_metabolites({leaf_model.metabolites.PALMITATE_p:-1})
    RXN2.add_metabolites({leaf_model.metabolites.Palmitoyl_ACPs_p:-1})
    tot = 1
RXN1.add_metabolites({leaf_model.metabolites.Fatty_Acids_p:tot})
RXN1.lower_bound = 0
RXN1.upper_bound = 1000
leaf_model.add_reaction(RXN1)

RXN2.add_metabolites({leaf_model.metabolites.Fatty_acyl_ACP_p:tot})
RXN2.lower_bound = 0
RXN2.upper_bound = 1000
leaf_model.add_reaction(RXN2)


generateMissingFormula(leaf_model)

rxn = Reaction("Biomass_leaf_tx")
for met in df["Unnamed: 0"]:
    if met in FA or float(df[df["Unnamed: 0"]==met][k])==0:
        continue
    rxn.add_metabolites({leaf_model.metabolites.get_by_id(met):-1*float(df[df["Unnamed: 0"]==met][k])})

rxn.lower_bound = 0
rxn.upper_bound = 1000
leaf_model.add_reaction(rxn)

rxn.objective_coefficient=1

################

#remove glutamine synthetase and glutamate synthase
#rxn = cobra_model.reactions.get_by_id("GLUTAMATE_SYNTHASE_FERREDOXIN_RXN_p")
#mets2remove=[cobra_model.metabolites.get_by_id("Reduced_ferredoxins_p"),cobra_model.metabolites.get_by_id("Oxidized_ferredoxins_p")]
#remove_metabolite_from_reaction(rxn,mets2remove)
#rxn = cobra_model.reactions.get_by_id("GLUTAMINESYN_RXN_p")
#mets2remove=[cobra_model.metabolites.get_by_id("ATP_p"),cobra_model.metabolites.get_by_id("aATP_p")]
#remove_metabolite_from_reaction(rxn,mets2remove)

#rxn = Reaction("NrefixationCostbypass")
#rxn.add_metabolites({cobra_model.metabolites.get_by_id("GLT_x"):1,cobra_model.metabolites.get_by_id("2_KETOGLUTARATE_x"):-1,cobra_model.metabolites.get_by_id("AMMONIUM_m"):-1})
#rxn.lower_bound = 0
#rxn.upper_bound = 1000
#cobra_model.add_reaction(rxn)

#provide energy for N fixation
#rxn = Reaction("NrefixationEnergy")
#rxn.add_metabolites({cobra_model.metabolites.get_by_id("ATP_p"):0.9,cobra_model.metabolites.get_by_id("aATP_p"):0.1,cobra_model.metabolites.get_by_id("ADP_p"):-0.8,cobra_model.metabolites.get_by_id("aADP_p"):-0.2,cobra_model.metabolites.get_by_id("Pi_p"):-1,cobra_model.metabolites.get_by_id("Reduced_ferredoxins_p"):2,cobra_model.metabolites.get_by_id("Oxidized_ferredoxins_p"):-2})
#rxn.lower_bound = 0
#rxn.upper_bound = 1000
#met = Metabolite("NrefixEnergyConstraint")
#rxn.add_metabolites({met:1})
#cobra_model.add_reaction(rxn)
#cobra_model.reactions.get_by_id("GCVMULTI_RXN_m").add_metabolites({met:-1})

#turn off phosphoserine transaminase
#cobra_model.reactions.get_by_id("PSERTRANSAM_RXN_p").lower_bound = 0
#cobra_model.reactions.get_by_id("PSERTRANSAM_RXN_p").upper_bound = 0

cobra_model = leaf_model.copy()


temp = cobra_model.copy()

PPFD = df2["Light intensity"][0]
#constrain maintenace
ATPase = (0.0049*PPFD) + 2.7851
temp.reactions.get_by_id("ATPase_tx").lower_bound = ATPase
temp.reactions.get_by_id("ATPase_tx").upper_bound = ATPase
temp.reactions.get_by_id("NADPHoxc_tx").lower_bound = ATPase/9
temp.reactions.get_by_id("NADPHoxc_tx").upper_bound = ATPase/9
temp.reactions.get_by_id("NADPHoxp_tx").lower_bound = ATPase/9
temp.reactions.get_by_id("NADPHoxp_tx").upper_bound = ATPase/9
temp.reactions.get_by_id("NADPHoxm_tx").lower_bound = ATPase/9
temp.reactions.get_by_id("NADPHoxm_tx").upper_bound = ATPase/9

#constraint Starch degradation rate
daylength = weather["daylength"]
StarchDegradationRate = df["Vstarch"][0]*daylength/(24-daylength)
temp.reactions.get_by_id("STARCH_p_accumulation").lower_bound = StarchDegradationRate
temp.reactions.get_by_id("STARCH_p_accumulation").upper_bound = StarchDegradationRate

temp.reactions.get_by_id("MAL_v_accumulation").lower_bound = -1*0.0698903487288*StarchDegradationRate
temp.reactions.get_by_id("MAL_v_accumulation").upper_bound = -1*0.0698903487288*StarchDegradationRate

temp.reactions.get_by_id("CIT_v_accumulation").lower_bound = -1*-0.056884259879*StarchDegradationRate
temp.reactions.get_by_id("CIT_v_accumulation").upper_bound = -1*-0.056884259879*StarchDegradationRate


#check if model works
temp.solver="glpk"
sol = flux_analysis.parsimonious.optimize_minimal_flux(temp)

rxn =  temp.reactions.get_by_id("Biomass_leaf_tx")
print("Biomass C accumulation rate ="+str(rxn.flux*595.7664799294965))
print("O2 uptake rate ="+str(sol.fluxes["O2_tx"]))

fout= open("./Nighttime_flux.csv","w")
for rxn in temp.reactions:
    fout.write(rxn.id+","+rxn.reaction+","+str(sol.fluxes[rxn.id])+"\n")
fout.close()
