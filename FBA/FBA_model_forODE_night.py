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
        df = pd.read_csv(BytesIO(msg))
        print(df.to_string())

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

#no external light energy
cobra_model.reactions.get_by_id("Photon_tx").lower_bound = 0
cobra_model.reactions.get_by_id("Photon_tx").upper_bound = 0

#set export of sugars as objective
cobra_model.reactions.get_by_id("Phloem_output_tx").objective_coefficient=1

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


temp = cobra_model.copy()

PPFD = df["Light intensity"][0]
#constrain maintenace
ATPase = (0.0049*PPFD) + 2.7851
ATPase = round(ATPase,3)
temp.reactions.get_by_id("ATPase_tx").lower_bound = ATPase
temp.reactions.get_by_id("ATPase_tx").upper_bound = ATPase

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
#temp.solver="glpk"
sol = flux_analysis.parsimonious.optimize_minimal_flux(temp)
rxn =  temp.reactions.get_by_id("Phloem_output_tx")
met = temp.metabolites.sSUCROSE_b
print("Sucrose export rate ="+str(rxn.metabolites[met]*rxn.flux))
print("O2 uptake rate ="+str(temp.reactions.O2_tx.flux))

fout= open("./Nighttime_flux.csv","w")
for rxn in temp.reactions:
    fout.write(rxn.id+","+rxn.reaction+","+str(rxn.flux)+"\n")
fout.close()
