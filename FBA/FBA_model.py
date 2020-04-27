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


from cobra import io,flux_analysis
from cobra.core import Reaction, Metabolite

#import model. Update file name and location in the next line
cobra_model = io.sbml.read_sbml_model("Data/PlantCoreMetabolism_v2_0_0.xml")

#Remove all metabolites except sucrose from Phloem
rxn = cobra_model.reactions.get_by_id("Phloem_output_tx")
mets2remove = list()

for met in rxn.metabolites.keys():
    #if "SUCROSE" in met.id:# or "GLC" in met.id or "FRU" in met.id:
    #    continue
    #else:
    #    mets2remove.append(met)
    mets2remove.append(met)

remove_metabolite_from_reaction(rxn,mets2remove)
rxn.add_metabolites({cobra_model.metabolites.get_by_id("sSUCROSE_b"):-1})
#rxn.add_metabolites({cobra_model.metabolites.get_by_id("GAP_c"):-1})
coeff = sum(rxn.metabolites.values())
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
rxn = Reaction("GAP_tx",name = "TP source")
rxn.add_metabolites({cobra_model.metabolites.get_by_id("GAP_c"):1})
rxn.upper_bound = 1000
rxn.lower_bound = 0
cobra_model.add_reaction(rxn)

#add source reaction for TP
rxn = Reaction("GLYCOLATE_tx",name = "Glycolate source")
rxn.add_metabolites({cobra_model.metabolites.get_by_id("GLYCOLLATE_c"):1})
rxn.upper_bound = 1000
rxn.lower_bound = 0
cobra_model.add_reaction(rxn)

#add source reaction for TP
rxn = Reaction("GLYCERATE_tx",name = "Glycerate sink")
rxn.add_metabolites({cobra_model.metabolites.get_by_id("GLYCERATE_c"):-1})
rxn.upper_bound = 1000
rxn.lower_bound = 0
cobra_model.add_reaction(rxn)

#remove mGS and cGS
cobra_model.reactions.get_by_id("GLUTAMINESYN_RXN_m").lower_bound = 0
cobra_model.reactions.get_by_id("GLUTAMINESYN_RXN_m").upper_bound = 0
cobra_model.reactions.get_by_id("GLUTAMINESYN_RXN_c").lower_bound = 0
cobra_model.reactions.get_by_id("GLUTAMINESYN_RXN_c").upper_bound = 0

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
rxn = Reaction("NrefixationEnergy")
rxn.add_metabolites({cobra_model.metabolites.get_by_id("ATP_p"):0.9,cobra_model.metabolites.get_by_id("aATP_p"):0.1,cobra_model.metabolites.get_by_id("ADP_p"):-0.8,cobra_model.metabolites.get_by_id("aADP_p"):-0.2,cobra_model.metabolites.get_by_id("Pi_p"):-1,cobra_model.metabolites.get_by_id("Reduced_ferredoxins_p"):2,cobra_model.metabolites.get_by_id("Oxidized_ferredoxins_p"):-2})
rxn.lower_bound = 0
rxn.upper_bound = 1000
met = Metabolite("NrefixEnergyConstraint")
rxn.add_metabolites({met:1})
cobra_model.add_reaction(rxn)
cobra_model.reactions.get_by_id("GCVMULTI_RXN_m").add_metabolites({met:-1})

#turn off phosphoserine transaminase
#cobra_model.reactions.get_by_id("PSERTRANSAM_RXN_p").lower_bound = 0
#cobra_model.reactions.get_by_id("PSERTRANSAM_RXN_p").upper_bound = 0


cobra_model.reactions.get_by_id("Pi_ec").lower_bound = -1000
cobra_model.reactions.get_by_id("Pi_ec").upper_bound = 1000

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


#check if model works
sol = flux_analysis.parsimonious.optimize_minimal_flux(cobra_model)


VSuc_list = list()
for i in range(0,len(df)):
    temp = cobra_model.copy()
    PPFD = df["Light intensity"][i]
    #constrain maintenace
    ATPase = (0.0049*PPFD) + 2.7851
    temp.reactions.get_by_id("ATPase_tx").lower_bound = ATPase
    temp.reactions.get_by_id("ATPase_tx").upper_bound = ATPase

    #constraint TP flux
    temp.reactions.get_by_id("GAP_tx").lower_bound = df["VT3P"][i]
    temp.reactions.get_by_id("GAP_tx").upper_bound = df["VT3P"][i]

    #constraint glycollate and glycerate fluxes flux
    temp.reactions.get_by_id("GLYCOLATE_tx").lower_bound = df["Vt_glycolate"][i]
    temp.reactions.get_by_id("GLYCOLATE_tx").upper_bound = df["Vt_glycolate"][i]
    temp.reactions.get_by_id("GLYCERATE_tx").lower_bound = df["Vt_glycerate"][i]
    temp.reactions.get_by_id("GLYCERATE_tx").upper_bound = df["Vt_glycerate"][i]

    #temp.reactions.get_by_id("NrefixationCostbypass").lower_bound = df270["Vt_glycolate"][i]
    #temp.reactions.get_by_id("NrefixationCostbypass").upper_bound = df270["Vt_glycolate"][i]

    #temp.reactions.get_by_id("NrefixationEnergy").lower_bound = df270["Vt_glycerate"][i]
    #temp.reactions.get_by_id("NrefixationEnergy").upper_bound = df270["Vt_glycerate"][i]

    temp.reactions.get_by_id("MAL_v_accumulation").lower_bound = 0.0698903487288*df["Vstarch"][i]
    temp.reactions.get_by_id("MAL_v_accumulation").upper_bound = 0.0698903487288*df["Vstarch"][i]

    temp.reactions.get_by_id("CIT_v_accumulation").lower_bound = -0.056884259879*df["Vstarch"][i]
    temp.reactions.get_by_id("CIT_v_accumulation").upper_bound = -0.056884259879*df["Vstarch"][i]

    #check if model works
    sol = flux_analysis.parsimonious.optimize_minimal_flux(temp)
    VSuc = temp.reactions.get_by_id("Phloem_output_tx").x
    #print(PPFD)
    #print(temp.reactions.get_by_id("Phloem_output_tx").x)
    #print(temp.reactions.get_by_id("GLYCOLATE_tx").x)
    #print(temp.reactions.get_by_id("GLYCERATE_tx").x)
    #print("------------")
    VSuc_list.append(VSuc)
df["VSuc"] = VSuc_list


import matplotlib.pyplot as plt
#plt.plot(df270["Light intensity"])
print(df.to_string())
plt.plot(df["Hour"],df["Vc"],label="Vc")
plt.plot(df["Hour"],df["Vo"],label="Vo")
plt.plot(df["Hour"],df["VT3P"],label="VT3P")
plt.plot(df["Hour"],df["VSuc"],label="VSuc")
plt.plot(df["Hour"],df["Vstarch"],label="Vstarch")
plt.xlabel("Hour")
plt.ylabel("Flux")
plt.legend(loc="center")
plt.show()
plt.close()