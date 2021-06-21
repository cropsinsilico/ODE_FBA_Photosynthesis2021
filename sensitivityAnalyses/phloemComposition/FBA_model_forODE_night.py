def addAllantoinMetabolism(backup):
    model = backup.copy()
    ############
    # Nino-Gonzalez et al 2019
    # Chen et al 2006
    # Ritzel et al 2001
    # Takagi et al 2018
    met1 = Metabolite("S_ALLANTOIN_c",name="S-ALLANTOIN:(S)-(+)-allantoin",
                      formula="C4H6N4O3",compartment="c",
                      charge=0)
    amet1 = Metabolite("aS_ALLANTOIN_c",name="S-ALLANTOIN:(S)-(+)-allantoin",
                      formula="C4H5N4O3",compartment="c",
                      charge=-1)

    rxn1 = Reaction("Allantoin_tx",name="Allantoin uptake")
    rxn1.add_metabolites({model.metabolites.PROTON_e:-1,
                          model.metabolites.PROTON_c:1,met1:0.78,amet1:0.22})
    rxn1.gene_reaction_rule='Glyma.15G066400'
    rxn1.lower_bound = 0
    rxn1.upper_bound = 1000
    ##########
    met2 = Metabolite("S_ALLANTOIN_r",name="(S)-(+)-allantoin",
                      formula="C4H6N4O3",compartment="r",
                      charge=0)
    amet2 = Metabolite("aS_ALLANTOIN_r",name="(S)-(+)-allantoin",
                      formula="C4H5N4O3",compartment="r",
                      charge=-1)


    proton_R = model.metabolites.PROTON_c.copy()
    proton_R.id = "PROTON_r"
    proton_R.compartment = "r"
    rxnProton = Reaction("PROTON_rc",name="ER proton exchange")
    rxnProton.add_metabolites({proton_R:-1,model.metabolites.PROTON_c:1})
    rxnProton.lower_bound = -1000
    rxnProton.upper_bound = 1000

    rxn2 = Reaction("Allantoin_rc",name="Allantoin peroxisome uptake/efflux")
    rxn2.add_metabolites({met1:-0.78,amet1:-0.22,proton_R:-0.12,met2:0.9,amet2:0.1})
    rxn2.lower_bound = -1000
    rxn2.upper_bound = 1000
    ###########
    met3 = Metabolite("ALLANTOATE_r",name="allantoate",
                      formula="C4H7N4O4",compartment="r",
                      charge=-1)

    water_R = model.metabolites.WATER_c.copy()
    water_R.id = "WATER_r"
    water_R.compartment = "r"
    rxnWater = Reaction("H2O_rc",name="ER water exchange")
    rxnWater.add_metabolites({water_R:-1,model.metabolites.WATER_c:1})
    rxnWater.lower_bound = -1000
    rxnWater.upper_bound = 1000


    rxn3 = Reaction("ALLANTOINASE_RXN_r",name="ALLANTOINASE-RXN")
    rxn3.add_metabolites({met2:-0.9,amet2:-0.1,water_R:-1,
                          met3:1,proton_R:0.9})
    rxn3.gene_reaction_rule='Glyma.15G073000 or Glyma.15G072900 or Glyma.13G240500 or Glyma.13G240600'
    #Glyma.15G073000 - cytoplasmic in Uniprot
    #Glyma.15G072900 - cytoplasmic in Uniprot
    #Glyma.13G240500 - cytoplasmic in Uniprot
    #Glyma.13G240600 - cytoplasmic in Uniprot
    #But according to Takagi et al 2018 ER

    rxn3.lower_bound = 0
    rxn3.upper_bound = 1000
    ############
    met4 = Metabolite("CPD0_2298_r",name="CPD0-2298:(S)-ureidoglycine",
                      formula="C3H7N3O3",compartment="r",
                      charge=0)
    amet4 = Metabolite("aCPD0_2298_r",name="CPD0-2298:(S)-ureidoglycine",
                      formula="C3H6N3O3",compartment="r",
                      charge=-1)

    met5 = model.metabolites.AMMONIUM_c.copy()
    met5.id = "AMMONIUM_r"
    met5.compartment = "r"
    rxnNH4 = Reaction("NH4_rc",name="Ammonium ER exchange")
    rxnNH4.add_metabolites({met5:-1,model.metabolites.AMMONIUM_c:1})
    rxnNH4.lower_bound = -1000
    rxnNH4.upper_bound = 1000

    co2_R = model.metabolites.CARBON_DIOXIDE_c.copy()
    co2_R.id = "CARBON_DIOXIDE_r"
    co2_R.compartment = "r"
    rxnCO2 = Reaction("CO2_rc",name="CO2 ER exchange")
    rxnCO2.add_metabolites({co2_R:-1,model.metabolites.CARBON_DIOXIDE_c:1})
    rxnCO2.lower_bound = -1000
    rxnCO2.upper_bound = 1000

    rxn4 = Reaction("ALLANTOATE_DEIMINASE_RXN_r",name="ALLANTOATE-DEIMINASE-RXN:allantoate deiminase")
    rxn4.gene_reaction_rule='Glyma.15G156900 or Glyma.09G050800'
    rxn4.add_metabolites({met3:-1,proton_R:-1.72,water_R:-1,
                          met4:0.72,amet4:0.28,met5:1,co2_R:1})
    #Glyma.15G156900 -ER in Uniprot
    #Glyma.09G050800 -ER in Uniprot
    rxn4.lower_bound = 0
    rxn4.upper_bound = 1000
    #############
    met6 = Metabolite("CPD_1091_r",name="CPD-1091:(S)-ureidoglycolate",
                      formula="C3H5N2O4",compartment="r",
                      charge=-1)

    water_R = model.metabolites.WATER_c.copy()
    water_R.id = "WATER_r"
    water_R.compartment = "r"
    rxnWATER = Reaction("H2O_rc",name="H2O ER exchange")
    rxnWATER.add_metabolites({water_R:-1,model.metabolites.WATER_c:1})
    rxnWATER.lower_bound = -1000
    rxnWATER.upper_bound = 1000

    rxn5 = Reaction("URUR_RXN_r",name="URUR-RXN:(S)-ureidoglycine aminohydrolase")
    rxn5.gene_reaction_rule='Glyma.17G148400 or Glyma.05G066500'
    #Glyma.17G148400
    #Glyma.05G066500
    rxn5.add_metabolites({met4:-0.72,amet4:-0.28,proton_R:-0.28,water_R:-1,
                          met6:1,met5:1})
    rxn5.lower_bound = 0
    rxn5.upper_bound = 1000
    #############
    # met7 = model.metabolites.UREA_c.copy()
    # met7.id = "UREA_r"
    # met7.compartment="r"

    # rxn6 = Reaction("ALLANTOICASE_RXN_r",name="ALLANTOICASE-RXN:allantoicase")
    # rxn6.add_metabolites({met3:-1,model.metabolites.WATER_r:-1,
    #                       met6:1,met7:1})
    # rxn6.lower_bound = 0
    # rxn6.lower_bound = 1000
    #############

    glyox_R = model.metabolites.GLYOX_x.copy()
    glyox_R.id = "GLYOX_r"
    glyox_R.compartment = "r"
    rxnGlyox = Reaction("glyox_rx",name="glyoxylate ER-peroxisome exchange")
    rxnGlyox.add_metabolites({glyox_R:-1,model.metabolites.GLYOX_x:1})

    rxn7 = Reaction("UREIDOGLYCOLATE_HYDROLASE_RXN_r",name="UREIDOGLYCOLATE-HYDROLASE-RXN:ureidoglycolate amidohydrolase")
    rxn7.gene_reaction_rule='Glyma.20G205500 or Glyma.10G184900'
    #Glyma.20G205500 - ER in Uniprot
    #Glyma.10G184900 - ER in Uniprot
    rxn7.add_metabolites({met6:-1,proton_R:-2,water_R:-1,
                          met5:2,co2_R:1,glyox_R:1})
    rxn7.lower_bound = 0
    rxn7.upper_bound = 1000
    #############
    # rxn8 = Reaction("UREIDOGLYCOLATE_LYASE_RXN_r",name="UREIDOGLYCOLATE-LYASE-RXN:ureidoglycolate lyase")
    # rxn8.add_metabolites({met6:-1,
    #                       met7:1,model.metabolites.GLYOX_r:1})
    # rxn8.lower_bound = 0
    # rxn8.lower_bound = 1000
    ############
    model.add_reactions([rxn1,rxn2,rxn3,rxn4,rxn5,rxn7,
                         rxnCO2,rxnGlyox,rxnNH4,rxnProton,rxnWater])
    return model

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
f1 = open("../../ePhotosynthesis/InputEvn.txt")
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
cobra_model = io.sbml.read_sbml_model("./model_phloemUpdated.sbml")

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

#dont allow NH4 uptake
cobra_model.reactions.get_by_id("NH4_tx").lower_bound = 0
cobra_model.reactions.get_by_id("NH4_tx").upper_bound = 0
cobra_model.reactions.get_by_id("Nitrate_tx").lower_bound = 0
cobra_model.reactions.get_by_id("Nitrate_tx").upper_bound = 0
cobra_model = addAllantoinMetabolism(cobra_model)

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
ATPase = round(ATPase,2)
temp.reactions.get_by_id("ATPase_tx").lower_bound = ATPase
temp.reactions.get_by_id("ATPase_tx").upper_bound = ATPase

#constraint Starch degradation rate
daylength = weather["daylength"]
StarchDegradationRate = df["Vstarch"][0]*daylength/(24-daylength)
temp.reactions.get_by_id("STARCH_p_accumulation").lower_bound = StarchDegradationRate
temp.reactions.get_by_id("STARCH_p_accumulation").upper_bound = StarchDegradationRate

temp.reactions.get_by_id("MAL_v_accumulation").lower_bound = -1*0.71*StarchDegradationRate
temp.reactions.get_by_id("MAL_v_accumulation").upper_bound = -1*0.71*StarchDegradationRate

temp.reactions.get_by_id("CIT_v_accumulation").lower_bound = -1*-0.56*StarchDegradationRate
temp.reactions.get_by_id("CIT_v_accumulation").upper_bound = -1*-0.56*StarchDegradationRate

#check if model works
temp.solver="glpk"
sol = flux_analysis.parsimonious.optimize_minimal_flux(temp)
rxn =  temp.reactions.get_by_id("Phloem_output_tx")
met = temp.metabolites.sSUCROSE_b
print("Sucrose export rate ="+str(rxn.metabolites[met]*sol.fluxes[rxn.id]))
print("O2 uptake rate ="+str(sol.fluxes["O2_tx"]))

fout= open("./Nighttime_flux.csv","w")
for rxn in temp.reactions:
    fout.write(rxn.id+","+rxn.reaction+","+str(sol.fluxes[rxn.id])+"\n")
fout.close()
