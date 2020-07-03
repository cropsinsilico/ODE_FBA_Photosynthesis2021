
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
cobra_model = io.sbml.read_sbml_model("./../Data/PlantCoreMetabolism_v2_0_0.xml")


for rxn in cobra_model.reactions:
    if rxn.lower_bound == -1000:
        rxn.lower_bound = -3000
    if rxn.upper_bound == 1000:
        rxn.upper_bound = 3000


from cobra import flux_analysis
leaf_model = cobra_model.copy()

leaf_model.reactions.GLC_tx.upper_bound = 0
leaf_model.reactions.GLC_tx.lower_bound = 0
leaf_model.reactions.Sucrose_tx.upper_bound = 1
leaf_model.reactions.Sucrose_tx.lower_bound = 1
leaf_model.reactions.Photon_tx.upper_bound = 0
leaf_model.reactions.Photon_tx.lower_bound = 0


from sweetlovegroup.transform import setupC3DielModel
leaf_model = setupC3DielModel(leaf_model)

leaf_model.reactions.Phloem_output_tx1.objective_coefficient = 0
leaf_model.reactions.Phloem_output_tx2.objective_coefficient = 0
PPFD = df["Light intensity"][0]
leaf_model.reactions.Photon_tx1.upper_bound = PPFD
leaf_model.reactions.Photon_tx1.lower_bound = 0
ATPase = (0.0049*PPFD) + 2.7851
leaf_model.reactions.ATPase_tx1.upper_bound = df["Light intensity"][0]
leaf_model.reactions.ATPase_tx1.lower_bound = df["Light intensity"][0]

fin = open("./../ePhotosynthesis/OutputFluxT.txt")
for line in fin:
    if "vz_1" in line:
        PSII = float(line.replace(" ","").replace("vz_1",""))
leaf_model.reactions.ATPase_tx1.upper_bound = PSII
leaf_model.reactions.ATPase_tx1.lower_bound = PSII

met = leaf_model.metabolites.rubisco_bal_p1
rxn = leaf_model.reactions.RXN_961_p1
rxn.add_metabolites({met:-3})
rxn.add_metabolites({met:df["Vc"][0]/df["Vo"][0]})

from cobra.flux_analysis import pfba
sol = pfba(leaf_model)
print("Phloem export rate ="+str(sol.fluxes["diel_biomass"]*4))


fout= open("./Daytime_flux.csv","w")
for rxn in temp.reactions:
    fout.write(rxn.id+","+rxn.reaction+","+str(sol.fluxes[rxn.id])+"\n")
fout.close()
