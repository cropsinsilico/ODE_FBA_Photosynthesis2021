
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

########################################################
#This function was used to set up a C3 leaf diel model #
########################################################
def setupC3DielModel(core_model,transferMets="",starch_sucrose_ratio=None):
    '''
    This function can be used to generate a fully constrained diel C3 leaf model
    from a core model.
    Inputs: 1) a cobra model 2) a list of metabolites allowed to accumulate 3)
    starch to sucrose accumulation rate ratio
    Outputs: a fully constrained diel C3 leaf model
    '''
    from cobra.core import Metabolite, Reaction
    import re

    #create two copies of model elements for day and night
    tempCompDict = dict()
    for comp in core_model.compartments:
        tempCompDict[comp+"1"] = core_model.compartments[comp]+" day"
        tempCompDict[comp+"2"] = core_model.compartments[comp]+" night"

    cobra_model2 = core_model.copy()
    for met in cobra_model2.metabolites:
        met.id = met.id+"1"
        met.compartment = met.compartment+"1"
    for rxn in cobra_model2.reactions:
        rxn.id = rxn.id+"1"

    cobra_model3 = core_model.copy()
    for met in cobra_model3.metabolites:
        met.id = met.id+"2"
        met.compartment = met.compartment+"2"
    for rxn in cobra_model3.reactions:
        rxn.id = rxn.id+"2"

    #merge the day and night model
    cobra_model = cobra_model2+cobra_model3
    for met in cobra_model3.metabolites:
        if not cobra_model.metabolites.__contains__(met.id):
            cobra_model.add_metabolites(met.copy())
    for comp in cobra_model.compartments:
        cobra_model.compartments[comp] = tempCompDict[comp]

    met1 = Metabolite("X_Phloem_contribution_t1",name="Phloem output during the day",compartment="b1")
    cobra_model.reactions.get_by_id("Phloem_output_tx1").add_metabolites({met1:1})
    met2 = Metabolite("X_Phloem_contribution_t2",name="Phloem output during at night",compartment="b1")
    cobra_model.reactions.get_by_id("Phloem_output_tx2").add_metabolites({met2:1})

    rxn = Reaction("diel_biomass")
    rxn.add_metabolites({met1:-3,met2:-1})
    rxn.lower_bound = 0
    rxn.upper_bound = 1000
    cobra_model.add_reaction(rxn)

    #Adding reactions to allow for day-night metabolite accumulations
    if transferMets!="":
        tmfile = open(transferMets,"r")
        tmset=set()
        for line in tmfile:
            tmset.add(line.replace("\n",""))
    else:
        tmset=set(["STARCH_p","SUCROSE_v","MAL_v","aMAL_v","NITRATE_v","CIT_v",
        "aCIT_v","GLN_v","ASN_v","SER_v","GLN_v","GLY_v","THR_v","L_ALPHA_ALANINE_v",
        "4_AMINO_BUTYRATE_v","VAL_v","ILE_v","PHE_v","LEU_v","LYS_v","ARG_v",
        "L_ASPARTATE_v","GLT_v","HIS_v","bHIS_v","MET_v","PRO_v","TRP_v","TYR_v",
        "CYS_v","FRUCTAN_v","AMMONIUM_v"])

    for met in tmset:
        if met == "AMMONIUM_v" or met=="FRUCTAN_v":
            continue
        tempRxn = Reaction(met+"_dielTransfer")
        tempRxn.add_metabolites({cobra_model.metabolites.get_by_id(met+"1"):-1,cobra_model.metabolites.get_by_id(met+"2"):1})
        tempRxn.lower_bound=-1000
        if not ((met == "STARCH_p") or (met == "SUCROSE_v") or (met == "MAL_v") or (met == "aMAL_v") or (met == "NITRATE_v") or (met == "CIT_v") or (met == "aCIT_v") or (met == "PROTON_v")):
            tempRxn.lower_bound=0
        tempRxn.upper_bound=1000
        cobra_model.add_reaction(tempRxn)

    fractionMets=dict()
    for rxn in cobra_model.reactions:
        for met in rxn.metabolites.keys():
            prefix=""
            a=re.search("^a{1,3}",met.id)
            anion=""
            if a:
                anion=a.group(0)
                prefix=anion
            b=re.search("^b{1,3}",met.id)
            basic=""
            if b:
                basic=b.group(0)
                prefix=basic
            if ((not prefix == "") and met.compartment == "v1"):
                fractionMets[met]=prefix

    temp=cobra_model.copy()
    for met in fractionMets.keys():
        for rxn in met.reactions:
            if rxn.id.__contains__("_dielTransfer"):
                continue
            else:
                mainMet = met.id[len(fractionMets[met]):]
                coeff1 = temp.reactions.get_by_id(rxn.id).metabolites.get(temp.metabolites.get_by_id(mainMet))
                coeff2 = temp.reactions.get_by_id(rxn.id).metabolites.get(temp.metabolites.get_by_id(met.id))
                if not coeff1:
                    coeff1=0
                if not coeff2:
                    coeff2=0
                total = coeff1 + coeff2
                coeff1 = float(coeff1)/total
                coeff2 = float(coeff2)/total
                if cobra_model.reactions.has_id(met.id[0:len(met.id)-1]+"_dielTransfer"):
                    ub = temp.reactions.get_by_id(met.id[0:len(met.id)-1]+"_dielTransfer").upper_bound
                    lb = temp.reactions.get_by_id(met.id[0:len(met.id)-1]+"_dielTransfer").lower_bound
                    temp.reactions.get_by_id(met.id[0:len(met.id)-1]+"_dielTransfer").remove_from_model()
                    temp.reactions.get_by_id(mainMet[0:len(mainMet)-1]+"_dielTransfer").remove_from_model()
                    Reac = Reaction(mainMet[0:len(mainMet)-1]+"_dielTransfer",name=mainMet+"_dielTransfer")
                    Reac.add_metabolites({temp.metabolites.get_by_id(met.id[0:len(met.id)-1]+"1"):-coeff2,temp.metabolites.get_by_id(met.id[0:len(met.id)-1]+"2"):coeff2,temp.metabolites.get_by_id(mainMet[0:len(mainMet)-1]+"1"):-coeff1,temp.metabolites.get_by_id(mainMet[0:len(mainMet)-1]+"2"):coeff1})
                    Reac.lower_bound=lb
                    Reac.upper_bound=ub
                    temp.add_reaction(Reac)
                    print(Reac.reaction)
                break
    ####ADD CONSTRAINTS TO MODEL####
    cobra_model = temp.copy()

    #objective function
    cobra_model.reactions.get_by_id("diel_biomass").objective_coefficient=1
    #Leaves - light
    cobra_model.reactions.get_by_id("Sucrose_tx1").lower_bound=0
    cobra_model.reactions.get_by_id("Sucrose_tx1").upper_bound=0
    cobra_model.reactions.get_by_id("GLC_tx1").lower_bound=0
    cobra_model.reactions.get_by_id("GLC_tx1").upper_bound=0
    cobra_model.reactions.get_by_id("CO2_tx1").lower_bound=0
    cobra_model.reactions.get_by_id("NH4_tx1").lower_bound=0
    cobra_model.reactions.get_by_id("NH4_tx1").upper_bound=0
    #Leaves - dark
    cobra_model.reactions.get_by_id("Sucrose_tx2").lower_bound=0
    cobra_model.reactions.get_by_id("Sucrose_tx2").upper_bound=0
    cobra_model.reactions.get_by_id("GLC_tx2").lower_bound=0
    cobra_model.reactions.get_by_id("GLC_tx2").upper_bound=0
    cobra_model.reactions.get_by_id("Photon_tx2").lower_bound=0
    cobra_model.reactions.get_by_id("Photon_tx2").upper_bound=0
    cobra_model.reactions.get_by_id("NH4_tx2").lower_bound=0
    cobra_model.reactions.get_by_id("NH4_tx2").upper_bound=0
    cobra_model.reactions.get_by_id("CO2_tx2").upper_bound=0

    #Set pG6P transporter to 0
    cobra_model.reactions.get_by_id("G6P_Pi_pc1").lower_bound=0
    cobra_model.reactions.get_by_id("G6P_Pi_pc1").upper_bound=0
    cobra_model.reactions.get_by_id("G6P_Pi_pc2").lower_bound=0
    cobra_model.reactions.get_by_id("G6P_Pi_pc2").upper_bound=0

    #Turn off PTOX
    cobra_model.reactions.get_by_id("Plastoquinol_Oxidase_p1").lower_bound=0
    cobra_model.reactions.get_by_id("Plastoquinol_Oxidase_p1").upper_bound=0

    #nitrate uptake constrain
    Nitrate_balance = Metabolite("Nitrate_bal_c", name = "Weights to balance nitrate uptake", compartment = "c1")
    cobra_model.reactions.get_by_id("Nitrate_ec1").add_metabolites({Nitrate_balance:-2})
    cobra_model.reactions.get_by_id("Nitrate_ec2").add_metabolites({Nitrate_balance:3})

    #Rubisco balance
    Rubisco_balance = Metabolite("rubisco_bal_p1", name = "Weights to balance RuBP carboxygenase oxygenase balance", compartment = "p1")
    cobra_model.reactions.get_by_id("RXN_961_p1").add_metabolites({Rubisco_balance:3})
    cobra_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p1").add_metabolites({Rubisco_balance:-1})

    #generic ATPase and NADPH oxidase
    Maintenance_constraint = Metabolite("ATPase_NADPHoxidase_constraint_c1",name =  "ATPase_NADPHoxidase_constraint_c1", compartment = "c1")
    Maintenance_constraint2 = Metabolite("ATPase_NADPHoxidase_constraint_c2",name =  "ATPase_NADPHoxidase_constraint_c2", compartment = "c2")
    Maintenance_constraint3 = Metabolite("Light_dark_maintainence_constraint",name =  "Light_dark_maintainence_constraint", compartment = "c1")
    cobra_model.reactions.get_by_id("ATPase_tx1").add_metabolites({Maintenance_constraint:1,Maintenance_constraint3:1})
    cobra_model.reactions.get_by_id("ATPase_tx2").add_metabolites({Maintenance_constraint2:1,Maintenance_constraint3:-1})
    cobra_model.reactions.get_by_id("NADPHoxc_tx1").add_metabolites({Maintenance_constraint:-3})
    cobra_model.reactions.get_by_id("NADPHoxc_tx2").add_metabolites({Maintenance_constraint2:-3})
    cobra_model.reactions.get_by_id("NADPHoxm_tx1").add_metabolites({Maintenance_constraint:-3})
    cobra_model.reactions.get_by_id("NADPHoxm_tx2").add_metabolites({Maintenance_constraint2:-3})
    cobra_model.reactions.get_by_id("NADPHoxp_tx1").add_metabolites({Maintenance_constraint:-3})
    cobra_model.reactions.get_by_id("NADPHoxp_tx2").add_metabolites({Maintenance_constraint2:-3})

    ##constrain sucrose and starch storage
    if starch_sucrose_ratio is not None:
        Sucorse_starch_balance = Metabolite("sucrose_starch_bal_c", name = "Weights to balance sucrose-starch uptake", compartment = "c1")
        cobra_model.reactions.get_by_id("SUCROSE_v_dielTransfer").add_metabolites({Sucorse_starch_balance:-1*starch_sucrose_ratio})
        cobra_model.reactions.get_by_id("STARCH_p_dielTransfer").add_metabolites({Sucorse_starch_balance:1})

    #Plastid enolase was not detected in Arabidopsis mesophyll tissue
    cobra_model.reactions.get_by_id("2PGADEHYDRAT_RXN_p1").lower_bound=0
    cobra_model.reactions.get_by_id("2PGADEHYDRAT_RXN_p1").upper_bound=0
    cobra_model.reactions.get_by_id("2PGADEHYDRAT_RXN_p2").lower_bound=0
    cobra_model.reactions.get_by_id("2PGADEHYDRAT_RXN_p2").upper_bound=0

    #Setting chloroplastic NADPH dehydrogenase to 0  ((Yamamoto et al., 2011)
    cobra_model.reactions.get_by_id("NADPH_Dehydrogenase_p1").lower_bound=0
    cobra_model.reactions.get_by_id("NADPH_Dehydrogenase_p1").upper_bound=0
    cobra_model.reactions.get_by_id("NADPH_Dehydrogenase_p2").lower_bound=0
    cobra_model.reactions.get_by_id("NADPH_Dehydrogenase_p2").upper_bound=0

    #ATP_ADP_Pi constrained to 0 because while there is evidence for its existance, it does not carry high flux
    cobra_model.reactions.get_by_id("ATP_ADP_Pi_pc1").lower_bound = 0
    cobra_model.reactions.get_by_id("ATP_ADP_Pi_pc1").upper_bound = 0
    cobra_model.reactions.get_by_id("ATP_ADP_Pi_pc2").lower_bound = 0
    cobra_model.reactions.get_by_id("ATP_ADP_Pi_pc2").upper_bound = 0

    return cobra_model

import pandas as pd
from io import BytesIO
import sys
import logging
logging.basicConfig()
logger = logging.getLogger('logger')

fin = open("./../ePhotosynthesis/InputEvn.txt")
for line in fin:
    if "PPFD" in line:
        PPFD = int(line.split(" ")[1])



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
leaf_model.reactions.Sucrose_tx.upper_bound = 0
leaf_model.reactions.Sucrose_tx.lower_bound = 0
leaf_model.reactions.Photon_tx.upper_bound = 0
leaf_model.reactions.Photon_tx.lower_bound = 0


leaf_model = setupC3DielModel(leaf_model)

leaf_model.reactions.Phloem_output_tx1.objective_coefficient = 0
leaf_model.reactions.Phloem_output_tx2.objective_coefficient = 0

leaf_model.reactions.Photon_tx1.upper_bound = PPFD
leaf_model.reactions.Photon_tx1.lower_bound = 0
ATPase = (0.0049*PPFD) + 2.7851
leaf_model.reactions.ATPase_tx1.upper_bound = ATPase
leaf_model.reactions.ATPase_tx1.lower_bound = ATPase


from cobra.flux_analysis import pfba
sol = pfba(leaf_model)
print("Phloem export rate ="+str(sol.fluxes["diel_biomass"]*4))


fout= open("./Diel_flux.csv","w")
for rxn in leaf_model.reactions:
    fout.write(rxn.id+","+rxn.reaction+","+str(sol.fluxes[rxn.id])+"\n")
fout.close()
