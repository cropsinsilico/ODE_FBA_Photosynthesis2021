
#################################################################################
# This function is a modified version of cobrapy pfba function			#
#										#
#################################################################################

import logging
from warnings import warn
from itertools import chain

from optlang.symbolics import Zero

from cobra.util import solver as sutil
from cobra.core.solution import get_solution

def pfba_Weighted(model, weightings, fraction_of_optimum=1.0, objective=None, reactions=None):
    """Perform basic pFBA (parsimonious Enzyme Usage Flux Balance Analysis)
    to minimize total flux.
    pFBA [1] adds the minimization of all fluxes the the objective of the
    model. This approach is motivated by the idea that high fluxes have a
    higher enzyme turn-over and that since producing enzymes is costly,
    the cell will try to minimize overall flux while still maximizing the
    original objective function, e.g. the growth rate.
    Parameters
    ----------
    model : cobra.Model
        The model
    fraction_of_optimum : float, optional
        Fraction of optimum which must be maintained. The original objective
        reaction is constrained to be greater than maximal_value *
        fraction_of_optimum.
    objective : dict or model.problem.Objective
        A desired objective to use during optimization in addition to the
        pFBA objective. Dictionaries (reaction as key, coefficient as value)
        can be used for linear objectives.
    reactions : iterable
        List of reactions or reaction identifiers. Implies `return_frame` to
        be true. Only return fluxes for the given reactions. Faster than
        fetching all fluxes if only a few are needed.
    Returns
    -------
    cobra.Solution
        The solution object to the optimized model with pFBA constraints added.
    References
    ----------
    .. [1] Lewis, N. E., Hixson, K. K., Conrad, T. M., Lerman, J. A.,
       Charusanti, P., Polpitiya, A. D., Palsson, B. O. (2010). Omic data
       from evolved E. coli are consistent with computed optimal growth from
       genome-scale models. Molecular Systems Biology, 6,
       390. doi:10.1038/msb.2010.47
    """
    reactions = model.reactions if reactions is None \
        else model.reactions.get_by_any(reactions)
    with model as m:
        add_pfba_Weighted(m, weightings, objective=objective,
                 fraction_of_optimum=fraction_of_optimum)
        m.slim_optimize(error_value=None)
        solution = get_solution(m, reactions=reactions)
    return solution


#################################################################################
# This function is a modified version of cobrapy add_pfba function			#
#										#
#################################################################################

def add_pfba_Weighted(model, weightings, objective=None, fraction_of_optimum=1.0):
    """Add pFBA objective
    Add objective to minimize the summed flux of all reactions to the
    current objective.
    See Also
    -------
    pfba
    Parameters
    ----------
    model : cobra.Model
        The model to add the objective to
    objective :
        An objective to set in combination with the pFBA objective.
    fraction_of_optimum : float
        Fraction of optimum which must be maintained. The original objective
        reaction is constrained to be greater than maximal_value *
        fraction_of_optimum.
    """
    if objective is not None:
        model.objective = objective
    if model.solver.objective.name == '_pfba_objective':
        raise ValueError('The model already has a pFBA objective.')
    sutil.fix_objective_as_constraint(model, fraction=fraction_of_optimum)
    reaction_variables = ((rxn.forward_variable, rxn.reverse_variable)
                          for rxn in model.reactions)
    variables = chain(*reaction_variables)
    model.objective = model.problem.Objective(
        Zero, direction='min', sloppy=True, name="_pfba_objective")
    #print([v for v in variables])
    tempDict = dict()
    for v in variables:
        w = str(v).split("=")[1].replace(" ","").replace("<","")
        found=False
        for rxn in weightings.keys():
            if w.__contains__(rxn):
                #print(v)
                #print(rxn)
                tempDict[v]=weightings[rxn]
                found=True
                break
        if not found:
            #print("Weightings for reaction "+w+" not found, so assuming weighting = 1")
            tempDict[v] = 1
    model.objective.set_linear_coefficients(tempDict)


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
cobra_model = io.sbml.read_sbml_model("./model_phloemUpdated.sbml")


#Remove all metabolites except sucrose from Phloem
rxn = cobra_model.reactions.get_by_id("Phloem_output_tx")
mets2remove = list()


fin= open("./../../ePhotosynthesis/InputATPCost.txt","r")
for line in fin:
    JATPase = float(line.replace("ATPCost ","").replace("ATPCost\t",""))
    break
fin.close()


fin= open("./../../ePhotosynthesis/InputNADPHCost.txt","r")
for line in fin:
    JNADPHox = float(line.replace("NADPHCost ","").replace("NADPHCost\t",""))
    break
fin.close()
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
rxn = Reaction("GAP_tx",name = "TP source")
rxn.add_metabolites({cobra_model.metabolites.get_by_id("GAP_c"):1})
rxn.upper_bound = 1000
rxn.lower_bound = 0
cobra_model.add_reaction(rxn)


rxn = Reaction("G3P_tx",name = "PGA source")
rxn.add_metabolites({cobra_model.metabolites.get_by_id("G3P_c"):1})
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


temp = cobra_model.copy()

PPFD = df["Light intensity"][0]
#constrain maintenace
ATPase = (0.0049*PPFD) + 2.7851
ATPase = round(ATPase,2)
temp.reactions.get_by_id("ATPase_tx").lower_bound = ATPase
temp.reactions.get_by_id("ATPase_tx").upper_bound = ATPase
temp.reactions.get_by_id("NADPHoxc_tx").lower_bound = ATPase/9
temp.reactions.get_by_id("NADPHoxc_tx").upper_bound = ATPase/9
temp.reactions.get_by_id("NADPHoxp_tx").lower_bound = ATPase/9
temp.reactions.get_by_id("NADPHoxp_tx").upper_bound = ATPase/9
temp.reactions.get_by_id("NADPHoxm_tx").lower_bound = ATPase/9
temp.reactions.get_by_id("NADPHoxm_tx").upper_bound = ATPase/9

#constraint TP flux
temp.reactions.get_by_id("GAP_tx").lower_bound = df["VT3P"][0]
temp.reactions.get_by_id("GAP_tx").upper_bound = df["VT3P"][0]
temp.reactions.get_by_id("G3P_tx").lower_bound = df["VPGA"][0]
temp.reactions.get_by_id("G3P_tx").upper_bound = df["VPGA"][0]

#constraint glycollate and glycerate fluxes flux
temp.reactions.get_by_id("GLYCOLATE_tx").lower_bound = df["Vt_glycolate"][0]
temp.reactions.get_by_id("GLYCOLATE_tx").upper_bound = df["Vt_glycolate"][0]
temp.reactions.get_by_id("GLYCERATE_tx").lower_bound = df["Vt_glycerate"][0]
temp.reactions.get_by_id("GLYCERATE_tx").upper_bound = df["Vt_glycerate"][0]

#temp.reactions.get_by_id("NrefixationCostbypass").lower_bound = df270["Vt_glycolate"][i]
#temp.reactions.get_by_id("NrefixationCostbypass").upper_bound = df270["Vt_glycolate"][i]
#temp.reactions.get_by_id("NrefixationEnergy").lower_bound = df270["Vt_glycerate"][i]
#temp.reactions.get_by_id("NrefixationEnergy").upper_bound = df270["Vt_glycerate"][i]

temp.reactions.get_by_id("MAL_v_accumulation").lower_bound = 0.0698903487288*df["Vstarch"][0]
temp.reactions.get_by_id("MAL_v_accumulation").upper_bound = 0.0698903487288*df["Vstarch"][0]

temp.reactions.get_by_id("CIT_v_accumulation").lower_bound = -0.056884259879*df["Vstarch"][0]
temp.reactions.get_by_id("CIT_v_accumulation").upper_bound = -0.056884259879*df["Vstarch"][0]



for rxn in cobra_model.reactions:
    if rxn.lower_bound == -1000:
        rxn.lower_boudn = -3000
    if rxn.upper_bound == 1000:
        rxn.upper_bound = 3000


#ADD ATP source reaction in FBA to represent ATP from JATPase
rxn = Reaction("ATP_source_from_ODE")
rxn.add_metabolites({temp.metabolites.get_by_id("ADP_p"):-0.8,
                     temp.metabolites.get_by_id("aADP_p"):-0.2,
                     temp.metabolites.get_by_id("Pi_p"):-1,
                     temp.metabolites.get_by_id("PROTON_p"):-0.9,
                     temp.metabolites.get_by_id("ATP_p"):0.9,
                     temp.metabolites.get_by_id("aATP_p"):0.1,
                     temp.metabolites.get_by_id("WATER_p"):1})
rxn.lower_bound = JATPase
rxn.upper_bound = JATPase
temp.add_reaction(rxn)

#ADD NADPH source reaction in FBA to represent NADPH from ODE
rxn = Reaction("NADPH_source_from_ODE")
rxn.add_metabolites({temp.metabolites.get_by_id("NADP_p"):-1,
                     temp.metabolites.get_by_id("WATER_p"):-1,
                     temp.metabolites.get_by_id("NADPH_p"):1,
                     temp.metabolites.get_by_id("OXYGEN_MOLECULE_p"):1,
                     temp.metabolites.get_by_id("PROTON_p"):1})
rxn.lower_bound = JNADPHox
rxn.upper_bound = JNADPHox
temp.add_reaction(rxn)

#check if model works
temp.solver="glpk"
#sol = custom_pFBA(temp)
weightings = dict()
for rxn in temp.reactions:
	if rxn.id == "ATP_ADP_Pi_pc":
		weightings[rxn.id]=0.5
	else:
		weightings[rxn.id]=1
sol = pfba_Weighted(temp,weightings=weightings)

rxn =  temp.reactions.get_by_id("Phloem_output_tx")
met = temp.metabolites.sSUCROSE_b
print("Sucrose export rate ="+str(rxn.metabolites[met]*sol.fluxes[rxn.id]))

total = JATPase
for rxn in temp.metabolites.ATP_p.reactions:
    if round(rxn.flux,3) != 0:
        coeff1 = rxn.metabolites[temp.metabolites.ATP_p]
        coeff2 = rxn.metabolites[temp.metabolites.aATP_p]
        ATPflux = sol.fluxes[rxn.id]*(coeff1+coeff2)
        print(rxn.id+"\t"+str(ATPflux)+"="+str(total))
        if rxn.id == "ATP_ADP_Pi_pc":
            total = total + ATPflux
            print(ATPflux)
print("Extra APTase flux ="+str(total))

fout= open("./../../ePhotosynthesis/InputATPCost.txt","w")
fout.write("ATPCost	"+str(round(total,4)))
fout.close()

total = JNADPHox
for rxn in temp.metabolites.NADPH_p.reactions:
    if round(rxn.flux,3) != 0:
        coeff1 = rxn.metabolites[temp.metabolites.NADPH_p]
        NADPHflux = sol.fluxes[rxn.id]*(coeff1)
        print(rxn.id+"\t"+str(NADPHflux)+"="+str(total))
        if rxn.id == "MALATE_DEH_RXN_p":
            total = total + NADPHflux
            print(NADPHflux)
for rxn in temp.metabolites.NADH_p.reactions:
    if round(rxn.flux,3) != 0:
        coeff1 = rxn.metabolites[temp.metabolites.NADH_p]
        NADPHflux = sol.fluxes[rxn.id]*(coeff1)
        print(rxn.id+"\t"+str(NADPHflux)+"="+str(total))
        if rxn.id == "MALATE_DEH_RXN_p":
            total = total + NADPHflux
            print(NADPHflux)
print("Extra NADPH flux ="+str(total))

fout= open("./../../ePhotosynthesis/InputNADPHCost.txt","w")
fout.write("NADPHCost	"+str(round(total,4)))
fout.close()

fout= open("./Daytime_flux.csv","w")
for rxn in temp.reactions:
    fout.write(rxn.id+","+rxn.reaction+","+str(sol.fluxes[rxn.id])+"\n")
fout.close()
