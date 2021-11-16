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
        df2 = pd.read_csv(BytesIO(msg))
        print(df2.to_string())


from cobra import io,flux_analysis
from cobra.core import Reaction, Metabolite

#import model. Update file name and location in the next line
cobra_model = io.sbml.read_sbml_model("./../Data/PlantCoreMetabolism_v2_0_0.xml")

#Remove all metabolites except sucrose from Phloem
rxn = cobra_model.reactions.get_by_id("Phloem_output_tx")
mets2remove = list()



fin= open("./../ePhotosynthesis/InputATPCost.txt","r")
for line in fin:
    JATPase = float(line.replace("ATPCost ","").replace("ATPCost\t",""))
    break
fin.close()


fin= open("./../ePhotosynthesis/InputNADPHCost.txt","r")
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
cobra_model.reactions.get_by_id("NH4_tx").lower_bound=0
cobra_model.reactions.get_by_id("NH4_tx").upper_bound=0
cobra_model.reactions.get_by_id("Nitrate_tx").lower_bound = 0
cobra_model.reactions.get_by_id("Nitrate_tx").upper_bound = 0
cobra_model = addAllantoinMetabolism(cobra_model)

#no external light energy
cobra_model.reactions.get_by_id("Photon_tx").lower_bound = 0
cobra_model.reactions.get_by_id("Photon_tx").upper_bound = 0

#no rubisco carboxylase or oxygenase
cobra_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p").lower_bound = 0
cobra_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p").upper_bound = 0
cobra_model.reactions.get_by_id("RXN_961_p").lower_bound = 0
cobra_model.reactions.get_by_id("RXN_961_p").upper_bound = 0

#no starch biosynthesis
#cobra_model.reactions.get_by_id("GLYCOGENSYN_RXN_p").lower_bound = 0
#cobra_model.reactions.get_by_id("GLYCOGENSYN_RXN_p").upper_bound = 0

#set export of sugars as objective
cobra_model.reactions.get_by_id("Phloem_output_tx").objective_coefficient=0

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

#constraint TP flux
temp.reactions.get_by_id("GAP_tx").lower_bound = df2["VT3P"][0]
temp.reactions.get_by_id("GAP_tx").upper_bound = df2["VT3P"][0]
temp.reactions.get_by_id("G3P_tx").lower_bound = df2["VPGA"][0]
temp.reactions.get_by_id("G3P_tx").upper_bound = df2["VPGA"][0]

#constraint glycollate and glycerate fluxes flux
temp.reactions.get_by_id("GLYCOLATE_tx").lower_bound = df2["Vt_glycolate"][0]
temp.reactions.get_by_id("GLYCOLATE_tx").upper_bound = df2["Vt_glycolate"][0]
temp.reactions.get_by_id("GLYCERATE_tx").lower_bound = df2["Vt_glycerate"][0]
temp.reactions.get_by_id("GLYCERATE_tx").upper_bound = df2["Vt_glycerate"][0]

#temp.reactions.get_by_id("NrefixationCostbypass").lower_bound = df270["Vt_glycolate"][i]
#temp.reactions.get_by_id("NrefixationCostbypass").upper_bound = df270["Vt_glycolate"][i]
temp.reactions.get_by_id("2PGADEHYDRAT_RXN_p").lower_bound = 0
temp.reactions.get_by_id("2PGADEHYDRAT_RXN_p").upper_bound = 0

temp.reactions.get_by_id("G6P_Pi_pc").lower_bound=0
temp.reactions.get_by_id("G6P_Pi_pc").upper_bound=0

temp.reactions.get_by_id("MAL_v_accumulation").lower_bound = 0.71*df["Vstarch"][0]
temp.reactions.get_by_id("MAL_v_accumulation").upper_bound = 0.71*df["Vstarch"][0]

temp.reactions.get_by_id("CIT_v_accumulation").lower_bound = -0.56*df["Vstarch"][0]
temp.reactions.get_by_id("CIT_v_accumulation").upper_bound = -0.56*df["Vstarch"][0]




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
rxn =  temp.reactions.get_by_id("Biomass_leaf_tx")
print("Biomass C accumulation rate ="+str(rxn.flux*595.7664799294965))

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

fout= open("./../ePhotosynthesis/InputATPCost.txt","w")
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

fout= open("./../ePhotosynthesis/InputNADPHCost.txt","w")
fout.write("NADPHCost	"+str(round(total,4)))
fout.close()

fout= open("./Daytime_flux.csv","w")
for rxn in temp.reactions:
    fout.write(rxn.id+","+rxn.reaction+","+str(sol.fluxes[rxn.id])+"\n")
fout.close()
