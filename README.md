# All code in this repository was used in the ODE-FBA leaf model study (paper details will be available here once the study is published)

## Index ##  
**Data** : This folder contains biomass compositions and metabolic models used in the study  
**Data/biomass_soy.csv** : Soy tissue specific biomass composition  
**Data/PlantCoreMetabolism_v2_0_0.xml** : The PlantCoreMetabolism v2.0.0 model  
**Data/Soy_core_model_GA.sbml** : PlantCoreMetabolism v2.0.0 with soy GRP associations  
**Data/GenerateSupplementalDataS2FromSoyModel.ipynb** : A jupyter notebook adding Allantoin degradation pathway into Soy core model  
**Data/SupplementalDataS2.xml** : SupplementalData S2  
**ePhotosynthesis** : This folder contains all code required to run the ePhotosynthesis model in coupled and uncoupled configuration including the yaml files to run it via yggdrasil
**ePhotosynthesis/InputEvn.txt** : A text file for specifying light intensity, [CO2], daylength and whether sucrose export needs to be considered in the ODE model (leave as 1 when attempting to replicate results in the study)  
**ePhotosynthesis/InputATPCost.txt** : A text file for specifying Vatp_ode  
**ePhotosynthesis/InputNADPHCost.txt** : A text file for specifying Vnadph_ode  
**ePhotosynthesis/Sim_Ephotosynthesis.m** : M file to run ephotosynthesis model  
**ePhotosynthesis/EphotosynthesisOnly.yml** : YAML file to run ephotosynthesis model via yggdrasil
**ePhotosynthesis/OutputFluxT.txt**: A text file containing results of an ODE model run  
**ePhotosynthesis/FluxAnnotation.xlsx**: An excel file containing annotations for the ODE outputs  
**FBA** : This folder contains all code required to run the FBA model in coupled and uncoupled configuration including the yaml files to run it via yggdrasil  
**FBA/dielFBA_model_only.py** : Python code to run the diel FBA model on its own (a.k.a FBA model in text)  
**FBA/yggrasil_ODE_FBA_dielFBA_only.yaml** : YAML file to run FBA model via yggdrasil
**FBA/dielFBA_model_forODE.py** : Python code to run diel FBA model in LC configuration  
**FBA/yggrasil_ODE_FBA_dielFBA.yml** : YAML file to run LC model via yggdrasil  
**FBA/FBA_model_for_ODE.py** : Python code to run FBA module representing day-time metabolism in TC configuration  
**FBA/yggrasil_ODE_FBA_testing.yaml** : YAML file to run day-time FBA module of TC via yggdrasil  
**FBA/FBA_model_for_ODE_night.py** : Python code to run FBA module representing night-time metabolism in TC configuration  
**FBA/yggrasil_ODE_FBA_night.yaml** : YAML file to run night-time FBA module of TC via yggdrasil  
**FBA/FBA_model_growing_forODE.py** : Python code to run FBA module representing growing leaf day-time metabolism in TC configuration  
**FBA/yggrasil_ODE_FBA_growing.yaml** : YAML file to run day-time growing leaf FBA module via yggdrasil  
**FBA/FBA_model_growing_forODE.py** : Python code to run FBA module representing growing leaf night-time metabolism in TC configuration  
**FBA/yggrasil_ODE_FBA_night_growing.yaml** : YAML file to run night-time growing leaf FBA module via yggdrasil  
**FluxMap** : This folder contains SVG/PNG files of flux maps used in the paper and code used to study fluxes  
**FVA** : Folder containing jupyter notebook demonstrating how we managed to ensure unique flux values through NTT during pFBA  
**sensitivityAnalyses** : Folder contains all code, results and final compilations associated with testing model sensitivity to FBA parameters  
**sensitivityAnalyses/NGAM** : Folder contains all code, results and final compilations associated with testing model sensitivity to non-growth associated maintenance cost  
**sensitivityAnalyses/phloemComposition** : Folder contains all code, results and final compilations associated with testing model sensitivity to phloem composition using data collected from Wilkinson and Douglas 2003  
**sensitivityAnalyses/starchOAAratios** : Folder contains all code, results and final compilations associated with testing model sensitivity to organic acid accumulation/remobilization constraint (equation 5 and 6)  
**Analysis** : Folder containing code, log, results and final compilation of results associated with studying metabolism in leaves under varying PPFD and [CO2] using TC  
**Validations** : Folder containing code, logs, results and final compilation of results associated with Figures 3,4,5 and 7, and Tables 1 and 2  
**Validations/ODE_Figure2A.py** : code used to run ODE model for Figure 3A  
**Validations/LC_Figure2A.py** : code used to run LC model for Figure 3A  
**Validations/TC_Figure2A.py** : code used to run TC model for Figure 3A  
**Validations/FBA_Figure2B.py** : code used to run TC model for Figure 3B
**Validations/ODE_Figure2B.py** : code used to run ODE model for Figure 3B  
**Validations/LC_Figure2B.py** : code used to run LC model for Figure 3B  
**Validations/TC_Figure2B.py** : code used to run TC model for Figure 3B  
**Validations/Script1.ipynb** : jupyter notebook code to generate Figure 3A  
**Validations/Script2.ipynb** : jupyter notebook code to generate Figure 3B  
**Validations/YggdrasilLooping_LooselyCoupledFACE.py** : Python code to model FACE conditions with LC model using yggdrasil  
**Validations/YggdrasilLooping.py** : Python code to model mature leaf metabolism using TC model via yggdrasil  
**Validations/YggdrasilLooping_growing.py** : Python code to model growing leaf metabolism using TC model via yggdrasil  
**Validations/YggdrasilLooping_FACE%_372.py** : Python code to model leaf metabolism on day % of FACE experiment under ambient CO2 conditions using TC model  
**Validations/YggdrasilLooping_FACE%_552.py** : Python code to model leaf metabolism on day % of FACE experiment under elevated CO2 conditions using TC model  
