Citation: A healthy but depleted herd: predators decrease prey disease and density

Authors: Laura K. Lopez 1, Michael H. Cortez 2*, Turner DeBlieux 3, Ilona A. Menel 4, Bruce O’Brien 1, Carla E. Cáceres 4, Spencer R. Hall 3, and Meghan A. Duffy 1
        1 Department of Ecology & Evolutionary Biology, University of Michigan, Ann Arbor, MI 48109, USA
        2 Department of Biological Science, Florida State University, Tallahassee, FL 32306, USA
        3 Department of Biology, Indiana University, Bloomington, IN 47405 USA
        4 School of Integrative Biology, University of Illinois Urbana-Champaign, Urbana, IL 61801 

Contact: * Author responsible for writing code: mcortez@fsu.edu
         
Date: June 2022

-----------------------------------------------------------------------------------

DATA SHEETS

Data are saved in three folders, corresponding to different experiments 

Folder: MesocosmData
        - genotype.frequencies_MHC_converted.csv
                Clone frequencies in each replicate mesocosm. Columns are
                Parasites: 1 is with parasite, 0 is without parasite
                Predation: density of predators in individuals/L
                Tank: Replicate number
                Week: Weeks were -1 to 8 converted to 0 to 9 in main text
                Day: Julian day
                no.hosts.sampled: Number of individuals sampled 
                Remaining columns: frequencies of the 9 clones
        - meso.expt.densities_MHC_converted.csv
                Prey densities in each replicate mesocosm.  Columns are
                Parasites: 1 is with parasite, 0 is without parasite
                Predation: density of predators in individuals/L
                Tank: Replicate number
                Week: Weeks were -1 to 8 converted to 0 to 9 in main text
                JulianDate: self explanatory
                host.density: Density of prey in individuals/L
                infected.host.density: density of infected prey in individuals/L
                Prevalence: Infected density/total density

Folder: PredationData
        - YY##-##.xlsx 
                Consumption data for clone YY##-##.  Columns are
                Density: Number of prey offered to predator
                number.eaten: Number of prey eaten

Folder: PredatorPreferenceData
        - predationtrialresultsclean_mhc.csv
                Consumption data for clone Mid-37.  Columns are
                treatment: control means no exposure, exposed means exposed to spores
                        but not visibly infected, infected means visibly infected
                block: relates to experimental set up
                day: days since exposure 
                data: self explanatory
                replicate: self explanatory
                daphnia.alive: number of Daphnia alive at end of predation trial
                daphnia.dead: number of Daphnia dead at end of predation trial
                day.block: Block1 is days 1, 4, 7 and Block2 is days 10, 13
                percent.dead: daphnia.dead/10 (all trials started with 10 Daphnia)

-----------------------------------------------------------------------------------

CODE

Code for parameter estimation
        - script_EstimateParameters.m 
                Matlab script that estimates predator attack rates on each clone 
        - Rscript_EstimatePredatorPreference.R
                R script that estimates predator attack rate on susceptible, exposed,
                and infected Mid-37 clones
        - script_EstimateParameters.m
                Matlab script that estimates per spore probability of infection and
                reduction in spore burst size for consumed infected individuals for each clone
        - script_EstimateParameters_NoSelectivity.m
                Same, but assumes predators have equal attack rates on susceptible and
                infected prey (omega =1)

Code for figures
        - R0_GeberMethod.m
                Matlab script that computes R0 & R values and applies Geber Method
                Generates Figures 2, S3
        - R0_GeberMethod_NoSelectivity.m
                Same, but assumes predators have equal attack rates on susceptible and
                infected prey (omega =1)
                Generates Figures S4, S5
        - script_HHHpredictions.m
                Matlab script that computes negative of numerator of equation (S26)
                Generates Figure S6
        - script_HHHpredictions_NoSelectivity.m
                Same, but assumes predators have equal attack rates on susceptible and
                infected prey (omega =1)
                Generates Figures S7

             
