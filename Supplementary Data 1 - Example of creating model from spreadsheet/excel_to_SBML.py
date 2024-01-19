#!/usr/bin/env python
# coding: utf-8


# import package
import pandas as pd
import numpy as np
import sys


from libsbml import *


def multistep_check(reactions, initial_state): # function for checking if any of the reactions have been flagged as having multistep
    drop_these = [] # list for catching duplicate reactions later on

    if sum(reactions["Multistep"] == True) > 0: # check if there are any multistep reactions

            multistep_rxns = reactions.loc[reactions["Multistep"] == True] # collect multistep reactions

            reactions = reactions.loc[reactions["Multistep"] == False] # collect reactions that are not multistep

            reactions = reactions.drop("Multistep", axis=1) # drop the axis with information on multistep reactions

            potential_dups = len(reactions) # count how many reactions are left that may need to be duplicated

            for i in multistep_rxns.index: # iterate through the multistep reactions

                reactants = multistep_rxns.loc[i]["Reactants"].replace(" ", "").split("+") # collect the reactants and remove white space. Separate based on the + symbol

                reaction_name = multistep_rxns.loc[i]["Reaction Name"] # get the name of the reactions

                k_enc = multistep_rxns.loc[i]["k"] # get the rate constant of the reactions

                multiplicity_reactant = [] # make placeholder lists for determining what will be binding in steps

                single_reactant = [] # make placeholder lists for determining what will be the substrate binding the molecule

                duplicate_rxns = [] # placeholder for generating the reactions

                for react in reactants: # iterate through the reactants

                    if "*" in react: # if there is a * sign in the reactant indicating that there are multiple participating in the reaction, it is one of the binding proteins

                        multiplicity_reactant.append(react.split("*")[1]) # split and collect the reactant name

                        num_react = int(react.split("*")[0]) # split and collect the number of the reactants
                    else:
                        single_reactant.append(react) # if it does not have a * in it then it is the substrate
                        
                for i in range(potential_dups): # check every reaction for if the substrate is participating in the reaction record it as a reaction that needs to be duplicated if it does

                    if (single_reactant[0] in reactions["Reactants"].iloc[i]) & (single_reactant[0] + "_" not in reactions["Reactants"].iloc[i]): 
                        duplicate_rxns.append(i)

                fx_react = multiplicity_reactant[0] # set the reactant that will be duplicated

                enc_step_reactant = single_reactant[0] + f" + {fx_react}" # set the binding first reaction reactants

                enc_step_product = single_reactant[0] + f"_{fx_react}1" # make the first binding reaction products


                new_reactant = pd.DataFrame([[enc_step_product, 0]] ,columns=["Molecule", "Amount"]) # make a data frame recording the new species added to the system

                initial_state = initial_state.append(new_reactant) # add it to the list of reactants

                new_reaction = pd.DataFrame([reaction_name,enc_step_reactant,enc_step_product,k_enc]).T#,list(reactions.columns)) # make the new reaction in the format of the reaction data frame

                new_reaction.columns = reactions.columns # rename the columns to line up with the reaction data frame

                reactions = reactions.append(new_reaction) # add it to the reaction data frame
                
                if len(duplicate_rxns) > 0: # if there are any reactions that need to be duplicated cycle through them and make a new version of them where the reactant and/or product is now the new reactant we made

                    for index in duplicate_rxns: 

                        duplicate_rxn_reactant = reactions["Reactants"].iloc[index].replace(single_reactant[0], enc_step_product)

                        duplicate_name = reactions["Reaction Name"].iloc[index] + f"_{index}0{i}" 
                        

                        if pd.isna(reactions["Immediate_Products"].iloc[index]): 

                            duplicate_rxn_product = reactions["Immediate_Products"].iloc[index]

                        else:

                            duplicate_rxn_product = reactions["Immediate_Products"].iloc[index].replace(single_reactant[0], enc_step_product)

                        duplicate_rxn_k = reactions["k"].iloc[index]

                        new_reaction = pd.DataFrame([duplicate_name,duplicate_rxn_reactant,duplicate_rxn_product,duplicate_rxn_k]).T

                        new_reaction.columns = reactions.columns

                        reactions = reactions.append(new_reaction)

                for j in range(1, num_react): # now do the same for every additional binding step
                    enc_step_reactant = single_reactant[0] + f"_{fx_react}{j}"+ f" + {fx_react}"
                    enc_step_product = single_reactant[0] + f"_{fx_react}{j+1}"
                    new_name = reaction_name + f"_{j}"
                    if j < num_react - 1:
                        new_reactant = pd.DataFrame([[enc_step_product, 0]] ,columns=["Molecule", "Amount"])
                        initial_state = initial_state.append(new_reactant)                        
                    new_reaction = pd.DataFrame([new_name,enc_step_reactant,enc_step_product,k_enc]).T#,list(reactions.columns))
                    new_reaction.columns = reactions.columns
                    reactions = reactions.append(new_reaction)

                    if len(duplicate_rxns) > 0:
                        for index in duplicate_rxns:
                            
                            duplicate_rxn_reactant = reactions["Reactants"].iloc[index].replace(single_reactant[0], enc_step_product)
                            duplicate_name = reactions["Reaction Name"].iloc[index] + f"_{index}{j}{i}"
                            
                            if pd.isna(reactions["Immediate_Products"].iloc[index]):
                                duplicate_rxn_product = reactions["Immediate_Products"].iloc[index]
                            else:
                                duplicate_rxn_product = reactions["Immediate_Products"].iloc[index].replace(single_reactant[0], enc_step_product)
                            duplicate_rxn_k = reactions["k"].iloc[index]
                            new_reaction = pd.DataFrame([duplicate_name,duplicate_rxn_reactant,duplicate_rxn_product,duplicate_rxn_k]).T
                            new_reaction.columns = reactions.columns
                            reactions = reactions.append(new_reaction)
                    drop_these = drop_these + duplicate_rxns 

    reactions.index = range(len(reactions))

    reactions = reactions.drop(np.unique(drop_these))

    return reactions, initial_state # return the initial species state and the reaction data frames



def check(value, message): # this function is taken from the libSBML tutorial
  """If 'value' is None, prints an error message constructed using
  'message' and then exits with status code 1.  If 'value' is an integer,
  it assumes it is a libSBML return status code.  If the code value is
  LIBSBML_OPERATION_SUCCESS, returns without further action; if it is not,
  prints an error message constructed using 'message' along with text from
  libSBML explaining the meaning of the code, and exits with status code 1.
  """
  if value == None:
    raise SystemExit('LibSBML returned a null value trying to ' + message + '.')
  elif type(value) is int:
    if value == LIBSBML_OPERATION_SUCCESS:
      return
    else:
      err_msg = 'Error encountered trying to ' + message + '.'                 + 'LibSBML returned error code ' + str(value) + ': "'                 + OperationReturnValue_toString(value).strip() + '"'
      raise SystemExit(err_msg)
  else:
    return



def create_model(reactions,initial_state): # now use the reaction and species state data frames to create an SBML model. The initialization steps for the model are taken from the libSBML tutorial

    try:
        document = SBMLDocument(3, 1) # make the documentation object
    except ValueError:
        raise SystemExit('Could not create SBMLDocumention object')

    model = document.createModel()
	
    check(model,                              'create model') # set units
    check(model.setTimeUnits("second"),       'set model-wide time units')
    check(model.setExtentUnits("item"),       'set model units of extent')
    check(model.setSubstanceUnits('item'),    'set model substance units')

    # Create the basic Model object inside the SBMLDocument object.  To
    # produce a model with complete units for the reaction rates, we need
    # to set the 'timeUnits' and 'extentUnits' attributes on Model.  We
    # set 'substanceUnits' too, for good measure, though it's not strictly
    # necessary here because we also set the units for invididual species
    # in their definitions.

    per_second = model.createUnitDefinition()
    check(per_second,                         'create unit definition')
    check(per_second.setId('per_second'),     'set unit definition id')
    unit = per_second.createUnit()
    check(unit,                               'create unit on per_second')
    check(unit.setKind(UNIT_KIND_SECOND),     'set unit kind')
    check(unit.setExponent(-1),               'set unit exponent')
    check(unit.setScale(0),                   'set unit scale')
    check(unit.setMultiplier(1),              'set unit multiplier')

    # Create a compartment inside this model, and set the required
    # attributes for an SBML compartment in SBML Level 3.

    c1 = model.createCompartment()
    check(c1,                                 'create compartment')
    check(c1.setId('c1'),                     'set compartment id')
    check(c1.setConstant(True),               'set compartment "constant"')
    check(c1.setSize(1),                      'set compartment "size"')
    check(c1.setSpatialDimensions(3),         'set compartment dimensions')
    check(c1.setUnits('litre'),               'set compartment size units')
    # Create a unit definition we will need later.  Note that SBML Unit
    # objects must have all four attributes 'kind', 'exponent', 'scale'
    # and 'multiplier' defined.
    species = []
    for i in range(len(initial_state)): # iterate through each species in our species data frame and add them to the model
        create_species = model.createSpecies()
        species_name = initial_state.iloc[i][0]
        starting_quantity = int(initial_state.iloc[i][1])
        check(create_species,                                 f'create species {species_name}')
        check(create_species.setId(species_name),                     f'set species {species_name} id')
        check(create_species.setCompartment('c1'),            f'set species {species_name} compartment')
        check(create_species.setConstant(False),              f'set "constant" attribute on {species_name}')
        check(create_species.setInitialAmount(starting_quantity),             f'set initial amount for {species_name}')
        check(create_species.setSubstanceUnits('mole'),       f'set substance units for {species_name}')
        check(create_species.setBoundaryCondition(False),     f'set "boundaryCondition" on {species_name}')
        check(create_species.setHasOnlySubstanceUnits(False), f'set "hasOnlySubstanceUnits" on {species_name}')

    for i in range(len(reactions)): # iterate through each reaction in our species data frame and add them to the model
        k = model.createParameter()
        check(k,                                  'create parameter k')
        check(k.setId(f'k{i}'),                       'set parameter k id')
        check(k.setConstant(True),                'set parameter k "constant"')
        check(k.setValue(float(reactions["k"].iloc[i])),                      'set parameter k value')
        check(k.setUnits('per_second'),           'set parameter k units')
        

         
            
            
        r1 = model.createReaction()
        check(r1,                                 'create reaction')
        check(r1.setId(reactions["Reaction Name"].iloc[i]),                     'set reaction id')
        check(r1.setReversible(False),            'set reaction reversibility flag')
        check(r1.setFast(False),                  'set reaction "fast" attribute')
        reactants = reactions["Reactants"].iloc[i].replace(" ", "").split("+")
        
        for react in reactants:
            species_ref1 = r1.createReactant()
            check(species_ref1,                       'create reactant')
            check(species_ref1.setSpecies(react),      'assign reactant species')
            check(species_ref1.setConstant(False),     'set "constant" on species ref 1')
        if isinstance(reactions["Immediate_Products"].iloc[i], str):
            products = reactions["Immediate_Products"].iloc[i].replace(" ", "").split("+")
            for prod in products:
                species_ref2 = r1.createProduct()
                check(species_ref2,                       'create product')
                check(species_ref2.setSpecies(prod),      'assign product species')
                check(species_ref2.setConstant(False),     'set "constant" on species ref 2')

        reactants = ' * '.join(reactants)
        math_ast = parseL3Formula(f'k{i} * {reactants} * c1')
        check(math_ast,                           'create AST for rate expression')
        kinetic_law = r1.createKineticLaw()
        check(kinetic_law,                        'create kinetic law')
        check(kinetic_law.setMath(math_ast),      'set math on kinetic law')
    return writeSBMLToString(document) # return the SBML model
    

def main(): # function to process and save the file
    reader = SBMLReader() # load SBML reader
    system_call = sys.argv # collect system calls
    
    reactions = pd.read_excel(system_call[1]) # get reaction data frame
    initial_state = pd.read_excel(system_call[2]) # get species data frame
    reactions, initial_state = multistep_check(reactions, initial_state) # process if they have multistep reactions
    mod = create_model(reactions,initial_state) # make the model
    mod = readSBMLFromString(mod) # get the SBML model
    writeSBMLToFile(mod, system_call[3]) # save the SBML Model.

if __name__ == "__main__":
    main()



