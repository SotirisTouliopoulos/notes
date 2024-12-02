
from cobra.io import load_json_model


bigg_cofactors = ['atp[c]', 'atp_c', 'adp_c', 'adp[c]',
                  'atp[c]', 'atp_c', 'adp_c', 'adp[c]',
                  'udp[c]', 'udp_c', 'ump[c]', 'ump_c',
                  'amp_c', 'amp[c]',
                  'gdp[c]', 'gdp_c', 'gtp[c]', 'gtp_c',
                  'accoa_c', 'accoa[c]', 'coa_c', 'coa[c]',  # acetyl-CoA
                  'q8[c]', 'q8_c', 'q8h2_c', 'q8h2[c]', 'mqn8_c', 'mqn8[c]', 'mql8_c', 'mql8[c]', 'q8h2_c', 'q8h2[c]',
                  'actp[c]', 'actp_c',
                  'h2o_c', 'h2o[c]', 'h2o_e', 'h2o[e]',
                  'pi_e', 'pi[e]', 'pi_c', 'pi[c]', 'ppi[c]', 'ppi_c',
                  'pep_c', 'pep[c]',
                  'h_c', 'h[c]', 'h_e', 'h[e]',
                  'o2_c', 'o2[c]', 'o2_e', 'o2[e]',
                  'co2_c', 'co2[c]', 'co2_e', 'co2[e]',
                  'nadp_c', 'nadp[c]', 'nadph_c', 'nadph[c]', 'nad_c', 'nad[c]', 'nadh_c', 'nadh[c]',
                  'nadp_e', 'nadp[e]', 'nadph_e', 'nadph[c]', 'nad_e', 'nad[e]', 'nadh_e', 'nadh[e]',
                  'fadh2_c', 'fadh2[c]', 'fad_c', 'fad[c]',
                  'nh4_c', 'nh4[c]', 'nh4_e', 'nh4[e]',
                  'pyr[c]', 'pyr_c'
                ]

bigg_building_blocks = ['ala_L[c]', 'asp_L[c]', ' gln_L[c]', 'glu_L[c]', 'glu_L[c]', 
                        'ser_L[c]', 'trp_L[c]', 'met_L[c]', 'lys_L[c]', 'cyst_L[c]'
                       ]


# this function may be used later if directionality is needed
"""
def arrow_index():
    # load cobra model
    cobra_model = load_json_model("../data/e_coli_core.json")
    # list with reactions names
    reactions_ids =  [ reaction.id for reaction in cobra_model.reactions ]
    
    # define arrow types
    arrow_types = ["-->", "<=>", "<--"]
    # define list to store indices
    arrow_indices = []
    
    # parse every reaction
    for reaction in reactions_ids:
        reaction_stoichiometry = str(cobra_model.reactions.get_by_id(reaction))
        print(reaction_stoichiometry.split())
        
        for arrow_type in arrow_types:
            # find index of arrow
            try:
                arrow_index = reaction_stoichiometry.split().index(arrow_type)
            except:
                arrow_index = None
            finally:
                if arrow_index != None:
                    arrow_indices.append(arrow_index)
"""         
    

def reactant_product_cofactor():
    # load cobra model
    cobra_model = load_json_model("../data/e_coli_core.json")
    
    # list with reactions names
    reactions_ids =  [ reaction.id for reaction in cobra_model.reactions ]
    
    # define variable found in the reaction stoichiometry str
    stoichiometry_elements = ["-->", "<=>", "<--", "+"]
    
    # define lists to store reactants, products and cofactors of all reactions
    reactant_product_all_reactions = []
    cofactor_all_reactions = []
    
    for reaction in reactions_ids:
        # define lists to store reactants, products and cofactors of one reaction
        reactant_product_one_reaction = []
        cofactor_one_reaction = []
        
        # convert stoichiometry to str
        reaction_stoichiometry = str(cobra_model.reactions.get_by_id(reaction))

        # iterate through every element of the stoichiometry str
        for i in range(1, len(reaction_stoichiometry.split()) , 1):
            element = reaction_stoichiometry.split()[i]
            
            # try convert stoichiometries from str to float to find coefficients
            try:
                element = float(element)
            except:
                pass
            
            # find reactants or products
            if  element not in bigg_cofactors and element not in bigg_building_blocks and \
                element not in stoichiometry_elements and type(element) != float :
                reactant_product_one_reaction.append(element)
            
            # find cofactors
            elif element in bigg_cofactors or element in bigg_building_blocks:
                cofactor_one_reaction.append(element)
                     
        reactant_product_all_reactions.append(reactant_product_one_reaction)
        cofactor_all_reactions.append(cofactor_one_reaction)
        
    return reactions_ids, reactant_product_all_reactions, cofactor_all_reactions
    

def cofactor_specificity():
    reactions_ids, reactant_product_all_reactions, cofactor_all_reactions = reactant_product_cofactor()
    
    # define list to store possible combinations of reactions
    cofactor_specificity = []

    for i in range(len(reactions_ids)):
        for j in range(len(reactions_ids)):
            
            # avoid comparison of the same reaction
            if reactions_ids[i] != reactions_ids[j]:
                # boolean to check if pairwise elements are the same
                identical_elements = (set(reactant_product_all_reactions[i]) == set(reactant_product_all_reactions[j]))
                
                if identical_elements == True:
                    # boolean to check if pairwise cofactors are the same
                    identical_cofactors = (set(cofactor_all_reactions[i]) != set(cofactor_all_reactions[j]))
                    
                    if identical_cofactors == False:
                        cofactor_specificity.append((reactions_ids[i], reactions_ids[j]))
                        
                        #print(reactions_ids[i], reactions_ids[j])
                        #print(cofactor_all_reactions[i], cofactor_all_reactions[j])

    unique_set = {tuple(sorted(pair)) for pair in cofactor_specificity}
    print("Combinations of reactions:" , list(unique_set))


cofactor_specificity()