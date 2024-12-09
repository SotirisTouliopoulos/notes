
from cobra.io import load_json_model
from collections import defaultdict


bigg_cofactors = ['atp[c]', 'atp_c', 'adp_c', 'adp[c]',
                  'atp[c]', 'atp_c', 'adp_c', 'adp[c]',
                  'udp[c]', 'udp_c', 'ump[c]', 'ump_c',
                  'amp_c', 'amp[c]',
                  'gdp[c]', 'gdp_c', 'gtp[c]', 'gtp_c',
                  'accoa_c', 'accoa[c]', 'coa_c', 'coa[c]',
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


def find_arrow_index():
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
        
        for arrow_type in arrow_types:
            # find index of arrow
            try:
                arrow_index = reaction_stoichiometry.split().index(arrow_type)
            except:
                arrow_index = None
            finally:
                if arrow_index != None:
                    arrow_indices.append(arrow_index)
                    
    return arrow_indices

                    
def find_reactants_products_cofactors():
    
    # get the arrow indices from the corresponding function
    arrow_indices = find_arrow_index()
    
    # load cobra model
    cobra_model = load_json_model("../data/e_coli_core.json")
    
    # list with reactions names
    reactions_ids =  [ reaction.id for reaction in cobra_model.reactions ]
    
    # define variable found in the reaction stoichiometry str
    stoichiometry_elements = ["-->", "<=>", "<--", "+"]
    
    # define lists to store reactants, products and cofactors of all reactions
    reactants_all_reactions = []
    products_all_reactions = []
    cofactors_all_reactions = []
    
    for i in range(len(reactions_ids)):
        
        reaction = reactions_ids[i]
        arrow_index = arrow_indices[i]
        
        # define lists to store reactants, products and cofactors of one reaction
        reactants_single_reaction = []
        products_single_reaction = []
        cofactors_single_reaction = []
        
        # convert stoichiometry to str
        reaction_stoichiometry = str(cobra_model.reactions.get_by_id(reaction))

        # iterate through every element of the stoichiometry str, avoiding the first which is the reaction id
        for j in range(1, len(reaction_stoichiometry.split()) , 1):
            element = reaction_stoichiometry.split()[j]
            
            # try convert stoichiometries from str to float to find coefficients
            try:
                element = float(element)
            except:
                pass
            
            # find reactants
            if  element not in bigg_cofactors and element not in bigg_building_blocks and \
                element not in stoichiometry_elements and type(element) != float and j <= arrow_index:
                reactants_single_reaction.append(element)
                
            # find products
            elif  element not in bigg_cofactors and element not in bigg_building_blocks and \
                element not in stoichiometry_elements and type(element) != float and j >= arrow_index:
                products_single_reaction.append(element)

            # find cofactors
            elif element in bigg_cofactors or element in bigg_building_blocks:
                cofactors_single_reaction.append(element)
                     
        reactants_all_reactions.append(reactants_single_reaction)
        products_all_reactions.append(products_single_reaction)
        cofactors_all_reactions.append(cofactors_single_reaction)
        
    return reactions_ids, reactants_all_reactions, products_all_reactions, cofactors_all_reactions


def find_reactions_combinations():
    
    reactions_ids, reactants_all_reactions, products_all_reactions, cofactors_all_reactions = find_reactants_products_cofactors()
    
    # store possible combinations of reactions with same reactants/products and different cofactors
    cofactor_specificity = []

    for i in range(len(reactions_ids)):
        for j in range(len(reactions_ids)):
            
            # avoid comparison of the same reaction
            if reactions_ids[i] != reactions_ids[j]:

                # avoid comparison with empty lists
                if len(reactants_all_reactions[i]) > 0 and len(reactants_all_reactions[j]) > 0 and \
                   len(products_all_reactions[i]) > 0 and len(products_all_reactions[j]) > 0:
                
                    # boolean to check if pairwise reactants/products are the same
                    identical_reactants = (set(reactants_all_reactions[i]) == set(reactants_all_reactions[j]))
                    identical_products = (set(products_all_reactions[i]) == set(products_all_reactions[j]))
                    
                    if identical_reactants == True and identical_products == True:
                        # boolean to check if pairwise cofactors are the same
                        identical_cofactors = (set(cofactors_all_reactions[i]) == set(cofactors_all_reactions[j]))
                        
                        if identical_cofactors == False:
                            cofactor_specificity.append((reactions_ids[i], reactions_ids[j]))
                            
                            # for debugging only
                            #print("Reaction 1:", reactions_ids[i], "| Reaction 2:", reactions_ids[j])
                            #print("Reaction 1 Reactants:", reactants_all_reactions[i], "| Reaction 2 Reactants:", reactants_all_reactions[j])
                            #print("Reaction 1 Products:", products_all_reactions[i], "| Reaction 2 Products:", products_all_reactions[j])                        
                            #print("Reaction 1 Co-factors:", cofactors_all_reactions[i], "| Reaction 2 Co-factors:", cofactors_all_reactions[j])

    unique_combinations = {tuple(sorted(pair)) for pair in cofactor_specificity}
    #print("Unique combinations of reactions:" , list(unique_combinations))
    
    return unique_combinations


def merge_connected_combinations():

    unique_combinations = find_reactions_combinations()
    
    # Step 1: Create a graph using an adjacency list representation
    # We'll use a defaultdict where each key is a node, and its value is a set of neighboring nodes
    graph = defaultdict(set)

    # Build the adjacency list by iterating through each pair in the input
    for x, y in unique_combinations:
        # Add y as a neighbor of x
        graph[x].add(y)
        # Add x as a neighbor of y (because the graph is undirected)
        graph[y].add(x)

    # At this point, for an input of this format: [('A', 'B'), ('A', 'C'), ('B', 'C'), ('D', 'K')]
    # the `graph` looks like:
    # {'A': {'B', 'C'}, 'B': {'A', 'C'}, 'C': {'A', 'B'}, 'D': {'K'}, 'K': {'D'}}

    # Step 2: Define a function to find connected components using Depth-First Search (DFS)
    def find_connected_components(graph):
        # A set to keep track of visited nodes
        visited = set()
        # A list to store the connected components
        components = []

        # Helper function to perform a recursive DFS
        def dfs(node, component):
            # Mark the current node as visited
            visited.add(node)
            # Add the node to the current component
            component.append(node)
            # Visit all unvisited neighbors of the current node
            for neighbor in graph[node]:
                if neighbor not in visited:
                    dfs(neighbor, component)

        # Iterate through all nodes in the graph
        for node in graph:
            # If the node has not been visited, it's the start of a new connected component
            if node not in visited:
                # Create a new component (list) to hold connected nodes
                component = []
                # Perform DFS starting from this node
                dfs(node, component)
                # Add the completed component to the list of components
                components.append(component)

        # Return the list of connected components
        return components

    # Step 3: Find connected components in the graph
    groups = find_connected_components(graph)

    # At this point, `groups` will contain:
    # [['A', 'B', 'C'], ['D', 'K']]

    # Print the result
    print(groups)
    
    
merge_connected_combinations()