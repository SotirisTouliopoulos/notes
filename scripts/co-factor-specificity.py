
from cobra.io import load_json_model
from collections import defaultdict
from collections import Counter


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


class Cofactor_Specificity:
    
    def __init__(self, model):
        
        self.model = model


    def find_reactants_products_cofactors(self):
        
        self.reactions_ids =  [ reaction.id for reaction in self.model.reactions ]
        
        self.reactants_list_all_reactions = []
        self.products_list_all_reactions = []
        self.cofactors_list_all_reactions = []
        self.reversibility_list_all_reactions = []
        
        for reaction in self.reactions_ids:
            reactants_list_single_reaction = []
            products_list_single_reaction = []
            cofactors_list_single_reaction = []

            reaction_information = cobra_model.reactions.get_by_id(reaction)
            
            reactants = reaction_information.reactants
            products = reaction_information.products
            
            reversibility = reaction_information.reversibility
            self.reversibility_list_all_reactions.append(reversibility)
            
            for reactant in reactants:
                reactant = str(reactant)
                if reactant in bigg_cofactors or reactant in bigg_building_blocks:
                    cofactors_list_single_reaction.append(reactant)
                else:
                    reactants_list_single_reaction.append(reactant)
                    
            for product in products:
                product = str(product)
                if product in bigg_cofactors or product in bigg_building_blocks:
                    cofactors_list_single_reaction.append(product)
                else:
                    products_list_single_reaction.append(product)         
                    
            #print(reaction_information, reactants_list_single_reaction, products_list_single_reaction, cofactors_list_single_reaction, end ="\n")
            self.reactants_list_all_reactions.append(reactants_list_single_reaction)
            self.products_list_all_reactions.append(products_list_single_reaction)
            self.cofactors_list_all_reactions.append(cofactors_list_single_reaction)


    def find_reactions_combinations(self):
        
        self.find_reactants_products_cofactors()
        
        combinations = []
        
        for i in range(len(self.reactions_ids)):
            for j in range(len(self.reactions_ids)):
                
                # avoid comparison of the same reaction
                if self.reactions_ids[i] != self.reactions_ids[j]:

                    # avoid comparison with empty lists
                    if len(self.reactants_list_all_reactions[i]) > 0 and len(self.reactants_list_all_reactions[j]) > 0 and \
                    len(self.products_list_all_reactions[i]) > 0 and len(self.products_list_all_reactions[j]) > 0:
                    
                        # boolean to check if pairwise reactants/products are the same
                        identical_reactants = (set(self.reactants_list_all_reactions[i]) == set(self.reactants_list_all_reactions[j]))
                        identical_products = (set(self.products_list_all_reactions[i]) == set(self.products_list_all_reactions[j]))
                        identical_cofactors = (set(self.cofactors_list_all_reactions[i]) == set(self.cofactors_list_all_reactions[j]))

                        if identical_reactants == True and identical_products == True and identical_cofactors == False:
                            # finish here and keep combination
                            combinations.append((self.reactions_ids[i], self.reactions_ids[j]))
                        else:
                            if self.reversibility_list_all_reactions[i] == True or self.reversibility_list_all_reactions[j] == True:
                                # switch reaction i and compare updates reactants-products
                                identical_reactants = (set(self.reactants_list_all_reactions[i]) == set(self.products_list_all_reactions[j]))
                                identical_products = (set(self.products_list_all_reactions[i]) == set(self.reactants_list_all_reactions[j]))

                                if identical_reactants == True and identical_products == True and identical_cofactors == False:
                                    # finish here and keep combination
                                    combinations.append((self.reactions_ids[i], self.reactions_ids[j]))

        self.unique_combinations = {tuple(sorted(pair)) for pair in combinations}


    def merge_connected_combinations(self):
        
        self.find_reactions_combinations()
        
        # Step 1: Create a graph using an adjacency list representation
        # We'll use a defaultdict where each key is a node, and its value is a set of neighboring nodes
        graph = defaultdict(set)

        # Build the adjacency list by iterating through each pair in the input
        for x, y in self.unique_combinations:
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
        self.groups = find_connected_components(graph)

        # At this point, `groups` will contain:
        # [['A', 'B', 'C'], ['D', 'K']]
        

    def cofactors_abundance(self):
        
        self.merge_connected_combinations()
                
        # Count occurrences of each cofactor
        flat_list = [item for sublist in self.cofactors_list_all_reactions for item in sublist]
        element_counts = Counter(flat_list)  
        self.sorted_counts = element_counts.most_common()
        
        print(self.sorted_counts)
                
        
    def filter_reactions(self):
        
        self.cofactors_abundance()
        self.knockout_indices_all = []
        
        for group in self.groups:
            keep_indices_group = []
            cofactors_count_min = float('+inf')
            
            print(group)
            
            # find minimum number of cofactors across reactions of a single group
            for reaction in group:
                reaction_index = self.reactions_ids.index(reaction)   
                reaction_cofactors = self.cofactors_list_all_reactions[reaction_index]
                cofactors_count = len(reaction_cofactors)
                
                if cofactors_count < cofactors_count_min:
                    cofactors_count_min = cofactors_count

                print(reaction, reaction_cofactors, cofactors_count, cofactors_count_min)
                
            # find which reactions match the minimum number of cofactors
            for reaction in group:
                # find the index from the reaction list (general)
                reaction_index = self.reactions_ids.index(reaction)                   
                reaction_cofactors = self.cofactors_list_all_reactions[reaction_index]
                cofactors_count = len(reaction_cofactors)
   
                # equal to minimum value
                if cofactors_count == cofactors_count_min:
                    keep_indices_group.append(reaction_index)
                    
                # greater than minimum value ==> knockout
                else:
                    self.knockout_indices_all.append(reaction_index)

            
            # if 2 or more reactions match the minimum number of cofactors
            if len(keep_indices_group) > 1:
                # find which reaction has the most abundant cofactor (general abundance from the model) 
                # and keep her if the most abundant cofactor appears in all reactions, 
                # then check the next most abundant cofactor
                
                
                def order_cofactors_by_abundance(sorted_counts_data, specific_reaction_cofactors):
                    
                    """
                    Function that takes the cofactors of a single reaction and sorts them 
                    based on their overall abundance in the model

                    Example input:
                    
                    sorted_counts_data = [('h_c', 35), ('h2o_c', 18), ('h_e', 17), ('atp_c', 13), 
                                        ('adp_c', 12), ('nad_c', 12), ('nadh_c', 12), ('pi_c', 12)]
                                        
                    specific_reaction_cofactors = ['h_e', 'h_c', 'adp_c', 'nadh_c']
                    
                    Example output:
                    
                    [35, 17, 12, 12]
                    """

                    # Convert the list of tuples to a dictionary for fast lookup
                    abundance_dict = dict(sorted_counts_data)
                    
                    # Sort the items based on their abundance in the dictionary
                    sorted_cofactors = sorted(specific_reaction_cofactors, key=lambda x: abundance_dict.get(x, 0), reverse=True)
                    sorted_abundances = [abundance_dict.get(item, 0) for item in sorted_cofactors]
                    
                    return sorted_abundances
                     
                      
                def find_first_winning_sublist(sublists):
                    
                    """
                    Function to find which reaction to keep from a certain group 

                    Example input:
                    
                    sublists = [
                        [9, 3, 5, 7],
                        [9, 3, 4, 6],
                        [0, 5, 3, 3],
                        [4, 2, 1, 6],
                        [7, 3, 2, 0]
                    ]
                    
                    Example output:
                    
                    [9, 3, 5, 7]
                    
                    """
                    
                    num_sublists = len(sublists)
                    max_index = 0  # Start with the first sublist as the default winner
                    losers_indices = []  # List to store indices of non-winning sublists
                    
                    # Loop through the sublists to compare each one
                    for i in range(1, num_sublists):
                        for j, (a, b) in enumerate(zip(sublists[max_index], sublists[i])):
                            if a > b:
                                break  # Current winner is better; move to the next sublist
                            elif b > a:
                                max_index = i  # Update the winner index if the current sublist is better
                                break
                        else:
                            # If all comparable elements are equal, assign the first as a winner
                            max_index = 0
                    
                    # Identify non-winning sublists and store their indices
                    for i in range(num_sublists):
                        if i != max_index:
                            losers_indices.append(i)
                    
                    
                    # group_reactions_knockout = group[losers_indices]
                    for loser_index in losers_indices:
                        reaction_general_index = self.reactions_ids.index(group[loser_index])
                        self.knockout_indices_all.append(reaction_general_index)
                    
                            
                
                # apply the 'order_cofactors_by_abundance' function to reactions in the 'keep_indices' list
                # of the current group. Add all sorted abundance of cofactors for each reactions in a list
                # then call the 'find_first_winning_sublist' and add indices of the loser sublists into a list
                                    
                sorted_abundance_all_reactions = [
                order_cofactors_by_abundance(self.sorted_counts, self.cofactors_list_all_reactions[idx])
                for idx in keep_indices_group ]
                
                find_first_winning_sublist(sorted_abundance_all_reactions)                                
                    
                


cobra_model = load_json_model("../data/e_coli_core.json")
cofactor_specificity = Cofactor_Specificity(cobra_model)
cofactor_specificity.filter_reactions()