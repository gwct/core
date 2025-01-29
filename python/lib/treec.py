#############################################################################
# Functions to read and handle phylogenetic trees
# Gregg Thomas
#############################################################################

import sys
import re
import copy
import random
import itertools

#############################################################################

class Tree:
# The treeParse function takes as input a rooted phylogenetic tree with branch lengths and returns the tree with node labels and a
# dictionary with usable info about the tree in the following format:
# node:[branch length (if present), ancestral node, node type, node label (if present)]

    def __init__(self, tree_string, get_subtrees=False, debug=False):

        tree = tree_string.strip();
        if tree[-1] != ";":
            tree += ";";
        # Some string handling

        self.orig_tree_str = tree;
        self.topo_str = "";
        self.labeled_topo_str = "";
        self.tree_str = "";
        # Tree string types

        self.type = {};
        self.bl = {};
        self.has_bl = "";
        self.label = {};
        self.has_label = "";
        self.anc = {};
        self.desc = {};
        self.sis = {};
        self.subtree = {};
        subtrees = {};
        # Node attributes

        self.nodes = [];
        self.num_nodes = 0;

        self.tips = [];
        self.num_tips = 0;

        self.internals = [];
        self.num_internals = 0;

        self.num_polytomies = 0;

        self.rooted = "NA";
        self.root = "";
        # Node lists and counts

        ## Class attributes
        #####

        self.topo_str = re.sub('[)][\d\w<>/.eE_:-]+', ')', tree);
        self.topo_str = re.sub(':[\d.eE-]+', '', self.topo_str);
        ## Remove the branch lengths and node labels from the input tree string
        #####     

        for tip_label in self.topo_str.replace("(","").replace(")","").replace(";","").split(","):
            self.tips.append(tip_label);
            self.type[tip_label] = 'tip';
        ## Retrieval of the tip labels
        #####
       
        if debug:
            print();
            print("ORIGINAL:", self.orig_tree_str);
            print("TOPOLOGY:", self.topo_str);
        ## Node counting and root checking
        #####

        # self.bl = { node : "NA" for node in self.nodes };
        # self.label = { node : "" for node in self.internals + [self.root] };
        # self.anc = { node : "NA" for node in self.nodes };
        ## With the node lists set, initialize the other node attributes
        #####

        # z = 0;
        # numnodes = 1;
        # while z < (len(self.topo_str)-1):
        #     self.labeled_topo_str += self.topo_str[z];
        #     if self.topo_str[z] == ")":
        #         node_label = "<" + str(numnodes) + ">";
        #         self.labeled_topo_str += node_label;
        #         numnodes += 1;
        #     z += 1;
        # if debug:
        #     print("LABELED TOPOLOGY:", self.labeled_topo_str);
        #     print("----------");
        ## Add the generated internal node labels onto the topology string
        #####

        subtrees = {};

        pairs = re.findall("\([\d\w/.,:_<>-]+\)", tree);
        node_count = 1;
        while pairs:
        # Loop over all pairs of opening an closing parens at each
        # tree level

            for pair in pairs:
            # Loop over every pair in the current set

                anc_label = "<" + str(node_count) + ">";
                self.desc[anc_label] = [];
                node_count += 1;
                # The ancestral label

                node_list = pair.replace("(", "").replace(")", "").split(",");
                num_nodes = len(node_list);
                if num_nodes > 2:
                    self.num_polytomies += 1;
                # In most cases, the number of branches in a pair will be 2, but
                # for polytomies it will be more

                cur_nodes = [];
                # A list of nodes in the current pair to add
                # as descendants for the ancestral node

                for node in node_list:
                # Loop over each node in the current pair list

                    if ":" in node:
                        node, bl = node.split(":");
                    else:
                        bl = "NA"
                    # Check if the current node has a branch length

                    if node in self.tips:
                        self.label[node] = "NA";
                        self.bl[node] = bl;
                        # Tip info

                        self.desc[node] = "NA";
                        # Tips have no descendants

                        if bl != "NA":
                            subtrees[node] = node + ":" + str(bl);
                        else:
                            subtrees[node] = node;
                        # For tips, the subtree is just that node and its branch length
                        # If the current node is a tip label, add no label
                    ## Tips

                    else:
                        cur_label = re.sub("<[\d]+>", "", node);
                        node = re.findall("<[\d]+>", node)[0];
                        # Parse out the provided labels in the tree from
                        # the ones generated here for all internal nodes

                        self.type[node] = "internal";
                        self.internals.append(node);
                        self.label[node] = cur_label;
                        self.bl[node] = bl;
                        # Node information

                        subtrees[node] = "(";
                        for d in self.desc[node]:
                            subtrees[node] += subtrees[d] + ",";
                        subtrees[node] = subtrees[node][:-1] + ")" + node;

                        if cur_label:
                            subtrees[node] += str(cur_label);
                        if bl != "NA":
                            subtrees[node] += ":" + str(bl);
                        # For internal nodes, the subtree is the combination of the subtrees of its
                        # descendants
                        # Add branch lengths and labels as needed
                    ## Internal nodes

                    cur_nodes.append(node);
                    self.nodes.append(node);
                    self.anc[node] = anc_label;
                    # Add the ancestral label for the current node
                ## End node loop

                self.desc[anc_label] = cur_nodes;
                # Set the descendants of the current ancestor as all nodes in the current pair

                for node in cur_nodes:
                    self.sis[node] = [ n for n in cur_nodes if n != node ];
                # Set the sister branches for each node in the current pair

                tree = tree.replace(pair, anc_label);
                # Replace the current pair with the ancestral label in the tree
            ## End pairs loop

            pairs = re.findall("\([\d\w/.,:_<>-]+\)", tree);
            # Find all new pairs in the tree
        # If there are no pairs left, we've reached the root of the tree
        ## End tree level loop

        ## Main tree parsing loop
        #####

        self.nodes.append(anc_label);
        self.internals.append(anc_label);
        self.type[anc_label] = "internal";
        self.root = anc_label;
        self.label[anc_label] = "NA";
        self.bl[anc_label] = "NA";
        self.anc[anc_label] = "NA";
        self.sis[anc_label] = "NA";
        # Add the final anc label as the root

        if num_nodes > 2:
            self.rooted = False;
        else:
            self.rooted = True;
        # Check if there is a polytomy at the root for rootedness

        subtrees[anc_label] = "(";
        for d in self.desc[anc_label]:
            subtrees[anc_label] += subtrees[d] + ",";
        subtrees[anc_label] = subtrees[anc_label][:-1] + ")" + anc_label;
        self.tree_str = subtrees[self.root];
        # Add the subtree for the root, which ends up being the fully labeled tree

        if debug:
            print();
            print("NODES:", self.nodes);
            print("ROOTED:", self.rooted);
            print("ROOTNODE:", self.root);
        ## Root block
        #####

        if get_subtrees:
            self.subtrees = subtrees;
        # Save the subtrees if flagged

        if all(self.bl[n] == "NA" for n in self.bl):
            self.has_bl = False;
        else:
            self.has_bl = True;
            self.bl = { n : self.bl[n] if self.bl[n] != "NA" else "0.0" for n in self.bl };

        if all(self.label[n] == "" for n in self.label):
            self.has_label = False;
        else:
            self.has_label = True;
            self.label = { n : self.label[n] if self.label[n] != "" else "NA" for n in self.label };
        ## Checks for labels and bls
        #####

        self.num_tips = len(self.tips);
        self.num_internals = len(self.internals);
        self.num_nodes = len(self.nodes);
        ## Counts
        #####        

        if debug:
            print();
            self.showAttrib("label","anc","desc","sis");
            print("FINAL: " + self.tree_str);
            
            if get_subtrees:
                print();
                for node in self.subtrees:
                    print(node);
                    print(self.subtrees[node]);
                    print();

#############################################################################

# class Tree:
# # The treeParse function takes as input a rooted phylogenetic tree with branch lengths and returns the tree with node labels and a
# # dictionary with usable info about the tree in the following format:
# # node:[branch length (if present), ancestral node, node type, node label (if present)]

#     def __init__(self, tree_string, debug=False):

#         tree = tree_string.strip();
#         if tree[-1] != ";":
#             tree += ";";
#         # Some string handling

#         self.orig_tree_str = tree;
#         self.topo_str = "";
#         self.labeled_topo_str = "";
#         self.tree_str = "";
#         # Tree string types

#         self.type = {};
#         self.bl = {};
#         self.has_bl = "";
#         self.label = {};
#         self.has_label = "";
#         self.anc = {};
#         self.desc = {};
#         self.sis = {};
#         self.subtree = {};
#         # Node attributes

#         self.nodes = [];
#         self.num_nodes = 0;

#         self.tips = [];
#         self.num_tips = 0;

#         self.internals = [];
#         self.num_internals = 0;

#         self.rooted = "NA";
#         self.root = "";
#         # Node lists and counts

#         ## Class attributes
#         #####

#         self.topo_str = remBranchLength(tree);  

#         ## Remove the branch lengths and node labels from the input tree string
#         #####        

#         for tip_label in self.topo_str.replace("(","").replace(")","").replace(";","").split(","):
#             self.nodes.append(tip_label);
#             self.tips.append(tip_label);
#             self.type[tip_label] = 'tip';
        
#         ## Retrieval of the tip labels
#         #####

#         labeled_tree = "";
#         z = 0;
#         numnodes = 1;
#         while z < (len(tree)-1):
#             labeled_tree += tree[z];
#             if tree[z] == ")":
#                 node_label = "<" + str(numnodes) + ">";
#                 labeled_tree += node_label;
#                 self.nodes.append(node_label);
#                 self.internals.append(node_label);
#                 self.type[node_label] = 'internal';
#                 numnodes += 1;
#             z += 1;

#         self.internals = self.internals[:-1];
#         self.root = node_label;
#         # Get the root node as the last node read

#         ## Generate internal node labels and add them to the original tree string
#         #####

#         self.num_tips = len(self.tips);
#         self.num_internals = len(self.internals);
#         self.num_nodes = len(self.nodes);

#         self.rooted = self.checkRooted();
        
#         if debug:
#             print();
#             print("TREE:", self.orig);
#             print("TOPOLOGY:", self.topo_str);
#             print("NODES:", self.nodes);
#             print("ROOTED:", self.rooted);
#             print("ROOTNODE:", self.root);
#         ## Node counting and root checking
#         #####

#         self.bl = { node : "NA" for node in self.nodes };
#         self.label = { node : "" for node in self.internals + [self.root] };
#         self.anc = { node : "NA" for node in self.nodes };

#         ## With the node lists set, initialize the other node attributes
#         #####

#         z = 0;
#         numnodes = 1;
#         while z < (len(self.topo_str)-1):
#             self.labeled_topo_str += self.topo_str[z];
#             if self.topo_str[z] == ")":
#                 node_label = "<" + str(numnodes) + ">";
#                 self.labeled_topo_str += node_label;
#                 numnodes += 1;
#             z += 1;

#         if debug:
#             print("LABELED TOPOLOGY:", self.labeled_topo_str);
#             print("----------");

#         ## Add the generated internal node labels onto the topology string
#         #####

#         for node in self.nodes:
#         # One loop through the nodes to retrieve all other info
#             if debug:
#                 print("NODE:", node);

#             if node in self.tips:
#                 if node + ":" in tree:
#                     cur_bl = re.findall(node + ":[\d.Ee-]+", labeled_tree);
#                     cur_bl = cur_bl[0].replace(node + ":", "");
#                     if debug:
#                         print("FOUND BL:", cur_bl);
#                     self.bl[node] = cur_bl;
#                 ## If there are branch lengths, parse with regex and add to bl dict
#             ## Parse tips
             
#             elif node in self.internals:
#                 if node + "(" in labeled_tree or node + "," in labeled_tree or node + ")" in labeled_tree:
#                     if debug:
#                         print("NO BL OR LABEL");
#                 ## If there's no labels or branch lengths, keep the NAs

#                 elif node + ":" in labeled_tree:
#                     cur_bl = re.findall(node + ":[\d.Ee-]+", labeled_tree);
#                     cur_bl = cur_bl[0].replace(node + ":", "");
#                     if debug:
#                         print("FOUND BL:", cur_bl);
#                     self.bl[node] = cur_bl;
#                 ## If there's only branch lengths, add them to the bl dict                               

#                 else:

#                     cur_label_str = re.findall(node + "[\d\w<>_*+.Ee/-]+:[\d.Ee-]+", labeled_tree);
#                     # If this pattern is found there is both a branch length and a label for the nod

#                     if cur_label_str:
#                         cur_label = cur_label_str[0].replace(node, "");
#                         cur_label = cur_label[:cur_label.index(":")];
#                         #cur_bl = cur_label_str[0].replace(node, "").replace(cur_label, "").replace(":", "");
#                         cur_bl = cur_label_str[0][cur_label_str[0].index(":")+1:]
#                         if debug:
#                             print("FOUND BL AND LABEL:", cur_bl, cur_label);
#                         self.label[node] = cur_label;
#                         self.bl[node] = cur_bl;
#                     ## Parse if there is both bl and label
                        
#                     else:
#                         cur_label = re.findall(node + "[\w*+.<> -]+", labeled_tree);
#                         cur_label = cur_label[0].replace(node, "");
#                         if debug:
#                             print("FOUND LABEL:", cur_label);
#                         self.label[node] = cur_label;
#                     ## Parse if there is only a label
#                 ## Check if there are branch lenghts and labels or just labels
#             ## Parse internal nodes

#             elif node == self.root:
#                 possible_label = labeled_tree[labeled_tree.index(node)+len(node):];
#                 if possible_label:
#                     self.label[node] = possible_label;
#                 continue;
#             ## Parse the root and continue since there is no ancestor to parse in the block below

#             #####

#             anc_tree = labeled_tree[labeled_tree.index(node):][1:];
#             # Get the tree string starting at the current node label
#             # Ancestral labels are always to the right of the node label in the text of the tree, 
#             # so we start our scan from the node label

#             if debug:
#                 print("NODE:", node);
#                 print("ANC_TREE:", anc_tree);
                
#             cpar_count = 0;
#             cpar_need = 1;

#             for i in range(len(anc_tree)):
#                 if anc_tree[i] == "(":
#                     cpar_need = cpar_need + 1;
#                 if anc_tree[i] == ")" and cpar_need != cpar_count:
#                     cpar_count = cpar_count + 1;
#                 if anc_tree[i] == ")" and cpar_need == cpar_count:
#                     anc_tree = anc_tree[i+1:];
#                     self.anc[node] = anc_tree[:anc_tree.index(">")+1];
#                     break;
#                 # When the parentheses counts match, the ancestor will start at the next position and
#                 # end at the >
#             # We find the ancestral label by finding the ) which matches the nesting of the number of ('s found

#             if debug:
#                 print("FOUND ANC:", self.anc[node]);
#                 print("---");

#             #####
#         ## End node loop

#         if debug:
#             for node in self.anc:
#                 print(node, self.anc[node]);

#         #####

#         if all(self.bl[n] == "NA" for n in self.bl):
#             self.has_bl = False;
#         else:
#             self.has_bl = True;
#             self.bl = { n : self.bl[n] if self.bl[n] != "NA" else "0.0" for n in self.bl };

#         if all(self.label[n] == "" for n in self.label):
#             self.has_label = False;
#         else:
#             self.has_label = True;
#             self.label = { n : self.label[n] if self.label[n] != "" else "NA" for n in self.label };

#         for node in self.nodes:
#             self.desc[node] = self.getDesc(node);
#             self.sis[node] = self.getSister(node);
#         #self.subtree = self.genSubtrees();
#         # A few more useful node attributes to store

#         if not self.rooted:
#             self.internals += [self.root];
#         # For unrooted trees, add the "root" node as an internal node
#         # In this context, the "root" node is simply the last node read

#         #self.tree_str = self.subtree[self.root];
#         # Full tree with node names, labels, and branch lengths

    ##########

    def checkRooted(self):
    # Checks if a tree object is rooted or not by counting the difference
    # in the number of tip and internal nodes

        if (self.num_internals + 1) != (self.num_tips - 1):
            return False;
        elif (self.num_internals + 1) == (self.num_tips - 1):
            return True;
        else:
            return -1;

    ##########

    def getDesc(self, node):
        # This function takes a node in the current tree object
        # and returns a list of the two direct descendant nodes of it.

        if node in self.tips:
            return [node];
        else:
            return [ n for n in self.nodes if self.anc[n] == node ];

    ##########

    def getSister(self, node):
        # This function takes a node in the current tree object
        # and returns the other direct descendant of its ancestral node

        if node == self.root:
            return "NA";
        anc_desc = self.getDesc(self.anc[node]);
        return [ n for n in anc_desc if n != node ];

    ##########

    def getClade(self, node, full=False):
    # This function takes a node in the current tree object
    # and finds all tip labels that are descendants of it.
    # This is done by getting the direct descendants of the node with getDesc and then
    # recursively calling itself on those descendants.

        clade = [];
        if node in self.tips:
            return [node];

        desc = self.desc[node];
        for d in desc:
            if d not in self.tips:
                clade += self.getClade(d, full);
                if full:
                    clade.append(d);
                # If full is true, the function will also return all internal nodes
                # descending from the given node
            else:
                clade.append(d);

        return clade;

    ##########

    def getClades(self, full=False):
    # Calls getClade() on every node in a tree object

        clades = {};
        for node in self.nodes:
            clades[node] = set(self.getClade(node, full=full));

        return clades;

    ##########

    def getSplit(self, node):
    # Returns the tips from a node that do not descend from it to define a split (with clade)

        if node == self.root:
            return "NA";

        return set(self.tips) - set(self.getClade(node));

    ##########

    def getSplits(self):
    # Calls getSplit on every node in a tree object

        splits = {};
        for node in self.nodes:
            splits[node] = self.getSplit(node);

        return splits;

    ##########

    def getQuartet(self, node):
    # Returns 4 sets of tips given an internal node:
    # 1: tips descending from one descendant node
    # 2: tips descending from the other descendant node
    # 3: tips descending from the sister node
    # 4: all other tips not in the previous 3 categories

        if self.rooted:
            if node == self.root or self.type[node] == "tip":
                return "NA";
            else:
                if self.anc[node] == self.root:
                    if self.type[self.sis[node][0]] == "tip":
                        return "NA";
                    # If the sister is a tip, quartets cannot be sampled from this branch    

                    d1 = self.desc[node][0];
                    d2 = self.desc[node][1];

                    q1 = set(self.getClade(d1));
                    q2 = set(self.getClade(d2));

                    sis_desc = self.desc[self.sis[node][0]];
                    q3 = set(self.getClade(sis_desc[0]));
                    q4 = set(self.getClade(sis_desc[1]));
                    # If the sister is not a tip, return both of its descendants as q3 and q4              
                ## For nodes that have the root as an ancestor, the sister node needs
                ## to be parsed differently

                else:
                    d1 = self.desc[node][0];
                    d2 = self.desc[node][1];
                    sis = self.sis[node][0];

                    q1 = set(self.getClade(d1));
                    q2 = set(self.getClade(d2));
                    q3 = set(self.getClade(sis));
                    q4 = set(self.tips) - q1 - q2 - q3;
                ## For most nodes, just look up the clades of the descendant and sister nodes
        ## For rooted trees

        else:
            if node == self.root or self.type[node] == "tip":
                return "NA";
            else:
                d1 = self.desc[node][0];
                d2 = self.desc[node][1];
                q1 = set(self.getClade(d1));
                q2 = set(self.getClade(d2));

                if len(self.desc[self.anc[node]]) > 2:
                    other = [ n for n in self.desc[self.anc[node]] if n != node ];
                    q3 = set(self.getClade(other[0]));
                    q4 = set(self.getClade(other[1]));    
                else:
                    q3 = set(self.getClade(self.sis[node][0]));
                    q4 = set(self.tips) - q1 - q2 - q3;
            ## If the node has more than 2 descendants, simply use the first two...
        ## For unrooted trees        

        return { 'd1' : q1, "d2" : q2, "s" : q3, "q4" : q4 };

    ##########

    def getQuartets(self, root=True):
    # Calls getQuartet on every internal node in a tree object

        quartets = {};
        if root: 
            nodes = self.internals;
        else:
            nodes = self.internals[:-1];
        for node in nodes:
            quartets[node] = self.getQuartet(node);

        return quartets;

    ##########

    def sampleQuartets(self, num_quartets=100):
    # For sCF, we treat the species tree as unrooted, so for each node(*)/branch, the possible clades
    # to sample quartets from are the two clades directly descendant from the node, the clade
    # descendant from the sister node, and all other species.
    #
    #         /\
    #        /  \
    #       /\   \
    #      /  \   \
    #     *    \   \
    #    /\     \   \
    #  d1  d2   sis  other
    # ROOTED TREE
    #
    #   sister           descendant 1 (left)
    #         \         /
    #           -------*
    #         /         \
    #    other           descendant 2 (right)
    # UNROOTED TREE
    #

        full_quartets = self.getQuartets();
        sampled_quartets = {};

        for node in self.nodes:
            if node in self.tips or node == self.root:
                continue;
            # Cannot calculate sCF for tips, the root, or node descendant from the root

            assert all(len(full_quartets[node][q]) > 0 for q in full_quartets[node]), \
                " * ERROR: quartet sampling failed for node: " + node + "\n" + \
                "\tleft:   " + len(full_quartets[node]['d1']) + "\n" + \
                "\tright:  " + len(full_quartets[node]['d2']) + "\n" + \
                "\tsister: " + len(full_quartets[node]['s']) + "\n" + \
                "\tother:  " + len(full_quartets[node]['q4']) + "\n"
            # Make sure each clade list has species or throw an error.. this shouldn't happen

            split1_pairs = list(itertools.product(full_quartets[node]['d1'], full_quartets[node]['d2']));
            split2_pairs = list(itertools.product(full_quartets[node]['s'], full_quartets[node]['q4']));
            # Get every possible pair of species from each clade side of the split
            # i.e. all pairs of left-right species (split1) and all pairs of sister-other species (split2)

            quartets = list(itertools.product(split1_pairs, split2_pairs));
            # Get all pairs from the pairs of species in split1 and split2 for all possible quartets at this node

            cur_num_quartets = len(quartets);
            # Count the total number of quartets at this node

            if cur_num_quartets > num_quartets:
                random.shuffle(quartets);
                quartets = quartets[:num_quartets];
            # If there are more quartets on the current node than the number to sample, sub-sample here
            # Otherwise, use all quartets
            ## SET A SEED FOR REPRODUCIBILITY

            sampled_quartets[node] = quartets;
            # Add the current set of quartets to the global dict

        return sampled_quartets;

    ##########

    def genSubtrees(self):
    # Generates sub-tree strings for every node in a tree object

        subtrees = {};
        # Subtree dict
        
        for node in self.tips:
            if self.has_bl:
                subtrees[node] = node + ":" + str(self.bl[node]);
            else:
                subtrees[node] = node;
        # For tips, the subtree is just that node and its branch length

        for node in self.internals + [self.root]:
            subtrees[node] = "(";
            for d in self.desc[node]:
                subtrees[node] += subtrees[d] + ",";
            subtrees[node] = subtrees[node][:-1] + ")" + node;

            if self.has_label:
                subtrees[node] += str(self.label[node]);
            if self.has_bl and node != self.root:
                subtrees[node] += ":" + str(self.bl[node]);
        # For internal nodes, the subtree is the combination of the subtrees of its
        # descendants

        return subtrees;

    ##########

    def findClades(self, tip_set):
    # A function that takes a set of tips and finds all the monophyletic
    # clades within that set, and returns the LCA of each clade found

        tip_set = set(tip_set);

        clade_nodes = [];
        for node in self.internals:
            if set(self.getClade(node)) <= tip_set:
                clade_nodes.append(node);
        # In the input tree, find all nodes that have clades that are
        # a subset of the given tip set

        nodes_to_rm = [];
        for node in clade_nodes:
            for node_check in clade_nodes:
                if node == node_check:
                    continue;

                if node in set(self.getClade(node_check, full=True)):
                #if node in self.desc[node_check]:
                    nodes_to_rm.append(node);
        clade_set = [ n for n in clade_nodes if n not in nodes_to_rm ];
        # We only want the deepest node from each possible clade, so remove
        # nodes from those found that are descendants of another node.

        nodes_to_rm = [];
        for node in tip_set:
            if any(node in set(self.getClade(n, full=True)) for n in clade_set):
                nodes_to_rm.append(node);
        clade_set += [ n for n in tip_set if n not in nodes_to_rm ];
        # From the original input set, remove any nodes that are now found
        # in a clade

        return set(clade_set);

    ##########

    def findSplits(self, tip_set, clades=False, splits=False):
    # A function that takes a set of tips and finds whether any branches
    # are defined on either side by them

        split_nodes = [];
        for node in self.nodes:
            if clades:
                clade = clades[node];
            else:
                clade = set(self.getClade(node));
            # Lookup the clade for the given node

            if splits:
                split = splits[node];
            else:
                split = self.getSplit(node);
            # Lookup the split for the given node

            if clade == tip_set or split == tip_set:
                split_nodes.append(node);
            # If the clade or the split matches the tip set, add it to the list
            # of splits

        return set(split_nodes);

    ##########

    def LCA(self, node_list):
    # Given a list of nodes, this function finds the least common ancestor of them,
    # and tells whether the nodes provided form a monophyletic clade.

        ancs = {};
        for node in node_list:
            ancs[node] = [node];
        # For each node in the input list, we make a list of the path from that node 
        # to the root of the tree, including that node

        for node in node_list:
            if node == self.root:
                continue;

            cur_anc = self.anc[node];
            ancs[node].append(cur_anc);

            while cur_anc != self.root:
                cur_anc = self.anc[cur_anc];
                ancs[node].append(cur_anc);
        # For each node, add every node between it and the root to its ancs list

        intersect_anc = set.intersection(*list(map(set, list(ancs.values()))));
        # Gets the intersect of all lists of paths to the root, unordered

        lcp = [t for t in list(ancs.values())[0] if t in intersect_anc];
        # Orders the nodes in the intersect of all paths based on their order in
        # the path of an arbitrary node (since it should be identical for all of them)

        return lcp[0];
        # Returns the first node in the common path as the LCA  

    ##########

    def Monophyletic(self, node_list):
    # Determines whether a set of nodes is within a monophyletic clade -- do
    # they all descend from the same ancestor with no other descendants?

        monophyletic = False;
        if set(self.getClade(self.LCA(node_list))) == set(node_list):
            monophyletic = True;
        return monophyletic;

    ##########

    def addBranchLength(self):
    # Re-writes the branch lengths onto a tree object topology

        tree = self.labeled_topo_str;

        for node in self.nodes:
            new_node_str = node;
            if self.has_label and node in self.internals:
                new_node_str += str(self.label[node]);

            if self.has_bl and node != self.root:
                new_node_str += ":" + str(self.bl[node]);

            tree = tree.replace(node, new_node_str);

        return tree;

    ##########

    def addLabel(self, label_dict, delim=""):
    # Given a dictionary with { node : label } format, adds those labels and the branch
    # lengths onto the given tree's topology
    # Returns: tree string

        new_tree_str = self.tree_str;
        # Extract labeled tree topology to add labels to

        for node in self.nodes:
            old_label = node;
            new_label = node;
            if self.type[node] == "tip":
                old_label = node;
                new_label = node;
            # If the node is a tip, always add the node as a label

            if node in label_dict:
                if self.has_label:
                    old_label += self.label[node];

                if delim and self.has_label:
                    new_label += self.label[node] + delim + str(label_dict[node]);
                else:
                    new_label += str(label_dict[node]);
            # If the node is in the given dictionary, add the new label

            # if self.has_bl and node != self.root:
            #     new_label += ":" + str(self.bl[node]);
            # If the tree has branch lengths, add the bl

            if new_label or node == self.root:
                new_tree_str = new_tree_str.replace(old_label, new_label);
            # If a new label was created, replace the old node label in the
            # tree with it

        return new_tree_str;

    ##########

    def Prune(self, node_list, debug=False):
    # Prunes a tree string by reconstructing it from sub-trees, skipping those
    # that are to be pruned and combining branch lengths of the unpruned neighbors

    ## PRUNING UNROOTED TREES IS BUGGYYYYYYYY

        # unrooted = False;
        # if not self.rooted:
        #     unrooted = True;
        #     self = self.Root(self.tips[0]);            

        subtrees = {};
        # A dict of subtrees for each node in the tree

        pruned_desc_dict = copy.deepcopy(self.desc);
        # A local copy of the descendants dict to adjust for pruning

        node_list = list(self.findClades(set(node_list)));
        # Make sure all the nodes to prune represent monophyletic groups (else this algorithm will error)

        for node in self.tips:
            if self.has_bl:
                subtrees[node] = node + ":" + str(self.bl[node]);
            else:
                subtrees[node] = node;
        ## Add tips to subtree dict
        #####

        for node in self.internals:
            
            # if not self.rooted and node == self.root:
            #     continue;

            desc_to_prune = [d for d in pruned_desc_dict[node] if d in node_list];
            desc_to_keep = [d for d in pruned_desc_dict[node] if d not in node_list];
            # For the current node, get which descendants to prune and which to keep

            if len(desc_to_prune) in [0, 2]:
                subtrees[node] = "(";
                for d in pruned_desc_dict[node]:
                    subtrees[node] += subtrees[d] + ",";
                subtrees[node] = subtrees[node][:-1] + ")" + node;

                if self.has_label:
                    subtrees[node] += str(self.label[node]);
                if self.has_bl and node != self.root:
                    subtrees[node] += ":" + str(self.bl[node]);
            # If neither or both descendants are pruned, build the subtree for the current node

            elif desc_to_prune:
            # If one descendant is pruned, make adjustments
                if node == self.root:
                    continue;

                anc = self.anc[node];
                desc_to_keep = desc_to_keep[0];
                # Unpack the node to keep and its ancestor

                if self.has_bl:
                    new_bl = str(float(self.bl[desc_to_keep]) + float(self.bl[node]));

                    self.bl[desc_to_keep] = new_bl;
                    # Update the branch length of the node in the main tree class so it is carried over

                    subtrees[desc_to_keep] = subtrees[desc_to_keep][:subtrees[desc_to_keep].rfind(":")] + ":" + new_bl;
                # If the tree has branch lengths, we need to add the bl of the current node with the bl
                # of its ancestor and adjust the subtree to match
                            
                anc_desc_ind = pruned_desc_dict[anc].index(node);
                pruned_desc_dict[anc][anc_desc_ind] = desc_to_keep;
                # In the ancestral node, replace this node as a descendant with the node to keep 
        ## Construct subtrees for each internal node, adjusting when a node needs to be pruned
        #####

        if any(d in node_list for d in self.desc[self.root]):
        # If either of the nodes descending from the last (root) node are to be pruned, then the final
        # tree is just the subtree from the other node.
            if self.rooted:
                desc_to_keep = [d for d in pruned_desc_dict[self.root] if d not in node_list][0];
                pruned_tree = subtrees[desc_to_keep] + ";";
                # Find the node that is to be kept and lookup its subtree
                
                if self.has_bl:
                    pruned_tree = pruned_tree[:pruned_tree.rfind(":")];
                # Remove the trailing branch length which led to the root

            else:
                desc_to_keep = [d for d in pruned_desc_dict[self.root] if d not in node_list];

                pruned_tree = "";
                if len(desc_to_keep) > 1:
                    pruned_tree = "(";

                for d in desc_to_keep:
                    pruned_tree += subtrees[d] + ",";
                pruned_tree = pruned_tree[:-1];

                if len(desc_to_keep) > 1:
                    pruned_tree += ")";

                pruned_tree += ";"
            # For unrooted trees, if there is more than one descendant to keep, combine the subtrees
            # Otherwise, the final tree is just the subtree from that node to keep

        else:
        # If both descendants from the last (root) node are to be kept, combine them
        # to the final tree
            if self.rooted:
                l = subtrees[pruned_desc_dict[self.root][0]];
                r = subtrees[pruned_desc_dict[self.root][1]];
                # Lookup both descendant subtrees

                pruned_tree = "(" + l + "," + r + ");";

            else:
                l = subtrees[pruned_desc_dict[self.root][0]];
                m = subtrees[pruned_desc_dict[self.root][1]];
                r = subtrees[pruned_desc_dict[self.root][2]];

                pruned_tree = "(" + l + "," + m + "," + r + ");";
            # Combine the subtrees and add the last branch length if necessary
        ## Make the final tree by combining the subtrees descending from the last (root) node
        #####

        pruned_tree = re.sub("<[\d]+>", "", pruned_tree);
        # Remove the added node labels

        # pruned_tree = Tree(pruned_tree);

        # if unrooted:
        #     pruned_tree = pruned_tree.Unroot;

        # if ",:" in pruned_tree:
        #     print()

        return pruned_tree;

    ##########

    def rmTips(self, debug=False):
    # Removes all tips from a tree, leaving internal branches as new tips

        subtrees = {};
        # A dict of subtrees for each node in the tree

        pruned_desc_dict = copy.deepcopy(self.desc);
        # A local copy of the descendants dict to adjust for pruning

        for node in self.internals:
            desc_to_prune = [d for d in pruned_desc_dict[node] if d in self.tips];
            desc_to_keep = [d for d in pruned_desc_dict[node] if d not in self.tips];
            # For the current node, get which descendants to prune and which to keep

            if len(desc_to_prune) == 2:
                subtrees[node] = node;
                if self.has_label:
                    subtrees[node] += str(self.label[node]);
                if self.has_bl and node != self.root:
                    subtrees[node] += ":" + str(self.bl[node]);
            # If both descendants are tips, then the subtree for this node is the node label

            elif len(desc_to_prune) == 0:
                subtrees[node] = "(";
                for d in pruned_desc_dict[node]:
                    subtrees[node] += subtrees[d] + ",";
                subtrees[node] = subtrees[node][:-1] + ")" + node;

                if self.has_label:
                    subtrees[node] += str(self.label[node]);
                if self.has_bl and node != self.root:
                    subtrees[node] += ":" + str(self.bl[node]);
            # If neither descendant is a tip, add the two descendant subtrees together as usual

            else:
            # If one descendant is a tip, make adjustments
                if not self.rooted and node == self.root:
                    continue;

                anc = self.anc[node];
                desc_to_keep = desc_to_keep[0];
                # Unpack the node to keep and its ancestor

                if self.has_bl:
                    new_bl = ":" + str(float(self.bl[desc_to_keep]) + float(self.bl[node]));
                    subtrees[desc_to_keep] = subtrees[desc_to_keep][:subtrees[desc_to_keep].rfind(":")] + new_bl;
                # If the tree has branch lengths, we need to add the bl of the current node with the bl
                # of its ancestor and adjust the subtree to match
                            
                anc_desc_ind = pruned_desc_dict[anc].index(node);
                pruned_desc_dict[anc][anc_desc_ind] = desc_to_keep;
                # In the ancestral node, replace this node as a descendant with the node to keep                
        ## Construct subtrees for each internal node, adjusting when a node needs to be pruned
        #####

        if any(d in self.tips for d in self.desc[self.root]):
        # If either of the nodes descending from the last (root) node is a tip, then the final
        # tree is just the subtree from the other node.
            if self.rooted:
                desc_to_keep = [d for d in pruned_desc_dict[self.root] if d not in self.tips][0];
                pruned_tree = subtrees[desc_to_keep] + ";";
                # Find the node that is to be kept and lookup its subtree
                
                if self.has_bl:
                    pruned_tree = pruned_tree[:pruned_tree.rfind(":")];
                # Remove the trailing branch length which led to the root

            else:
                desc_to_keep = [d for d in pruned_desc_dict[self.root] if d not in self.tips];
                pruned_tree = "(";
                for d in desc_to_keep:
                    pruned_tree += subtrees[d] + ",";
                pruned_tree = pruned_tree[:-1] + ");";

        else:
        # If both descendants from the last (root) node are to be kept, combine them
        # to the final tree
            if self.rooted:
                l = subtrees[pruned_desc_dict[self.root][0]];
                r = subtrees[pruned_desc_dict[self.root][1]];
                # Lookup both descendant subtrees

                pruned_tree = "(" + l + "," + r + ");";

            else:
                l = subtrees[pruned_desc_dict[self.root][0]];
                m = subtrees[pruned_desc_dict[self.root][1]];
                r = subtrees[pruned_desc_dict[self.root][2]];

                pruned_tree = "(" + l + "," + m + "," + r + ");";
            # Combine the subtrees and add the last branch length if necessary
        ## Make the final tree by combining the subtrees descending from the last (root) node
        #####

        pruned_tree = re.sub("<[\d]+>", "", pruned_tree);
        # Remove the added node labels

        return pruned_tree;

    ##########

    def Unroot(self):
    # A function to convert a rooted tree to an unrooted tree

        if not self.rooted:
            print("tree is already unrooted");
            return self;
        # If the tree is already unrooted, just return it

        subtrees = self.genSubtrees();
        # Parse the subtrees at each node in the tree

        root_desc = self.desc[self.root];
        # Get the branches descending directly from the root

        split_node = root_desc[1];
        sis_node = root_desc[0];
        # We will take one of the clades descending from the root and split it to unroot
        # while the other (its sister) will be intact
        
        if split_node in self.tips:
            split_node = root_desc[0];
            sis_node = root_desc[1];
        # We can't split a descendant clade if it is just a single tip, so
        # check for that and switch here

        if split_node in self.tips:
            return False;
        # If both descendants from the root are tips, then this is a simple, unrooted tree (a,b);

        split_desc = self.desc[split_node];
        # Get the descendants of the clade to split

        unrooted_tree_str = "(" + subtrees[sis_node] + "," + subtrees[split_desc[0]] + "," + subtrees[split_desc[1]] + ");";
        # Concatenate the intact sister descendant from the root, and both clades from the other descendant into a single unrooted
        # tree string

        unrooted_tree = Tree(unrooted_tree_str);
        # Parse the tree string into a tree class

        return unrooted_tree;

    ##########

    def Root(self, node_list):
    # A function to root or re-root a tree
    # Node list can either be multiple tips or a single internal node

        if self.rooted:
            self = self.Unroot();
        # If the tree is rooted, unroot it first

        if all(node in self.tips for node in node_list):
        # If a list of tips is given as input
        
            if self.Monophyletic(node_list):
                new_root = self.LCA(node_list);
            # If the tips form a monophyletic group, then their LCA is the new root node

            else:
                if self.LCA(node_list) == self.root:
                # If the LCA of the given tips is the root (of an unrooted tree), then we
                # try to root on the remaining tips

                    possible_split = [ tip for tip in self.tips if tip not in node_list ];
                    # Get all the other tips to check if they form a monophyletic group

                    if self.Monophyletic(possible_split):
                        print("rooting on ingroup")
                        new_root = self.LCA(possible_split);
                    # If the remaining tips form a monophyletic group, the new root node is their LCA

                else:
                    print("cannot root");
                    return False;
                # If we can't root on either the provided tips or their complement, we cannot root this tree because
                # the given tips are paraphyletic
            ##

        elif len(node_list) > 1:
            print("can't root on both tips and internal nodes or multiple internal nodes")
            return False;
        # If the list isn't all tips and is longer than 1 entry, we can't root it

        else:
            new_root = node_list[0];
        # If the given node list is 1 entry long and it is an internal node, root on that node

        subtrees = self.genSubtrees();
        # Parse the subtrees at each node in the tree

        anc = self.anc[new_root];
        sis = self.sis[new_root];
        # Get the ancestral node and sister node of the new root

        if self.has_bl:
            root_bl = float(self.bl[new_root]);
            new_root_bl = ":" + str(root_bl / 2);
            # For the new branch length, just divide the old one by 2... probably not the best...

            tree_str = "(" + subtrees[new_root].replace(":" + self.bl[new_root], new_root_bl);
            # Start by adding the subtree at the new root to the new tree string

        else:
            new_root_bl = "";
            tree_str = "(" + subtrees[new_root];
            # Start by adding the subtree at the new root to the new tree string
        ## Initialize the new rooted tree with the root and its branch length (if the input tree had them)


        closings = [];
        # Closings is a list of closing parentheses (and branch lengths)
        # As each subtree is added to the tree string, an opening parentheses is added and will
        # need to be closed

        while anc != self.root:
        # This loops through every node, starting from the new root node

            tree_str +=  ",(" + subtrees[sis[0]];
            # Add the subtree from the SISTER node to the tree string, again as a sister

            if self.has_bl:
                closings.append("):" + self.bl[anc]);
            else:
                closings.append(")");
            # Add a closing parentheses (and branch length) to the closings list for this subtree

            sis = self.sis[anc];
            anc = self.anc[sis[0]];
            # Get the next SISTER branch, which is sister of the ancestor to the current ancestor
            # and the next ancestor
        ## End subtree rooting loop            

        tree_str += ",(" + subtrees[sis[0]] + "," + subtrees[sis[1]] + "".join(reversed(closings)) + ")" + new_root_bl + ");";
        # Add both other branches that are sister to the current branch descending from the root of the tree, then add
        # the closing parentheses and the remainder of the original branch length

        rooted_tree = Tree(tree_str);

        return rooted_tree;

    ##########

    def showAttrib(self, *attrib_list):
    # This function simply prints out the given tree attributes to the screen
    # Valid attributes are: "type", "desc", "anc", "sis", "clade", "split", "quartet"

        print("-" * 60);
        ## Seperator

        pad = 40;
        ## Width of columns

        valid_attribs = ["type", "length", "label", "desc", "anc", "sis", "clade", "split", "quartet"];
        ## List of valid attributes

        attrib_list = [ attrib for attrib in attrib_list if attrib in valid_attribs ];
        attrib_rm = [ attrib for attrib in attrib_list if attrib not in valid_attribs ];
        if attrib_rm:
            print("WARNING: The following are not valid Tree attributes to display: " + ",".join(attrib_rm));
        ## Parse passed attributes and warn if any are invalid

        headers = ["NODE"] + [ attrib.upper() for attrib in valid_attribs if attrib in attrib_list ];
        outline = [ spacedOut(header, pad) for header in headers ];
        print("".join(outline).strip());
        ## Display the attributes as headers

        header_len = (pad * len(attrib_list)) + len(attrib_list[-1]);
        print("-" * header_len);
        ## Seperator between headers and rows

        for node in self.nodes:
            outline = spacedOut(node, pad);
            for attrib in valid_attribs:
                if attrib in attrib_list:
                    if attrib == "type":
                        if node == self.root:
                            outline += spacedOut(self.type[node] + ",root", pad);
                        else:
                            outline += spacedOut(self.type[node], pad);
                    ## type

                    if attrib == "label":
                        if self.type[node] == "tip":
                            outline += spacedOut("", pad);
                        else:
                            outline += spacedOut(self.label[node], pad);
                    # label

                    if attrib == "length":
                        outline += spacedOut(self.bl[node], pad);
                    # length

                    if attrib == "desc":
                        if node in self.tips:
                            outline += spacedOut(self.desc[node], pad);
                        else:
                            outline += spacedOut(",".join(self.desc[node]), pad);
                    ## desc

                    if attrib == "anc":
                        outline += spacedOut(self.anc[node], pad);
                    ## anc

                    if attrib == "sis":
                        if node == self.root:
                            outline += spacedOut(self.sis[node], pad);
                        else:
                            outline += spacedOut(",".join(self.sis[node]), pad);
                    ## sis

                    if attrib == "clade":
                        clade = str(set(self.getClade(node)));
                        outline += spacedOut(clade, pad);
                        if len(clade) > pad:
                            outline += "\t";
                    ## clade
                        
                    if attrib == "split":
                        split = str(self.getSplit(node))
                        outline += spacedOut(split, pad);
                        if len(split) > pad:
                            outline += "\t";
                    ## split

                    if attrib == "quartet":
                        if self.type[node] == "tip" or node == self.root:
                            outline += "NA";
                        else:
                            #outline += spacedOut("", pad);
                            quartet = self.getQuartet(node);

                            outline += str(quartet["d1"]);

                            for key in ['d2', 's', 'q4']:
                                outline += "\t" + str(quartet[key]);
                    ## quartet

                ## End attrib loop
            ## End valid loop
        ## End node loop
            
            print(outline.strip());
        print("-" * 60);
        ## Seperator      

## END TREE CLASS
#############################################################################
## BEGIN TREE STRING FUNCTIONS

def remBranchLength(tree_str):
# Removes branch lengths and labels from a tree string

    tree_str = re.sub('[)][\d\w<>/.eE_:-]+', ')', tree_str);
    tree_str = re.sub(':[\d.eE-]+', '', tree_str);

    return tree_str;

#############################################################################

def getSubtree(node, tree_str):
# Gets the subtree string at a given node from a labeled tree string from treeParse
# Much slower than genSubtrees method

    subtree = "";
    # Initialize the subtree string

    partree = tree_str[:tree_str.index(node)][::-1]
    # Slice the tree string at the index of the node label and reverse it

    cp = 0;
    op = 0;
    # Counts of closing an opening parentheses. The subtree will be complete when they
    # are equal

    for c in partree:
        if c == ")":
            cp = cp + 1;
        if c == "(":
            op = op + 1;
        subtree = subtree + c;
        if cp == op:
            break;
    # Loop through every character in the sliced and reversed tree, counting parentheses and adding
    # charcters to the subtree string one at a time. Stop when the number of closing parentheses equals
    # the number of opening

    return subtree[::-1];
    # Return the reverse of the subtree, which is the original orientation

#############################################################################

# def adjustTreeDict(tree_dict, root):
# # Gets some info about branch lengths from the tree dict and retrieves support
# # label

#     num_no_supp = 0;
#     # Number of nodes with no support

#     bls = [];
#     # A list of all branch lengths in the tree

#     na_bl_spec = [];
#     # A list of species that have no branch length (i.e. tips from Astral)

#     for n in tree_dict:
#         if n == root:
#             continue;
#         # Skip the root since it won't have a length or support

#         if tree_dict[n][0] == "NA":
#             tree_dict[n][0] = 0.0;
#             na_bl_spec.append(n);
#         # If the current node has no branch length, set it as 0.0 in the dictionary for later and
#         # add to the na_bl_spec dict to add the NA back in later. Also don't add this to the list
#         # of bls to calculate the percentile from
#         else:
#             tree_dict[n][0] = float(tree_dict[n][0]);
#             bls.append(tree_dict[n][0]);
#             # Add the floated bl to the bl list to calculate percentile. The NA branch lengths should
#             # be left out
#         # Convert the branch lengths in the tree dict to floats

#         if tree_dict[n][3] == 'NA':
#             num_no_supp += 1;
#             continue;
#         cur_supp = tree_dict[n][3];
#         if "/" in cur_supp:
#             cur_supp = cur_supp.split("/")[1];
#         cur_supp = cur_supp.replace("_", "");
#         tree_dict[n][3] = float(cur_supp);
#         # Parse out the gCF from the supports, skip if there is no support
#         # TODO: This is still very specific to my input format

#     return tree_dict, bls, na_bl_spec;

#############################################################################

def mapNodes(query_tree, query_clades, query_splits, ref_tree, ref_clades, ref_splits):
# Species must have same tip labels, except for missing/pruned species
# Smaller tree must be query tree

    node_map = {};
    # Node map dict

    for query_node in query_tree.internals:
    # We map every internal query node in the query tree

        max_node = "";
        max_node_set = {};
        # We map nodes by finding the clade in reference tree that has 
        # that the query clade is a subset of and that maximizes number of 
        # nodes in the reference tree

        query_clade = query_clades[query_node];
        query_split = query_splits[query_node];
        # Look up the query clade and split for the current node

        for ref_node in ref_tree.internals:
        # Check every internal ref node in the ref tree

            check_node = False;
            # Flag if the current query clade is contained in the current ref clade

            ref_clade = ref_clades[ref_node];
            ref_split = ref_splits[ref_node];
            # Look up the ref clade and split

            if not query_tree.rooted and query_node == query_tree.root:
                if query_clade <= ref_clade or query_clade <= ref_split:
                    check_node = True;
            # If the tree is NOT rooted and the query node is the "root", check if
            # the query node is a subset ofeither side of the ref node

            elif query_clade <= ref_clade and query_split <= ref_split:
                check_node = True;
            # If the tree is rooted, check whether the query split is a subset on both
            # sides of the ref split

            # elif query_tree.split[query_node] <= ref_tree.clade[ref_node] and query_tree.clade[query_node] <= ref_tree.split[ref_node]:
            #     check_node = True;                    

            if check_node:
                node_set = ref_clade | ref_split;
                if len(node_set) > len(max_node_set):
                    max_node_set = node_set;
                    max_node = ref_node;
            # If we find the query node matches the ref node, here we check if it is the deepest
            # possible ref node by counting the tips and comparing it to the previous match
        ## End ref node loop

        node_map[query_node] = max_node;
        # Save the deepest matching node as the map
    ## End query node loop

    rev_map = { v : k for k,v in node_map.items() if v != "" };
    # In some cases we may want the map from the query tree to the ref tree, though
    # not all nodes will map

    return node_map, rev_map;

#############################################################################

def spacedOut(string, totlen, sep=" "):
# Properly adds spaces to the end of a message to make it a given length
    spaces = sep * (totlen - len(string));
    return string + spaces;  

#############################################################################

def debugTree(globs):

    import lib.core as CORE

    step = "Reading input species tree";
    step_start_time = CORE.report_step(globs, step, False, "In progress...");
    # Status update

    if globs['st-input-type'] == "file":
        globs['orig-st-str'] = open(globs['st-input'], "r").read().strip();
    else:
        globs['orig-st-str'] = globs['st-input']
    # If the input type is a file, read the file here, otherwise take the string as the tree input

    print("\n\n");

    tree_str = globs['orig-st-str']
    #tree_str = "((a,b),c,(d,e));" # Unrooted
    #tree_str = "(((a,b),c),(d,e));" # Rooted
    #tree_str = ("((a,b),c)")
    #tree_str = "(((a,b),c),(((d,e),f),h),(i,j),(k,l));";
    #tree_str = "((((a,b),c),(((d,e),f),h),i),j);";
    #tree_str = "(mspr:0.0018190933,((mwsb:0.0018060852,((mcas:0.0005908566,mpwk:0.0000009867):0.0000076480,mmus:0.0000009867):0.0000009867):0.0011628549,(((((((((((cgri:0.0152686696,((pcam:0.0012137705,psun:0.0006517732):0.0083093402,prob:0.0087755222):0.0112747594):0.0027806008,maur:0.0194531500):0.0098942322,moch:0.0330046438):0.0022180414,pman:0.0173033334):0.0096202861,((jjac:0.0617468174,itri:0.0714744806):0.0197552895,ngal:0.0425089249):0.0299773246):0.0088542882,mung:0.0286172457):0.0205036203,rnor:0.0311801045):0.0011873972,((rdil:0.0199545510,gdol:0.0117952834):0.0120430301,rsor:0.0202951611):0.0008975499):0.0043168585,(hall:0.0119687298,(pdel:0.0121964374,mnat:0.0120591970):0.0031186198):0.0070437051):0.0091272024,mpah:0.0130568501):0.0093423702,mcar:0.0073123519):0.0001651381):0.0030773597,mspi:0.0005888274);"
    #tree_str = "((jjac:0.1003165447,itri:0.1044182968)0.891:0.010544,(ngal:0.069324448,((pman:0.0406478328,(moch:0.0497980904,((cgri:0.0250402305,maur:0.0328120476)0.668:0.00352051,(prob:0.0144107446,(pcam:0.0028058377,psun:0.0031617491)0.977:0.0129734762)0.984:0.0237499666)0.93:0.0123999599)0.395:0.0020879791)0.901:0.0112615439,(mung:0.0594804999,(rnor:0.0396699258,(rsor:0.0265720262,((rdil:0.019967991,gdol:0.021105056)0.769:0.0085977828,((hall:0.0161568992,(pdel:0.0154147548,mnat:0.0156744369)0.42:0.0026038866)0.86:0.0102277821,(mpah:0.0196602678,(mcar:0.0092755425,((mmus:0.0019110393,mpwk:0.0051314407999999995)0.629:0.002015403,mspi:0.0049513531)0.808:0.0044021912)0.839:0.0074024138)0.821:0.0100953924)0.576:0.0041320338)0.387:0.0023496865)0.671:0.0051006134)0.982:0.0242160502)0.775:0.0069231995)0.98:0.0418120466)0.0:0.010544);";
    tree_str = "(mm10:0.0003564481,((((((((((((((((((Anonymomys_mindorensis_FMNH2222010:0.0050386474,(Niviventer_lepturus_JAE4746:0.0091125409,(Niviventer_rapit_LSUMZ38196:0.0000010072,Niviventer_bukit_LSUMZ38136:0.0000010072)100:0.0043913495)100:0.0006369911)100:0.0005408369,Leopoldomys_sabanus_Z25191:0.0079023684)100:0.0016674827,(Margaretamys_parvus_NMVZ21833:0.0006902620,Margaretamys_christinae_TMar1:0.0018652480)100:0.0021109997)100:0.0045438766,Saxatilomys_paulinae_ABTC116266:0.0068792519)100:0.0028417948,((((((((Taeromys_celebensis_NMVZ20637:0.0013401157,Paruromys_dominator_JAE4870:0.0035925299)100:0.0019976354,Eropeplus_canus_NMVZ21733:0.0016548458)93:0.0005572723,(((Sundamys_muelleri_JAE1283:0.0000010072,Sundamys_annandalei_YS118:0.0012186300)100:0.0011068134,Lenomys_meyeri_LSUMZ39079:0.0016578859)95:0.0005564736,Genus_sp-nov_NMVZ21816:0.0011053067)30:0.0000028202)33:0.0000010071,((Halmaheramys_bokimekot_33262-WD225:0.0027704015,Bullimus_bagobus_JAE1712:0.0050198446)35:0.0000022713,((Bunomys_torajae_NMVZ21736:0.0005438317,Bunomys_andrewsi_sample1:0.0012083604)96:0.0000075197,Bunomys_chrysocomus_JAE4867:0.0005520733)100:0.0011027590)9:0.0000025470)8:0.0000010071,((Taeromys_taerae_LSUMZ39224:0.0000010072,Taeromys_hamatus_NMVZ63317:0.0000010072)99:0.0011566003,(Taeromys_punicans_MVZ225815:0.0000010072,Taeromys_callitrichus_LSUMZ37153:0.0000010072)100:0.0022146614)99:0.0005564230)19:0.0000022656,Papagomys_armandvillei_WAMM32592:0.0055870600)98:0.0016861359,(((((rnor6:0.0028298450,Limnomys_sibuanus_FMNH147948:0.0038941816)92:0.0005610837,(Rattus_xanthurus_LSUMZ39196:0.0018552956,((((Rattus_tunneyi_RAT130:0.0000010072,Rattus_praetor_sample1:0.0012265319)58:0.0000027008,(((Rattus_lutreolus_sample2:0.0011004214,Rattus_colletti_RAT25:0.0000010072)95:0.0000010071,Rattus_villosissimus_sample1:0.0017137023)100:0.0005542109,Rattus_leucopus_sample1:0.0022808135)56:0.0000028774)99:0.0005632869,Kadarsanomys_sodyi_MZB-Sample:0.0024854938)95:0.0005371273,Rattus_morotaiensis_ASAM29:0.0005541689)91:0.0005510229)88:0.0000017353)99:0.0005508264,(Rattus_exulans_sample1:0.0283793809,Rattus_argentiventer_LSUMZ38199:0.0017084937)75:0.0000010071)74:0.0000010071,Rattus_hoffmanni_LSUMZ39187:0.0000010072)98:0.0011051989,(Nesokia_indica_ABTC117074:0.0004951904,Bandicota_indica_ABTC119185:0.0019190720)100:0.0039724363)99:0.0044922891)99:0.0047931170,Berylmys_bowersi_ABTC116431:0.0054971741)98:0.0019235154)86:0.0005644883,(((((Waiomys_mamasae_Z20624:0.0010861691,Gracilimus_radix_Z08508:0.0029571126)100:0.0063028248,Sommeromys_macrorhinos_Z21740:0.0068438549)71:0.0004861625,(Tateomys_macrocercus_Z21834:0.0036288262,Melasmothrix_naso_Z21808:0.0008285755)98:0.0008431259)73:0.0005450849,((Paucidentomys_vermidax_Z21914:0.0039756312,Echiothrix_centrosa_Z27993:0.0011952311)100:0.0006532961,Hyorhinomys_stuempkei_JAE4661:0.0031011021)68:0.0000022774)86:0.0000025866,Tateomys_rhinogradoides_Z21785:0.0012033978)100:0.0027059259)98:0.0040580267,(((Crunomys_melanius_KU168558:0.0011275489,Crunomys_celebensis_Z21894:0.0016961576)100:0.0055038648,Maxomys_musschenbroekii_sample1:0.0059793357)100:0.0025137005,Maxomys_surifer_JAE3632:0.0104312938)99:0.0023021403)98:0.0066143541,((((((((((Thamnomys_venustus_LSUMZ37935:0.0096054165,Thamnomys_kempi_TCD4437:0.0029137067)100:0.0033212501,(Grammomys_dryas_LSUMZ37716:0.0125539005,(Aethomys_nyikae_FMNH196708:0.0024408518,Aethomys_hindei_MVZ196257:0.0039796475)100:0.0023094270)93:0.0006307089)98:0.0010633738,Thallomys_paedulcus_FMNH158089:0.0093795603)89:0.0005667842,((((((Pelomys_minor_FMNH219213:0.0015499803,Pelomys_isseli_FMNH137886:0.0025902198)100:0.0031706829,(Lemniscomys_striatus_TCD4711:0.0024670057,(Arvicanthus_niloticus_MNHN1999201:0.0032505861,Arvicanthus_neumanni_FMNH158037:0.0000010072)100:0.0044798003)100:0.0043925011)100:0.0005730607,Mylomys_dybowskii_MNHN1997072:0.0025355870)99:0.0000629964,Rhabdomys_pumilio_MNHN1995065:0.0050591120)100:0.0040865401,(Dasymys_montanus_FMNH144389:0.0000010072,Dasymys_incomtus_LSUMZ37702:0.0000010072)100:0.0078803041)99:0.0005555579,Micaelemys_namaquensis_MNHN1987009:0.0057417291)74:0.0000025068)100:0.0011266698,Oenomys_hypoxanthus_LSUMZ37838:0.0123346009)74:0.0000020698,(((Hybomys_univittatus_TCD4026:0.0006980667,Hybomys_lunaris_LSUMZ37734:0.0000010072)100:0.0066401405,(Stochomys_longicaudatus_MNHN1999459:0.0034141260,Dephomys_defua_MNHNR24062:0.0036277447)99:0.0020274150)100:0.0040952634,Hybomys_planifrons_R24117:0.0085817477)92:0.0000028470)99:0.0019142499,(Otomys_tropicalis_LSUMZ37841:0.0029323676,Otomys_denti_LSUMZ37840:0.0050084528)100:0.0074511163)98:0.0024409990,Pithecheir_melanurus_LSUMZ38198:0.0169721886)93:0.0019329066,Vandeleuria_oleracea_ABTC116404:0.0206277771)95:0.0014798667,(((Phloeomys_pallidus_JEC43:0.0004916650,Phloeomys_cumingi_FMNH198472:0.0031448938)100:0.0143388898,((Musseromys_anacuao_FMNH209523:0.0053179885,Carpomys_melanurus_FMNH198723:0.0086958966)100:0.0081237068,(Crateromys_heaneyi_KU165014:0.0052816947,Batomys_russatus_KU165776:0.0059547788)100:0.0033983305)99:0.0018819004)100:0.0124097643,(Lophuromys_woosnami_LSUMZ37793:0.0195090284,Lophiomys_imhausi_UM5152:0.0361739705)100:0.0122182286)96:0.0041971200)53:0.0010239875)47:0.0000026730,((((((((((((((Mesembriomys_macrurus_ABTC7337:0.0013103799,Mesembriomys_gouldii_ABTC7449:0.0016918315)100:0.0026186657,Conilurus_albipes_NMV-c7585:0.0017047680)59:0.0000025682,Conilurus_penicillatus_M14056:0.0022276491)90:0.0005527487,(((((((Pseudomys_albocinereus_WAMM18855:0.0016726293,Pseudomys_apodemoides_M16672:0.0028962151)86:0.0005451980,(((Pseudomys_oralis_RAT7:0.0011091976,(((Pseudomys_novaehollandiae_ABTC08140:0.0022157251,Pseudomys_australis_M25541:0.0011282052)78:0.0000020886,(Pseudomys_bolami_ABTC107389:0.0027984403,(Pseudomys_higginsi_UC1421:0.0024726625,Pseudomys_hermannsburgensis_ABTC08067:0.0017279500)55:0.0000010071)72:0.0000024498)79:0.0000023374,Pseudomys_delicatulus_U1509:0.0005973140)89:0.0005550486)92:0.0000559453,Pseudomys_shortridgei_Z25113:0.0000010072)30:0.0000023211,Pseudomys_desertor_ABTC41464:0.0005809169)58:0.0004886312)29:0.0000010071,((Pseudomys_nanus_U1448:0.0005621579,Pseudomys_gracilicaudatus_ABTC08163:0.0016960157)100:0.0028127011,Pseudomys_fumeus_ABTC08168:0.0016621072)23:0.0000010071)96:0.0006559211,Pseudomys_occidentalis_ABTC08043:0.0028828871)97:0.0013487763,((Pseudomys_chapmani_WAMM29440:0.0022631655,((Pseudomys_patrius_JM10865:0.0000021113,Pseudomys_johnsoni_ABTC08053:0.0005876467)100:0.0018678885,Pseudomys_calabyi_U4324:0.0045191259)95:0.0011167157)99:0.0016467705,((((Notomys_mitchellii_M21518:0.0000010072,Notomys_fuscus_M22830:0.0005952391)76:0.0000010071,Notomys_aquilo_M-32455:0.0000010072)76:0.0000010071,Notomys_alexis_U1308:0.0006158271)100:0.0033214780,Notomys_cervinus_M22857:0.0033607483)54:0.0006014685)58:0.0013560963)98:0.0013447466,((((Zyzomys_woodwardi_WAMM21724:0.0101461668,Zyzomys_maini_ABTC07900:0.0035720137)87:0.0010824090,(Zyzomys_palatalis_ABTC30744:0.0030892092,Zyzomys_argurus_ABTC61630:0.0011312085)100:0.0036424235)99:0.0014075178,Zyzomys_pedunculatus_Z34925:0.0033740124)79:0.0000025723,(Leggadina_lakedownensis_WAMM21623:0.0006421711,Leggadina_forresti_WAMM62323:0.0011683609)100:0.0083238138)99:0.0011159087)100:0.0027003267,(Uromys_sp-nov_YS380:0.0016565831,Uromys_caudimaculatus_DT10Uc01:0.0011211213)100:0.0044874370)92:0.0006702667)68:0.0000027029,(((Xeromys_myoides_ABTC30709:0.0050329415,((Parahydromys_asper_ABTC45798:0.0000010072,(Hydromys_sp-nov_YS391:0.0000010072,Hydromys_chrysogaster_Z5248:0.0000010072)100:0.0031033852)99:0.0024837544,Crossomys_moncktoni_ABTC46614:0.0028039460)99:0.0027774643)62:0.0000020861,(Pseudohydromys_pumehanae_ABTC43920:0.0182044961,Microhydromys_richardsoni_ANWCM35851:0.0016972492)90:0.0033776105)87:0.0005576494,(Paraleptomys_wilhelmina_H126-ABTC133716:0.0011063434,(Mirzamys_louiseae_ABTC133754:0.0038995095,Leptomys_elegans_ABTC45741:0.0027746374)50:0.0000010071)89:0.0005484154)100:0.0027783678)71:0.0006450135,(Solomys_salebrosus_ABTC64864:0.0037950591,Melomys_rufescens_ABTC44798:0.0028961927)100:0.0042379254)100:0.0032513703,Coccymys_ruemmleri_ABTC49489:0.0226810555)67:0.0000023526,(((Xenuromys_barbatus_ABTC129204:0.0169019693,Mammelomys_lanosus_ABTC47208:0.0098483198)100:0.0028056679,(Pogonomelomys_mayeri_AMM21849:0.0056417434,Abeomelomys_sevia_KUM161018:0.0016624479)100:0.0016560878)100:0.0007312521,Mallomys_rothschildi_ABTC47402:0.0066653187)98:0.0015213963)79:0.0011205207,(Lorentzimys_nouhuysi_ABTC-46314:0.0098385566,Anisomys_imitator_ABTC45107:0.0039652549)52:0.0011010092)37:0.0000027587,Macruromys_major_ABTC43909:0.0091206240)41:0.0002938785,((Pogonomys_macrourus_ABTC43144:0.0123576946,Hyomys_goliath_ABTC42697:0.0079268166)68:0.0013830432,Chiruromys_vates_ABTC43096:0.0051157773)47:0.0006752480)97:0.0012539332,(((Rhynchomys_isarogensis_JAE2195:0.0033859572,Archboldomys_luzonensis_EAR1826:0.0035896145)63:0.0000021200,(Chrotomys_whiteheadi_KU168311:0.0006468897,Chrotomys_mindorensis_JAE520:0.0000010072)100:0.0059196381)95:0.0033774371,Apomys_lubangensis_JAE169:0.0064808591)100:0.0056982870)100:0.0005224163,Haeromys_minahassae_JAE4996:0.0120257274)99:0.0014312385,Chiropodomys_gliroides_Z25153:0.0168322213)95:0.0025791671)76:0.0006383061,Malacomys_longipes_TCD4253:0.0172680532)72:0.0005485985,Apodemus_sylvaticus_SMG3823:0.0230003368)93:0.0010422110,(((((Zelotomys_woosnami_TTU-M115260:0.0021539614,Colomys_goslingi_TCD4175:0.0017886760)100:0.0008993289,((Tarsomys_apoensis_KU167148:0.0183939279,Berylmys_mackenziei_ABTC116402:0.0046769407)99:0.0250562005,(Praomys_degraaffi_TCD4312:0.0048070690,Myomyscus_verreauxii_MNHN1999102:0.0051788691)68:0.0000021224)100:0.0502468276)57:0.0000020137,Praomys_jacksoni_REC31:0.0035174507)63:0.0000027654,Mastomys_natalensis_T4162:0.0067168805)94:0.0016381565,((Hylomyscus_vulcanorum_TCD4382:0.0051106389,Hylomyscus_stella_LSUMZ37742:0.0039853414)100:0.0022108609,Heimyscus_fumosus_FMNH167810:0.0045000126)99:0.0023108517)98:0.0111111505)100:0.0068495664,Mus_pahari_10460X1:0.0045598593)94:0.0011295370,(Mus_platythrix_10460X13:0.0059807461,(Mus_musculoides_LSUMZ37820:0.0030117047,Mus_minutoides_10460X2:0.0039949805)100:0.0060807088)52:0.0000026844)100:0.0016712566,((Mus_cookii_10460X3:0.0027891368,Mus_cervicolor_10460X10:0.0033382703)100:0.0005599928,Mus_caroli_10460X9:0.0036535119)99:0.0016565961)99:0.0016811367,(Mus_spretus_10460X4:0.0016939387,Mus_macedonicus_10460X12:0.0011108994)100:0.0005590220)99:0.0011118513,(Mus_musculus-reference_10460X5:0.0000010072,Mus_musculus-castaneus_10460X7:0.0005654919)84:0.0000020567)74:0.0001810016,Mus_musculus-musculus_10460X14:0.0000010072)85:0.0000135183,Mus_musculus-domesticus_10460X15:0.0005527500);";
    st = Tree(tree_str, get_subtrees=False, debug=False);
    # Some test trees

    #gt_str = "(mspr:0.0016324816,((mwsb:0.0000029568,mcar:0.0153129287):0.0021281260,(mcas:0.0019633923,(mmus:0.0000020851,((((hall:0.0174240612,(pdel:0.0130029362,mnat:0.0129445878):0.0029314306):0.0103158970,((((mung:0.0649440245,(((((((pcam:0.0031385524,psun:0.0033197989):0.0260127153,prob:0.0098921329):0.0196507658,cgri:0.0281001191):0.0019130533,maur:0.0355556774):0.0047041281,moch:0.0402020094):0.0057618080,pman:0.0379711174):0.0163580030,((jjac:0.0776442621,itri:0.1479591099):0.0249499942,ngal:0.0723486821):0.0508154463):0.0069335015):0.0123544338,rnor:0.0413244014):0.0021945352,(rdil:0.0203480143,gdol:0.0196458955):0.0075988551):0.0020342947,rsor:0.0270726716):0.0054596486):0.0105131536,mpah:0.0220127233):0.0095631120,mpwk:0.0009405335):0.0014971022):0.0021557143):0.0016426560):0.0028428833,mspi:0.0008185431);";

    #st = Tree(globs['orig-st-str'], debug=False);
    #gt = TREE.Tree(gt_str, debug=False);
    #st = TREE.Tree(tree_str);
    
    # st.subtree = st.genSubtrees();
    # print(st.subtree["<176>"]);
    # sys.exit();
    #st.tree_str = st.subtree[st.root];
    print(st.tree_str);

    # #print(globs['orig-st-str']);

    prune_list = ['rnor6', 'Genus_sp-nov_NMVZ21816', 'Zyzomys_palatalis_ABTC30744', 'Zyzomys_maini_ABTC07900', 'Pseudomys_patrius_JM10865', 'Pseudomys_johnsoni_ABTC08053', 'Pseudomys_chapmani_WAMM29440', 'Pseudomys_calabyi_U4324', 'Leggadina_lakedownensis_WAMM21623', 'Uromys_sp-nov_YS380', 'Leporillus_apicalis_C90', 'Hydromys_sp-nov_YS391', 'Arvicanthus_niloticus_MNHN1999201', 'Hybomys_lunaris_LSUMZ37734', 'Pithecheir_melanurus_LSUMZ38198', 'Praomys_degraaffi_TCD4312', 'Mus_spretus_10460X4', 'Mus_spicilegus_10460X11', 'Mus_macedonicus_10460X12', 'Mus_musculus-musculus_10460X14', 'Mus_musculus-castaneus_10460X7', 'mm10', 'Mus_musculus-domesticus_10460X15']

    # # print(st.labeled_topo_str);
    # # st.showSplit();
    st.showAttrib("length", "anc", "desc");
    print(st.rooted);
    print(st.root);
    #print(st.tree_str);
    #print(st.tips);
    #pruned_tree = st.rmTips()
    #print(pruned_tree);

    #prune_list = ["a", "b", "c"]
    pruned_tree = st.Prune(prune_list);
    print(pruned_tree);

    #st_quartets = st.getQuartets(root=False);


    # text = open("check.txt", "r").read().split("\n")[:-1];
    # check = [];
    # for line in text:
    #     #print(line);
    #     spec = line.split("\t");
    #     if spec[1] == "1":
    #         check.append(spec[0]);

    # for tip in st.tips:
    #     #print(tip);
    #     if tip in ["mm10", "rnor6"]:
    #         continue;
    #     tip = tip.split("_");
    #     tip = tip[0] + " " + tip[1];
    #     if tip not in check:
    #         print(tip);

    # print(st.num_tips)

    # 
    # print(st.root);
    # print(st.num_tips);

#############################################################################
