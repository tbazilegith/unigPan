#!/usr/bin/env python3

# Script to find unique genes of a target strain(subspecies) from a pangenome (ex. roary output)

# Author: Tassy J-S. Bazile, Bioinformatician, CDC-APHL Bioinformatics Fellow 2021-2023

gg = """
    Business Logic:
    The progrom finding_unique_genesFromPan03_extract.py uses (1)core genes for all (strains)subspecies and  the (2)core genes within a subspecies (or strain) to identify unique genes for a target strain (subspecies). 
    """ 

"""   
Excution:
    python find_unique_genes.py --dir_allCgenes PATH-TO-CORE-GENES --dir_targetSubsp PATH-TO-GENE-PRES-ABS --uniqueGenes UNIQUEGENES

"""

import os
import argparse
import pandas as pd
import json
import itertools
import fileinput
import glob
import re
# Parsing argument
parser = argparse.ArgumentParser(prog = 'finding_unique_genesFromPan_extract.py [Options]', description=gg)

#parser.add_argument('--interSubsp_cgenes', type = argparse.FileType('r'), help =' provide the core gene file of all strains or its path', required=True)
#parser.add_argument('--intraSubsp_cgenes', type = argparse.FileType('r'), help =' provide the core gene file of the subspecies or its path', required=True)
parser.add_argument('--dir_allCgenes', type = str, help =' provide the directory path to the core gene files from pangenome between target and other strains(subspecies)', required=True)
parser.add_argument('--dir_targetSubsp', type = str, help =' provide the directory path to the  genes_presence_absence file for the target strains(subspecies)', required=True)
#parser.add_argument('--dir_targetcgenes', type = str, help =' provide the directory path to the core genes file for the target strain(subspecies)', required=True)
parser.add_argument('--uniqueGenes', type=str,help= ' provide a file name such as <<yourTargetStrain>> and its path to write the unique genes',required=True )
args = parser.parse_args()

#inter_cgenes = args.interSubsp_cgenes

# Input all core genes from all (other) strains
p_cgenes =args.dir_allCgenes # directory path to all core genes files
list_cgfle = os.listdir(p_cgenes) # list of files
#ll=[p_cgenes + i for i in list_cgfle]
#tup_file= tuple(ll)

# Input core genes for target strain or subspecies
genes_target = args.dir_targetSubsp # core genes for target subspecies or strain
#intra_cgenes = args.dir_targetcgenes # dir path core genes for target subspecies or strain
#list_intrafile = os.listdir(intra_cgenes)

# output files
un_gen_st = args.uniqueGenes 
un_gen_stdout = un_gen_st.split('/')[-1] # to print on screen 
un_gen = args.uniqueGenes + '_unique_genes.txt'
un_genlist=open(un_gen, 'w')
#print(ll)

#print(list_cgfle) 

# All core genes
# open gene_presence_absence to get list all the genes
def open_gene_pres(path_to_csv):
    gene_pres_df = pd.pandas.read_csv(path_to_csv,header=0)
    target_genes_list = gene_pres_df.Gene.values.tolist()
    annotation_list = gene_pres_df.Annotation.values.tolist()
    return target_genes_list, annotation_list
# All core genes


#Opening core_genes_result files (all files in a directory)
def open_cgenes(p_cgenes,list_cgfle):
    #filpath = [os.path.join(p_cgenes , cgr) for cgr in list_cgfle]
    #all_cg_files =[]
    ll=[p_cgenes + i for i in list_cgfle] # list of files and their paths
    tup_file= tuple(ll)# fileinput.input() takes tuple of the list of paths
    #for path in filpath:
        #with open(path, 'r') as cg_fh : 
    with fileinput.input(tup_file) as f: 
        inter_cgenesList= [line for line in f]   
            #all_cg_files.append(inter_cgenesList)       
        clean_list= [sub.replace("\n", "") for sub in inter_cgenesList]
        core_genes_l=[g.split(":")[0] for g in clean_list] # a list for each file
        #llist_cg.append(core_genes_l) # list of list
        #all_merg = list(itertools.chain.from_iterable(llist_cg)) # merging all inner lists
        #return all_merg # was one return
        return ll, core_genes_l
#lll =open_cgenes(ll)
#print(lll) 

# Getting unique elements in list
def core_genes(core_genes_l):
    #ll=[g.split(":") for g in core_genes1] # each gene group as a list
    #tpgenes=[tuple(gn) for gn in all_merg] # each list as a tuple
    #dd=dict(tpgenes)  # all tuples as a dictionary
    #genes_key=[k for k in dd.keys()]
    set_coreg = set(core_genes_l) # removing duplication
    list_intercoreg = list(set_coreg)
    number_cgenes = len(list_intercoreg)
    return  number_cgenes, list_intercoreg

# List of unique genes  for each subspecies
#Difference between genes in target train and intersubspecies core genes (intersection)
def unique_genes(target_genes_list,list_intercoreg):
    uniqueG_target = list(set(target_genes_list).difference(list_intercoreg))
    number_uniG = len(uniqueG_target)
    return uniqueG_target,number_uniG 

# looking for uniques or complement genes between list
#unique_raph = list(set(raph_geneslist).difference(coregenes_list))

# filter defined genes from list
def gene_defined(uniqueG_target):
   start_name='group'
   gene_with_group = [x for x in uniqueG_target if re.match(r'^{}'.format(start_name), x)]
   gene_without_group = [x for x in uniqueG_target if x not in gene_with_group]
   num_gene_def = len(gene_without_group)
   return gene_without_group, num_gene_def

# Genes and Annotation
# using target list genes, annotation, and unique genes
def gene_N_annotation(target_genes_list, annotation_list, gene_without_group):
    gene_annotation = dict(zip(target_genes_list, annotation_list))
    uniqueG_annotat = [gene_annotation[g] for g in gene_without_group]
    uniqueG_annotation_df = pd.DataFrame({"Unique_Gene": gene_without_group,
                                          "Annotation" : uniqueG_annotat})
    return uniqueG_annotation_df

#Creating csv file
def gene_annotation_to_csv(uniqueG_annotation_df):
    uniqueG_annotation_df.to_csv(args.uniqueGenes + '_defined_unique_genes.csv', index=False)

# Execution
def main():
    # Intersubspecies core genes
    listoffile, intercg = open_cgenes(p_cgenes,list_cgfle)
    number_cg1,list_intercoreg = core_genes(intercg)

    # Intrasubspecies core genes
    #intracg = open_cgenes(list_intrafile,intra_cgenes)
    #listoffile2, intracg = open_cgenes(intra_cgenes,list_intrafile)
    genes_list, annotat = open_gene_pres(genes_target)
    number_cg2,list_targgene = core_genes(genes_list)
    
    # unique genes for target
    unique_genesTg, number_intracg = unique_genes(genes_list, list_intercoreg)
    # defined genes
    def_genes, num_gene_def = gene_defined(unique_genesTg)
    # Unique Genes and Annotation
    uniqueG_annotation_df  = gene_N_annotation(genes_list, annotat, def_genes)
    # csv file
    gene_annotation_to_csv(uniqueG_annotation_df) 
    #print("Number of unique genes for strain / subspecies ", un_gen_stdout, ":", number_intracg, "\n")
    print(listoffile) #
    print(uniqueG_annotation_df)
    #print(intercg)
    print("number of genes shared by all strains:", number_cg1)
    #print(genes_list)
    #print(intracg)
    print("numbre of genes for target strain:", number_cg2)
    print("Number of unique genes including group_XXX for strain / subspecies ", un_gen_stdout, ":", number_intracg, "\n")
    print("Number of defined unique genes for strain / subspecies ", un_gen_stdout, ":", num_gene_def,"\n")
    print("defined unique genes: ", def_genes)
    #print("unique genes: ", unique_genesTg)
    un_genlist.write('Number of unique genes: ' + str(number_intracg) +'\n') 
    un_genlist.write('\n'.join(unique_genesTg))

if __name__ =="__main__":
    main()
# closing the file
un_genlist.close()

