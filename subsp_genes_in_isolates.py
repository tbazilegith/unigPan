#!/usr/bin/env python3

# Script to find unique genes of a target strain(subspecies) from a pangenome (ex. roary output)

# Author: Tassy J-S. Bazile, Bioinformatician, CDC-APHL Bioinformatics Fellow 2021-2023

gg = """
    Business Logic:
    The program subsp_genes_in_isolates.py take (1) the genes  of an isolate from the genes_presence_absence file of all (strains)subspecies, the (2) defined unique genes (file) of a subspecies (or strain) to identify them in an isolate. 
    """ 

"""   
Excution:
    python subsp_genes_in_isolates.py --dir_gene_pres PATH-TO-CORE-GENES_inter --isolates' ISOLATE_to_compare --uniqueGenes FILENAME.csv  
"""

import os
import sys
import argparse
import pandas as pd
import json
import itertools
import fileinput
import glob
import re

#sys.path.append("/blue/bphl-florida/t.bazile1/aphl_bi_fellow/data/legionella/roary_dir/unique_genes_programs/uni_gene_module")
#from uni_gene_module.finding_unique_genesFromPan03_extract import gene_N_annotation
#from finding_unique_genesFromPan03_extract import gene_N_annotation  # to generate table genes - function
# Parsing argument
parser = argparse.ArgumentParser(prog = 'subsp_genes_in_isolates.py [Options]', description=gg)

#parser.add_argument('--interSubsp_cgenes', type = argparse.FileType('r'), help =' provide the core gene file of all strains or its path', required=True)
#parser.add_argument('--intraSubsp_cgenes', type = argparse.FileType('r'), help =' provide the core gene file of the subspecies or its path', required=True)
#parser.add_argument('--dir_allCgenes', type = str, help =' provide the directory path to the core gene files for all strains(subspecies)', required=True)
parser.add_argument('--dir_gene_pres', type = str, help =' provide the directory path to the genes_pres_abs file for the all isolates and subspecies', required=True)
parser.add_argument('--isolates', type = str, help =' provide the directory the isolate or sample ID to compare', required=True)
parser.add_argument('--uniqueGenes', type=str,help= ' provide a file name such as <<*_unique_genes.sv>> or its path',required=True )
args = parser.parse_args()


# Input file gene_pres_abs for genes and sample
genes_samp = args.dir_gene_pres
 # directory path to all core genes files
#list_cgfle = os.listdir(p_cgenes)
#ll=[p_cgenes + i for i in list_cgfle]
#tup_file= tuple(ll)

# Input for isolate to compare
isolate_column = args.isolates # core genes for target subspecies or strain
#intra_cgenes = args.dir_targetcgenes # dir path core genes for target subspecies or strain
#list_intrafile = os.listdir(intra_cgenes)

#input defined unique genes
def_uniq_genes = args.uniqueGenes
# output files
#un_gen_st = args.uniqueGenes 
#un_gen_stdout = un_gen_st.split('/')[-1] # to print on screen 
#un_gen = args.uniqueGenes + '_unique_genes.txt'
#un_genlist=open(un_gen, 'w')
#print(ll)

#print(list_cgfle) 

# All core genes
# open gene_presence_absence to get list all the genes
def open_gene_pres(path_to_csv):
    gene_pres_df = pd.pandas.read_csv(path_to_csv)
    #selecting column
    gene_pres_df2 = gene_pres_df[["Gene", "Annotation", str(isolate_column)]]
    gene_pres_df3 = gene_pres_df2.dropna(subset = [str(isolate_column)]) # drop missing values
    isol_genes_list = gene_pres_df3.Gene.values.tolist()
    isol_annot_list = gene_pres_df3.Annotation.values.tolist()    
    #genes_list = gene_pres_df.Gene.values.tolist()
    #isolate_gene_protein = gene_pres_df.isolate_column.values
    return isol_genes_list, isol_annot_list
 
# Filter isolate's genes
def iso_gene_defined(isol_genes_list):
   start_name='group'
   gene_with_group = [x for x in isol_genes_list if re.match(r'^{}'.format(start_name), x)]
   gene_without_group = [x for x in isol_genes_list if x not in gene_with_group]
   num_gene_def = len(gene_without_group)
   return gene_without_group, num_gene_def

# open subsp unique genes
def open_unique_genes(path_uniq_genes):
    uniq_df = pd.pandas.read_csv(path_uniq_genes)
    genes_list = uniq_df.Unique_Gene.values.tolist()
    return genes_list

# Genes check in isolate
def genes_in_isolates(genes_list, isol_genes_list):
    subs_isol_genes = [i for i in genes_list if i in isol_genes_list]
    return subs_isol_genes
 
# Generate table gene-annotation
def iso_gene_annotation(isol_genes_list, isol_annot_list,subs_isol_genes):
    isol_gene_annot = dict(zip(isol_genes_list, isol_annot_list))
    isol_gene_annotation = [isol_gene_annot[g] for g in subs_isol_genes]
    iso_gene_annot_def = pd.DataFrame({"Unique_Gene": subs_isol_genes,
                                      "Annotation" : isol_gene_annotation})
    return iso_gene_annot_def

#Creating csv file
def iso_gene_annotation_to_csv(iso_gene_annot_def):
    iso_gene_annot_def.to_csv(isolate_column +'_'+ def_uniq_genes.split('_')[0]+'_genes.csv', sep=',', index=False, header=True)



"""
#Opening core_genes_result files (all files in a directory)
def open_cgenes(p_cgenes,list_cgfle):
    #filpath = [os.path.join(p_cgenes , cgr) for cgr in list_cgfle]
    #all_cg_files =[]
    ll=[p_cgenes + i for i in list_cgfle] # list of files and their paths
    tup_file= tuple(ll)# fileinput.input() takes tuple
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

# converting each list into tuplee and dictionary
def core_genes(core_genes_l):
    #ll=[g.split(":") for g in core_genes1] # each gene group as a list
    #tpgenes=[tuple(gn) for gn in all_merg] # each list as a tuple
    #dd=dict(tpgenes)  # all tuples as a dictionary
    #genes_key=[k for k in dd.keys()]
    set_coreg=set(core_genes_l)
    list_setcoreg=list(set_coreg)
    number_cgenes = len(list_setcoreg)
    return  number_cgenes

# list of core genes for each subspecies

def unique_genes(intra_cg, inter_cg):
    unique_target = list(set(intra_cg).difference(inter_cg))
    number_intracg = len(unique_target)
    return unique_target, number_intracg 

# looking for uniques or complement genes between list
#unique_raph = list(set(raph_geneslist).difference(coregenes_list))
"""
def main():
    # Intersubspecies core genes
    isol_genes, isol_annot_list = open_gene_pres(genes_samp)
    #print(isol_genes) #  contains undined genes
    
    # filtered genes in isolates
    isol_filtered_genes, number_genes = iso_gene_defined(isol_genes)
    #print(isol_filtered_genes)
    print("Isolate filtered genes for ",isolate_column, ": ", number_genes)
    #number_cg1 = core_genes(intercg)
    # open defined unique genes in suspecies
    genes_list = open_unique_genes(def_uniq_genes)
    print("defined unique genes in subspecies: ",'\n', genes_list)
    # subspecies genes in isolate
    subs_isol_genes = genes_in_isolates(genes_list,isol_filtered_genes)
    print("suspecies genes in isolate: ",'\n', subs_isol_genes)
    # Isolate_genes_annotation
    iso_gene_annot_def = iso_gene_annotation(isol_genes, isol_annot_list,subs_isol_genes)
    print("subspecies genes in isolates and function: ",'\n', iso_gene_annot_def)
    # csv file
    iso_gene_annotation_to_csv(iso_gene_annot_def) 
    # Intrasubspecies core genes
    #intracg = open_cgenes(list_intrafile,intra_cgenes)
    #listoffile2, intracg = open_cgenes(intra_cgenes,list_intrafile)
    #genes_list = open_gene_pres(genes_target)
    #number_cg2 = core_genes(genes_list)
    
    # unique genes for target
    #unique_genesTg, number_intracg = unique_genes(genes_list, intercg)

    #print("Number of unique genes for strain / subspecies ", un_gen_stdout, ":", number_intracg, "\n")
    #print(listoffile) #
    #print(intercg)
   # print("numbre of genes shared by all strains:",number_cg1)
   # print(genes_list)
    #print(intracg)
    #print("numbre of genes for target strain:",number_cg2)
    #print("Number of unique genes for strain / subspecies ", un_gen_stdout, ":", number_intracg, "\n")
    #print("unique genes: ", unique_genesTg)
    #un_genlist.write('Number of unique genes: ' + str(number_intracg) +'\n') 
    #un_genlist.write('\n'.join(unique_genesTg))

if __name__ == "__main__":
    main()
# closing the file
#un_genlist.close()

