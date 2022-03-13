#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import sys
import yaml 
import os
import pandas as pd
import xlsxwriter
from datetime import datetime, date
from connection import get_connection
from query_patients import query4_patients
import numpy as np
from scipy.stats import fisher_exact
import scipy.stats as stats
import datetime
import time
import argparse  


# In[ ]:


"""
python FE_test.py /shares/archive/develdata/ajames/GD_PD/results/lancet/pats/ /shares/archive/develdata/ajames/GD_PD/results/lancet/control/ 
"""
__author__ = 'Alva James'

parser = argparse.ArgumentParser()

parser.add_argument('input',
                        help='input PD patient path from unidb'
                        )
parser.add_argument('control',
                        help='control, can be restricted or clear control from unidb: just the path'
                        )

parser.add_argument('-output',
                        help='output path')
parser.add_argument('-control_type',
                        help='output path')
args = parser.parse_args()

if len(sys.argv) < 2:
    parser.print_help()
    sys.exit(0)
        
pd_pats       = args.input
contrl_pats   = args.control
outpath       = args.output
control       = args.control_type


# In[ ]:


#pd_pats      = "/shares/archive/develdata/ajames/GD_PD/results/lancet/pats/"
#contrl_pats  = "/shares/archive/develdata/ajames/GD_PD/results/lancet/control/"
files_pats   = os.listdir(pd_pats)
files_ctrl   = os.listdir(contrl_pats)
files_pat    = [os.path.join(pd_pats,i) for i in files_pats]
files_ctrl   = [os.path.join(contrl_pats,i) for i in files_ctrl]
len(files_pat)


# In[ ]:


total_control_pd = 8971
clear_control    = 22918
total_pats_pd    = 3272


# In[ ]:


cols =['chrom','vcf_pos','filter_type','gene_name','patientid','country_name','analysis_tag','dob','gnomad_genomes_af','gnomad_exomes_af','cadd_raw']
def merger(files):
    """
    Merge all inputs from unidb
    """
    data_frames = []
    for filename in files:
        df    = pd.read_csv(filename,header=0,sep='\t')
        name  = os.path.basename(os.path.normpath(filename))
        genename    = name.split('_')[2]
        filter_type = name.split('_')[3]
        filter_type = filter_type.split('.')[0]
        #df.rename(columns={'id':'varid'},inplace=True)
        df['gene_name']  =genename
        df['filter_type']=filter_type
        df               =df[cols]
        data_frames.append(df)       
    merged_df  = pd.concat(data_frames)
    return merged_df

def get_unique_patients(files_pat,string1):
    """
    Get the number of unique patients from control and cases
    """
    gene_dfs =[]
    dcts = {}
    for filename in files_pat:
        df      = pd.read_csv(filename, sep='\t', header=0)
        if df.shape[0] >0:            
            name = os.path.basename(os.path.normpath(filename))
            genename    = name.split('_')[2]
            filter_type = name.split('_')[3]
            filter_type = filter_type.split('.')[0]
            #print(filename,df.head(),filter_type,genename)
            #print(filter_type)
            patient_id=df.patientid.drop_duplicates().shape[0]
            print(filename,filter_type,genename,patient_id)
            #print(genename,filter_type,patient_id)
            print("Preparing dataframe for",genename,"with",patient_id, "patients", "for filter type", filter_type)
            dct = {   
            'gene_name' : genename,
            string1: patient_id, 'filter_type':filter_type}
        #dcts     = {k:[v] for k,v in dct.items()}
            gene_dfs.append(dct) 
    gene_df = pd.DataFrame.from_dict(gene_dfs)
    return gene_df

def get_final_df(df,total,string1,newstring):
    """
    Make the final table with unique number of patients from three different filter modes
    The filter modes are class, filter and cadd
    """
    class_df= df.query('filter_type=="class"')
    class_df[newstring]=class_df[string1].apply(lambda x:total-x)
    filter_df =df.query('filter_type=="filter"')
    filter_df[newstring]=filter_df[string1].apply(lambda x:total-x)
    cadd_df=df.query('filter_type=="cadd"')
    cadd_df[newstring]=cadd_df[string1].apply(lambda x:total-x)
    return class_df,filter_df,cadd_df

#alternative=greater: means then the "greater" direction is to the right 
#alternative= less, and the "lesser" direction is to the left from that table.
# is important to put the group expected to have higher odds of the event in the first column, here we are expecting to have higher OR in PD cohort
#An OR higher than 1 means that the first group (in this case, control) was more likely to experience the enrichment than the second group
def fisher_test(control_withoutvars, control_withvars, PD_withoutvar, PD_withvars):
    """
    FE test for Pval
    """
    contingency_table = [[PD_withvars,PD_withoutvar],[control_withvars,control_withoutvars]]
    return fisher_exact(contingency_table, alternative='greater')[1]
def odds_ratio(control_withoutvars, control_withvars, PD_withoutvar, PD_withvars):
    """
    FE test for odds ratio
    """
    contingency_table = [[PD_withvars,PD_withoutvar],[control_withvars,control_withoutvars]]
    oddsratio, pvalue = stats.fisher_exact(contingency_table,alternative='greater')
    return oddsratio


# In[ ]:


def gene_num(cadd,class_f,filter_f,string):
    """
    Gives the number of genes captured in the final output for each filters
    """
    print("cadd res",string,len(cadd_f),"\n",
      "class res",string,len(class_f),
      "\n","filter res",string,len(filter_f))


# In[ ]:


removelist  = pd.date_range(start="1776-09-07",end="1901-01-01")
pats        = merger(files_pat)
pats        = pats[pats['analysis_tag'].notna()]
pats        = pats.query('cadd_raw != -1.000000 or gnomad_exomes_af != -1.000000 or gnomad_genomes_af != -1.000000')
pats        = pats[~pats['dob'].isin(removelist)]


# In[ ]:


controls = merger(files_ctrl)
controls = controls.query('country_name !="United Arab Emirates" and country_name !="Egypt" and country_name !="Canada"')
controls = controls[controls['analysis_tag'].notna()]
controls     = controls.query('cadd_raw != -1.000000 or gnomad_exomes_af != -1.000000 or gnomad_genomes_af != -1.000000 ')
controls     = controls[~controls['dob'].isin(removelist)]


# In[ ]:


class_f    = list(controls.query('filter_type=="class"')['gene_name'].drop_duplicates())
cadd_f     = list(controls.query('filter_type=="cadd"')['gene_name'].drop_duplicates())
filter_f   = list(controls.query('filter_type=="filter"')['gene_name'].drop_duplicates())
class_pat  = list(pats.query('filter_type=="class"')['gene_name'].drop_duplicates())
cadd_pat   = list(pats.query('filter_type=="cadd"')['gene_name'].drop_duplicates())
filter_pat = list(pats.query('filter_type=="filter"')['gene_name'].drop_duplicates())


# In[ ]:


controls_genes = gene_num(cadd_f,class_f,filter_f,"control")
patient_genes = gene_num(cadd_pat,class_pat,filter_pat,"PD_pats")


# In[ ]:


gene_df_pat  = get_unique_patients(files_pat,"pd_with_var")
gene_df_ctrl = get_unique_patients(files_ctrl,"ctrl_with_var")
class_df,filter_df,cadd_df = get_final_df(gene_df_pat,total_pats_pd,'pd_with_var','pd_without_var')
pd_final                   = pd.concat([filter_df,class_df,cadd_df],axis=0)
class_df_ctrl,filter_df_ctrl,cadd_df_ctrl = get_final_df(gene_df_ctrl,total_control_pd,'ctrl_with_var','ctrl_without_var')
control_final                            = pd.concat([filter_df_ctrl,class_df_ctrl,cadd_df_ctrl],axis=0)


# In[ ]:


contrl_pd_final        = pd.merge(control_final,pd_final, on=['gene_name','filter_type'])
contrl_pd_final["pval"]= contrl_pd_final.apply(lambda r: fisher_test(r.pd_with_var,r.pd_without_var,r.ctrl_with_var,r.ctrl_without_var),axis=1)
contrl_pd_final["OR"]  = contrl_pd_final.apply(lambda r: odds_ratio(r.pd_with_var,r.pd_without_var,r.ctrl_with_var,r.ctrl_without_var),axis=1)


# In[ ]:


TM          = time.localtime()
if control=="restricted":
    outfile     = 'FE_OR_{}.csv'.format(time.strftime('%b-%d-%Y_%H%M', TM))
    ctrl_genes  = 'control_genes_{}.csv'.format(time.strftime('%b-%d-%Y_%H%M', TM))
    pats_genes  = 'pats_genes_{}.csv'.format(time.strftime('%b-%d-%Y_%H%M', TM))
    controls_genes.to_csv(ctrl_genes,sep=',',index=False)
    patient_genes.to_csv(ctrl_genes,sep=',',index=False)
    contrl_pd_final.to_csv(outfile,sep=',',index=False)

else:
    outfile     = 'FE_OR_clear_{}.csv'.format(time.strftime('%b-%d-%Y_%H%M', TM))
    ctrl_genes  = 'control_genes_clear_{}.csv'.format(time.strftime('%b-%d-%Y_%H%M', TM))
    pats_genes  = 'pats_genes_cear_{}.csv'.format(time.strftime('%b-%d-%Y_%H%M', TM))
    controls_genes.to_csv(ctrl_genes,sep=',',index=False)
    patient_genes.to_csv(ctrl_genes,sep=',',index=False)
    contrl_pd_final.to_csv(outfile,sep=',',index=False)

