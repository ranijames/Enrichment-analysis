#!/usr/bin/env python
# coding: utf-8

import sys
import yaml 
import os
import pandas as pd
import xlsxwriter
from datetime import datetime, date
from connection import get_connection
import argparse   
from utility_filters import *

__author__ = 'Alva James'

path2config ='/data/ajames/genetic_modifiers/scripts/pipeline/jupyter-notebooks/config.yml'
with open(path2config, 'r') as file:
     settings = yaml.safe_load(file)
cursor, conn = get_connection(settings)

parser = argparse.ArgumentParser()

parser.add_argument('input',
                        help='input file with variants from unidb'
                        )
parser.add_argument('-genename',
                        help='mode, can be patientid, orderid, familyid'
                        )

parser.add_argument('-mode',
                        help='mode filter')
parser.add_argument('output',
                        help='output path')
args = parser.parse_args()

if len(sys.argv) < 3:
    parser.print_help()
    sys.exit(0)


genename      = args.genename
mode          = args.mode
outpath       = args.output
gene_varaints = args.input

ropad                  = "2021 week 27_ROPAD order IDs_ank.xlsx"
ropad_df               = pd.read_excel(ropad,na_values = "Missing", 
                                       sheet_name='ROPAD',engine='openpyxl') 
ropad_WGS              = ropad_df["ID order step 3"]
ropad_WGS              = pd.DataFrame(ropad_WGS.dropna())
ropad_WGS.rename(columns={'ID order step 3':'orderid'},inplace=True)
ropad_WGS['orderid']  = ropad_WGS['orderid'].astype(str)
ropad_WGS['orderid']  = ropad_WGS['orderid'].str.strip()
orderid_wgs             = list(ropad_WGS['orderid'])
orderid_wgs_ropad       = ",".join(map(lambda x: "'"+x+"'", map(str, orderid_wgs))) 


def all_pats_test(var,order):
    query = """
         SELECT
            v.id AS varid,
            v.chrom AS chrom,
            v.vcf_pos AS vcf_pos,
            v.vcf_ref AS vcf_ref,
            v.vcf_alt AS vcf_alt,
            group_concat(distinct term.term) as HPO_terms,
            group_concat(distinct term.name) as HPO_names,
            group_concat(distinct pp.patientid,"-",pp.person_status,"-",pp.affected_status) as family_label,
            familyid as familyid,
            p.patientid as patientid,
            orderid as orderid,
            analysis_tag as analysis_tag,
            zygosity as zygosity,
            variant_caller_info as variant_caller_info,
            coverage as coverage,
            frequency as frequency,
            d.quality as quality,
            p.person_status as person_status,
            p.affected_status as affected_status,
            p.first_name as first_name,
            p.surname as surname,
            p.gender as gender,
            p.dob as dob,
            c.name as country_name,
            geo.name as region_name,
            is_consanguineous as is_consanguineous,
            g.gene_name as gene_name,
            t.tx_name as tx_name,
            ta.*,
            va.*,
            vc.*,
            ga.*,
            group_concat(
                concat_ws(':', ifnull(g.gene_name,'.'), ifnull(t.tx_name,'.'), ifnull(ta.hgvsc,'.'), ifnull(ta.hgvsp,'.'))
                SEPARATOR '|'
            ) as All_Transcript_Annotations
        FROM variant AS v
            LEFT JOIN variant_annotation va on v.id=va.variant_id and va.status='active'
            LEFT JOIN variant_classification vc on v.id=vc.variant_id
            LEFT JOIN transcript_annotation ta on v.id=ta.variant_id and ta.status='active'
            LEFT JOIN transcript t on t.id=ta.transcript_id and t.status='active'
            LEFT JOIN gene g on g.id=t.gene_id and g.status='active'
            LEFT JOIN gene_annotation ga on g.id=ga.gene_id and ga.status='active'
            JOIN detection d on v.id=d.variant_id
            JOIN variantcaller_genotype vcg on vcg.id=d.variantcaller_genotype_id
            JOIN sample s on s.id=d.sample_id and s.status='Normal'
            LEFT JOIN analysistype `at` on `at`.id=s.analysistype_id
            JOIN `order` o on o.id=s.order_id
            JOIN patient p on p.id=o.patient_id and p.status='active'
            LEFT JOIN country c on c.id=p.country_id
            LEFT JOIN georegion geo on geo.id=c.georegion_id
            JOIN family f on p.family_id=f.id
            JOIN patient pp on pp.family_id=f.id and pp.status='active'
            LEFT JOIN patient2term pt on pt.patient_id=p.id
            LEFT JOIN term on pt.term_id=term.id
        WHERE v.id in ({0}) and orderid in ({1}) 
        GROUP BY v.id,s.id
        HAVING 1
        """.format(var,order)
    return query

# WHERE d.coverage >=20 and d.frequency >=20 and d.quality >=220 and `at`.analysis_tag in ('ILLUDRAGEN_WGS', 'ILLUWGS') and `at`.analysis_tag not like ('%%_SOMATIC') or ('%%_PICOP_TUMOR')
# 

print(genename)

if gene_varaints.endswith('xls') or gene_varaints.endswith('csv') :
    df                      =  pd.read_csv(gene_varaints, sep='\t', header=0)
    df                      =  df[df['gene_name'] == genename]
    #df                      = df[pd.to_numeric(df.gnomad_genomes_af, errors='coerce').notnull()]
else:
    df       = pd.read_excel(gene_varaints,sheet_name ='Z22HC1A_all_variants',engine='openpyxl')
    df       =  df[df['gene_name'] == genename]
    #df       = df[pd.to_numeric(df.gnomad_genomes_af, errors='coerce').notnull()]

var_1    = list(general_filters.filter4_gemod_exac_class_cadd(df).varid.drop_duplicates())
var_2    = list(general_filters.filter4_class(df).varid.drop_duplicates())
var_3    = list(general_filters.filter4_gemod_exac_class(df).varid.drop_duplicates())
print("for class:",len(var_2))

if mode =='cadd' and len(var_1) >=1:
    my_variants     = ",".join(map(lambda x: "'"+x+"'", map(str, var_1)))
    pd_pats         = pd.read_sql(all_pats_test(my_variants,orderid_wgs_ropad), conn)
elif mode =='class' and len(var_2) >=1:
    my_variants     = ",".join(map(lambda x: "'"+x+"'", map(str, var_2)))
    pd_pats = pd.read_sql(all_pats_test(my_variants,orderid_wgs_ropad), conn)
elif mode =='filter' and len(var_3) >=1:
    my_variants     = ",".join(map(lambda x: "'"+x+"'", map(str, var_3)))
    pd_pats = pd.read_sql(all_pats_test(my_variants,orderid_wgs_ropad), conn)
else:
    print(genename,"has got no variants after", mode ,"filtering!") 
    sys.exit(0)   

outfile       = 'pats_pd_{}_{}.csv'.format(genename,mode)
pd_pats_final = general_filters.filter4_QC_frq_excl_somatic(pd_pats)
print("before filter:",pd_pats.shape,"after filter :",pd_pats_final.shape)
pd_pats_final.to_csv(outpath+outfile,sep='\t',index=False)
print ('Successfully saved the results as', outpath+outfile, "with",pd_pats_final.patientid.drop_duplicates().shape[0],"patients")
