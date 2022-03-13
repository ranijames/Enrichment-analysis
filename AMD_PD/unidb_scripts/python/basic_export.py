#!/usr/bin/env python2

"""
Usage: program  [options] [--help]

Options:
    --help                      displays Usage
    --mode=[variants|patients]  defines output type: "variant export" (default) or "patient export"
    --bedfile=<bedfile>         bedfile for variant export
    --genelist=<genelist>       comma-separated gene list for variant/patient export (may be given as file or string)
    --varlist=<varids>          comma-separated list of variant ids for patient export (may be given as file or string)
    --filters=<filters>         file format: "column"<tab>"operator"<tab>"value" (may be given as file or string)
    --sql-filter=<sqlfilter>    allows to extend the HAVING statement of the sql query directly (may be given as file or string)
    --outfile=<outfile>         defines output file

Note:
    filer operators include "<=", ">=", "=", "!=", "like", "not like", "in", "not in"

"""
# from __future__ import print_function
from docopt import docopt
import os
import common
import re
import sys
import traceback

TMPTIME = re.sub(r"\W+", "", "{}".format(common.getTimestamp()))

LOGFILE = '/dev/null'
logger = common.getLogger(LOGFILE)

cfg = common.getConfig()

TMP_DIR = cfg["paths"]["devel_tmp"]

cur = common.getUnidbCur()

MAX_ROW_LIMIT = 10000000

MYHEADER = ["varid", "chrom", "vcf_pos", "vcf_ref", "vcf_alt", "HPO_terms", "HPO_names", "family_label", "familyid",
    "patientid", "orderid", "runname", "barcodename", "analysis_tag", "zygosity", "vc_info", "coverage", "frequency",
    "quality", "person_status", "affected_status", "first_name", "last_name", "gender", "dob", "country_name",
    "region_name", "is_consanguineous", "HomPatCount", "HetPatCount", "HomAffCount", "HetAffCount", "FamCount",
    "TotalPatCount", "TotalAffCount", "TotalFamCount", "AF_Pat", "AF_Aff", "AF_Fam", "analysistypelist",
    "HomAssayCount", "HetAssayCount", "TotalAssayCount", "AF_Assay",
    "HomControlCount", "HetControlCount", "TotalControlCount", "AF_Control", "gene_name", "tx_name", "hgvsc", "hgvsp",
    "effect", "effectclass", "class", "comment", "AllTranscriptAnnotations", "sift_pred", "rsid", "provean_pred",
    "polyphen2_hvar_pred", "polyphen2_hdiv_pred", "mutationtaster_pred", "lrt_pred", "fathmm_pred", "caddgt20",
    "cadd_raw", "ada_score", "pgmd_accession", "hgmd_variantclass", "hgmd_accession", "cosmic70", "clinvar",
    "centomd_patho", "exac_all", "PopFreqMax"]
MYOPERATORS = ["<=",">=","=","!=","like","not like","in","not in"]

def getVariantsQuery(my_variants, filters=[], sql_filter=""):

    filter_str = sql_filter
    for col,op,val in filters:
        if op == "<=":
            filter_str += " AND ({} IS NULL OR {} {} {})".format(col,col,op,val)
        elif op == ">=":
            filter_str += " AND ({} IS NULL OR {} {} {})".format(col,col,op,val)
        elif op == "=":
            filter_str += " AND ({} {} '{}')".format(col,op,val)
        elif op == "!=":
            filter_str += " AND ({} {} '{}')".format(col,op,val)
        elif op == "like":
            filter_str += " AND ({} {} '%{}%')".format(col,op,val)
        elif op == "not like":
            filter_str += " AND ({} {} '%{}%')".format(col,op,val)
        elif op == "in":
            filter_str += " AND ({} {} ('{}'))".format(col,op,val.replace(",","','"))
        elif op == "not in":
            filter_str += " AND ({} {} ('{}'))".format(col,op,val.replace(",","','"))

    logger.info("Filter string is: {}".format(filter_str))

    ### list all variants ###
    query = """
        SELECT
            v.id AS varid,
            v.chrom AS chrom,
            v.vcf_pos AS vcf_pos,
            v.vcf_ref AS vcf_ref,
            v.vcf_alt AS vcf_alt,
            vcc.HomPatCount,
            vcc.HetPatCount,
            vcc.HomAffCount,
            vcc.HetAffCount,
            vcc.FamCount,
            vcc.TotalPatCount,
            vcc.TotalAffCount,
            vcc.TotalFamCount,
            if(vcc.AF_Pat between 0 AND 1, vcc.AF_Pat, NULL) AS AF_Pat,
            if(vcc.AF_Aff between 0 AND 1, vcc.AF_Aff, NULL) AS AF_Aff,
            if(vcc.AF_Fam between 0 AND 1, vcc.AF_Fam, NULL) AS AF_Fam,
            replace(vcc.analysistypelist,',',';') AS analysistypelist,
            vcc2.HomUnaffCount AS HomControlCount,
            vcc2.HetUnaffCount AS HetControlCount,
            vcc2.TotalUnaffCount AS TotalControlCount,
            if(vcc2.AF_healthy between 0 AND 1, vcc2.AF_healthy, NULL) AS AF_Control,
            g.gene_name as gene_name,
            t.tx_name as tx_name,
            ta.*,
            va.*,
            vc.*,
            ga.*,
            group_concat(
                concat_ws(':', ifnull(g.gene_name,'.'), ifnull(t.tx_name,'.'), ifnull(ta.hgvsc,'.'), ifnull(ta.hgvsp,'.'))
                SEPARATOR '|'
            ) as `AllTranscriptAnnotations`
        FROM {} AS v
            LEFT JOIN view_variant_count_static vcc ON vcc.variant_id=v.id
            LEFT JOIN view_variant_count_healthy vcc2 ON vcc2.variant_id=v.id
            LEFT JOIN variant_annotation va on v.id=va.variant_id and va.status='active'
            LEFT JOIN variant_classification vc on v.id=vc.variant_id
            LEFT JOIN transcript_annotation ta on v.id=ta.variant_id and ta.status='active'
            LEFT JOIN transcript t on t.id=ta.transcript_id and t.status='active'
            LEFT JOIN gene g on g.id=t.gene_id and g.status='active'
            LEFT JOIN gene_annotation ga on g.id=ga.gene_id and ga.status='active'
        GROUP BY v.id
        HAVING 1 {}
        """.format(my_variants,filter_str)
    return query

def getPatientsQuery(my_variants, filters=[], sql_filter=""):

    filter_str = sql_filter
    for col,op,val in filters:
        if op == "<=":
            filter_str += " AND ({} IS NULL OR {} {} {})".format(col,col,op,val)
        elif op == ">=":
            filter_str += " AND ({} IS NULL OR {} {} {})".format(col,col,op,val)
        elif op == "=":
            filter_str += " AND ({} {} '{}')".format(col,op,val)
        elif op == "!=":
            filter_str += " AND ({} {} '{}')".format(col,op,val)
        elif op == "like":
            filter_str += " AND ({} {} '%{}%')".format(col,op,val)
        elif op == "not like":
            filter_str += " AND ({} {} '%{}%')".format(col,op,val)
        elif op == "in":
            filter_str += " AND ({} {} ('{}'))".format(col,op,val.replace(",","','"))
        elif op == "not in":
            filter_str += " AND ({} {} ('{}'))".format(col,op,val.replace(",","','"))

    logger.info("Filter string is: {}".format(filter_str))

    ### list all detections ###
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
            r.name as runname,
            b.name as barcodename,
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
            vcc.HomPatCount,
            vcc.HetPatCount,
            vcc.HomAffCount,
            vcc.HetAffCount,
            vcc.FamCount,
            vcc.TotalPatCount,
            vcc.TotalAffCount,
            vcc.TotalFamCount,
            if(vcc.AF_Pat between 0 AND 1, vcc.AF_Pat, NULL) AS AF_Pat,
            if(vcc.AF_Aff between 0 AND 1, vcc.AF_Aff, NULL) AS AF_Aff,
            if(vcc.AF_Fam between 0 AND 1, vcc.AF_Fam, NULL) AS AF_Fam,
            replace(vcc.analysistypelist,',',';') AS analysistypelist,
            vcc2.HomPatCount AS HomAssayCount,
            vcc2.HetPatCount AS HetAssayCount,
            vcc2.TotalPatCount AS TotalAssayCount,
            if(vcc2.AF_Pat between 0 AND 1, vcc2.AF_Pat, NULL) AS AF_Assay,
            vcc3.HomUnaffCount AS HomControlCount,
            vcc3.HetUnaffCount AS HetControlCount,
            vcc3.TotalUnaffCount AS TotalControlCount,
            if(vcc3.AF_healthy between 0 AND 1, vcc3.AF_healthy, NULL) AS AF_Control,
            g.gene_name as gene_name,
            t.tx_name as tx_name,
            ta.*,
            va.*,
            vc.*,
            ga.*,
            group_concat(
                concat_ws(':', ifnull(g.gene_name,'.'), ifnull(t.tx_name,'.'), ifnull(ta.hgvsc,'.'), ifnull(ta.hgvsp,'.'))
                SEPARATOR '|'
            ) as `AllTranscriptAnnotations`
        FROM {} AS v
            LEFT JOIN view_variant_count_static vcc ON vcc.variant_id=v.id
            LEFT JOIN view_variant_count_healthy vcc3 ON vcc3.variant_id=v.id
            LEFT JOIN variant_annotation va on v.id=va.variant_id and va.status='active'
            LEFT JOIN variant_classification vc on v.id=vc.variant_id
            LEFT JOIN transcript_annotation ta on v.id=ta.variant_id and ta.status='active'
            LEFT JOIN transcript t on t.id=ta.transcript_id and t.status='active'
            LEFT JOIN gene g on g.id=t.gene_id and g.status='active'
            LEFT JOIN gene_annotation ga on g.id=ga.gene_id and ga.status='active'
            JOIN detection d on v.id=d.variant_id
            JOIN variantcaller_genotype vcg on vcg.id=d.variantcaller_genotype_id
            JOIN sample s on s.id=d.sample_id and s.status='Normal'
            LEFT JOIN run r on r.id=run_id
            LEFT JOIN barcode b on b.id=barcode_id
            LEFT JOIN analysistype `at` on `at`.id=s.analysistype_id
            LEFT JOIN view_variant_count_assay vcc2 ON vcc2.variant_id=v.id and vcc2.analysistype_id=`at`.id
            JOIN `order` o on o.id=s.order_id
            JOIN patient p on p.id=o.patient_id
            LEFT JOIN country c on c.id=p.country_id
            LEFT JOIN georegion geo on geo.id=c.georegion_id
            JOIN family f on p.family_id=f.id
            JOIN patient pp on pp.family_id=f.id and pp.status='active'
            LEFT JOIN patient2term pt on pt.patient_id=p.id
            LEFT JOIN term on pt.term_id=term.id
        GROUP BY v.id,s.id
        HAVING 1 {}
        """.format(my_variants,filter_str)

    return query

def parseFiltersFile(filters_raw):

    filters = []

    in_string = filters_raw.replace(';','\n').replace("::",'\t')
    if os.path.isfile(filters_raw):
        with open(filters_raw,'r') as f:
            in_string = f.read().replace('\r','')

    for line in in_string.split('\n'):
        if line.strip()=="":
            continue
        if line.count('\t')!=2:
            logger.error("Column count is not %d (but %d)! Filter entry will be skipped ... "%(2,line.count('\t')))
            continue
        column,operator,value = [x.strip() for x in line.split('\t')]
        if column not in MYHEADER:
            logger.error("Column %s not in MYHEADER! Filter entry will be skipped ... "%column)
            continue
        if operator not in MYOPERATORS:
            logger.error("Operator %s not in MYOPERATORS! Filter entry will be skipped ... "%operator)
            continue
        filters.append([column,operator,value])

    return filters


def main(argv=None):

    args = docopt(__doc__, None)

    if args['--mode'] is None:
        mode = "variants"
        logger.info("using default mode %s"%mode)
    elif args['--mode'] in ["variant","variants"]:
        mode= "variants"
        logger.info("mode was specified as input: %s"%mode)
    elif args['--mode'] in ["patient","patients"]:
        mode="patients"
        logger.info("mode was specified as input: %s"%mode)
    else:
        raise NameError("mode not recognized!")

    gene_bed = []
    gene_list = []
    gene_list_raw = []
    varid_list = []

    if args['--bedfile']:
        for line in open(args['--bedfile']).readlines():
            row = line.split('\t')
            gene_bed.append({'chrom':row[0],'start':row[1],'end':row[2]})
        logger.info("bedfile was specified as input: %s"%args['--bedfile'])
    if args['--genelist']:
        if os.path.isfile(args['--genelist']):
            gene_list_raw = [x.strip() for x in open(args['--genelist']).read().replace(",","\n").split('\n') if x.strip()!=""]
        else:
            gene_list_raw = [x.strip() for x in args['--genelist'].split(',') if x.strip()!=""]
        logger.info("genelist was specified as input: %s"%args['--genelist'])
    if args['--varlist']:
        if os.path.isfile(args['--varlist']):
            varid_list = [x.strip() for x in open(args['--varlist']).read().replace(",", "\n").split('\n') if x.strip()!=""]
        else:
            varid_list = [x.strip() for x in args['--varlist'].split(',') if x.strip()!=""]
        logger.info("varlist was specified as input: %s"%args['--varlist'])

    filters = []
    if args['--filters']:
        filters = parseFiltersFile(args['--filters'])

    sql_filter = ""
    if args['--sql-filter']:
        if os.path.isfile(args['--sql-filter']):
            sql_filter = " AND ({})".format(open(args['--sql-filter']).read())
        else:
            sql_filter = " AND ({})".format(args['--sql-filter'])

    if args['--outfile']:
        outfile = args['--outfile']
        logger.info("outfile was specified via cmd: %s" % outfile)
    else:
        outfile = "./unidb_export_{}.xls".format(TMPTIME)
        # outfile = "/data/NGS-SeqData/mweiss/unidb_export_{}.xls".format(TMP_TIMESTAMP)
        logger.info("using default outfile %s" % outfile)

    if len(gene_list_raw)>0:

        cur.execute("SELECT gene,group_concat(name) AS tx_names from dbtx.refSeqGenes_UCSC WHERE gene in ({}) group by gene".format(
            ','.join('"%s"'%(x.strip()) for x in gene_list_raw)
        ))
        for row in cur.fetchall():
            logger.info("Using the following transcripts for '%s': %s"%(row['gene'],row['tx_names']))
            gene_list.append(row['gene'].upper())

        if len(gene_list)!=len(gene_list_raw):
            logger.error("The following genes will be ignored because no transcripts were found: %s"%(', '.join(x for x in gene_list_raw if x.upper() not in gene_list)))

        logger.info("fetching variants for %d genes"%len(gene_list))
        gene_bed = common.getBedRegionsForGeneList_UCSC(gene_list)

    if len(gene_bed)>0:
        for bed in gene_bed:
            cur.execute("SELECT id FROM variant WHERE chrom=%s AND vcf_pos BETWEEN %s AND %s",(bed['chrom'],bed['start'],bed['end']))
            varid_list += [x['id'] for x in cur.fetchall()]
            
    if len(varid_list) > MAX_ROW_LIMIT:

        logger.error("Number of selected variants exceeds the MAX_ROW_LIMIT ({}>{})! Please use a divide-and-conquer approach!".format(len(varid_list),MAX_ROW_LIMIT))
        raise NameError("Number of selected variants exceeds the MAX_ROW_LIMIT ({}>{})! Please use a divide-and-conquer approach!".format(len(varid_list),MAX_ROW_LIMIT))

    elif len(varid_list)>0:

        logger.info("exporting annotations for %d variants" % len(varid_list))
        my_variants = "(select * from variant where id in ({}))".format(','.join(str(x) for x in varid_list))

        if mode=="variants":
            query = getVariantsQuery(my_variants, filters, sql_filter)
        elif mode=="patients":
            logger.info("performing safety check ... ")
            cur.execute("SELECT count(*) N FROM {} AS v join detection d on v.id=d.variant_id".format(my_variants))
            dcount = cur.fetchone()['N']
            if dcount > MAX_ROW_LIMIT:
                logger.error("Number of selected detections exceeds the MAX_ROW_LIMIT ({}>{})! Please use a divide-and-conquer approach!".format(dcount,MAX_ROW_LIMIT))
                raise NameError("Number of selected detections exceeds the MAX_ROW_LIMIT ({}>{})! Please use a divide-and-conquer approach!".format(dcount,MAX_ROW_LIMIT))
            else:
                logger.info("safety check passed ... without additional filters, {} detections will be considered".format(dcount))
            print(my_variants)
            print("filter is al: ",filters)
            print("SQlfilter is al: ",sql_filter)
            query = getPatientsQuery(my_variants, filters, sql_filter)
        else:
            raise NameError("Uknown mode '%s'"%mode)

        try:
            logger.info("running query ... ")
            counter = cur.execute(query)
            logger.info("query returned %d rows"%counter)
        except Exception,e:
            logger.error("query failed: "+query)
            raise

        logger.info("writing into outfile %s"%outfile)
        with open(outfile,'w') as f:
            header = None
            for row in cur.fetchall():
                if header is None:
                    header = [x for x in MYHEADER if x in row.keys()] + [x for x in row.keys() if x not in MYHEADER]
                    f.write('\t'.join(k for k in header) +'\n')
                f.write('\t'.join(str(row[k]).replace('None','-1')[:10000].replace('\n', ' ').replace('\r', ' ') for k in header) +'\n')


if __name__=="__main__":

    LOGFILE = '{}/basic_export_{}.log'.format(cfg['paths']['log_dir'],TMPTIME)
    logger = common.getLogger(LOGFILE)

    task_id = common.register_task(TMPTIME,logger=logger)
    try:
        common.activate_task(task_id,logger=logger)
        main()
        common.update_task(task_id,'finished',logger=logger)
    except:
        common.update_task(task_id,'crashed',logger=logger)
        mail_adresses = cfg['mysql_unidb']['emergency_contacts']
        mail_prefix = "PROBLEM"
        mail_subject = os.path.abspath(sys.argv[0])
        mail_msg = "ERROR: Exception caught\n{}".format(traceback.format_exc())
        if not mail_msg.strip().endswith("SystemExit"):
            logger.error(mail_msg)
            mail_msg += '\n\n' + open(LOGFILE).read()
            if len(mail_adresses)>0:
                logger.warn("This incidence will be reported!")
                common.mail(mail_adresses, mail_msg, mail_prefix, mail_subject)
        raise
