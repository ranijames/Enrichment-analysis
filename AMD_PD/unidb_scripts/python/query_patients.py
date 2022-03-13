#!/usr/bin/env python
import pymysql
__author__ = 'Alva James'
def query4_patients(my_variants):
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
            ) as All_Transcript_Annotations
        FROM variant AS v
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
        WHERE v.id in ({0})
        GROUP BY v.id,s.id
        HAVING 1
        """.format(my_variants)
    return query
