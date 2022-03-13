#!/usr/bin/python
# -*- coding: utf-8 -*-

# In[17]:
__author__ = 'Alva James'

# Comment Robert: Can you write a small comment for each query what it is doing and what is its purpose. ONe/two sentences is enought
class control_queries:

    def __init__(self, my_variants, orderid):
        df = self.my_variants
        orderid = self.orderid

    def query4_park_control_r(my_variants, orderid):
        queryps = \
            """
             SELECT DISTINCT p.patientid, v.id,o.orderid,d.coverage,
             d.frequency,d.quality,`at`.analysis_tag,
    c.`name` as country_name, geo.`name` as region_name FROM
    unidb.variant AS v JOIN unidb.detection AS d ON d.variant_id = v.id
    JOIN unidb.variant_classification as vc on vc.variant_id = d.variant_id
    JOIN unidb.sample AS s ON d.sample_id = s.id
    JOIN unidb.`order` AS o ON o.id = s.order_id
    JOIN unidb.patient p on p.id = o.patient_id
    JOIN unidb.family f on f.id = p.family_id
    JOIN unidb.variant_annotation AS va ON va.variant_id = v.id
    JOIN unidb.transcript_annotation as ta on ta.variant_id=v.id
    JOIN unidb.country c on c.id=p.country_id
    JOIN unidb.georegion geo on geo.id=c.georegion_id
    JOIN unidb.analysistype AS `at` on `at`.id=s.analysistype_id where 
    p.affected_status IN ('unaffected', 'unaffected (healthy)') 
    AND DATE_ADD(p.dob, INTERVAL 18 YEAR) < NOW() 
    and v.id in ({0}) and o.orderid not in ({1}) 
    and `at`.analysis_tag in ('ILLUWES','ILLUDRAGEN_WGS','ILLUWGS','ILLUWESTWIST_V2','ILLUWESTWIST','ILLUWESAGI') 
    and d.frequency >=20 and d.quality >=220 and d.coverage >=20
    and p.`status`='active' """.format(my_variants,
                orderid)
        return queryps

    # Comment Robert: 'illudragne_wgs' should be probably 'illudragen_wgs'
    def quey4_all_control_pd(my_variants):
        query = \
            """
            SELECT *
            FROM view_family_to_sample AS vfts
        JOIN patient AS p ON p.id = vfts.patient_id
        JOIN view_family_to_variant AS vftv ON vftv.patientid = vfts.patientid
        WHERE vftv.variant_id in ({0}) and vfts.analysis_tag IN ('ILLUWES','ILLUWESAGI','illuwestwist','illuwestwist_v2','illuwgs','illudragne_wgs')
        AND p.affected_status IN ('unaffected', 'unaffected (healthy)')
        AND vfts.georegion_name IN ('North America' , 'Europe')
        OR vfts.country_name IN ('Israel', 'Brazil')
        AND p.status = ('active') AND DATE_ADD(vfts.dob, INTERVAL 18 YEAR) < NOW()
        """.format(my_variants)
        return query

    def quey4_control_exl_ropad(my_variants):
        query = \
            """
 SELECT v.id,p.dob,c.`name` as country_name,geo.`name` as region_name,p.patientid,`at`.analysis_tag,p.affected_status,d.coverage,d.frequency,d.quality
FROM unidb.variant AS v JOIN unidb.detection AS d ON d.variant_id = v.id
JOIN unidb.variant_classification as vc on vc.variant_id = d.variant_id
JOIN unidb.sample AS s ON d.sample_id = s.id
JOIN unidb.`order` AS o ON o.id = s.order_id
JOIN unidb.patient p on p.id = o.patient_id
JOIN unidb.analysistype AS `at` on `at`.id=s.analysistype_id
JOIN unidb.variant_annotation as va on va.variant_id=v.id
JOIN country c on c.id=p.country_id
JOIN georegion geo on geo.id=c.georegion_id
WHERE v.id in ({0}) and va.`status`=('active') 
AND p.affected_status IN ('unaffected','unaffected (healthy)')
AND geo.`name` IN ('North America', 'Europe') OR c.`name` IN ('Israel', 'Brazil')
AND DATE_ADD(p.dob, INTERVAL 18 YEAR) < NOW();""".format(my_variants)
        return query
