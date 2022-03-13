#!/usr/bin/env python
import pandas as pd
from datetime import datetime, date

__author__ = 'Alva James'
class general_filters:
    def __init__(self, df):
        df = self.df
    """
   Filter for the variants for PD patients and control
   """
    def filter4_gemod_exac_class_cadd(df):
        filtered_df = df[(df["exac_afr"] <= 0.0001)]
        filtered_df = filtered_df[filtered_df["gnomad_genomes_af"] <= 0.0001]
        filtered_df = filtered_df[filtered_df["gnomad_exomes_af"] <= 0.0001]
        filtered_df = filtered_df[filtered_df["cadd_raw"] >= 4]
        gemod_exac_class_cadd = filtered_df[
            (filtered_df["effectclass"] == "HIGH")
            | (filtered_df["effectclass"] == "MODERATE")
        ]
        return gemod_exac_class_cadd

    def filter4_gemod_exac_class(df):
        filtered_df = df[(df["exac_afr"] <= 0.0001)]
        filtered_df = filtered_df[filtered_df["gnomad_genomes_af"] <= 0.0001]
        filtered_df = filtered_df[filtered_df["gnomad_exomes_af"] <= 0.0001]
        gemod_exac_class_fil = filtered_df[
            (filtered_df["effectclass"] == "HIGH")
            | (filtered_df["effectclass"] == "MODERATE")
        ]
        return gemod_exac_class_fil

    def filter4_class(df):
        class_filtered = df[
            (df["effectclass"] == "HIGH") | (df["effectclass"] == "MODERATE")
        ]
        return class_filtered
    

    def filter4_control_geo_tags(QC_somatic_filtered):
        var_pat_f2_unaffected_healthy = QC_somatic_filtered[
            (QC_somatic_filtered["analysis_tag"]   == "ILLUWES")
            | (QC_somatic_filtered["analysis_tag"] == "ILLUWESAGI")
            | (QC_somatic_filtered["analysis_tag"] == "ILLUWESTWIST")
            | (QC_somatic_filtered["analysis_tag"] == "ILLUWESTWIST_V2")
            | (QC_somatic_filtered["analysis_tag"] == "ILLUWGS")
            | (QC_somatic_filtered["analysis_tag"] == "ILLUDRAGEN_WGS")
        ]
        var_pat_f2_unaffected_healthy_eur_na = var_pat_f2_unaffected_healthy[
            (var_pat_f2_unaffected_healthy["region_name"]   == "North America")
            | (var_pat_f2_unaffected_healthy["region_name"] == "Europe")
        ]
        var_pat_f2_unaffected_healthy_bra_isra = var_pat_f2_unaffected_healthy[
            (var_pat_f2_unaffected_healthy["country_name"]   == "Israel")
            | (var_pat_f2_unaffected_healthy["country_name"] == "Brazil")
        ]
        control_18_geo = pd.concat(
            [
                var_pat_f2_unaffected_healthy_bra_isra,
                var_pat_f2_unaffected_healthy_eur_na,
            ]
        )
        control_18_geo = control_18_geo.query(
            'country_name !="United Arab Emirates" and country_name !="Egypt" and country_name !="Canada"'
        )
        control_18_geo = control_18_geo.query(
            'affected_status =="unaffected"  or affected_status == "unaffected (healthy)"'
        )
        return control_18_geo

    def filter4_QC_frq_excl_somatic(var_pat):
        var_pat_f1 = var_pat[(var_pat["coverage"] >= 20)]
        var_pat_f1 = var_pat_f1[var_pat_f1["frequency"] >= 20]
        var_pat_f1 = var_pat_f1[var_pat_f1["quality"] >= 220]
        var_pat_f1 = var_pat_f1[
            var_pat_f1["analysis_tag"].str.contains("_SOMATIC|_PICOP_TUMOR") == False
        ]
        return var_pat_f1
