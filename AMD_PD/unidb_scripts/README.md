# Existing workflow**

The existing workflow is depicted in Fig1. The current pipeline directly fetches data from UniDB and operates on disease cohort (ROPAD PD cohort ca. 3300 patients). Two control groups are used;

restricted cohort  (based on age, affected status, country, exclusion of analysis types and including certain analysis types) and unrestricted group which includes all uni DB patients except disease cohort or somatic (Q: Does it mean somatic symptoms? A: No, the somatic term here refers to cancer patients who carry somatic mutations). For a given set of genes, in this case, study 27 PD-Gaucher Disease-related GWAS genes (Q: Can we link the accession IDs of these genes and their GWAS scores here?), these genes and their associated variants are filtered/selected. On these selected genes, the variant filter is applied; these filtering criteria are H; H, AF frequency; H, M, AF, Cadd >= 4. For each cohort, the number of unique patients with variants per GWAS genes is reported. The existing WF performs Fischer's exact test/hypergeometric test and multiple hypothesis tests which enables us to estimate the statistical significance of the variants between control and disease cohorts and enrichment/depletion level. 

# Existing pipeline 
Data set used for developing the above pipeline

The existing pipeline is designed to identify the enriched genes within the ROPAD Pd patients within the unidb database versus two types of controls. The first type of control is restricted control. The second type of control is clear control.

The restricted control is classified based on the following criteria:
- Age above 18
- Affected status: unaffected or unaffected (healthy)
- Only include Country/Georegion: EU, Israel, Brasil, USA
-  Exclude analyses types RODPAD
- Exclude all related panel and WGS used in disease cohort, and SOMATIC (SOM, PiCOP) 
- Include ILLUWES, ILLUWESAGI, ILLUWESTWIST, ILLUWESTWIST_V2,  ILLUWGS, ILLUDRAGEN_WGS
- The clear control conditions are:
- All uniDB patients, except disease cohort or somatic
- The disease group 
- ROPAD PD patients




Besides the above conditions for the control and patients. The gene's variants have to undergo certain filters.

The variant level filters

The variant level filtered is termed as, "CADD", "CLASS", and "FILTER".  The conditions applied to generate the output for each filter type is described in the following section.

CADD filter	FILTER	CLASS


exac_afr <= 0.0001

gnomad_genomes_af <= 0.0001

gnomad_exomes_af <= 0.0001

cadd_raw >= 4

effectclass  == "HIGH"| "effectclass" =="MODERATE"

	

exac_afr <= 0.0001

gnomad_genomes_af <= 0.0001

gnomad_exomes_af <= 0.0001

effectclass == "HIGH"| "effectclass" =="MODERATE"

	

effectclass == "HIGH"| "effectclass" =="MODERATE"

Patient-level  or control individual filters

The above filters are used to fetch each gene's variants from UNIDB. Then after that for the extracted variants we then need to extract the controls and ROPAD patients. After that, there are some patient or control level filters applied to fetch the patients and both controls from the unidb.  The patient-level/control individual level  filters are,

On patient or control individual


coverage >= 20

frequency >= 20

quality >= 220




The data which we have in-house

Total restricted control	All control      	ROPAD PD patients      
8971         	   22918	3272 
Modules in the pipeline

The pipeline is developed mainly in python, MSQL queries with helper scripts written in shell.  

Within the GitHub page, we have the following python tool to derive all the restricted control from the unidb

The filters are defined as functions within class general filters:

The pipeline is designed to do the following tasks,

Fetch all variants
Filter the variants for 3 filter types
Fetch all ROPAD patients and healthy control individuals for the filtered variants 
 Apply the patient level filters
Apply the age filter
Save the results in the output matrix for three filter types

The scripts used to accomplish the above tasks are described in the below section

 Utility_filters.py

Within the Utility_filters.py there are different functions for each filter type.

1. filter4_gemod_exac_class_cadd : CADD filter

2. filter4_gemod_exac_class : Filter

3. filter4_class: Class filter

4. filter4_control_geo_tags : For country and georegion

5. filter4_QC_frq_excl_somatic : for somatic plus patient level filters for coverage, frequency and quality

       2.  Basic export.py

     Basic export helps to fetch all variants within UNIDB for a given list of genes as a new line separated input file. The basic export outputs the variants for the given list of genes with a lot of other information about the variants. The output is sufficient to apply the CADD, CLASS, and FILTER filter types. The filtered variant table is then used for further fetching of healthy individuals and ROPAD patients from UNIDB

    4.  controls_restricted_Pd.py

The controls_restricted_Pd.py  is a small python wrapper tool that first establishes a connection to the unidb database  via a function called get_connection. After that, it calls the variant filter functions from the class utility_filters.py and apply the user mentioned filter on the input variants of the gene of interest. 

The inputs for the tool are,

Defined as parameter "–mode" one of the three filters for variants
The input file with genes variants
output folder's path
The name of gene of interest as parameter "-genename"

Where "-mode" is one of the variant filter modes (CLASS,  CADD, or FILTER). The "-genename"  is the gene_symbol of the gene of interest. output path is the path where you need your output and the input is the variants of the gene of interest, this is the output from the basic_export.py script

The main part of the tool is mYSQL query to fetch the individuals with the given list of variants of the gene of interest. The result or the output contains the restricted control for the enrichment analysis on gene level. The SQl query contains almost all filters defined in the section "The restricted control is classified based on the following criteria:" above. The database Unidb has got only data of birth(dob) for each records. Therefore, a dob to Age converter function is also defined within the tool. The function is named as, from_dob_to_age within the script. From the age then, >18 is filtered. Besides  that, all oddly indicated ages are also removed (<120 years for example).

The tool throws an error or exit when there are no records for a given gene after filtering for its variants. The output is always for a given gene, given mode of filter a CSV file with patientid and all corresponding information. The timestamp is attached to the output filename for simplification, and for future reference.

 
# Final Fisher extact test or enrichment anaysis
```
python FE_test.py /shares/archive/develdata/ajames/GD_PD/results/lancet/pats/ /shares/archive/develdata/ajames/GD_PD/results/lancet/control/
```
