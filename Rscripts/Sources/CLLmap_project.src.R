### Main CLLmap config and functions used for managing data, analysis, plotting and manuscript prep
### Author: Binyamin A. Knisbacher
### Importantly, this script is provided for review and educational purposes only.

# Usage:
# suppressPackageStartupMessages(source('~/github/Rscripts/Sources/CLL1085/CLLmap_project.src.R', proj<-new.env()))

suppressPackageStartupMessages(library(data.table))
GITHUB_DIR = '~/github'
RSCRIPTS_DIR = paste0(GITHUB_DIR, '/', 'Rscripts')
source(paste0(RSCRIPTS_DIR, '/', 'Sources/stats.src.R'))
source(paste0(RSCRIPTS_DIR, '/', 'Sources/RNAseq.src.R'), rnaseq <- new.env())
source('~/github/Rscripts/plotting/comut.utils.src.R', comut.utils <- new.env())
source('~/github/Rscripts/Sources/RNAseq.src.R', rnaseq.utils <- new.env())

ANALYSIS_SETS = c('all', 'mcll', 'ucll')
ANALYSIS_LABELS = list(all='All', mcll='M-CLL', ucll='U-CLL')
ANALYSIS_SETS_N = length(ANALYSIS_SETS)

MASTER_METADTA_TABLE = '/xchip/cga_home/bknisbac/CLL/results_cll1085/tables/cll1085_master_table_primary_dataset.20210702.txt'
RNA_WORKSPACES_SAMPLE_TO_PARTICIPANT_FILE = '/xchip/cga_home/bknisbac/CLL/metadata/terra_workspaces_data/all_samps_to_pats.with_dietrich.20210208.tsv'
KIPPS_COHORT_LONGITUDINAL_MD = '/xchip/cga_home/bknisbac/CLL/metadata/fludarabine_resistance_data_Kipps_and_IBM_2019_withNormals.txt'
#poorly formatted / problematic ID conversion table
CW_SAMPLE_RENAME_FILE = '/xchip/cga_home/bknisbac/CLL/metadata/rename_cw_bad_id_format_conversion_table.txt'

get_longitudinal_args <- function(opt='full'){
  long_args <- new.env()
  if(opt=='full'){
    long_args$run_mode = 'full'
    long_args$PRE_ONLY = F
  } else if (opt=='minimal_pre'){
    long_args$run_mode = 'full'
    long_args$PRE_ONLY = F
  }
}

rename_samples_by_file <- function(dt, cols=NULL, convert_file=CW_SAMPLE_RENAME_FILE, from_col='original_id', to_col='new_id'){
  # get converter
  converter_dt = fread(convert_file)
  # convert specified columns (or any match if cols not specified)
  if(is.null(cols))
    cols = colnames(dt)
  for (col in cols){ #col = 'WES_TUMOR_SAMPLE_ID'
    tmp_arr = converter_dt[['new_id']][match(dt[[col]], converter_dt[['original_id']])]
    na_positions = which(is.na(tmp_arr))
    tmp_arr[na_positions] = dt[[col]][na_positions]
    dt[[col]] = tmp_arr
  }
  return(dt)
}


#### Misc
DATA_SHARING_ROOT = '/xchip/cga_home/bknisbac/CLL/results_cll1085/sharing'

#### Genomic annotations
GENE_BED_FILE = '/xchip/cga_home/bknisbac/resources/gencode19_noChrPrefix_mitoMT.withNewHugo.bed.gz'
GENES_GTF_FILE = '/Volumes/xchip_cga_home/bknisbac/resources/gencode19_noChrPrefix_mitoMT.gtf.gz'

### Project plot dirs
PROJECT_PLOT_DIR = '/xchip/cga_home/bknisbac/CLL/results_cll1085/project_plots'
PROJECT_PLOTS_TIMESTAMP = '20210707'
plot_files = list()
plot_files[['de_volcanoplot']] = paste0(PROJECT_PLOT_DIR, '/', 'volcanoplot_de_with_markers_with_acetylation.',PROJECT_PLOTS_TIMESTAMP,'.pdf')
plot_files[['marker_umaps']] = paste0(PROJECT_PLOT_DIR, '/umap_ec_all_markers.', PROJECT_PLOTS_TIMESTAMP,'.pdf')

RESULT_TABLES_DIR = '/xchip/cga_home/bknisbac/CLL/results_cll1085/tables'
#result_table specifics
N_GENES_PER_FILE_FOR_PREDICTOR_TABLE = paste0(proj$RESULT_TABLES_DIR, '/', 'n_genes_per_file_for_predictor_table.20210707.tsv')

### tri(12) phenocopy analysis
REVISION_ANALYSIS_OUTDIR = '/xchip/cga_home/bknisbac/CLL/results_cll1085/revision_natgen'
PHENOCOPY_TRI12_ROOT_DIR = '/Volumes/xchip_cga_home/bknisbac/CLL/results_cll1085/sample_sets/core2_tumors/degs_bc9_v1.5_clean1_original_counts_phenocopy_tri12'
PHENOCOPY_TRI12_ANALYSIS_OUTDIR = paste0(PHENOCOPY_TRI12_ROOT_DIR, '/analyses')
PHENOCOPY_TRI12_VAR_LIST = list()
PHENOCOPY_TRI12_VAR_LIST[['phenocopy_runs']] = c('mcll_tri12', 'mcll_non_tri12', 'ucll_tri12', 'ucll_non_tri12')
PHENOCOPY_TRI12_VAR_LIST[['ec_per_pheno_run']] = c(mcll_tri12='EC3', mcll_non_tri12='EC3', ucll_tri12='EC5', ucll_non_tri12='EC5')
PHENOCOPY_TRI12_VAR_LIST[['comp_per_pheno_run']] = c(mcll_tri12='EC3vsAll', mcll_non_tri12='EC3vsAll', ucll_tri12='EC5vsAll', ucll_non_tri12='EC5vsAll')
PHENOCOPY_TRI12_VAR_LIST[['fdr_cutoff']] = 0.1

SUPP_TABLE_ROOT_DIR = '/xchip/cga_home/bknisbac/CLL/results_cll1085/tables/supp_tables'
SUPP_TABLE_TIMESTAMP = '20211016'
SUPP_TABLE_DIR = paste0(SUPP_TABLE_ROOT_DIR, '/', 'supp_tables_', SUPP_TABLE_TIMESTAMP)

### For Supplementary table assembly
NOTCH1_TARGET_SEQ_TERRA_MODEL_SAMPLES = '/xchip/cga_home/bknisbac/CLL/results_cll1085/tables/supp_tables/for_supp/Broad_CLL_NOTCH1_targeted_seq_rerun20201013_samples.tsv'
NOTCH1_TARGET_SEQ_WELL_TO_BARCODE = '/xchip/cga_home/bknisbac/data/cll/targeted_notch1_seq_broad/array_positional_index_mapping.txt'
WES_SAMPLE_COVERAGE = '/xchip/cga_home/bknisbac/CLL/results_cll1085/tables/supp_tables/for_supp/wes_samples_coverage_20210326.txt'
WGS_SAMPLE_COVERAGE = '/xchip/cga_home/bknisbac/CLL/results_cll1085/tables/supp_tables/for_supp/wgs_samples_coverage_20210326.txt'
RNA_SEQUENCING_METRICS = '/xchip/cga_home/bknisbac/CLL/results_cll1085/sample_sets/all_tumors_20201216_with_dietrich/all_tumors_20201216_with_dietrich.metrics.txt.gz'

SURFACEOME_GENES_FILE = '/xchip/cga_home/bknisbac/resources/proteomics/surfaceome/surfaceome_proteins_genes_insilico_Bausch-Fluck_2018.tsv'

get_surfaceome_genes <- function(drop_dups=T){
  surf_dt = fread(SURFACEOME_GENES_FILE)
  setnames(surf_dt, 'Surfaceom_Label', 'Surfaceome_Label')
  if(drop_dups){
    surf_dt = surf_dt[!duplicated(UniProt_gene)]
  }
  return(surf_dt)
}

#### IGHV
IGHV_FREEZE_FILE = '/xchip/cga_home/bknisbac/CLL/results_cll1085/tables/cll1100_ighv_status_extended_20201223.tsv'

##### DNA data
MUTSIG_PATIENT_COUNTS_AND_RATES_WES = '/xchip/cga_home/bknisbac/CLL/results_cll1085/dna/mutsig/mutsig_20200703/patient_counts_and_rates.txt'
MUTSIG_PATIENT_COUNTS_AND_RATES_WES_PROCESSED = '/xchip/cga_home/bknisbac/CLL/results_cll1085/dna/mutsig/mutsig_20200703/patient_counts_and_rates.processed.txt'
MUTSIG_PATIENT_COUNTS_AND_RATES_WGS = '/xchip/cga_home/bknisbac/CLL/results_cll1085/dna/nc_drivers/patient_counts_and_rates_pseudo-WES-177.txt'
MUTSIG_PATIENT_COUNTS_AND_RATES_WES_BLACKLIST = c('804TD') #to prefer WGS for this one
MUTSIG_PATIENT_COUNTS_AND_RATES_WES_AND_WGS = '/xchip/cga_home/bknisbac/CLL/results_cll1085/tables/patient_counts_and_rates_wes984_wgs177.txt'

#### CNA data
SEG_FILES = list()
SEG_FILES[['wes_all']] = '/xchip/cga_home/bknisbac/CLL/results_cll1085/dna/gistic/gistic_20201110_cindy/raw/seg_files/WES_pairs_rerun_final_20200703.aggregated.seg'
SEG_FILES[['wgs_all']] = '/xchip/cga_home/bknisbac/CLL/results_cll1085/dna/wgs/seg_files/WGS177_aggregated_20201223.seg' #I gave timestamp by date in google drive

### CLUMPS files
CLUMPS_ROOT_DIR = '/xchip/cga_home/bknisbac/CLL/results_cll1085/dna/clumps/clumps_run_100220'
CLUMPS_DIRS = list()
CLUMPS_DIRS[['root']] = CLUMPS_ROOT_DIR
CLUMPS_DIRS[['all']] = paste0(CLUMPS_ROOT_DIR, '/', 'wes984')
CLUMPS_DIRS[['mcll']] = paste0(CLUMPS_ROOT_DIR, '/', 'mcll512')
CLUMPS_DIRS[['ucll']] = paste0(CLUMPS_ROOT_DIR, '/', 'ucll459')

get_clumps_files <- function(ftype='clumps', s='all', ret='dt'){
  croot = CLUMPS_DIRS[[s]]
  # Significance
  if(ftype=='clumps'){
    cf = paste0(croot, '/', 'clumps_ouput.tsv')
    if(ret=='file')
      return(cf)
    c_dt = fread(cf)[,-c('V1')]
  }
  # Mutations per patient
  if(ftype=='muts.gp'){
    cf = paste0(croot, '/', 'muts.gp')
    if(ret=='file')
      return(cf)
    c_dt = fread(cf)
    c_dt[['uniprot']] = sub(':.*', '', c_dt$uniprot_change)
    setnames(c_dt, 'patient', 'tumor_sample_id')
  }
  return(c_dt)
}

### Drivers
DRIVER_MASTER_TABLE = '/xchip/cga_home/bknisbac/CLL/results_cll1085/tables/drivers_merged_past_and_present_20201109_concise.edited_20201229.txt'
DRIVER_BLACKLIST_FILE = '/xchip/cga_home/bknisbac/CLL/results_cll1085/tables/driver_gene_blacklist_manual_20201130.txt'
DRIVER_BLACKLIST_BY_EXPRESSION_FILE = '/xchip/cga_home/bknisbac/CLL/results_cll1085/dna/drivers/drivers_20201018/cll_driver_gene_expression_blacklist.txt'
DRIVER_SYMBOL_UPDATE_PAIRS = list(c('MLL2', 'KMT2D'))

CLUMPS_DRIVERS_MULTILIST = list(all=c('RPS23', 'NCAPG', 'MAP2K2', 'DICER1', 'RAF1', 'DIS3'), mcll=c('MYLK4', 'BRAF', 'DICER1'), ucll=c('RPS23','RRM1', 'RAF1', 'MAP2K2'))
CLUMPS_DRIVERS_MULTILIST_NO_RHT = list(all=c('RPS23', 'NCAPG', 'MAP2K2'), mcll=c('MYLK4'), ucll=c('RPS23','RRM1'))
CLUMPS_DRIVERS_NOVELTY_TABLE = '/xchip/cga_home/bknisbac/CLL/results_cll1085/tables/clumps_drivers_all_mcll_ucll_novelty_20201229.tsv'
DRIVERS_BED_FILE = '/xchip/cga_home/bknisbac/CLL/results_cll1085/tables/drivers_snv_all_mcll_ucll_union_mutsig_clumps_20201211.bed'
DRIVER_GSEA_TOP5BEST = '/xchip/cga_home/bknisbac/CLL/results_cll1085/dna/drivers/drivers_20210107/pathways/results_20210107/gprofiler_All_M-CLL_U-CLL_mutsig_clumps_withClusters_C2-CP-REACTOME_C5-BP_H_top5byFDR.tsv'

#Drivers used in timing analysis 2021/04/05
DRIVERS_LIST_IN_TIMING_ANALYSIS = '/xchip/cga_home/bknisbac/CLL/results_cll1085/tables/driver_genes_included_in_timing_20210405.txt'

#MAFs
#Note that curveball MAF has samples not in the final curveball analysis (IGHV NAs)
DRIVER_CURVEBALL_MAF = '/xchip/cga_home/bknisbac/CLL/results_cll1085/dna/mafs/MAF_clinical_assembly_20201105/wes_wgs_broad_rfcaller_notch1_merged.withCLUMPS.for_curveball.20201209.maf'
DRIVER_MAF = DRIVER_CURVEBALL_MAF
#The difference from ".for_curveball." maf is that here MBLs are dropped - BETTER TO REFERENCE ONLY THE CURVEBALL COMUT WHEN POSSIBLE
DRIVER_CLINICAL_MAF = '/xchip/cga_home/bknisbac/CLL/results_cll1085/dna/mafs/MAF_clinical_assembly_20201105/wes_wgs_broad_rfcaller_notch1_merged.withCLUMPS.for_clinical.20201209.maf'

DRIVER_MERGE_UNFILTERED_MAF = '/xchip/cga_home/bknisbac/CLL/results_cll1085/dna/mafs/MAF_clinical_assembly_20201105/wes_wgs_broad_rfcaller_notch1_merged.20201209.UNFILTERED_EXTENDED.maf'

BROAD_MAF_LIST = list()
BROAD_MAF_LIST[['wes']] = '/xchip/cga_home/bknisbac/CLL/results_cll1085/dna/mafs/MAF_20200914/WES_pairs_rerun_final_20200703.aggregated.exclude.bt.select.columns.tsv'
BROAD_MAF_LIST[['wes_mcll']] = '/xchip/cga_home/bknisbac/CLL/results_cll1085/dna/mafs/MAF_20200914/mcll512_with_sanger_notch1_coding_20200924.tsv.gz'
BROAD_MAF_LIST[['wes_ucll']] = '/xchip/cga_home/bknisbac/CLL/results_cll1085/dna/mafs/MAF_20200914/ucll459_with_sanger_notch1_coding_20200924.tsv.gz'
BROAD_MAF_LIST[['wgs']] = '/xchip/cga_home/bknisbac/CLL/results_cll1085/dna/wgs/wgs_20200701/concat_maf.aggregated.WGS-177_20200701.mafcols.tsv.gz'

ABSOLUTE_WES_CCF = '/xchip/cga_home/bknisbac/CLL/results_cll1085/dna/absolute/purity_ploidy_ccf/wes_maf_ccf_20210406.txt'
ABSOLUTE_WES_PURITY_PLOIDY = '/xchip/cga_home/bknisbac/CLL/results_cll1085/dna/absolute/purity_ploidy_ccf/wes_pairs_purity_ploidy_20210405.txt'

ABSOLUTE_WGS_PURITY_PLOIDY = '/xchip/cga_home/bknisbac/CLL/results_cll1085/dna/wgs/wgs_pairs_purity_20210714.tsv'
MUTATIONAL_SIGNATURES_H_MATRIX_WGS = '/xchip/cga_home/bknisbac/CLL/results_cll1085/dna/wgs/mutational_signatures_H_matrix_20210812.tsv'


#Driver-related tables
DRIVER_CCF_TABLE_WES = '/xchip/cga_home/bknisbac/CLL/results_cll1085/dna/drivers/driver_genes_ccf_table.20210324.txt'

#Tables
U1_MUT_STATUS_FILE_ICGC_RAW = '/xchip/cga_home/bknisbac/CLL/results_cll1085/sharing/ferran_u1_mut/U1_status_ICGC_cohort_20200723.tsv'
U1_MUT_STATUS_FILE_BROAD_RAW = '/xchip/cga_home/bknisbac/CLL/results_cll1085/sharing/ferran_u1_mut/U1_status_predicted_based_on_splicing_Broad_cohort_20200723.tsv'
U1_MUT_STATUS_FILE = '/xchip/cga_home/bknisbac/CLL/results_cll1085/tables/ICGC_DFCI_U1_mut_status_20200723.tsv'
U1_MUT_WITH_BROAD_MODELS_RAW = '/xchip/cga_home/bknisbac/CLL/results_cll1085/sharing/ferran_u1_mut/ICGC_DFCI_U1_prediction_RNAseqSplicing_20200724.added_ICGC_pat.tsv'

IGLV3_21_R110_MUTATION_TABLE = '/xchip/cga_home/bknisbac/CLL/results_cll1085/tables/IGLV3-21_R110_mutations_per_participant.tsv'
BCR_STEREOTYPE_FILE = '/xchip/cga_home/bknisbac/CLL/results_cll1085/tables/cll1100_bcr_stereotype.extracted_from_sup_table.20211019.tsv'

EVENT_NOVELTY_TABLE = '/xchip/cga_home/bknisbac/CLL/results_cll1085/tables/event_novelty_table.20210120.tsv'


CLINICAL_ANALYSIS_DATA_DIR = '/xchip/cga_home/bknisbac/CLL/results_cll1085/clinical_analysis/data'
CLINICAL_ANALYSIS_DATA_TIMESTAMP = '20201211'
CLINICAL_COMUT_TIMESTAMP = '20210702'
CLINICAL_COMUT_FILE_LIST = list()
for(s in ANALYSIS_SETS){
  CLINICAL_COMUT_FILE_LIST[[s]] = paste0(CLINICAL_ANALYSIS_DATA_DIR, '/', s, '_comut_', CLINICAL_COMUT_TIMESTAMP, '.tsv')
}
CLINICAL_MULTIOMIC_TIMESTAMP = '20210702'
CLINICAL_COMUT_FILE_LIST[['all_multiomic']] = paste0(CLINICAL_ANALYSIS_DATA_DIR, '/', 'all_multiomic', '_comut_', CLINICAL_MULTIOMIC_TIMESTAMP, '.tsv')

CLINICAL_DATA_ROOTDIR = '/xchip/cga_home/bknisbac/CLL/clinical'
CLINICAL_DATA_DIR = '/xchip/cga_home/bknisbac/CLL/clinical/clinical_data_20210222'
CLINICAL_DATA_FULL_FILE = '/xchip/cga_home/bknisbac/CLL/clinical/clinical_data_20210222/output_cluster_fullset.20210222.csv'
CLINICAL_RNA_SEX_PREDICTIONS_FILE = '/xchip/cga_home/bknisbac/CLL/clinical/clinical_sex_predicted_KNN_all_rnaseqc_20200323.txt'
#NOTE: uneven column length
CLINICAL_ANALYSIS_SETS = '/xchip/cga_home/bknisbac/CLL/results_cll1085/tables/clinical_analysis_sets_fromKristen.20210701.tsv'

CURVEBALL_ROOT_DIR = '/xchip/cga_home/bknisbac/CLL/results_cll1085/curveball'
CURVEBALL_COMUT_DIR = '/xchip/cga_home/bknisbac/CLL/results_cll1085/curveball/comuts'
CURVEBALL_COMUT_GENOMIC_TIMESTAMP = '20210702' #not used in manuscript
CURVEBALL_COMUT_MULTIOMIC_TIMESTAMP = '20210702'
CURVEBALL_COMUT_FILE_LIST = list()

CURVEBALL_COMUT_FILE_LIST[['mvu_union_genomic']] = paste0(CURVEBALL_COMUT_DIR, '/', 'comut_for_curveball_wes_wgs_noContam_withMBL_withIGHV_', CURVEBALL_COMUT_GENOMIC_TIMESTAMP, '.tsv')
CURVEBALL_COMUT_FILE_LIST[['mvu_union_multiomic']] = paste0(CURVEBALL_COMUT_DIR, '/', 'comut_for_curveball_wes_wgs_noContam_withMBL_withIGHV_multiomic_', CURVEBALL_COMUT_MULTIOMIC_TIMESTAMP, '.tsv')

CURVEBALL_RESULT_FILES = list()
CURVEBALL_RESULT_FILES[['multiomic_odds_ratio']] = '/xchip/cga_home/bknisbac/CLL/results_cll1085/curveball/results/multiomics_curveball/revision_20210706/multiomics_odds_ratio_iter_0.txt'
CURVEBALL_RESULT_FILES[['multiomic_co_p']] = '/xchip/cga_home/bknisbac/CLL/results_cll1085/curveball/results/multiomics_curveball/revision_20210706/multiomics_positive_p_values_df_iter_0.txt'
CURVEBALL_RESULT_FILES[['multiomic_me_p']] = '/xchip/cga_home/bknisbac/CLL/results_cll1085/curveball/results/multiomics_curveball/revision_20210706/multiomics_negative_p_values_df_iter_0.txt'

STATS_COMUT_DIR = '/xchip/cga_home/bknisbac/CLL/results_cll1085/curveball/comuts'
STATS_COMUT_GENOMIC_TIMESTAMP = CURVEBALL_COMUT_GENOMIC_TIMESTAMP
STAT_COMUT_FILE_LIST = list()
#this has 1064 patients as in Fig. 1 comut (CNAs are like fig. 1, but note it has union SNV)
STAT_COMUT_FILE_LIST[['all_cna_union_snv_genomic']] = paste0(STATS_COMUT_DIR, '/', 'comut_for_stats_wes_wgs_noContam_withMBL_withIGHV_', STATS_COMUT_GENOMIC_TIMESTAMP, '.tsv')

### CNAs
CNA_ROOT_DIR = '/xchip/cga_home/bknisbac/CLL/results_cll1085/dna/gistic/gistic_20201120'
CNA_TIMESTAMP = '20201229'
CNA_ALL_DRIVER_TABLE = paste0(CNA_ROOT_DIR, '/', 'gistic_master_table_All_20201227.tsv')
CNA_ALL_DRIVER_TABLE_ANNOTATED = sub('.tsv$', '.annotated.tsv', CNA_ALL_DRIVER_TABLE)
CNA_ALL_WIDE_PEAK_GENES = paste0(CNA_ROOT_DIR, '/', 'genes_by_gistic_All.tsv')
CNA_MVU_DRIVER_COMPARISON_TABLE = paste0(CNA_ROOT_DIR, '/', 'gistic_master_table_MvU_mcll_vs_ucll_20201206.tsv')
CNA_MVU_DRIVER_COMPARISON_TABLE_ANNOTATED = sub('.tsv$', '.annotated.tsv', CNA_MVU_DRIVER_COMPARISON_TABLE)
CNA_MVU_WIDE_PEAK_GENES = paste0(CNA_ROOT_DIR, '/', 'genes_by_gistic_MvU.tsv')
CNA_MVU_BED = paste0(CNA_ROOT_DIR, '/', 'gistic_master_table_MvU_mcll_vs_ucll_20201206.bed')

### CNA comuts
CNA_COMUT_FILE_LIST = list()
for (d in c('wes', 'wgs', 'wes_wgs', 'wes_wgs_withContam')){
  for (s in c('all','mcll', 'ucll', 'mcll_mvu_union', 'ucll_mvu_union')){
    ds = paste0(d, '_', s)
    CNA_COMUT_FILE_LIST[[ds]] = paste0(CNA_ROOT_DIR, '/', 'assembled_cna_comuts', '/', 'comut_cna_arm_focal_', ds, '_', CNA_TIMESTAMP, '.tsv')
  }
}

CNA_BED_FILE_LISTS = list()
CNA_BED_FILE_LISTS[['wes_wgs_withContam_all_withMvUspecific']] = paste0('/xchip/cga_home/bknisbac/CLL/results_cll1085/tables/cna_focalOnly_wes_wgs_withContam_all_withMvUspecific.', CNA_TIMESTAMP, '.bed')
CNA_BED_FILE_LISTS[['mvu_focals']] = sub('\\.[a-z]*$', '.focal.bed', CNA_MVU_DRIVER_COMPARISON_TABLE)
CNA_BED_FILE_LISTS[['mvu_focals_noHeader']] = sub('.bed', '.noHeader.bed', CNA_BED_FILE_LISTS[['mvu_table']])
CNA_BED_FILE_LISTS[['mvu_focals_extended']] = sub('\\.[a-z]*$', '.focal.extended.bed', CNA_MVU_DRIVER_COMPARISON_TABLE)


#### COMUTS
VARIANT_CLASS_REPLACE_LIST = list(c("_Mutation", ""), c("Silent", "Synonymous"), c("Splice_Site.*|Splice_Region", "Splice_site"),
                                  c("Nonstop|De_novo_Start.*|Start_.*|Translation_.*|Read\\-through.*|Stop_Codon_Del", "Other_non_syn."), c("In_frame.*", "In_frame_indel"), c("Frame_Shift.*", "Frame_shift"),
                                  c("3'\\-?UTR|5'\\-?UTR|3'\\-?Flank|5'\\-?Flank|IGR|Intron|RNA", "Non-coding"))
VARIANT_CLASSES_ORDERED = c("Frame_shift", "Nonsense", "Splice_site", "In_frame_indel", "Missense", "Other_non_syn.", "Non-coding", "Synonymous")


### Structural variants (SVs)
SV_FILES = list()
SV_FILES[['sv_annot']] = '/xchip/cga_home/bknisbac/CLL/results_cll1085/dna/SV/SV_merge_Broad_Puente2015_IgCaller.20210308.annot.20210308.tsv'
SV_FILES[['sv_annot_ccf']] = '/xchip/cga_home/bknisbac/CLL/results_cll1085/dna/SV/SV_merge_Broad_Puente2015_IgCaller.20210308.annot.20210308.withCCF.tsv'
SV_FILES[['sv_annot_wgsClustFN']] = '/xchip/cga_home/bknisbac/CLL/results_cll1085/dna/SV/from_ferran/sv_clustering_ferran/translocations_20210321.tsv'
SV_FILES[['sv_ccf']] = '/xchip/cga_home/bknisbac/CLL/results_cll1085/tables/sv_ccf_v3_wgs_20210930.tsv' #concat of raw results from Terra



### Epigenetic data
METHYLATION_MASTER_TABLE_RAW = '/xchip/cga_home/bknisbac/CLL/results_cll1085/tables/CLL.meth.complete.freeze.20210205.tsv'
METHYLATION_CONCISE_FILE = '/xchip/cga_home/bknisbac/CLL/results_cll1085/tables/CLL.meth.complete.freeze.20210205.soft_filtered.concise.tsv'

METHYLATION_SAMPLE_BLACKLIST = c() #list of IDs removed due to sharing restrictions in public repo


H3K27AC_DIR = '/xchip/cga_home/bknisbac/CLL/results_cll1085/epigenetics/h3k27ac'
H3K27AC_MARKERS_SUBDIR = 'EC_genes_enhancers_H3K27ac_differential_20210707'
H3K27AC_PER_MARKERS_UP_FILE = paste0(H3K27AC_DIR, '/', H3K27AC_MARKERS_SUBDIR, '/' , 'ALL_DESeq_regions_H3K27ac_DESeq_hg38_FC_Up_FDR_0.05_.tsv')
H3K27AC_PER_MARKERS_DOWN_FILE = paste0(H3K27AC_DIR, '/', H3K27AC_MARKERS_SUBDIR, '/', 'ALL_DESeq_regions_H3K27ac_DESeq_hg38_FC_Down_FDR_0.05_.tsv')
#Saved to compare to new results for rebuttal
H3K27AC_MARKERS_SUBDIR_FIRST_SUBMISSION = 'EC_leading_genes_20210219' #OLD
H3K27AC_PER_MARKERS_UP_FILE_FIRST_SUBMISSION = paste0(H3K27AC_DIR, '/', H3K27AC_MARKERS_SUBDIR_FIRST_SUBMISSION, '/' , 'ALL_DESeq_regions_H3K27ac_DESeq_hg38_FC_Up_FDR_0.05_.tsv') #OLD
H3K27AC_PER_MARKERS_DOWN_FILE_FIRST_SUBMISSION = paste0(H3K27AC_DIR, '/', H3K27AC_MARKERS_SUBDIR_FIRST_SUBMISSION, '/', 'ALL_DESeq_regions_H3K27ac_DESeq_hg38_FC_Down_FDR_0.05_.tsv') #OLD

H3K27AC_SAMPLE_FILE = '/xchip/cga_home/bknisbac/CLL/results_cll1085/epigenetics/h3k27ac/EC_genes_enhancers_H3K27ac_differential_20210707/SCLL_H3K27ac_RNA_seq_20210707.tsv'

H3K27AC_GENOMEWIDE_SUBDIR = 'genomewide_20210707'
H3K27AC_PER_GENOMEWIDE_UP_FILE_HG38 = paste0(H3K27AC_DIR, '/', H3K27AC_GENOMEWIDE_SUBDIR, '/' , 'Genome-wide.DESeq.hg38.EC.Up.tsv')
H3K27AC_PER_GENOMEWIDE_UP_FILE = paste0(H3K27AC_DIR, '/', H3K27AC_GENOMEWIDE_SUBDIR, '/' , 'Genome-wide.DESeq.hg19.EC.Up.tsv')
H3K27AC_PER_GENOMEWIDE_DOWN_FILE_HG38 = paste0(H3K27AC_DIR, '/', H3K27AC_GENOMEWIDE_SUBDIR, '/', 'Genome-wide.DESeq.hg38.EC.Down.tsv')
H3K27AC_PER_GENOMEWIDE_DOWN_FILE = paste0(H3K27AC_DIR, '/', H3K27AC_GENOMEWIDE_SUBDIR, '/', 'Genome-wide.DESeq.hg19.EC.Down.tsv')


### utilities
get_meth_dt <- function(opt='default', cols='default', fill.na='NA', add_cols=c()){
  meth_dt = fread(METHYLATION_MASTER_TABLE_RAW)
  meth_pat_convert = list(c('MBL1', 'SCLL-###')) # SCLL patient ID masked in public repo due to sharing restrictions

  if(opt=='default'){
    meth_dt = meth_dt[meth_dt$Blacklist=='FALSE' & meth_dt$Greylist_soft=='FALSE',]
  }

  cols_from = c('Participant_id', 'Sequencing_set', 'CLL_epitype', 'Cohort', 'Platform')
  cols_to = c('participant_id', 'sequencing_set', 'epitype', 'cohort_epi', 'platform')
  setnames(meth_dt, cols_from, cols_to)

  if(!is.null(meth_pat_convert)){
    for (pair in meth_pat_convert){
      meth_dt[meth_dt$participant_id==pair[1]][['participant_id']] = pair[2]
    }
  }

  if(cols=='default'){
    meth_dt = meth_dt[,unique(c(c('participant_id', 'epitype', 'epiCMIT.clinics'), add_cols)), with=F]
  } else if(cols=='extended_basic'){

    meth_dt = meth_dt[,unique(c(c('participant_id', 'cohort_epi', 'platform', 'sequencing_set', 'epitype', 'epiCMIT.clinics'), add_cols)), with=F]
  }

  if(!is.null(fill.na)){
    meth_dt[is.na(meth_dt)] = fill.na
  }

  meth_dt = meth_dt[order(participant_id),]

  return(meth_dt)
}

GENERATE_METH_TABLES = F
if(GENERATE_METH_TABLES){
  fwrite(get_meth_dt(cols='extended_basic'), METHYLATION_CONCISE_FILE, sep='\t')
}

#Note: By default (when opt=='pat_to_binary_mut', fill.na==NULL), then NA get a value of 0; use fill.na=NA to return NAs
get_iglv321_mut <- function(opt='pat_to_binary_mut', fill.na=NULL){
  if (opt=='pat_comut_col_dt'){
    r110_mut_dt = fread(IGLV3_21_R110_MUTATION_TABLE)
    return(r110_mut_dt)
  }

  if (opt=='pat_to_binary_mut'){
    r110_mut_dt = fread(IGLV3_21_R110_MUTATION_TABLE)
    if(!is.null(fill.na)){
      na_pats = (is.na(r110_mut_dt$IGLV321_R110_comut) | r110_mut_dt$IGLV321_R110_comut=='NA' | r110_mut_dt$IGLV321_R110_comut=='')
      r110_mut_dt[na_pats, 'IGLV321_R110_binary'] = fill.na
    }
    r110_mut_dt = r110_mut_dt[,c('participant_id', 'IGLV321_R110_binary')]
    setnames(r110_mut_dt, 'IGLV321_R110_binary', 'IGLV321_R110')
    return(r110_mut_dt)
  }

  if (opt=='unknown_as_binary'){
    r110_unk_dt = fread(IGLV3_21_R110_MUTATION_TABLE)
    r110_unk_dt$IGLV321_R110_unk = as.integer(is.na(r110_unk_dt$IGLV321_R110_comut))
    return(r110_unk_dt[,c('participant_id', 'IGLV321_R110_unk'),with=F])
  }
}

get_master_dt <- function(full=F, drop_dups=T, drop_redund_of_crc=T,
                          cohort_include = 'all',
                          cohort_drop = c('Beekman', 'DKFZ'),
                          clinical_pats = F,
                          tumors_only = F,
                          ighv = NULL,
                          data='any'){

  master_dt = fread(MASTER_METADTA_TABLE)
  if(full){ #overrides all options
    return(master_dt)
  }
  if(drop_dups)
    master_dt = master_dt[is_duplicate != 1]
  if(drop_redund_of_crc)
    master_dt = master_dt[is_redund_of_crc != 1]
  if(!is.null(cohort_drop) && length(cohort_drop)>0)
    master_dt = master_dt[!cohort %in% cohort_drop]
  if(cohort_include != 'all')
    master_dt = master_dt[cohort %in% cohort_include]
  if(clinical_pats){
    master_dt = master_dt[clinical_blacklist_pat != 1]
    ### TO DO: drop MBLs on demand
  }
  if(tumors_only)
    master_dt = master_dt[has_tumor == 1]
  if(!is.null(ighv)){
    master_dt = master_dt[IGHV_mut_freeze %in% ighv]
  }

  if(data=='rna'){
    master_dt = master_dt[(!is.na(rna_tumor_sample_id) & !rna_tumor_sample_id=='') | (!is.na(rna_normal_sample_id) & !rna_normal_sample_id==''),]
  } else if(data=='wes'){
    master_dt = master_dt[!is.na(wes_tumor_sample_id) & !wes_tumor_sample_id=='',]
  } else if(data=='wgs'){
    master_dt = master_dt[!is.na(wgs_tumor_sample_id) & !wgs_tumor_sample_id=='',]
  } else if(data=='dna'){
    master_dt = master_dt[(!is.na(wes_tumor_sample_id) & !wes_tumor_sample_id=='') | (!is.na(wgs_tumor_sample_id) & !wgs_tumor_sample_id==''),]
  }

  return(master_dt)
}


get_driver_blacklist <- function(file=DRIVER_BLACKLIST_FILE){
  driver_blacklist = sort(unique(fread(file, header=F)[['V1']]))
  return(driver_blacklist)
}

update_driver_symbols <- function(glist, update_pairs=DRIVER_SYMBOL_UPDATE_PAIRS){
  for (pair in update_pairs){
    if(pair[1] %in% glist)
      glist[glist==pair[1]] = pair[2]
  }
  return(glist)
}

get_genome_annot <- function(opt, subset=NULL, bed_col='gene_name'){
  if(opt=='gene_bed'){
    gene_bed = fread(GENE_BED_FILE)
    if(!is.null(subset)){
      gene_bed = gene_bed[gene_bed[[bed_col]] %in% subset,]
    }
    return(gene_bed)
  }
}

get_drivers <- function(sets=ANALYSIS_SETS,
                        driver_master_dt = fread(DRIVER_MASTER_TABLE),
                        driver_blacklist = get_driver_blacklist(),
                        symbol_update_pairs = DRIVER_SYMBOL_UPDATE_PAIRS,
                        with_clumps=T,
                        clumps_no_rht=F,
                        clumps_drivers_multilist=CLUMPS_DRIVERS_MULTILIST,
                        clumps_drivers_multilist_no_rht=CLUMPS_DRIVERS_MULTILIST_NO_RHT,
                        ret='multilist'){

  drivers_list = list()
  for (s in sets){
    drivers_list[[s]] = driver_master_dt[[paste0(s, '_sig')]][driver_master_dt[[paste0(s, '_sig')]] != ""]
    if(with_clumps){
      if(clumps_no_rht){
        drivers_list[[s]] = sort(unique(c(drivers_list[[s]], clumps_drivers_multilist_no_rht[[s]])))
      } else {
        drivers_list[[s]] = sort(unique(c(drivers_list[[s]], clumps_drivers_multilist[[s]])))
      }
    }
    drivers_list[[s]] = drivers_list[[s]][!drivers_list[[s]] %in% driver_blacklist]
    if(!is.null(symbol_update_pairs))
      drivers_list[[s]] = update_driver_symbols(drivers_list[[s]]) #update names (needed because mutsig had old names)
  }

  if (ret=='multilist'){
    return(drivers_list)
  } else if(ret=='union') {
    drivers_union = sort(unique(unlist(drivers_list)))
    return(drivers_union)
  }
}

get_drivers_union <- function(with_clumps=T){
  get_drivers(ret='union', with_clumps=with_clumps)
}

get_cna_comut <- function(dtype='wes', s='all', cna_comut_files=CNA_COMUT_FILE_LIST){
    if(s=='mvu_union'){ #special handling - merge
      cm_mcll = fread(cna_comut_files[[paste(dtype, 'mcll_mvu_union', sep='_')]])
      cm_ucll = fread(cna_comut_files[[paste(dtype, 'ucll_mvu_union', sep='_')]])
      cm = rbindlist(list(cm_mcll, cm_ucll))
      setorderv(cm, 'participant_id')
    } else {
      if(is.list(dtype))
        d = tolower(paste(sort(dtype), collapse='_'))
      else
        d = dtype

      ds = paste0(d, '_', s)
      cm = fread(cna_comut_files[[ds]])
    }
  return(cm)
}

#Note that MvU option has some All-only events, use sets=c('mcll', 'ucll') for those only
get_cna_events <- function(s='all', sets=NULL, all_analysis_sets=ANALYSIS_SETS){ #no diff between WGS and WES so no need for dtype
  cna_events = c()
  if(s=='union' || !is.null(sets)){ #get multiple sets
    if(s=='union' || (!is.null(sets) && sets=='union')){
      sets = all_analysis_sets
    }
    for (s in sets){
      cm = get_cna_comut(s=s)
      cna_events0 = colnames(cm)[colnames(cm) != 'participant_id']
      cna_events = sort(unique(c(cna_events, cna_events0)))
    }
  } else if(s=='MvU' || s=='mvu') { # use s value for single set
    cna_events = fread(CNA_MVU_DRIVER_COMPARISON_TABLE)$event_id
  } else { # use s value for single set
      cm = get_cna_comut(s=s)
      cna_events = colnames(cm)[colnames(cm) != 'participant_id']
  }
  return(cna_events)
}

get_U1_muts <- function(infile_raw=U1_MUT_STATUS_FILE_ICGC_RAW, infile_raw2=U1_MUT_STATUS_FILE_BROAD_RAW, file=U1_MUT_STATUS_FILE, id='participant_id', mdt=master_dt_full){
  ### For clinical data rows
  if(file.exists(file)){
    u1_muts = fread(file)
  }
  else {
    u1_muts_icgc = fread(infile_raw)
    u1_muts_icgc$participant_id = sub('^(\\d+)$', 'SCLL-\\1', sub('^(\\d\\d)$', '0\\1', sub('^(\\d)$', '0\\1', u1_muts_icgc$Case)))
    u1_muts_icgc[,Case:=NULL]
    setnames(u1_muts_icgc, c('U1_status'), c('U1_status'))

    #Note: this drops the Normals, which were just for control on the prediction process
    u1_muts_broad = fread(infile_raw2)
    setnames(u1_muts_broad, c('submitted_sample_id', 'U1_status_splicing_prediction'), c('rna_tumor_sample_id', 'U1_status'))
    u1_muts_broad = merge(u1_muts_broad, mdt[mdt$rna_tumor_sample_id!="", c('rna_tumor_sample_id', 'participant_id'),with=F], by='rna_tumor_sample_id')[,-c('rna_tumor_sample_id'),with=F]

    u1_muts = rbind(u1_muts_icgc, u1_muts_broad)
    setcolorder(u1_muts, rev(colnames(u1_muts)))
    u1_muts$U1 = ifelse(u1_muts$U1_status=='MUT', 1, 0)
    setorder(u1_muts, participant_id)

    if(!is.null(file))
      fwrite(u1_muts, file, sep='\t')
  }

  return(u1_muts)
}

u1_muts = get_U1_muts()

## Options:
# analysis = c('clinical', 'curveball')
# s = c('all', 'mcll', 'ucll')
# driver_set = c('all', 'mcll', 'ucll', 'union')
assemble_analysis_comut <- function(s, driver_set, timestamp, drop_pat_blacklist=T, drop_mbl=F, add_iglv321_mut=T, add_u1=T, outfile=NULL,
                                    mdt=get_master_dt(full=T), driver_maf=DRIVER_MAF, add_ighv_cols=F, ighv_cols=c('MCLL', 'UCLL'), fill.na=NULL, retval=NULL){

  #subset and merge with driver comut
  if(driver_set=='union'){
    drivers_list = get_drivers(ret='union')
  } else {
    drivers_list = get_drivers()[[s]]
  }

  #Get driver gene comut (mutsig+clumps)
  driver_comut = comut.utils$parse_maf(driver_maf, id_map_file = NULL,
                                       classes_to_drop = c(),
                                       classes_to_drop_renamed = c(),
                                       class_replace_pairs= NULL,
                                       drop_silent=F,
                                       sample_by_gene=F,
                                       maf_sample_col='participant_id',
                                       variant_classes_ordered=VARIANT_CLASSES_ORDERED,
                                       retval='comut', verbose=T)

  #get driver CNA comut
  cna_cm = get_cna_comut(dtype='wes_wgs', s=s)

  #merge CNA + SNV (CNAs define the set, hence all.x=T)
  cm = merge(cna_cm, driver_comut[, c('participant_id', drivers_list), with=F], by='participant_id', all.x=T)
  cm[is.na(cm)] = 0

  #Add IGLV32-1
  if(add_iglv321_mut){ #even if the function level fill.na is not NULL, this is handled later
    cm = merge(cm, get_iglv321_mut('pat_to_binary_mut', fill.na=NULL), by='participant_id', all.x=T)
  }

  #Add U1
  if(add_u1){ #even if the function level fill.na is not NULL, this is handled later
    cm = merge(cm, get_U1_muts()[,c('participant_id', 'U1'),with=F], by='participant_id', all.x=T)
  }

  if(!is.null(fill.na)){
    cm[is.na(cm)] = fill.na
    cm[cm=='NA'] = fill.na
    cm[cm==''] = fill.na
  }

  #convert numeric cols to numeric
  numeric_cols = colnames(cm)[2:ncol(cm)]
  cm[, (numeric_cols) := lapply(.SD, as.numeric), .SDcols=numeric_cols]

  if(drop_mbl){
    cm = cm[!cm$participant_id %in% mdt[mdt$tumor_type=='MBL']$participant_id]
  }

  if(drop_pat_blacklist){
    cm = cm[!cm$participant_id %in% mdt[mdt$clinical_blacklist_pat==1]$participant_id]
  }

  if(add_ighv_cols){
      cm[[ighv_cols[1]]] = 0
      cm[cm$participant_id %in% mdt[mdt$IGHV_mut_freeze=='mutated']$participant_id, ighv_cols[1]] = 1
      cm[[ighv_cols[2]]] = 0
      cm[cm$participant_id %in% mdt[mdt$IGHV_mut_freeze=='unmutated']$participant_id, ighv_cols[2]] = 1
  }

  if(!is.null(outfile))
    fwrite(cm, outfile, sep='\t')

  if(!is.null(retval)){
    if(retval=='data.table')
      return(cm)
    if(retval=='file')
      return(outfile)
  }

}


coords_to_bed <- function(coords, convert_to_broad_chr=F, decrement_start=F){
  bed = data.table()
  bed$chr = sub(':.*', '', coords)
  if(convert_to_broad_chr){
    bed$chr = sub('chrM', 'MT', bed$chr)
    bed$chr = sub('^chr', '', bed$chr)
  }
  bed$start = sub('-.*','', sub('.*:', '', coords))
  bed$end = sub('.*-','', sub('.*:', '', coords))

  if(decrement_start)
    bed$start = bed$start - 1

  return(bed)
}


assemble_cna_comut_for_sv_intersect <- function(outfile=CNA_BED_FILE_LISTS[['wes_wgs_withContam_all_withMvUspecific']]){
  #Get all event coords
  all_event_dt = fread(CNA_ALL_DRIVER_TABLE)[!grepl('arm', Copy_number)][,c('event_id', 'All_focal_wide_peak')]
  setnames(all_event_dt, c('event_id', 'wide_peak'))
  all_event_bed = cbind(coords_to_bed(all_event_dt$wide_peak, convert_to_broad_chr = T), all_event_dt)

  cna_all = fread(CNA_COMUT_FILE_LIST[['wes_wgs_withContam_all']])
  cna_all = melt(cna_all, id.vars = 'participant_id', variable.name = 'event_id', value.name='has_cna')[has_cna==1]
  cna_all = merge(all_event_bed, cna_all, by='event_id') #drops Arm events


  mu_only_events = fread(CNA_MVU_DRIVER_COMPARISON_TABLE)[!grepl('All|ALL', Category) & !grepl('arm', Copy_number)][,c('event_id', 'Wide_peak')]
  setnames(mu_only_events, c('event_id', 'wide_peak'))
  mu_event_bed = cbind(coords_to_bed(mu_only_events$wide_peak, convert_to_broad_chr = T), mu_only_events)

  cna_m = fread(CNA_COMUT_FILE_LIST[['wes_wgs_withContam_mcll_mvu_union']])[,c('participant_id', mu_only_events$event_id), with=F]
  cna_u = fread(CNA_COMUT_FILE_LIST[['wes_wgs_withContam_ucll_mvu_union']])[,c('participant_id', mu_only_events$event_id), with=F]
  cna_m_u = rbindlist(list(cna_m, cna_u))

  cna_m_u = melt(cna_m_u, id.vars = 'participant_id', variable.name = 'event_id', value.name='has_cna')[has_cna==1]
  cna_m_u = merge(mu_event_bed, cna_m_u, by='event_id') #drops Arm events

  cna_union = rbindlist(list(cna_all, cna_m_u))
  cna_union$score = '.'
  cna_union$strand = ifelse(grepl('[Dd]el', cna_union$event_id), '-', '+')
  setcolorder(cna_union, c('chr', 'start', 'end', 'event_id', 'score', 'strand', 'wide_peak', 'participant_id', 'has_cna'))
  setorderv(cna_union, c('participant_id', 'event_id'))

  fwrite(cna_union, outfile, sep='\t')
}

cna_table_to_bed <- function(table_file=CNA_MVU_DRIVER_COMPARISON_TABLE,
                             coord_col='Wide_peak',
                             out_bedfile=CNA_BED_FILE_LISTS[['mvu_table']], with_header=T, chr_order=c(1:22,'X','Y','MT'),
                             other_cols=NULL){

  bed = fread(table_file)[!grepl('[Aa]rm', Copy_number),c('event_id', coord_col, other_cols),with=F]
  bed = cbind(coords_to_bed(bed$Wide_peak, convert_to_broad_chr = T), bed)
  bed$strand = ifelse(grepl('[Dd]el', bed$event_id), '-', '+')
  bed$score = '.'
  setcolorder(bed, c('chr', 'start', 'end', 'event_id', 'score', 'strand'))
  if(!is.null(chr_order))
    bed$chr = factor(bed$chr, levels=chr_order, ordered = T)
  setorderv(bed, c('chr', 'start'))

  fwrite(bed, out_bedfile, sep='\t', col.names=with_header)
}

ASSEMBLE_CNA_COMUT_FOR_SV_INTERSECT = F
if(ASSEMBLE_CNA_COMUT_FOR_SV_INTERSECT){
  assemble_cna_comut_for_sv_intersect(outfile=CNA_BED_FILE_LISTS[['wes_wgs_withContam_all_withMvUspecific']])

  cna_table_to_bed(table_file=CNA_MVU_DRIVER_COMPARISON_TABLE, coord_col='Wide_peak', out_bedfile=CNA_BED_FILE_LISTS[['mvu_focals']])
  cna_table_to_bed(table_file=CNA_MVU_DRIVER_COMPARISON_TABLE, coord_col='Wide_peak', out_bedfile=CNA_BED_FILE_LISTS[['mvu_focals_noHeader']], with_header = F)

  #CNAs of MvU but include M and U frequency info
  cna_table_to_bed(table_file=CNA_MVU_DRIVER_COMPARISON_TABLE, coord_col='Wide_peak', other_cols=c('MvU_MCLL_npat', 'MvU_UCLL_npat'), out_bedfile=CNA_BED_FILE_LISTS[['mvu_focals_extended']])
}


#### Derived ####
master_dt_full = fread(MASTER_METADTA_TABLE) #full master dt
master_dt = get_master_dt() #get master dt with default filtering

get_pat_set <- function(pset, ret_uniq=T){
  if(pset=='clin_blacklist_pat'){
    ret_set = master_dt_full[clinical_blacklist_pat==1]$participant_id
  } else if(pset=='nhlbi_ibrutinib_trial'){ #this was TN in ibrutinib trial (should include also NHLBI-0029 for WES, but had no RNA)
    ret_set = c('NHLBI-0007', 'NHLBI-0013', 'NHLBI-0014', 'NHLBI-0015', 'NHLBI-0016', 'NHLBI-0022', 'NHLBI-0023', 'NHLBI-0024', 'NHLBI-0025', 'NHLBI-0026', 'NHLBI-0030', 'NHLBI-0031', 'NHLBI-0032', 'NHLBI-0034', 'NHLBI-0035', 'NHLBI-0038', 'NHLBI-0039', 'NHLBI-0040', 'NHLBI-0043', 'NHLBI-0044', 'NHLBI-0045', 'NHLBI-0046', 'NHLBI-0047', 'NHLBI-0048', 'NHLBI-0049', 'NHLBI-0051')
  } else if(pset=='master'){
    ret_set = get_master_dt(drop_redund_of_crc = F)$participant_id
  } else if(pset=='master_tumors'){
    ret_set = get_master_dt(tumors_only = T)$participant_id
  } else if(pset=='master_with_redund'){
    ret_set = get_master_dt(drop_redund_of_crc = F)$participant_id
  } else if(pset=='master_with_redund_tumors'){
    ret_set = get_master_dt(drop_redund_of_crc = F, tumors_only=T)$participant_id
  } else if(pset=='mcll'){
    ret_set = get_ighv()[IGHV_mut=='mutated']$participant_id
  } else if(pset=='ucll'){
    ret_set = get_ighv()[IGHV_mut=='unmutated']$participant_id
  } else if(pset=='mcll_wes' || pset=='wes_mcll'){
    ret_set = get_master_dt()[IGHV_mut_freeze=='mutated' & wes_tumor_sample_id!='']$participant_id
  } else if(pset=='ucll_wes' || pset=='wes_ucll'){
    ret_set = get_master_dt()[IGHV_mut_freeze=='unmutated' & wes_tumor_sample_id!='']$participant_id
  } else if(pset=='wes'){
    ret_set = get_master_dt()[wes_tumor_sample_id!='']$participant_id
  }


  if(ret_uniq)
    ret_set = sort(unique(ret_set))

  return(ret_set)
}

get_sample_set <- function(sset, ret_uniq=T){
  if(sset=='rna_discovery'){
    ret_set = master_dt_full[expression_cluster!='']$rna_tumor_sample_id
  } else if (sset=='ec_train_and_validation'){
    ret_set = fread(EC_TRAIN_SET_SAMPLE_LIST, header=F)$V1
  } else if (sset=='ec_test'){ #leave-out set
    ec_train_set = fread(EC_TRAIN_SET_SAMPLE_LIST, header=F)$V1
    ec_discovery_set = master_dt_full[expression_cluster!='']$rna_tumor_sample_id
    ret_set = ec_discovery_set[!ec_discovery_set %in% ec_train_set]
  }

  if(ret_uniq)
    ret_set = unique(ret_set)

  return(ret_set)
}

get_driver_master_dt = function(file = DRIVER_MASTER_TABLE,
                                symbol_update_pairs = DRIVER_SYMBOL_UPDATE_PAIRS){
  dm_dt = fread(file)
  if(!is.null(symbol_update_pairs)){
    dm_dt$gene = update_driver_symbols(dm_dt$gene, update_pairs = symbol_update_pairs)
  }
  return(dm_dt)
}

driver_master_dt = fread(DRIVER_MASTER_TABLE)
driver_blacklist = sort(unique(fread(DRIVER_BLACKLIST_FILE, header=F)[['V1']]))

clumps_drivers = sort(unique(unlist(CLUMPS_DRIVERS_MULTILIST)))

cna_master_all_dt = fread(CNA_ALL_DRIVER_TABLE)
cna_master_mvu_dt = fread(CNA_MVU_DRIVER_COMPARISON_TABLE)


GENERATE_GENOMIC_ANALYSIS_COMUTS = F
if(GENERATE_GENOMIC_ANALYSIS_COMUTS){

  ## Curveball comut for actuall curveball analysis
  #curveball MvU-union-CNA / Union-SNV (this version is restricted to contain samples among union of (450 UCLL, 512 MCLL + 177 WGS) that have IGHV status, no contamination or GCLL-0136 (MCL))
  outfile = assemble_analysis_comut('mvu_union', driver_set='union', timestamp=CURVEBALL_COMUT_GENOMIC_TIMESTAMP,
                                    drop_mbl=F, add_iglv321_mut=T, add_u1=F, fill.na=0, outfile=CURVEBALL_COMUT_FILE_LIST[['mvu_union_genomic']], mdt=get_master_dt(full=T), driver_maf=DRIVER_MAF, add_ighv_cols=T, ighv_cols=c('MCLL', 'UCLL'), retval='file')
  print(outfile)

  ## Curveball-like table for stats
  #curveball All-CNA / Union-SNV (no contamination or GCLL-0136 (MCL))
  outfile = assemble_analysis_comut('all', driver_set='union', timestamp=STATS_COMUT_GENOMIC_TIMESTAMP,
                                    drop_mbl=F, add_iglv321_mut=T, add_u1=F, fill.na=0, outfile=STAT_COMUT_FILE_LIST[['all_cna_union_snv_genomic']], mdt=get_master_dt(full=T), driver_maf=DRIVER_MAF, add_ighv_cols=T, ighv_cols=c('MCLL', 'UCLL'), retval='file')

  print(outfile)

  #clinical analysis
  for (s in ANALYSIS_SETS){ #s='all'
    outfile = assemble_analysis_comut(s, driver_set=s, timestamp=CLINICAL_COMUT_TIMESTAMP,
                                      drop_mbl=T, add_iglv321_mut=T, add_u1=F, fill.na='NA',
                                      outfile=CLINICAL_COMUT_FILE_LIST[[s]], mdt=get_master_dt(full=T), driver_maf=DRIVER_MAF, add_ighv_cols=F, retval='file')
    print(outfile)
  }

}


GENERATE_DRIVER_BED = F
if(GENERATE_DRIVER_BED){
  drivers_dt <- as.data.table(get_drivers(ret='union'))
  setnames(drivers_dt, 'V1', 'gene')
  drivers_bed = get_genome_annot('gene_bed', subset=drivers_dt$gene)
  drivers_bed$chr = factor(drivers_bed$chr, levels = c(1:10, 11:22, 'X', 'Y', 'M'), ordered=T)
  setorderv(drivers_bed, c('chr', 'start'))

  fwrite(drivers_bed, DRIVERS_BED_FILE, sep='\t')
}


PREPROCESS_COUNTS_AND_RATES=F
if(PREPROCESS_COUNTS_AND_RATES){
  comut.utils$get_patient_counts_and_rates(wes_counts_file=MUTSIG_PATIENT_COUNTS_AND_RATES_WES,
                                           wgs_counts_file=MUTSIG_PATIENT_COUNTS_AND_RATES_WGS,
                                           wes_blacklist=MUTSIG_PATIENT_COUNTS_AND_RATES_WES_BLACKLIST,
                                           id_map_file = MASTER_METADTA_TABLE,
                                           outfile=MUTSIG_PATIENT_COUNTS_AND_RATES_WES_AND_WGS, retval = NULL)

  comut.utils$get_patient_counts_and_rates(wes_counts_file=MUTSIG_PATIENT_COUNTS_AND_RATES_WES,
                                           wgs_counts_file=NULL,
                                           wes_blacklist=NULL,
                                           id_map_file = MASTER_METADTA_TABLE,
                                           outfile=MUTSIG_PATIENT_COUNTS_AND_RATES_WES_PROCESSED, retval = NULL)
}

get_mutsig_counts_and_rates <- function(opt='merged', pats=NULL){
  if(opt=='merged'){
    rates_dt = fread(MUTSIG_PATIENT_COUNTS_AND_RATES_WES_AND_WGS)
  }
  if(opt=='wes'){
    rates_dt = fread(MUTSIG_PATIENT_COUNTS_AND_RATES_WES_PROCESSED)
  }

  if(opt=='wgs'){
    rates_dt = fread(MUTSIG_PATIENT_COUNTS_AND_RATES_WGS)
    rates_dt = merge(rates_dt, master_dt_full[master_dt_full$wgs_pair_id!="", c('participant_id', 'wgs_tumor_sample_id'), with=F], by.x='name', by.y='wgs_tumor_sample_id')
  }

  if(!is.null(pats))
    rates_dt = rates_dt[rates_dt$participant_id %in% pats,]

  return(rates_dt)
}

get_epitypes <- function(fill.na='NA', spain_format=F, oakes_labs=c('LP', 'IP', 'HP'), spain_labs=c('n-CLL', 'i-CLL', 'm-CLL'), unclassified_relabel='NA', drop_unclassified=F){
  # old_annot_file = '/xchip/cga_home/bknisbac/CLL/resources/Puente_Nature_2015_epitype_as_LPIPHP.tsv'
  # epitype_dt = fread(old_annot_file)[,c('participant_id', 'epitype')]
  epitype_dt = get_meth_dt()[,c('participant_id', 'epitype')]

  if(!is.null(fill.na)){
    epitype_dt[is.na(epitype_dt$epitype) | epitype_dt$epitype=='' | epitype_dt$epitype=='-']$epitype = fill.na
    epitype_dt[is.na(epitype_dt$epitype)]$epitype = fill.na
  }

  if(!is.null(unclassified_relabel)){
    epitype_dt[epitype_dt$epitype=='unclassified']$epitype = unclassified_relabel
  }

  if(spain_format){
    for (i in 1:length(oakes_labs)){
      epitype_dt$epitype[epitype_dt$epitype==oakes_labs[i]] = spain_labs[i]
    }
  }

  if(drop_unclassified){
    epitype_dt = epitype_dt[!is.na(epitype) & !epitype%in%c('unclassified', '-', '', 'NA')]
  }

  return(epitype_dt)
}
epitypes_dt = get_epitypes()


#### utilities ####
focal_regex = '(del|amp)_([0-9]+)[pq][0-9]+'

### more functions
get_novel_drivers <- function(driver_sets='union', cna_sets='union', with_clumps=T, DRIVER_SYMBOL_UPDATE_PAIRS, ret='novelty_table', ret_colnames=c('event_id', 'is_novel')){

  novelty_dt = data.table()

  if(!is.null(driver_sets)){
    if(is.character(driver_sets) && driver_sets=='union'){
      driver_list = get_drivers_union(with_clumps=F) #clumps is added later, if applicable
    } else {
      driver_list = get_drivers(sets=driver_sets, with_clumps=F, ret='union') #clumps is added later, if applicable
    }
    dm_dt = get_driver_master_dt()
    novelty_dt = rbindlist(list(novelty_dt, dm_dt[dm_dt$gene %in% driver_list, c('gene', 'is_novel'),with=F]))

    if(with_clumps){
      clumps_novelty_dt = fread(CLUMPS_DRIVERS_NOVELTY_TABLE)
      novelty_dt = rbindlist(list(novelty_dt, clumps_novelty_dt))
    }

    novelty_dt = unique(novelty_dt)
    setnames(novelty_dt, ret_colnames)
  }

  if(!is.null(cna_sets)){
    if(cna_sets=='union' || 'union'%in%c(cna_sets)){
      cna_nov_dt = rbindlist(list(cna_master_all_dt[,c('event_id', 'Reference'),with=F], cna_master_mvu_dt[,c('event_id', 'Reference'),with=F]))
    } else if (tolower(cna_sets)=='MvU' || 'mvu'%in%tolower(c(cna_sets))){
      cna_nov_dt = cna_master_mvu_dt[,c('event_id', 'Reference'),with=F]
    } else if (tolower(cna_sets)=='all' || 'all'%in%tolower(c(cna_sets))){
      cna_nov_dt = cna_master_all_dt[,c('event_id', 'Reference'),with=F]
    } else {
      print('Bad get_cna_events() option, returning union!')
      cna_nov_dt = rbindlist(list(cna_master_all_dt[,c('event_id', 'Reference'),with=F], cna_master_mvu_dt[,c('event_id', 'Reference'),with=F]))
    }
    cna_nov_dt[['is_novel']] = ifelse(cna_nov_dt[['Reference']]=='', 1, 0)
    cna_nov_dt = cna_nov_dt[,c('event_id', 'is_novel'),with=F]
    cna_nov_dt = unique(cna_nov_dt)
    setnames(cna_nov_dt, ret_colnames)

    #retain events of interest
    cna_events = get_cna_events(sets=cna_sets)
    cna_nov_dt = cna_nov_dt[event_id %in% cna_events]

    novelty_dt = rbindlist(list(novelty_dt, cna_nov_dt))
  }

  ## prevent duplicates
  novelty_dt = unique(novelty_dt)

  if(ret=='novelty_table'){
    setorderv(novelty_dt, ret_colnames[1])
    return(novelty_dt)
  } else if(ret=='novel'){
    return(novelty_dt[novelty_dt[[ret_colnames[2]]]==1][[ret_colnames[1]]])
  } else if(ret=='known'){
    return(novelty_dt[novelty_dt[[ret_colnames[2]]]==0][[ret_colnames[1]]])
  } else {
    return(novelty_dt)
  }
}


get_pats <- function(d='any', s='any'){ #tn='any'
  mdtf = copy(master_dt_full)

  if(d=='wes' && s=='mcll'){
    mdtf = mdtf[set_mcll_wes==1,]
  }
  if(d=='wes' && s=='ucll'){
    mdtf = mdtf[set_ucll_wes==1,]
  }

  if(d=='wgs' && s=='mcll'){
    mdtf = mdtf[wgs_tumor_sample_id!='' & IGHV_mut_freeze=='mutated',]
  }
  if(d=='wgs' && s=='ucll'){
    mdtf = mdtf[wgs_tumor_sample_id!='' & IGHV_mut_freeze=='unmutated',]
  }

  return(mdtf$participant_id)
}

get_ighv <- function(id_col='participant_id', fill.na=NULL, fill.na2=NULL, fill.na.all=NULL, drop.na.mut=F, drop.na.pc=F, ret='all'){
  ighv_master_cols = c('IGHV_mut_freeze', 'IGHV_identity_freeze')
  ighv_ret_cols = c('IGHV_mut', 'IGHV_identity')
  ret_dt = copy(master_dt_full[,c(id_col, ighv_master_cols),with=F])
  ret_dt = ret_dt[!is.na(ret_dt[[id_col]]) & ret_dt[[id_col]]!='']
  setnames(ret_dt, ighv_master_cols, ighv_ret_cols)

  if(!is.null(fill.na.all)){
    fill.na=fill.na.all
    fill.na2=fill.na.all
  }

  if(drop.na.mut){
    ret_dt = ret_dt[!(is.na(ret_dt[[ighv_ret_cols[1]]]) | ret_dt[[ighv_ret_cols[1]]]=='')]
  }

  if(drop.na.pc){
    ret_dt = ret_dt[!(is.na(ret_dt[[ighv_ret_cols[1]]]) | ret_dt[[ighv_ret_cols[1]]]=='')]
  }

  if(!is.null(fill.na)){
    ret_dt[[ighv_ret_cols[1]]][is.na(ret_dt[[ighv_ret_cols[1]]]) | ret_dt[[ighv_ret_cols[1]]]==''] = fill.na
  }

  if(!is.null(fill.na2)){
    ret_dt[[ighv_ret_cols[2]]][is.na(ret_dt[[ighv_ret_cols[2]]]) | ret_dt[[ighv_ret_cols[2]]]==''] = fill.na
  }

  if(ret=='mut'){
    return(ret_dt[,c('participant_id', 'IGHV_mut'),with=F])
  } else {
    return(ret_dt)
  }
}



####### RNA and expression clusters ######
EC_SAMPLE_SET = 'core2_tumors'
EC_RUN_SUFFIX = 'bc9_v1.5'
EC_RUN_SUBSET_ROOTDIR = 'subsets'
EC_RUN_SUBSET = 'clean1'
EC_N = 8
EC_LIST = paste0('EC', 1:EC_N) #Original labels in input files are EC1, EC2, ..., EC8 respective to EC_NAMES
EC_NAMES = c('EC-m1', 'EC-u1', 'EC-m2', 'EC-o', 'EC-u2', 'EC-m3', 'EC-m4', 'EC-i')

EC_RUN_SUBSET_DIR_STR = ifelse(is.null(EC_RUN_SUBSET), '', paste0('/', EC_RUN_SUBSET_ROOTDIR, '/', EC_RUN_SUBSET)) #DERIVED
EC_RUN_SUBSET_SUFFIX = ifelse(is.null(EC_RUN_SUBSET), '', paste0('_', EC_RUN_SUBSET))
EC_BNMF_DIR = paste0('/Volumes/xchip_cga_home/bknisbac/CLL/results_cll1085/sample_sets/', EC_SAMPLE_SET, '/', 'bayesnmf_tpm_', EC_RUN_SUFFIX, EC_RUN_SUBSET_DIR_STR)
EC_BNMF_FILE = paste0(EC_BNMF_DIR, '/', 'g.Bayes.txt')

EC_BNMF_DOWNSAMP_SIZES = seq(300,600,20)
EC_BNMF_DOWNSAMP_RUNS_PER_SIZE = 10

MARKER_GENES_DIR = paste0('/xchip/cga_home/bknisbac/CLL/results_cll1085/sample_sets/',EC_SAMPLE_SET,'/bayesnmf_tpm_',EC_RUN_SUFFIX, EC_RUN_SUBSET_DIR_STR,'/marker_genes_tpms_', EC_RUN_SUFFIX)
MARKER_GENES_DOWNSAMP_DIR = paste0(MARKER_GENES_DIR, '/', 'parallel_marker_bnmf_runs_ds0.8')
EC_N_SAMPS_FILE = '/xchip/cga_home/bknisbac/CLL/results_cll1085/tables/n_samps_per_EC_H3K27ac_bc9_v1.5_20210707.tsv'
GENE_COUNTS_FILE_FOR_MARKERS = paste0('/Volumes/xchip_cga_home/bknisbac/CLL/results_cll1085/sample_sets/core2_tumors/core2_tumors.rnaseqc_tpm_deseqLog10_', EC_RUN_SUFFIX,'.txt.gz')
GENE_TYPES_TO_RETAIN_FOR_MARKER_BNMF = c('protein_coding')


PEER_FACTORS_FILE = '/xchip/cga_home/bknisbac/CLL/results_cll1085/sample_sets/core2_tumors/peers/core2_tumors_peer25.PEER_covariates.txt'
BATCH_AND_COVARIATE_ANNOTATION_FILE = paste0('/Volumes/xchip_cga_home/bknisbac/CLL/results_cll1085/tables/batch_and_covariate_annotation.', EC_SAMPLE_SET, '.', EC_RUN_SUFFIX, '.forClustering.txt')
CUSTOM_DE_COMPARISONS_FILE = paste0(GITHUB_DIR, '/', 'Rscripts/Sources/CLL1085/custom_de_comparisons_',EC_RUN_SUFFIX,'.txt')

RAW_GENE_COUNTS_FILE = '/Volumes/xchip_cga_home/bknisbac/CLL/results_cll1085/sample_sets/all_rnaseqc_20200506/all_rnaseqc_20200506.rnaseqc_counts.gct.gz'
RAW_GENE_TPMS_FILE = '/Volumes/xchip_cga_home/bknisbac/CLL/results_cll1085/sample_sets/all_rnaseqc_20200506/all_rnaseqc_20200506.rnaseqc_tpm.gct.gz'
RAW_EXON_COUNTS_PARQUET = '/xchip/cga_home/bknisbac/CLL/results_cll1085/sample_sets/all_tumors_20201216/all_tumors_20201216.all_tumors_20201216.exon_reads.parquet'

#Use all samples - primary DE analysis
EC_DE_GENES_DIR = paste0('/xchip/cga_home/bknisbac/CLL/results_cll1085/sample_sets/',EC_SAMPLE_SET,'/degs_',EC_RUN_SUFFIX, ifelse(is.null(EC_RUN_SUBSET), '', paste0('_', EC_RUN_SUBSET)), '_original_counts')
EC_DE_GENES_FILE = paste0(EC_DE_GENES_DIR, '/eachvsall_w_covariates/de_genes_limma_voom_results.txt')
#Use only train samples (for classifier)
EC_DE_GENES_TRAIN_DIR = paste0('/xchip/cga_home/bknisbac/CLL/results_cll1085/sample_sets/',EC_SAMPLE_SET,'/degs_',EC_RUN_SUFFIX, ifelse(is.null(EC_RUN_SUBSET), '', paste0('_', EC_RUN_SUBSET)), '_original_counts')
EC_DE_GENES_TRAIN_FILE = paste0(EC_DE_GENES_TRAIN_DIR, '/custom_comparisons/de_genes_limma_voom_results.txt')

EC_PLOTTING = list()
EC_PLOTTING[['marker_plot_priority']] = '/xchip/cga_home/bknisbac/CLL/results_cll1085/tables/ec_markers_heatmap_priority_bc9_v1.5.tsv'
EC_PLOTTING[['pathway_plot_priority']] = '/xchip/cga_home/bknisbac/CLL/results_cll1085/tables/ec_pathways_heatmap_priority_bc9_v1.5_clean1.20210707.tsv'
EC_PLOTTING[['pathway_fgsea_ec']] = paste0(MARKER_GENES_DIR, '/', 'fgsea_res_per_ec.tsv')
EC_PLOTTING[['pathway_fgsea_ec_sig']] = paste0(MARKER_GENES_DIR, '/', 'fgsea_res_per_ec_sig.tsv')
EC_MARKERS = list()
EC_MARKERS[['markers_dir']] = MARKER_GENES_DIR
EC_MARKERS_UP_FILE = paste0(MARKER_GENES_DIR,'/markers_up_as_geneset_list.txt')
EC_MARKERS_UP_FILE_DS_FILTERED = paste0(MARKER_GENES_DIR,'/markers_up_as_geneset_list.dsFiltered.txt')
EC_MARKERS_UP_FILE_DS_BLACKLIST = paste0(MARKER_GENES_DIR,'/markers_up.ds_filter_blacklist.txt')
EC_MARKERS[['up_list']] = rnaseq.utils$read_gmt(EC_MARKERS_UP_FILE_DS_FILTERED, file_format = 'flat')
EC_MARKERS[['up_list2']] = EC_MARKERS[['up_list']]
names(EC_MARKERS[['up_list2']]) = sapply(names(EC_MARKERS[['up_list']]), function(x){str_extract(x, '(EC\\d+)')})
EC_MARKERS_DN_FILE = paste0(MARKER_GENES_DIR,'/markers_dn_as_geneset_list.txt')
EC_MARKERS_DN_FILE_DS_FILTERED = paste0(MARKER_GENES_DIR,'/markers_dn_as_geneset_list.dsFiltered.txt')
EC_MARKERS_DN_FILE_DS_BLACKLIST = paste0(MARKER_GENES_DIR,'/markers_dn.ds_filter_blacklist.txt')
EC_MARKERS[['down_list']] = rnaseq.utils$read_gmt(EC_MARKERS_DN_FILE_DS_FILTERED, file_format = 'flat')
EC_MARKERS[['down_list2']] = EC_MARKERS[['down_list']]
names(EC_MARKERS[['down_list2']]) = sapply(names(EC_MARKERS[['down_list']]), function(x){str_extract(x, '(EC\\d+)')})
EC_MARKERS[['all']] = sort(unique(c(as.vector(unlist(EC_MARKERS[['up_list']])), as.vector(unlist(EC_MARKERS[['down_list']])))))

EC_MARKER_DOWNASMP_N_ITER_CUTOFF = 80
if(file.exists(EC_MARKERS_UP_FILE_DS_BLACKLIST)){
  EC_MARKERS[['up_blacklist']] = fread(EC_MARKERS_UP_FILE_DS_BLACKLIST, header=F)$V1
}
if(file.exists(EC_MARKERS_DN_FILE_DS_BLACKLIST)){
  EC_MARKERS[['down_blacklist']] = fread(EC_MARKERS_DN_FILE_DS_BLACKLIST, header=F)$V1
}

EC_MARKERS[['all_annot_table']] = paste0(EC_MARKERS[['markers_dir']], '/', 'ec_marker_annotations_', EC_RUN_SUFFIX, '.txt')

#major GSEA file
EC_GENESET_ANNOTATION_FILE = paste0(EC_MARKERS[['markers_dir']], '/', 'ec_geneset_annotations_', EC_RUN_SUFFIX, EC_RUN_SUBSET_SUFFIX, '.tsv')


EC_MODEL_SELECTION_TIMESTAMP = '20210712'
EC_TRAIN_SET_TIMESTAMP = '20210707'
EC_TRAIN_SET_SAMPLE_LIST = paste0(MARKER_GENES_DIR,'/predictor/train_sample_ids.',EC_TRAIN_SET_TIMESTAMP,'.txt')
EC_TEST_SET_SAMPLE_LIST = paste0(MARKER_GENES_DIR,'/predictor/test_sample_ids.',EC_TRAIN_SET_TIMESTAMP,'.txt')

get_ecs <- function(ret='dt', retcols='default', samp_id_ret='sample_id', as_factor=T){
  if(ret=='dt'){
    if(retcols=='default')
      retcols = c('participant_id', 'rna_tumor_sample_id', 'expression_cluster')
    ec_dt = copy(master_dt_full[master_dt_full$expression_cluster!="", retcols, with=F])
    ec_dt$ec_label = paste0('EC', ec_dt$expression_cluster)
    if(!is.null(samp_id_ret))
      setnames(ec_dt, 'rna_tumor_sample_id', 'sample_id')
    if(as_factor)
      ec_dt$ec_label = factor(ec_dt$ec_label)
    return(ec_dt)
  } else if(ret=='g.bayes.dt'){
    ec_dt = fread(EC_BNMF_FILE)
    ec_dt$ec_label = paste0('EC', ec_dt$bnmf_expression_cluster)
    ec_dt$ec_label = factor(ec_dt$ec_label)
    if(ret=='dt')
      return(ec_dt)
  }
}

get_ec_markers <- function(opt=NULL, ec_set=EC_RUN_SUFFIX, s='all', v=NULL){
  if(ec_set==EC_RUN_SUFFIX && (is.null(opt) || opt=='umap') ){
    marker_list = unique(fread(MARKER_GENE_FILES[['all_umap']], header=F)[[1]])
  } else if(grepl('^EC\\d+$', opt)){ #e.g. EC1
    marker_list = c(EC_MARKERS[['up_list2']][[opt]], EC_MARKERS[['down_list2']][[opt]])
  } else if(grepl('^EC\\d+_up$', opt)){ #e.g. EC1_up
    marker_list = EC_MARKERS[['up_list2']][[str_extract(opt, 'EC\\d+')]]
  } else if(grepl('^EC\\d+_(down|dn)$', opt)){ #e.g. EC1_down or EC1_dn
    marker_list = EC_MARKERS[['down_list2']][[str_extract(opt, 'EC\\d+')]]
  }
  return(marker_list)
}

ec_lab_convert <- function(x, from=EC_LIST, to=EC_NAMES, na.warn=T, by_regex=F){
  res = to[match(x, from)]
  if(by_regex){
    for (i in 1:length(from)){
      x = gsub(from[i], to[i], x)
    }
    res = x
  }

  if(any(is.na(res)) & na.warn)
    warning('ec_lab_to_name function produced NAs')
  return(res)
}

get_ec_preds_args <- function(opt='tpm_genes',cv_folds=5){
  #This version with "normal"
  if(opt=='tpm_genes'){ #TPMs
    #initial setup
    ec_pred <- new.env()
    ec_pred$pred_rootdir = '/xchip/cga_home/bknisbac/CLL/results_cll1085/ec_predictor/alex_ng_revision_20210712/basic_tpms'
    ec_pred$pred_outdir = '/xchip/cga_home/bknisbac/CLL/results_cll1085/ec_predictor/alex_ng_revision_20210712/basic_tpms_bk_results'
    ec_pred$output_timestamp = EC_MODEL_SELECTION_TIMESTAMP
    ec_pred$gs_sizes = c(1:10, 20, 50)
    ec_pred$geneset_subdir_suffix = 'basic'
    #derived based on best model selection
    ec_pred$pred_dir = '/xchip/cga_home/bknisbac/CLL/results_cll1085/ec_predictor/alex_ng_revision_20210712/basic_tpms/results/topBot_7_a_basic'
    ec_pred$pred_full_file_regex = paste0(ec_pred$pred_dir, '/', 'rf2_ADASYN_split_[0-9]_ADASYN_full_predict_fold_[0-9].csv')
    ec_pred$pred_soft_clust_file_regex = paste0(ec_pred$pred_dir, '/soft_clust/', 'rf2_ADASYN_split_[0-9]_ADASYN_*.csv')
    ec_pred$pred_model = 'rf2_ADASYN'
    ec_pred$gene_list = paste0(ec_pred$pred_rootdir, '/', 'genes_for_predictor_topBot_7_a.20210707.txt')
  } else if(opt=='bc_tpm_genes'){ #Batch-corrected TPMs
    #initial setup
    ec_pred <- new.env()
    ec_pred$pred_rootdir = '/xchip/cga_home/bknisbac/CLL/results_cll1085/ec_predictor/alex_ng_revision_20210712/batchcorrected'
    ec_pred$pred_outdir = '/xchip/cga_home/bknisbac/CLL/results_cll1085/ec_predictor/alex_ng_revision_20210712/batchcorrected_bk_results'
    ec_pred$output_timestamp = EC_MODEL_SELECTION_TIMESTAMP
    ec_pred$gs_sizes = c(1:10)
    ec_pred$geneset_subdir_suffix = 'batchcorrected'

    #derived based on best model selection
    ec_pred$pred_dir = '/xchip/cga_home/bknisbac/CLL/results_cll1085/ec_predictor/alex_ng_revision_20210712/batchcorrected/results/topBot_9_a_batchcorrected'
    ec_pred$pred_full_file_regex = paste0(ec_pred$pred_dir, '/', 'rf_SMOTE_split_[0-9]_SMOTE_full_predict_fold_[0-9].csv')
    ec_pred$pred_soft_clust_file_regex = paste0(ec_pred$pred_dir, '/soft_clust/', 'rf_SMOTE_split_[0-9]_SMOTE_*.csv')
    ec_pred$pred_model = 'rf_SMOTE'
    ec_pred$gene_list = paste0(ec_pred$pred_rootdir, '/', 'genes_for_predictor_topBot_9_a.20210707.txt')
  } else {
   warning('Bad option for opt in get_ec_preds_args(), returning NULL')
    return(NULL)
  }

  ## overridable defaults
  ec_pred$cv_folds = cv_folds

  ## Derived outputs
  #The for-entropy file is changed in script that generates it due to multiple options/file-suffixes
  #New format:  paste0(ec_pred$pred_outdir, '/', ec_pred$pred_model, '_longitudinal_sample_table_for_entropy.',run_mode,'.', ec_pred$output_timestamp,'.csv')
  # ec_pred$output_for_entropy = paste0(ec_pred$pred_outdir, '/', ec_pred$pred_model, '_longitudinal_sample_table_for_entropy',ec_pred$output_timestamp,'.csv')

  ec_pred$output_preds_full = paste0(ec_pred$pred_outdir, '/', ec_pred$pred_model, '.preds_full.',ec_pred$output_timestamp,'.csv')
  ec_pred$output_preds_softclust = paste0(ec_pred$pred_outdir, '/', ec_pred$pred_model, '.preds_softclust.',ec_pred$output_timestamp,'.csv')
  dir.create(ec_pred$pred_outdir, showWarnings=F)
  return(ec_pred)
}

parse_ec_pred_files <- function(files_regex=NULL, file=NULL, file_list=NULL, parse_opt='full'){
  default_cols = c('sample_id', EC_LIST, 'ec_pred')
  if(parse_opt=='full'){
    my_colnames = default_cols
  } else if(parse_opt=='train_test'){
    my_colnames = c('sample_id', EC_LIST, 'diff_top2','ec_pred', 'ec_real', 'test_sample')
  } else {
    my_colnames = default_cols
  }

  if(!is.null(files_regex)){ #regex mode
    file_list = Sys.glob(files_regex)
  } else if (!is.null(file)){ #file list mode
    file_list = list(file)
  } else { #for completeness
    file_list = file_list
  }

  if(length(file_list)==1){ #one file
    pred_dt = fread(file_list[[1]])
    setnames(pred_dt, my_colnames)
    setorder(pred_dt, sample_id)
  } else { #multiple files
    pred_dt_list = list()
    for (i in 1:length(file_list)){
      pred_dt0 = fread(file_list[i])
      setnames(pred_dt0, my_colnames)
      setorder(pred_dt0, sample_id)
      pred_dt0$cv_fold = paste0('fold', i)
      pred_dt_list[[length(pred_dt_list)+1]] = pred_dt0
    }
    pred_dt = rbindlist(pred_dt_list)
  }

  return(pred_dt)
}

get_ec_preds <- function(opt='tpm_genes', file_type='full_soft', sample_regex = '.*', samp_blacklist=c(), filter_opt='master_tumors_mainCohorts', include_dietrich=T, ret_mode='fold_sum_max_dt'){
  ## Get config
  ec_pred = get_ec_preds_args(opt=opt)

  ## Get soft clustering,
  if(ret_mode=='soft_clust_train_test'){
    pred_dt = parse_ec_pred_files(files_regex=ec_pred$pred_soft_clust_file_regex, parse_opt='train_test')
    return(pred_dt)
  }

  ## Compute EC preds
  pred_dt = parse_ec_pred_files(files_regex=ec_pred$pred_full_file_regex)
  pred_dt = pred_dt[!sample_id %in% samp_blacklist,]
  pred_dt = pred_dt[grep(sample_regex, pred_dt$sample_id),]

  pred_l = melt.data.table(pred_dt, id.vars = c('sample_id', 'cv_fold'), measure.vars = EC_LIST, variable.name='EC', value.name = 'P')

  pred_l_sum = pred_l[,sum(P)/ec_pred$cv_folds, by=c('EC', 'sample_id')]
  setnames(pred_l_sum, c('EC', 'V1'), c('ec_pred', 'ec_pred_p'))

  pred_summax = pred_l_sum[order(sample_id, -ec_pred_p)][!duplicated(sample_id)][,c('sample_id', 'ec_pred', 'ec_pred_p')]

  if(ret_mode=='fold_sum_max_dt'){
    pred_dt = pred_summax
  } else if(ret_mode=='fold_sum_dt'){
    pred_dt = pred_l_sum
  } else if(ret_mode=='fold_dt_long'){
    setnames(pred_l, c('EC', 'P'), c('ec_pred', 'ec_pred_p'))
    pred_dt = pred_l
  }

  if(filter_opt=='master_tumors_mainCohorts'){
    #the mode of get_master_dt is to include all 610 ECs but not cohorts like Beekman and duplicates
    pred_dt = merge(pred_dt, get_master_dt(drop_redund_of_crc = F, tumors_only = T)[,c('participant_id', 'rna_tumor_sample_id')], by.x='sample_id', by.y='rna_tumor_sample_id', all.x=T)
    if(include_dietrich){
      pred_dt = pred_dt[!is.na(participant_id) | grepl('^H\\d+',sample_id)]
      pred_dt[grepl('^H\\d+',sample_id)]$participant_id = pred_dt[grepl('^H\\d+',sample_id)]$sample_id
    } else {
      pred_dt = pred_dt[!is.na(participant_id)]
    }
  } else { #don't filter anything - will complete participant_id annotation from file
    #Note: can't just use the samp-pat file because in master_dt I made some manual changess for samp-pat associations
    pred_dt = merge(pred_dt, master_dt_full[,c('participant_id', 'rna_tumor_sample_id')], by.x='sample_id', by.y='rna_tumor_sample_id', all.x=T)
    tmp_has = pred_dt[!is.na(participant_id)]
    tmp_missing = pred_dt[is.na(participant_id)]
    dim_pre = dim(tmp_missing)[1]
    tmp_missing = merge(tmp_missing[,-c('participant_id')], fread(RNA_WORKSPACES_SAMPLE_TO_PARTICIPANT_FILE), by='sample_id')
    dim_post = dim(tmp_missing)[1]
    if(dim_pre!=dim_post)
      warning(paste0('You are missing',dim_pre-dim_post,' samples in RNA_WORKSPACES_SAMPLE_TO_PARTICIPANT_FILE that are being dropped'))
    pred_dt = rbindlist(list(tmp_has, tmp_missing))
  }
  return(pred_dt)
}

EC_PREDS = list()
if(!is.null(get_ec_preds_args(opt='tpm_genes')$pred_full_file_regex)){ #Only if model already selected and variables set up for getting results from best result
  EC_PREDS[['tpm_genes_best']] = get_ec_preds(opt='tpm_genes', filter_opt='master_tumors_mainCohorts', ret_mode='fold_sum_max_dt')
  EC_PREDS[['tpm_genes_full']] = get_ec_preds(opt='tpm_genes', filter_opt='', ret_mode='fold_sum_dt')
  EC_PREDS[['tpm_genes_full_best']] = get_ec_preds(opt='tpm_genes', filter_opt='', ret_mode='fold_sum_max_dt')
  EC_PREDS[['tpm_genes_full_per_fold']] = get_ec_preds(opt='tpm_genes', filter_opt='', ret_mode='fold_dt_long')
  EC_PREDS[['tpm_genes_soft_clust_train_test']] = get_ec_preds(opt='tpm_genes', ret_mode='soft_clust_train_test') #no filtering applied regardless
}

EC_MCLL_DISCORDANT=c('EC2', 'EC5') #excludes EC-o and EC-i
EC_UCLL_DISCORDANT=c('EC1', 'EC3', 'EC6', 'EC7') #excludes EC-o and EC-i

get_ec_data <- function(marker_dir=MARKER_GENES_DIR, bnmf_dir=EC_BNMF_DIR, k_for_get_bnmf=length(EC_LIST),
                        bnmf_rdata_file=NULL,
                        transpose_H=T, lab_to_name=F, data_type='W_norm_up',
                        blacklist_for_W=F, ret='dt'){

  if(data_type=='best_k'){
    bnmf_all_file = paste0(bnmf_dir, '/res.L1EU.Bayes.all.txt')
    best_k = fread(bnmf_all_file)[is_max_freq==1]$K[1]
    return(best_k)
  }

  if(data_type=='best_full' || data_type=='best_res_run' || data_type=='best_rdata_file'){
    best.run = rnaseq$get_bnmf_results(bnmf_dir, K=k_for_get_bnmf, retval = 'best_res_run')
    best.rdata.file = paste(bnmf_dir, '/', paste("res.L1EU.Bayes",best.run,"RData",sep="."),sep="")

    if(data_type=='best_res_run'){
      best.run = rnaseq$get_bnmf_results(bnmf_dir, K=k_for_get_bnmf, retval = 'best_res_run')
      return(best.run)
    }
    if(data_type=='best_rdata_file'){
      return(best.rdata.file)
    }
    if(data_type=='best_full'){
      load(file=best.rdata.file)
      return(res.Bayes)
    }

  }

  if(data_type=='samps_to_ec'){
    samp_to_ec_dt = fread(paste0(bnmf_dir, '/g.Bayes.txt'))
    return(samp_to_ec_dt)
  }

  if(data_type=='H' || data_type=='H_norm'){
    if(is.null(bnmf_rdata_file)){ # Get BNMF data from best run (default mode)
      best.run = rnaseq$get_bnmf_results(bnmf_dir, K=k_for_get_bnmf, retval = 'best_res_run')
      load(file=paste(bnmf_dir, '/', paste("res.L1EU.Bayes",best.run,"RData",sep="."),sep=""))
    } else { # Get specific BNMF run data, as specified
      load(file=bnmf_rdata_file)
    }
    H = res.Bayes[[2]]
    H <- H[rowSums(H)!=0,]
    if(data_type=='H_norm'){
      H.norm <- apply(H,2,function(x) x/sum(x)) #H is now like H.norm in BNMF script
      H = H.norm
    }

    if(transpose_H || ret=='dt'){
      H = t(H)
    }
    #set if to use labels or names
    if(lab_to_name){
      ec_l = EC_NAMES
    } else {
      ec_l = EC_LIST
    }

    if(ret=='dt'){ #return data.table
      tmp_rownames = row.names(H)
      max_ec_arr = ec_l[apply(H, 1, which.max)]
      h_dt = cbind(data.table(tmp_rownames), data.table(H), max_ec_arr)
      setnames(h_dt, c('sample_id', ec_l[1:ncol(H)], ifelse(lab_to_name, 'ec_name', 'ec_label'))) #ec_l[1:ncol(H)] added for downsampling result fetching
      return(h_dt)
    } else { #return matrix
      if(transpose_H){
        names(H) = ec_l
      } else {
        row.names(H) = ec_l
      }
      return(H)
    }
  }

  if(data_type=='W_norm_up' || data_type=='W_norm_dn'){
    w_marker_norm_file = paste0(marker_dir, '/', 'W_markers_norm',ifelse(data_type=='W_norm_dn', '_down', ''),'.txt') #file without suffix is just up

    W_marker_norm = suppressWarnings(fread(w_marker_norm_file))
    ### Filter by blacklists
    if(blacklist_for_W){
      if(grepl('up', data_type, ignore.case = T)){
        mrk_blacklist = fread(EC_MARKERS_UP_FILE_DS_BLACKLIST, header=F)$V1
      } else {
        mrk_blacklist = fread(EC_MARKERS_DN_FILE_DS_BLACKLIST, header=F)$V1
      }
      W_marker_norm = W_marker_norm[!V1 %in% mrk_blacklist]
    }

    if(ret=='dt'){ #return dt
      setnames(W_marker_norm, 'V1', 'gene')
      return(W_marker_norm)
    } else { #return matrix
      W = as.matrix(W_marker_norm[,2:ncol(W_marker_norm),with=F])
      row.names(W) = W_marker_norm$V1
      return(W)
    }
  }
  if(data_type=='N_pats'){
    return(fread(EC_N_SAMPS_FILE))
  }

  if(data_type=='full_preds'){
    return(fread(EC_PREDICTIONS_TPM_RAW_FULL_FILE))
  }

}

generate_ec_and_h3k27ac_n_samps_file <- function(outfile=EC_N_SAMPS_FILE, dry_run=F){
  ac_samp_dt = fread(H3K27AC_SAMPLE_FILE)
  ac_N_dt = fread(H3K27AC_SAMPLE_FILE)[,c('EC'),with=F][,.N,by='EC']
  setnames(ac_N_dt, c('ec_label', 'H3K27ac_N'))
  ec_ac_samp_dt = merge(get_ecs()[,.N,by='ec_label'], ac_N_dt, by='ec_label')
  if(dry_run==T){
    return(ec_ac_samp_dt)
  }
  fwrite(ec_ac_samp_dt, EC_N_SAMPS_FILE, sep='\t')
}

GENERATE_EC_N_SAMPS_FILE = F
if(GENERATE_EC_N_SAMPS_FILE){
  generate_ec_and_h3k27ac_n_samps_file()
}

#Used originally just for Fig 3b (ECs, markers, H3k27ac support), but moved here for use in stats (checking effect of FC cutoff)
get_acetylation_data_per_markers <- function(ac_up_file=H3K27AC_PER_MARKERS_UP_FILE,
                                             ac_dn_file=H3K27AC_PER_MARKERS_DOWN_FILE,
                                             fc_cutoff=NA,
                                             ecs_exclude=NA,
                                             fill.na='NA', reverse_test=F, ret='default'){
  ## This rev test is to see if H3K27ac signal exists in reverse direction of expected
  #Answer was that there were 0 zero events!
  up_down_labs = c('UP', 'DOWN')
  if(reverse_test){
    up_down_labs = rev(up_down_labs)
    #Note: no need to change fc_cutoff to 1/fc_cutoff because the filter is on the original EC and not the H3K27ac result
  }

  ac_up = fread(ac_up_file)
  #Rename markers so Up-markers are identifiable
  ac_up$Gene.annot.Gene.EC.Up.Down.Broad = sub('MARKERS$', 'MARKERS_UP', ac_up$Gene.annot.Gene.EC.Up.Down.Broad)
  if(!is.na(fc_cutoff)){ #not used; added for checking difference in revision
    ac_up = ac_up[Peak.annot.log2FoldChange >= log2(fc_cutoff)]
  }
  ac_up = ac_up[grepl(up_down_labs[1], ac_up[['Gene.annot.Gene.EC.Up.Down.Broad']])] #drop originally down markers
  ac_up$orig_ec = sub('.*_(EC\\d+)_.*', '\\1', ac_up[['Gene.annot.Gene.EC.Up.Down.Broad']])
  setnames(ac_up, c('Peak.annot.EC.specific', 'Gene.annot.Gene.hgnc_symbol'), c('ac_ec', 'gene'))
  ac_up = unique(ac_up[,c('gene', 'ac_ec', 'orig_ec')][ac_ec==orig_ec])
  ac_up[['H3K27ac']] = 'Up'

  ac_dn = fread(ac_dn_file)

  #Rename markers so Up-markers are identifiable
  ac_dn$Gene.annot.Gene.EC.Up.Down.Broad = sub('MARKERS$', 'MARKERS_UP', ac_dn$Gene.annot.Gene.EC.Up.Down.Broad)
  if(!is.na(fc_cutoff)){ #not used; added for checking difference in revision
    ac_dn = ac_dn[Peak.annot.log2FoldChange <= log2(1/fc_cutoff)]
  }
  ac_dn = ac_dn[grepl(up_down_labs[2], ac_dn[['Gene.annot.Gene.EC.Up.Down.Broad']])] #drop originally up markers
  ac_dn$orig_ec = sub('.*_(EC\\d+)_.*', '\\1', ac_dn[['Gene.annot.Gene.EC.Up.Down.Broad']])
  setnames(ac_dn, c('Peak.annot.EC.specific', 'Gene.annot.Gene.hgnc_symbol'), c('ac_ec', 'gene'))
  ac_dn = ac_dn[,c('gene', 'ac_ec', 'orig_ec')]
  ac_dn = unique(ac_dn[,c('gene', 'ac_ec', 'orig_ec')][ac_ec==orig_ec])
  ac_dn[['H3K27ac']] = 'Down'

  #Use W to get order correct as in figure
  w_up_dt = get_ec_data(data_type='W_norm_up', blacklist_for_W=T, ret='dt')
  w_up_dt[['orig_ec']] = paste0('EC', apply(w_up_dt[,grep('EC', colnames(w_up_dt)),with=F], 1, which.max))
  w_up_dt = merge(w_up_dt[,c('gene', 'orig_ec')], ac_up, all.x=T, by=c('gene', 'orig_ec'), sort=F)

  w_dn_dt = get_ec_data(data_type='W_norm_dn', blacklist_for_W=T, ret='dt')
  w_dn_dt[['orig_ec']] = paste0('EC', apply(w_dn_dt[,grep('EC', colnames(w_dn_dt)),with=F], 1, which.max))
  w_dn_dt = merge(w_dn_dt[,c('gene', 'orig_ec')], ac_dn, all.x=T, by=c('gene', 'orig_ec'), sort=F)

  if(!any(is.na(ecs_exclude))){
    w_up_dt = w_up_dt[!orig_ec %in% ecs_exclude]
    w_dn_dt = w_up_dt[!orig_ec %in% ecs_exclude]
  }

  ac_dt = rbindlist(list(w_up_dt, w_dn_dt))[,c('gene', 'H3K27ac')]

  if(ret=='ac_dt'){
    return(ac_dt)
  }
  return(ac_dt[['H3K27ac']])
}


get_ec_discordance <- function(disc_by='ighv', clean=T, ec_drop=c('EC4', 'EC8')){ #disc_by='ighv'
  if (!disc_by %in% c('ighv', 'epitype', 'both')) { #check input
    warning('bad disc_by input, returning NULL')
    return(NULL)
  }

  #get ECs (drop those not of interest)
  ec_dt = get_ecs()
  ec_dt = ec_dt[!ec_label %in% ec_drop]

  #get ighv
  if(disc_by=='ighv' || disc_by=='both'){
    ec_dt = merge(ec_dt, get_ighv(drop.na.mut=T)[,-c('IGHV_identity')])
  }

  #get epitype
  if(disc_by=='epitype' || disc_by=='both')
    ec_dt = merge(ec_dt, get_epitypes(drop_unclassified = T))

  if(disc_by=='both')
    ec_dt = ec_dt[(epitype=='m-CLL' & IGHV_mut=='mutated') | (epitype=='n-CLL' & IGHV_mut=='unmutated')]

  if(disc_by=='ighv' || disc_by=='both'){
    ec_dt$discordant = 'Concordant'
    ec_dt[IGHV_mut=='mutated' & ec_label %in% EC_MCLL_DISCORDANT]$discordant = 'Discordant'
    ec_dt[IGHV_mut=='unmutated' & ec_label %in% EC_UCLL_DISCORDANT]$discordant = 'Discordant'
  } else if(disc_by=='epitype') {
    ec_dt$discordant = 'Concordant'
    ec_dt[epitype=='m-CLL' & ec_label %in% EC_MCLL_DISCORDANT]$discordant = 'Discordant'
    ec_dt[epitype=='n-CLL' & ec_label %in% EC_UCLL_DISCORDANT]$discordant = 'Discordant'
  }

  if(clean){
    ec_dt = ec_dt[!participant_id %in% master_dt_full[clinical_blacklist_pat==1]$participant_id]
  }

  return(ec_dt)
}


get_clinical_data <- function(opt='default', sset='all', clincols_save=CLINDATA_COLS, clin_blacklist_pat=get_pat_set('clin_blacklist_pat')){

  clindata = fread(CLINICAL_DATA_FULL_FILE)[,clincols_save,with=F]

  if(sum(table(clindata$patid)>1)){
    print(paste('Warning: found ', sum(table(clindata$patid)>1), 'duplicates in clinical data:') )
    print(names(table(clindata$patid)[table(clindata$patid)>1]))
    print('Selecting max years_ddx_to_status, if applicable')
    clindata = clindata[order(-years_ddx_to_status)]
    clindata = unique(clindata, by='patid')
    clindata = clindata[order(patid)]
  }
  if('patid' %in% colnames(clindata) && !'participant_id' %in% colnames(clindata))
    setnames(clindata, 'patid', 'participant_id')
  clindata$status_text = clindata$status #save original column
  #Be careful: need to deal with NAs and also make sure there's no 'Dead' or 'Alive', only lower case
  clindata$status_text = tolower(clindata$status_text)
  clindata$status = ifelse(clindata$status_text=='dead', 1, 0) #prepare status col as in example
  clindata$status[clindata$status_text==''] = NA #very important

  ## Add variables
  clindata$years_dsamp_to_status = clindata$time_dsamp_to_status / 365
  clindata$years_ddx_to_status = clindata$time_ddx_to_status / 365
  clindata$years_dsamp_to_d1trt = clindata$time_dsamp_to_d1trt / 365

  if(!is.null(clin_blacklist_pat) && !is.na(clin_blacklist_pat)){
    clindata = clindata[!participant_id %in% clin_blacklist_pat]
  }

  if(sset=='rna_discovery'){
    clindata = merge(clindata, get_ecs(), by='participant_id')
  }

  return(clindata)
}


if(F){
  ### Driver novelty table - union SNV(All/M/U) / union CNA(All/MvU)
  fwrite(get_novel_drivers(), EVENT_NOVELTY_TABLE, sep='\t')
}

#### Multiomics ###
## Multiomic curveball (EC-centric, including only samples with DNA, IGHV, ECs & epitype)
assemble_multiomic_comut <- function(genomic_comut_file=CURVEBALL_COMUT_FILE_LIST[['mvu_union_genomic']],
                                     ec_mode='discovery_cohort',
                                     epitype_mode='classified',
                                     outfile=NULL,
                                     add_u1=T,
                                     add_iglv321_mut=F,
                                     iglv321_na_mode=0,
                                     drop_zero_sum_cols=F,
                                     dash_replace=NULL,
                                     ret='file'){
  library(fastDummies)
  #get genonomics
  mo_cm_dt = fread(genomic_comut_file)

  #get ECs
  if(ec_mode=='discovery_cohort'){
    my_ecs_dt = get_ecs()[,c('participant_id', 'expression_cluster')]
    setnames(my_ecs_dt, 'expression_cluster', 'EC')
    my_ecs_dt = fastDummies::dummy_cols(my_ecs_dt, select_columns = "EC")[,-c('EC')]
    setnames(my_ecs_dt, sub('EC_', 'EC', colnames(my_ecs_dt)))
    mo_cm_dt = merge(mo_cm_dt, my_ecs_dt, by='participant_id')
  }

  #get U1
  if(add_u1)
    mo_cm_dt = merge(mo_cm_dt, get_U1_muts()[,c('participant_id', 'U1')], by='participant_id')

  #Add IGLV3-21
  if(add_iglv321_mut){ #even if the function level fill.na is not NULL, this is handled later
    if('IGLV321_R110' %in% colnames(mo_cm_dt)){
      print('Dropping existing IGLV321_R110 from comut')
      mo_cm_dt = mo_cm_dt[,-c('IGLV321_R110'),with=F]
    }
    if(iglv321_na_mode=='unk_column'){ #add binary column for Unknown
      mo_cm_dt = merge(mo_cm_dt, get_iglv321_mut(opt='pat_to_binary_mut', fill.na=0), by='participant_id', all.x=T)
      mo_cm_dt = merge(mo_cm_dt, get_iglv321_mut(opt='unknown_as_binary'), by='participant_id', all.x=T)
    } else {
      mo_cm_dt = merge(mo_cm_dt, get_iglv321_mut(opt='pat_to_binary_mut', fill.na=iglv321_na_mode), by='participant_id', all.x=T)
    }
  }

  #get epitypes
  if(epitype_mode=='classified'){
    my_epi_dt = get_epitypes(drop_unclassified=T)
    my_epi_dt = fastDummies::dummy_cols(my_epi_dt, select_columns = "epitype")[,-c('epitype')]
    setnames(my_epi_dt, sub('epitype_', '', colnames(my_epi_dt)))
    mo_cm_dt = merge(mo_cm_dt, my_epi_dt, by='participant_id')
  }

  if(drop_zero_sum_cols){
    cols_to_drop = colnames(mo_cm_dt[,-c('participant_id')])[colSums(mo_cm_dt[,-c('participant_id')])==0]
    print(paste('Dropping zero-sum columns:', paste(cols_to_drop, collapse = ' ')))
    mo_cm_dt = mo_cm_dt[, -cols_to_drop, with=F]
  }

  if(!is.null(dash_replace)){
    colnames(mo_cm_dt) = sub('-', dash_replace, colnames(mo_cm_dt))
  }

  if(!is.null(outfile)){
    fwrite(mo_cm_dt, outfile, sep='\t')
  }

  if(!is.null(ret)){
    if(ret=='file')
      return(outfile)
    if(ret=='dt')
      return(mo_cm_dt)
  }
}

GENERATE_MULTIOMIC_ANALYSIS_COMUTS = F
if(GENERATE_MULTIOMIC_ANALYSIS_COMUTS){
  #multiomic for curveball
  outfile = assemble_multiomic_comut(genomic_comut_file=CURVEBALL_COMUT_FILE_LIST[['mvu_union_genomic']],
                                     ec_mode='discovery_cohort',
                                     epitype_mode='classified',
                                     add_u1=T,
                                     # dash_replace='_', #Note: no dash replace here. Different from curveball comut
                                     outfile=CURVEBALL_COMUT_FILE_LIST[['mvu_union_multiomic']],
                                     ret='file')
  print(outfile)

  #multiomic for clinical analysis
  outfile = assemble_multiomic_comut(genomic_comut_file = CLINICAL_COMUT_FILE_LIST[['all']],
                                     ec_mode = 'discovery_cohort',
                                     epitype_mode = 'classified',
                                     add_u1 = T,
                                     add_iglv321_mut = T,
                                     iglv321_na_mode = 'unk_column',
                                     outfile = CLINICAL_COMUT_FILE_LIST[['all_multiomic']],
                                     drop_zero_sum_cols = T,
                                     dash_replace='_',
                                     ret = 'file')
  print(outfile)

}




process_multiomic_curveball_results <- function(multiomic_co_p=CURVEBALL_RESULT_FILES[['multiomic_co_p']],
                                                multiomic_comut_file=CURVEBALL_COMUT_FILE_LIST[['mvu_union_multiomic']],
                                                include_ighv=T, include_epitype=T,
                                                fdr_cutoff=0.1,
                                                min_n_cutoff=5,
                                                ret='mo_sig_dt'){

  ## helper function
  get_ec_val <- function(x,y,regex='^EC\\d+$', invert=F){
    if(grepl('^EC\\d+$', x)){
      return(ifelse(!invert, x, y))
    } else {
      return(ifelse(!invert, y, x))
    }
  }

  ### Get N events per pair ###
  cmo_full = fread(multiomic_comut_file)
  cmo_full = melt.data.table(cmo_full, id.vars = c('participant_id', grep('^EC\\d+$', colnames(cmo_full), value=T)), variable.name = 'event', value.name = 'event_binary', value.factor = F)[event_binary==1]
  cmo_full = melt.data.table(cmo_full[,-c('event_binary'),with=F], id.vars = c('participant_id', 'event'), variable.name='ec_label', value.name='ec_binary')[ec_binary==1]
  cmo_N = cmo_full[,.N,by=c('event', 'ec_label')]

  if(!is.null(min_n_cutoff)){
    cmo_N = cmo_N[N >= min_n_cutoff]
  }

  ### Get multiomic curveball p.values ###
  mo_cb_p = fread(multiomic_co_p)
  mo_p_long = melt.data.table(mo_cb_p, id.vars = 'V1')
  setnames(mo_p_long, c('entity1', 'entity2', 'p.value'))
  mo_p_long$entity2 = as.character(mo_p_long$entity2)
  mo_p_long = mo_p_long[entity1!=entity2]
  make_key <- function(a,b){paste0(sort(c(a,b)), collapse = '__')}
  mo_p_long$key = mo_p_long[, mapply( make_key, entity1, entity2 ) ]
  mo_p_long = mo_p_long[duplicated(key)][,-c('key'),with=F] #drop lower triangle (all are 1)
  mo_p_long = mo_p_long[(grepl('^EC\\d+$', mo_p_long$entity1) | grepl('^EC\\d+$', mo_p_long$entity2)) &
                          !(grepl('^EC\\d+$', mo_p_long$entity1) & grepl('^EC\\d+$', mo_p_long$entity2))] #retain only EC vs non-EC pairs
  if(!include_ighv){
    mo_p_long = mo_p_long[!(grepl('^(MCLL|UCLL)$', mo_p_long$entity1) | grepl('^(MCLL|UCLL)$', mo_p_long$entity2))]
  }
  if(!include_epitype){
    mo_p_long = mo_p_long[!(grepl('^[imn]-CLL$', mo_p_long$entity1) | grepl('^[imn]-CLL$', mo_p_long$entity2))]
  }

  mo_p_long$ec_label = mo_p_long[,mapply(get_ec_val, entity1, entity2)]
  mo_p_long$event = mo_p_long[,mapply(get_ec_val, entity1, entity2, invert=T)]
  mo_p_long = mo_p_long[,-c('entity1','entity2'), with=F]
  mo_p_long = merge(mo_p_long, cmo_N, by=c('ec_label', 'event'))
  if(!is.null(min_n_cutoff)){  #filter by a minimum of N events to include in FDR (already applied above to cmo_N, so merge actually did this already)
    mo_p_long = mo_p_long[N >= min_n_cutoff]
  }
  mo_p_long$p.adj = p.adjust(mo_p_long$p.value)
  mo_sig_dt = mo_p_long[p.adj < fdr_cutoff]

  if(ret=='mo_pair_n'){
    return(cmo_N)
  }
  if(ret=='mo_sig_dt'){
    return(mo_sig_dt)
  }
  if(ret=='mo_sig_list'){
    return(mo_sig_event_list)
  }
}


comut_to_group_freq <- function(cm, group_col='participant_id'){
  cm_freq = melt(cm, id.vars = group_col)[,sum(value), by=group_col]
  setnames(cm_freq, 'V1', 'freq')
  return(cm_freq)
}



##### project-wide COLORS utility
get_colors <- function(x){
  color_list = list(mcll='#361379', ucll='#e65100',
                    LP='#006E93', IP='#FDC010', HP='#963736',
                    oxphos='#00A651',
                    novel_driver='#cc1d7c')
  color_list[['n-CLL']] = color_list[['LP']]
  color_list[['i-CLL']] = color_list[['IP']]
  color_list[['m-CLL']] = color_list[['HP']]

  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#CC79A7", "#924900",
                 "#490092", "#920000", "#24ff24", "#000000", "#db6d00", "#006ddb")

  color_list[['ecs']] = cbPalette[1:EC_N]

  safe_cbpal <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
  color_list[['safe_cbpal']] = safe_cbpal

  if(x[1]=='mu')#shorthand
    x=c('mcll', 'ucll')

  if(is.null(x) || x=='all'){
    return(color_list)
  } else {
    color_arr = c()
    for(myx in x){
      color_arr = c(color_arr, color_list[[myx]])
      if(any(is.null(color_arr)))
        warning('Warning: error fetching project colors!')
    }
    return(color_arr)
  }
}
mu_colors = get_colors('mu')


#### Mutation categories (used for icomut plot in CLLmap data portal)
MUTCATEG_FILES = list()
MUTCATEG_FILES[['wgs_maf']] = '/xchip/cga_home/bknisbac/CLL/results_cll1085/dna/wgs/mutsig2cv_final_set_20200630/call-tool_mutsig2cv_hg19/final_set_20200630.final_analysis_set.maf'
MUTCATEG_FILES[['wgs_mutcateg']] = '/xchip/cga_home/bknisbac/CLL/results_cll1085/dna/wgs/mutsig2cv_final_set_20200630/call-tool_mutsig2cv_hg19/mutcategs.txt'
MUTCATEG_FILES[['wes_maf']] = '/xchip/cga_home/bknisbac/CLL/results_cll1085/dna/mutsig/mutsig_20200703/all_files_for_mutcateg/WES_pairs_rerun_final_20200703.final_analysis_set.maf'
MUTCATEG_FILES[['wes_mutcateg']] = '/xchip/cga_home/bknisbac/CLL/results_cll1085/dna/mutsig/mutsig_20200703/all_files_for_mutcateg/mutcategs.txt'

get_mut_categs <- function(s='wgs', mutcateg_file=NULL, maf_file=NULL, prefer='wes'){ #s='wes'; #'wes_wgs'
  if(!is.null(mutcateg_file) && !is.null(maf_file)){
    ms_maf = fread(maf_file)[,c('Tumor_Sample_Barcode', 'categ_idx'),with=F][categ_idx!=0]

    ms_maf = merge(ms_maf, master_dt_full[,c('participant_id', paste0(s, '_tumor_sample_id')) ,with=F],
                    by.x='Tumor_Sample_Barcode', by.y=paste0(s, '_tumor_sample_id'))
    mut_categs = fread(mutcateg_file)
    mut_categs$categ_idx = 1:nrow(mut_categs)
    categ_dt = merge(ms_maf, mut_categs[,c('categ_idx', 'name'),with=F], by='categ_idx')
    categ_n_dt = categ_dt[,.N,by=c('participant_id', 'name')]
    categ_n_dt[,tot := sum(N),by='participant_id']
    categ_n_dt[,pc := round(N/tot*100,1)]
    categ_pc_dt = dcast.data.table(categ_n_dt, formula = 'participant_id ~ name', value.var = 'pc', fill=0)
    setnames(categ_pc_dt, paste0('MutsigCateg_', colnames(categ_pc_dt)))
    setnames(categ_pc_dt, 'MutsigCateg_participant_id', 'participant_id')
    return(categ_pc_dt)
  }
  if(s=='wes'){
    wes_mutcateg_dt = get_mut_categs(s='wes', mutcateg_file = MUTCATEG_FILES[['wes_mutcateg']], maf_file = MUTCATEG_FILES[['wes_maf']])
    return(wes_mutcateg_dt)
  }
  if(s=='wgs'){
    wgs_mutcateg_dt = get_mut_categs(s='wgs', mutcateg_file = MUTCATEG_FILES[['wgs_mutcateg']], maf_file = MUTCATEG_FILES[['wgs_maf']])
    return(wgs_mutcateg_dt)
  }
  if(s=='wes_wgs'){
    wes_mutcateg_dt = get_mut_categs(s='wes', mutcateg_file = MUTCATEG_FILES[['wes_mutcateg']], maf_file = MUTCATEG_FILES[['wes_maf']])
    wgs_mutcateg_dt = get_mut_categs(s='wgs', mutcateg_file = MUTCATEG_FILES[['wgs_mutcateg']], maf_file = MUTCATEG_FILES[['wgs_maf']])
    if(!is.null(prefer) && prefer=='wes'){ #prefer WES over WGS entries if two for same patient exist
      wgs_mutcateg_dt = wgs_mutcateg_dt[!participant_id %in% wes_mutcateg_dt$participant_id]
    } else { #return both
      wes_mutcateg_dt$dtype = 'WES'
      wgs_mutcateg_dt$dtype = 'WGS'
    }
    wes_wgs_mutcateg_dt = rbindlist(list(wes_mutcateg_dt, wgs_mutcateg_dt), fill=T)
    return(wes_wgs_mutcateg_dt)
  }
}
