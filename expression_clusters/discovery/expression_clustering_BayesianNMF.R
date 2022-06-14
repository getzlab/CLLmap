### Comments
# This script is the main script for Expression cluster discovery
# Bayesian NMF core developed by Jaegil Kim for "Robertson et al. Comprehensive Molecular Characterization of Muscle-Invasive Bladder Cancer. Cell 171, 540–556.e25 (2017)."

### Usage from commandline
### Specific downsampling (enables commandline parallel)
#Rscript ~/github/Rscripts/rnaseq/expression_clustering_BayesianNMF.R --run bc9_v1.5 --subset clean1 --downsamp_size 300 --downsamp_threads 30 > ~/clean1_and_downsamp.300.log 2>&1 &
#Rscript ~/github/Rscripts/rnaseq/expression_clustering_BayesianNMF.R --run bc9_v1.5 --subset clean1 --downsamp_size 300 --downsamp_iters 30 --downsamp_threads 10 > ~/clean1_and_downsamp.300.full.log 2>&1 &
#for i in 300 320 340 360 380 400 420 440 460 480 500 520 540 560 580 600 ; do Rscript ~/github/Rscripts/rnaseq/expression_clustering_BayesianNMF.R --run bc9_v1.5 --subset clean1 --downsamp_size $i --downsamp_threads 10 --downsamp_runs 10 --downsamp_iters 30 > ~/clean1_and_downsamp.$i.log & done > ~/downsamp.parallel.loop.log 2>&1 &

#### Load packages ####
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(clusterSim))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(reshape))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(FactoMineR)) #Added for PCA analysis
suppressPackageStartupMessages(library(explor))
suppressPackageStartupMessages(library(Hmisc)) # for capitalize()

#### Project sources ####
#Scripts and source files
GITHUB_DIR = '/Volumes/xchip_cga_home/bknisbac/github'
RSCRIPTS_DIR = paste0(GITHUB_DIR, '/', 'Rscripts')
suppressPackageStartupMessages(source(paste0(RSCRIPTS_DIR, '/', 'Sources/graphics.src.R')))
suppressPackageStartupMessages(source(paste0(RSCRIPTS_DIR, '/', 'Sources/parallel.src.R')))
<<<<<<< HEAD
suppressPackageStartupMessages(source(paste0(RSCRIPTS_DIR, '/', 'Sources//CLL1085/CLL1100_current.src.R'), proj<-new.env()))
suppressPackageStartupMessages(source(paste0(GITHUB_DIR, '/Rscripts/Sources/RNAseq.src.R'))) #For get_gene_counts()
=======
suppressPackageStartupMessages(source(paste0(RSCRIPTS_DIR, '/Sources/RNAseq.src.R'))) #For get_gene_counts()
>>>>>>> 5892dc2afe20ea72e74026403ee0dddf16754e3a

##### CONSTS ####
### Defining run - option exists to specify from command line
SAMPLE_SET = 'core2_tumors'

RUN_SUFFIX = 'bc9_v1.5' #cleaned by removing redundant samples
SUBSET_SET = 'clean1' #NULL

## BNMF
NUM_BAYESNMF_ITERS_SAMP_BY_SAMP_CLUSTERING = 1000 #1000 #100 only
NUM_PARALLEL_BNMF_ITERATIONS = 100 #95 used on big machine

#overridable downsamp args
DOWNSAMP_SIZE_SINGLE = NULL
DOWNSAMP_N_PARALLEL_BNMF_ITERS = 10
DOWNSAMP_RUNS_PER_SIZE = 10
DOWNSAMP_ITERS_PER_RUN = 30

## other downsamp args
DOWNSAMP_SIZES = seq(300,600,20)

##### CONSTS ####
#Parse commandline args
get_option_list <- function(config=NULL){
  option_list <- list(
    make_option(c("-r", "--run"), action="store", type="character",help="", metavar = 'RUN', default=RUN_SUFFIX),
    make_option(c("-s", "--subset"), action="store", type="character",help="", metavar = 'SUBSET_SET', default=SUBSET_SET),
    make_option(c("-d", "--downsamp_size"), action="store", type="numeric",help="", metavar = 'DOWNSAMP_SIZE_SINGLE', default=DOWNSAMP_SIZE_SINGLE),
    make_option(c("--downsamp_runs"), action="store", type="numeric",help="", metavar = 'DOWNSAMP_RUNS_PER_SIZE', default=DOWNSAMP_RUNS_PER_SIZE),
    make_option(c("--downsamp_iters"), action="store", type="numeric",help="", metavar = 'DOWNSAMP_ITERS_PER_RUN', default=DOWNSAMP_ITERS_PER_RUN),
    make_option(c("--downsamp_threads"), action="store", type="numeric",help="", metavar = 'DOWNSAMP_N_PARALLEL_BNMF_ITERS', default=DOWNSAMP_N_PARALLEL_BNMF_ITERS)
  )
  return(option_list)
}
option_list = get_option_list()
args <- parse_args(OptionParser(option_list=option_list, usage = "Rscript %prog [options]"), print_help_and_exit=F)
print(args)

#override consts by args, if applicable
RUN_SUFFIX = args$run
SUBSET_SET = args$subset
DOWNSAMP_SIZE_SINGLE = args$downsamp_size
DOWNSAMP_RUNS_PER_SIZE = args$downsamp_runs
DOWNSAMP_ITERS_PER_RUN = args$downsamp_iters
DOWNSAMP_N_PARALLEL_BNMF_ITERS = args$downsamp_threads

# Note path consts
subset_samples_for_bnmf <- function(ge, sample_set, run_suffix, subset_set='', gene_cols=c('gene'), downsamp_n=NULL, downsamp_seed=NULL, make_copy=T){
  if(make_copy){
    ge = copy(ge)
  }
  #sample whitelist and blacklist
  include_list = NULL
  exclude_list = NULL
  if(run_suffix=='bc9_v1.5' && subset_set=='clean1'){
    blacklist_file = paste0('/xchip/cga_home/bknisbac/CLL/results_cll1085/sample_sets/',sample_set, '/bayesnmf_tpm_', run_suffix,'/sample_blacklists/sample_blacklist_', run_suffix, '_', subset_set,'.tsv')
    exclude_list = fread(blacklist_file, header=F)[['V1']]
  }

  if(!is.null(include_list)){
    ge = ge[,c(gene_cols, include_list),with=F]
  }
  if(!is.null(exclude_list)){
    exclude_list = exclude_list[exclude_list %in% colnames(ge)]
    if(length(exclude_list)>0)
      ge = ge[,-exclude_list,with=F]
  }

  ### Down-sampling
  if(!is.null(downsamp_n)){
    if(!is.null(downsamp_seed)){
      set.seed(downsamp_seed)
    }
    downsamp_samp_list = sample(colnames(ge)[!colnames(ge) %in% gene_cols], downsamp_n)
    ge = ge[,c(gene_cols, downsamp_samp_list),with=F]
  }

  return(ge)
}

##################
#####  MAIN  #####
##################
set.seed(1)

### Plotting vars and options ####
MATRIX_PLOT_SIZE = 48
H_MATRIX_PLOT_WIDTH = 48

color.axis="grey" # color of sample names from H martix plot
scale <- 0.8
.theme_ss <- theme_bw(base_size=14) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size=10*scale, family="mono"),
        axis.text.y = element_text(hjust = 0.5,size=10*scale, family="mono"),
        axis.title.x = element_text(face="bold",colour="black",size=14*scale),
        axis.title.y = element_text(face="bold",colour="black",size=14*scale),
        axis.text = element_text(size = 16*scale, family = "mono"),
        strip.text = element_text(lineheight=0.5),
        strip.text.x = element_text(size=10*scale,face='bold',angle=00),
        strip.text.y = element_text(size=10*scale,face="bold"),
        strip.background = element_rect(colour="black",fill="gray85"),
        panel.margin = unit(0.20,"lines"),
        plot.title=element_text(lineheight=1.0,face="bold",size=12*scale))

#colorblind palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#CC79A7", "#924900", "#490092", "#920000", "#24ff24", "#000000", "#db6d00", "#006ddb")

#### User-specific setup and input files ####
#Goal: have a TPM file ready for ananlysis
TUMOR_TYPE <- 'CLL'

#Input and output files
RESDIR_SAMPLE_SETS = '/Volumes/xchip_cga_home/bknisbac/CLL/results_cll1085/sample_sets'
GENES_GTF_FILE = '/Volumes/xchip_cga_home/bknisbac/resources/gencode19_noChrPrefix_mitoMT.gtf.gz'

GTF_COLNAMES = c('seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute')

GENE_COUNTS_FILE = paste0('/Volumes/xchip_cga_home/bknisbac/CLL/results_cll1085/sample_sets/', SAMPLE_SET, '/', SAMPLE_SET,'.rnaseqc_tpm_deseqLog10_', RUN_SUFFIX, '.txt.gz')

### Important filters and outdir specification
PARTICIPANT_INCLUSION_FILE1 = NULL
PARTICIPANT_INCLUSION_FILE1_MODE = 'participant_no_posttreatment'
SAMPLE_INCLUSION_REGEX = '.*' #'.*' to disable; Previously used temporarily to restrict analysis

SUBSETS_ROOT_SUBDIR = 'subsets'
BNMF_EXPRESSION_CLUSTERING_SUBDIR = paste0('bayesnmf_tpm_', RUN_SUFFIX,
                                           ifelse(is.null(SUBSET_SET), '', paste0('/', SUBSETS_ROOT_SUBDIR, '/', SUBSET_SET)))

MASTER_ID_MAPPING_TABLE = '/Volumes/xchip_cga_home/bknisbac/CLL/results_cll1085/tables/cll1085_master_table_primary_dataset.20200223.txt'
id_map = fread(MASTER_ID_MAPPING_TABLE)

TITLE_SUFFIX_PER_ANALYSIS = NULL #'GCLL'

### Run consts
PAIRWISE_CORRELATION_METHOD = 'pearson'
GENE_TYPES_TO_RETAIN = c('protein_coding', 'lincRNA')

DELETE_SUBOPTIMAL_L1EU_DATA = F #reduce storage footprint of runs

SKIP_BAYESNMF_RUN = F
RUN_MAIN_BNMF = T
RUN_DOWNSAMPLING = F #enables down-sampling for robustness analysis
SAVE_R_IMAGE_BEFORE_CLUSTERING = F
SAVE_R_IMAGE_AFTER_CLUSTERING = F

### Derived consts ###
#Gene counts file
if(is.null(GENE_COUNTS_FILE)){
  GENE_COUNTS_FILE = paste0('/Volumes/xchip_cga_home/bknisbac/CLL/results_cll1085/sample_sets', '/', SAMPLE_SET, '/', SAMPLE_SET, '.rnaseqc_tpm.gct.gz')
}

## Get gene counts
gene_counts = get_gene_counts(counts_file = GENE_COUNTS_FILE, genes_gtf_file=GENES_GTF_FILE, gtf_feature='gene', gene_types_to_retain=GENE_TYPES_TO_RETAIN,
                              return_gene_col='Description', unlog_base = 10, dup_gene_mode='sum', negatives_assign_val=0, deseq_norm=F,
                              log_norm_base=NULL, log_increment_val=1, return_class='data.table', rename_gene_col='gene')

if(!is.null(PARTICIPANT_INCLUSION_FILE1)){
  if(PARTICIPANT_INCLUSION_FILE1_MODE=='participant_no_posttreatment'){
    pt_inclusion_dt = fread(PARTICIPANT_INCLUSION_FILE1)
    pt_inclusion_dt = merge(pt_inclusion_dt, id_map[,c('participant_id', 'rna_tumor_sample_id')])
    rna_samp_id_to_include = pt_inclusion_dt[samp_post_treatment==0 & rna_tumor_sample_id != "",]$rna_tumor_sample_id
    gene_counts = gene_counts[,colnames(gene_counts)[colnames(gene_counts) %in% c('gene', rna_samp_id_to_include)],with=F]
  }
}

run_gexp_bnmf <- function(gcts, outdir,
                          max.K = 20, ### Hierarchical clustering iterations with varying K (# clusters to produce by cutting tree) up to max.K
                          pItem = 0.8, ### re-sampling rate = 80% = randomly choose 80% samples in each hclust iteration
                          innerLinkage = "average", #
                          finalLinkage = "ward.D", #Not used in current version of function
                          n.base = 250, # num of basic bootstrap iterations for given K in consensus matrix calculation
                          pairwise.cor.method = PAIRWISE_CORRELATION_METHOD,
                          n.iter = NUM_BAYESNMF_ITERS_SAMP_BY_SAMP_CLUSTERING,
                          num_par = NUM_PARALLEL_BNMF_ITERATIONS,
                          skip_bnmf_run = SKIP_BAYESNMF_RUN,
                          force_bnmf_override = F,
                          l1eu.num.iters=200000, #original defaults
                          l1eu.tol=1e-7, #original defaults
                          del_suboptimal_L1EU = DELETE_SUBOPTIMAL_L1EU_DATA,
                          matrix_plot_size = MATRIX_PLOT_SIZE,
                          save_image_before_clustering=SAVE_R_IMAGE_BEFORE_CLUSTERING,
                          save_image_after_clustering=SAVE_R_IMAGE_AFTER_CLUSTERING,
                          consensus_clustering_only=F){

  ### Functions (needed here for parallel to work)
  get_correlation_matrix <- function(counts, corr_method='spearman', outdir=".", plot_title=NULL, plot_title_suffix=NULL, width=24, height=24, file_suffix=NULL){
    #Set outfile names
    pdf_name = paste0("correlation_", corr_method, "_ordered_by_HC",ifelse(is.null(file_suffix), '', paste0('_', file_suffix)) ,".pdf")
    rdata_name = paste0("corr.all.", corr_method, ifelse(is.null(file_suffix), '', paste0('_', file_suffix)), ".RData")

    #Set plot title
    if(is.null(plot_title)){
      plot_title = paste0(capitalize(corr_method), ' correlation matrix', ifelse(is.null(plot_title_suffix), '', paste0(' - ', plot_title_suffix)))
    }

    if(!grepl("/$", outdir)){
      outdir = paste0(outdir,"/")
    }

    #Compute correlation
    corr_mat <- cor(counts, method=corr_method)

    #Generate heatmap PDF
    pdf(file=paste(outdir, pdf_name, sep=""), width=width, height=height)
    plot.heatmap.3(corr_mat,T,T,plot_title)
    dev.off()

    #Save RData
    save(corr_mat,file=paste(outdir, rdata_name, sep=""))

    return(corr_mat)
  }


  #Function: Select most probable cluster assignment per sample (uses the H matrix)
  get.consensus.clustering <- function(res,consensus.norm) {
    W <- res[[1]] #Not used in this function
    H <- res[[2]]
    W <- W[,colSums(W)!=0]
    H <- H[rowSums(H)!=0,]
    rownames(H) <- paste("G",seq(1:nrow(H)),sep="")
    H.norm <- apply(H,2,function(x) x/sum(x))
    g.Bayes <- apply(H.norm,2,function(x) which.max(x)) # select most probable cluster per sample
    K0 <- nrow(H.norm) #Get K* (probable number of clusters)

    x <- list(g.Bayes)
    return(x)
  }

  #Function: Bayesian NMF method using L1 regularization and eucledian distance metric (ref. PMID 23681989)
  #Arguments:
  # V0 - correlation-based consensus matrix
  # num.iters - Bayesian-NMF iters
  # a0 - hyperparameter
  BayesNMF.L1EU <- function(V0,num.iters=200000,a0=10,tol=1.e-07,K=20,K0=20,phi=1.0, verbose=TRUE) {
    eps <- 1.e-50
    del <- 1.0
    active_nodes <- colSums(V0) != 0 #only use pairs with non-0 association
    V0 <- V0[,active_nodes]
    V <- V0-min(V0) #transform to non-negative
    Vmin <- min(V)
    Vmax <- max(V)
    N <- dim(V)[1]
    M <- dim(V)[2]

    W <- matrix(runif(N * K)*Vmax,ncol=K)
    H <- matrix(runif(M * K)*Vmax,ncol=M)
    V.ap <- W%*%H+eps
    I <- array(1,dim=c(N,M))

    phi <- sd(V)^2*phi
    C <- N+M+a0+1
    b0 <- sqrt((a0-1)*(a0-2)*mean(V,na.rm=T)/K0)
    lambda.bound <- b0/C
    lambda <- (colSums(W)+rowSums(H)+b0)/C
    lambda.cut <- 1.5*lambda.bound

    n.like <- list()
    n.evid <- list()
    n.error <- list()
    n.lambda <- list()
    n.lambda[[1]] <- lambda
    iter <- 2
    count <- 1
    while (del >= tol & iter < num.iters) {
      #Update H and W seperately and approximate V (==V.ap)
      H <- H * (t(W) %*% V) / (t(W) %*% V.ap + phi * matrix(rep(1/lambda,M),ncol=M) + eps)
      V.ap <- W %*% H + eps
      W <- W * (V %*% t(H)) / (V.ap %*% t(H) + phi * t(matrix(rep(1/lambda,N),ncol=N)) + eps)
      V.ap <- W %*% H + eps
      #Caculate learning rate, delta and likelihood
      lambda <- (colSums(W)+rowSums(H)+b0)/C
      del <- max(abs(lambda-n.lambda[[iter-1]])/n.lambda[[iter-1]])
      like <- sum((V-V.ap)^2)/2
      #Store iteration results
      n.like[[iter]] <- like
      n.evid[[iter]] <- like + phi*sum((colSums(W)+rowSums(H)+b0)/lambda+C*log(lambda)) #evidence per iteration
      n.lambda[[iter]] <- lambda
      n.error[[iter]] <- sum((V-V.ap)^2)
      if (iter %% 100 == 0 && verbose==TRUE) {
        cat(iter,n.evid[[iter]],n.like[[iter]],n.error[[iter]],del,sum(colSums(W)!=0),sum(lambda>=lambda.cut),'\n')
      }
      iter <- iter+1
    }
    return(list(W,H,n.like,n.evid,n.lambda,n.error))
  }

  plot.heatmap.3 <- function(x,rowTF,colTF,main, ColSideColors=NULL, color_gradient=greenred(40)) {
    s1 <- 0.75
    s2 <- 1.0
    s3 <- 1.5
    mydist <- function(c) {dist(c,method="euclidean")}
    myclust <- function(c) {hclust(c,method="ward.D")}
    if(is.null(ColSideColors))
      heatmap.2(as.matrix(x), hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", dendrogram="both",margins=c(5,5),
                Rowv=rowTF, Colv=colTF, symbreaks=F, key=TRUE, symkey=F,main=main,
                density.info="none", trace="none",labCol=colnames(x),labRow=rownames(x),col=color_gradient,cex.lab=s1,cexRow=0.50,cexCol=0.50,keysize=s1)
    else
      heatmap.2(as.matrix(x), hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", dendrogram="both",margins=c(5,5),
                Rowv=rowTF, Colv=colTF, symbreaks=F, key=TRUE, symkey=F,main=main,
                density.info="none", trace="none",labCol=colnames(x),labRow=rownames(x),col=color_gradient,cex.lab=s1,cexRow=0.50,cexCol=0.50,keysize=s1,
                ColSideColors=ColSideColors)
  }

  #***To do: parallelize the bootstrapping process (array of K of len K*n.base)
  get.consensus.mat <- function(corr,max.K,pItem,innerLinkage,finalLinkage,n.base, verbose=TRUE) {
    min.K = 2 # const
    sample <- rownames(corr)
    n.sample <- ncol(corr)
    n.resample <- round(pItem*ncol(corr))
    consensus <- array(0,dim=c(n.sample,n.sample))
    n.count <- 0
    for (K in min.K:max.K) { #Bootstrap per range of K
      n.iter <- K*n.base
      for (i in 1:n.iter) { #num bootstraps per K
        index <- sample.int(n.sample,n.resample)
        corr.resample <- corr[index,index]
        dt <- as.dist(1-corr.resample)
        h.inner <- hclust(dt,method=innerLinkage) #hierarchical cluster
        g <- cutree(h.inner,K) #get clusters
        for (j in 1:K) { #increment all pairs per cluster by 1
          id <- match(names(g)[g==j],sample,nomatch=0)
          consensus[id,id] <- consensus[id,id]+1
        }
        n.count <- n.count+1
      }
      if (verbose==TRUE){
        cat(K,n.count,'\n')
      }
    }
    colnames(consensus) <- rownames(corr)
    rownames(consensus) <- rownames(corr)
    consensus.norm <- consensus/n.count # norm by number of total i-loop iterations: sum(min.K:max.K) * n.base
    return(list(consensus,consensus.norm))
  }

  get.sample.association.heatmap <- function(H,g.Bayes,scale) {
    g.ordering <- paste0("G", sort(unique(g.Bayes), decreasing = T))
    sample.ordering <- colnames(H)[order(g.Bayes,decreasing=F)]
    df <- t(H)
    df <- melt(df)
    colnames(df) <- c("sample","cluster","value")
    df$sample <- factor(df$sample,sample.ordering)
    df$cluster <- factor(df$cluster,g.ordering)
    df1 <- df
    df1[,"type"] <- "H matrix"

    H.norm <- apply(H,2,function(x) x/sum(x))
    df <- t(H.norm)
    df <- melt(df)
    colnames(df) <- c("sample","cluster","value")
    df$sample <- factor(df$sample,sample.ordering)
    df$cluster <- factor(df$cluster,g.ordering)
    df2 <- df
    df2[,"type"] <- "Normalized H"

    p = ggplot(df1,aes(x=sample,y=cluster,fill=value))+geom_tile() #geom_tile(colour="yellow")
    #p = p + facet_grid(. ~ type, scale = "free_y")
    p = p + scale_fill_gradient2(low="blue",mid="white",high ="darkblue",name=paste("Activity",sep=""))
    p = p + .theme_ss
    p = p + ggtitle("H matrix")
    p = p + xlab("Sample") + ylab("mRNA Clusters")
    p = p + theme(axis.title.x = element_text(face="bold",colour="black",size=14*scale))
    p = p + theme(axis.title.y = element_text(face="bold",colour="black",size=14*scale))
    p = p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size=5*scale, family="mono",face='bold',color=color.axis))
    p1 = p + theme(legend.position="top")

    p = ggplot(df2,aes(x=sample,y=cluster,fill=value))+geom_tile() #geom_tile(colour="yellow")
    #p = p + facet_grid(. ~ type, scale = "free_y")
    p = p + scale_fill_gradient2(low="blue",mid="white",high ="darkblue",name=paste("Activity",sep=""))
    p = p + .theme_ss
    p = p + ggtitle("Normalized H matrix")
    p = p + xlab("Sample") + ylab("mRNA Clusters")
    p = p + theme(axis.title.x = element_text(face="bold",colour="black",size=14*scale))
    p = p + theme(axis.title.y = element_text(face="bold",colour="black",size=14*scale))
    p = p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size=5*scale, family="mono",face='bold',color=color.axis))
    p2 = p + theme(legend.position="top")

    df <- data.frame(t(H))
    colnames(df) <- paste("H",seq(1:nrow(H)),sep="")
    df[,"mRNA"] <- g.Bayes
    df$mRNA <- paste("G",df$mRNA,sep="")
    df <- melt(df,id="mRNA")
    colnames(df) <- c("Subtype","cluster","Association")
    p = ggplot(df,aes(x=Association))
    p = p + geom_histogram(color="black",fill="gray75",binwidth=0.05)
    p = p + facet_grid(Subtype ~ cluster,scale='free_y')
    p = p + .theme_ss
    p = p + ggtitle("Sample Associations")
    p = p + xlab("Association") + ylab("Sample Counts")
    p = p + theme(axis.title.x = element_text(face="bold",colour="black",size=14*scale))
    p = p + theme(axis.title.y = element_text(face="bold",colour="black",size=14*scale))
    p = p + theme(axis.text.x = element_text(angle=0,vjust=0.0, size=12*scale, family="mono",face='bold',color=color.axis))
    p3 = p + theme(legend.position="top")
    return(list(p1,p2,p3))
  }


  ### Init
  gcts = copy(gcts)
  dir.create(outdir, recursive = T, showWarnings = F)
  cat('printing outdir in function', outdir, '\n')

  ### Feature (gene) selection based on expression centrality and dispersion across samples
  gexp = subset_gene_expression(gcts,
                                min_quantile_expression_mean=0, #previously 0.2
                                min_quantile_expression_sd=0,
                                min_quantile_expression_cv=0,
                                max_na_or_zero_samples_frac=0.95, #previously 0.1 (currently more permissive)
                                top_n_mode='cv',
                                top_n=2500)

  #This drops the gene columns (but they are not needed for the consensus clustering)
  #Placed after subsetting gene selection, to allow filtering based on raw values
  gexp = gexp[,lapply(.SD, function(x){log2(x+1)}) ,.SDcols=colnames(gexp)[sapply(gexp, is.numeric)]]

  ####### transform to fold-changed expression (in relation to median per gene)
  counts.fold <- t(apply(gexp[,sapply(gexp,is.numeric),with=F],1,function(x) x-median(x,na.rm=T)))

  ####### compute the spearman correlations among samples with selected genes above
  corr.all = get_correlation_matrix(counts.fold, corr_method=pairwise.cor.method, outdir=outdir, plot_title_suffix=TITLE_SUFFIX_PER_ANALYSIS)

  ###### Consensus matrix by bootstrapping ######
  consensus.mat.file = paste0(outdir, "consensus.mat.rds")
  consensus.norm.file = paste0(outdir, "consensus.norm.rds")
  if(!file.exists(consensus.mat.file) || !file.exists(consensus.norm.file)){
    res.consensus <- get.consensus.mat(corr.all,max.K,pItem,innerLinkage,finalLinkage,n.base)
    consensus.mat <- res.consensus[[1]]
    consensus.norm <- res.consensus[[2]]
    saveRDS(consensus.mat, consensus.mat.file)
    saveRDS(consensus.norm, consensus.norm.file)
  } else {
    consensus.mat = readRDS(consensus.mat.file)
    consensus.norm = readRDS(consensus.norm.file)
  }

  x <- consensus.mat
  pdf(file=paste(outdir,"consensus.ordered_by_HC.pdf",sep=""),width=matrix_plot_size, height= matrix_plot_size)
  plot.heatmap.3(x,T,T,"Consensus matrix by hierarchical clustering")
  dev.off()

  if(save_image_before_clustering){
    Rimage_before_bnmf = paste0(outdir, 'Rsession.image.beforeBNMF.RData')
    save.image(Rimage_before_bnmf)
  } #load(Rimage_before_bnmf)

  if(consensus_clustering_only){
    print('Exiting after consensus clustering.')
    return()
  }

  ############################################################
  ################## NMF - sample discovery
  ################## run Bayesian n.iter times
  ############################################################

  if(skip_bnmf_run){
    print('Skipping BayesNMF run for clustering, this setting should be used if all bnmf output files exist.')
  } else if(is.null(num_par) || num_par==1){
    for (i in 1:n.iter) {
      iter.outfile = paste(outdir,paste("res.L1EU.Bayes",i,"RData",sep="."),sep="")
      if(!file.exists(iter.outfile) || force_bnmf_override){
        print(paste('Starting iteration', i, 'of', n.iter))
        set.seed(i)
        res.Bayes <- BayesNMF.L1EU(as.matrix(consensus.norm),l1eu.num.iters,10,l1eu.tol,20,20,1.0)
        save(res.Bayes,file=iter.outfile)
      } else {
        print(paste('Skipping iteration', i, 'of', n.iter, 'for existing output', iter.outfile))
      }
    }
  } else { #Run in parallel
    parallel.registerCluster(num_par) #create and register parallel backend
    foreach(i=1:n.iter) %dopar% {
      iter.outfile = paste(outdir,paste("res.L1EU.Bayes",i,"RData",sep="."),sep="")
      if(!file.exists(iter.outfile) || force_bnmf_override){
        print(paste('Starting iteration', i, 'of', n.iter))
        set.seed(i)
        res.Bayes <- BayesNMF.L1EU(as.matrix(consensus.norm),l1eu.num.iters,10,l1eu.tol,20,20,1.0)
        save(res.Bayes,file=iter.outfile)
      } else {
        print(paste('Skipping iteration', i, 'of', n.iter, 'for existing output', iter.outfile))
      }
    }
    parallel.unregisterClusters() #unregister clusters upon activation
  }
  #n.iter=2
  tmpK <- rep(0,n.iter) #K per run
  tmpE <- rep(0,n.iter) #evid per run
  for (i in 1:n.iter) {
    bayesnmf_run_file = paste(outdir,paste("res.L1EU.Bayes",i,"RData",sep="."),sep="")
    if(file.exists(bayesnmf_run_file)){
      load(file=bayesnmf_run_file)
      lambda <- res.Bayes[[5]]
      lambda <- unlist(lambda[length(lambda)])
      lambda <- lambda-min(lambda)
      cat(lambda,sum(lambda!=0),'\n')
      tmpK[i] <- sum(lambda > 0)
      tmpE[i] <- res.Bayes[[4]][length(res.Bayes[[4]])]
    } else { #Missing (I previously had a failure on a specific iteration 1/1000, so added this)
      print(paste0('File does not exist on iteration ', i, ': ', bayesnmf_run_file))
      tmpK[i] <- NA
      tmpE[i] <- NA
    }

  }

  ##################################
  #### summary of BayesNMF runs ####
  ##################################
  ## Frequency of K values determined per run
  x <- table(tmpK)
  pdf(file=paste(outdir,paste("BayesNMF.freq.pdf",sep="."),sep=""),width=4,height=5)
  s1 <- 1.5
  s2 <- 2.0
  par(mfrow=c(1,1))
  par(mar=c(5,5,2,1))
  barplot(x,cex=s1,cex.axis=s1,cex.main=s1,cex.names=s1,cex.lab=s1,xlab="# of clusters",ylab="Freq.",main=paste(TUMOR_TYPE,sep="."))
  dev.off()

  #########################################
  ################ df has all info for BayesNMF runs and run # with the lowest "evid" for given K (== runK) must be selected
  #########################################
  ## Select run with lowest evid among runs with n_clusts equal to most frequent K (sort df decreasing by evid and select last)
  #to df
  df <- data.frame(seq(1:n.iter),tmpK,unlist(tmpE))
  df[df$tmpK,]
  colnames(df) <- c("run","K","evid")
  #get most frequent K
  k_freq = table(df$K)
  df$is_max_freq = ifelse(df$K %in% names(k_freq[k_freq==max(k_freq)]), 1,0)
  #get run with lowest evidence from most frequent K
  df <- df[order(df$is_max_freq, df$evid, decreasing=c(F,T)),]
  run.K = df[nrow(df),'run'] #automatic selection for the lowest scoring run #Note: sometimes there's a tie (thus selecting last 'randomly')
  write.table(df, paste0(outdir, '/', 'res.L1EU.Bayes.all.txt'), quote=F, row.names=F, sep='\t')


  ## Delete files from non-best runs
  if(del_suboptimal_L1EU){
    for (i in 1:n.iter){
      if (i != run.K){ #del all but best run
        file.remove(paste0(outdir,paste("res.L1EU.Bayes",i,"RData",sep=".")))
      }
    }
  }


  ## Reload data for best run
  load(file=paste(outdir,paste("res.L1EU.Bayes",run.K,"RData",sep="."),sep=""))
  res.clust <- get.consensus.clustering(res.Bayes,consensus.norm)
  g.Bayes <- res.clust[[1]]
  save(g.Bayes,file=paste0(outdir,"g.Bayes.RData")) ### cluster membership is saved here. The membership is determined by the maximum association of H
  write.table(data.table(sample_id = names(g.Bayes), bnmf_expression_cluster = g.Bayes), paste0(outdir,"g.Bayes.txt"), sep='\t', row.names=F, quote=F) #save as table
  ordering <- order(g.Bayes,decreasing=F)
  x <- corr.all
  rownames(x) <- names(g.Bayes)[ordering]
  colnames(x) <- names(g.Bayes)[ordering]
  pdf(file=paste(outdir,"correlation.ordered_by_bayesNMF.pdf",sep=""),width=matrix_plot_size,height= matrix_plot_size)
  plot.heatmap.3(x[ordering,ordering],F,F,paste(capitalize(pairwise.cor.method), "correlation matrix"))
  dev.off()

  ### Plot consensus matrix
  x <- consensus.mat
  ### By HC
  pdf(file=paste(outdir,"consensus.ordered_by_HC_withBNMFcolors.pdf",sep=""),width=matrix_plot_size,height= matrix_plot_size)
  plot.heatmap.3(x,T,T,"Consensus matrix by hierarchical clustering", ColSideColors=cbPalette[g.Bayes[match(names(g.Bayes), colnames(x))]])
  dev.off()

  ### By BNMF
  pdf(file=paste(outdir,"consensus.ordered_by_bayesNMF.pdf",sep=""),width=matrix_plot_size,height= matrix_plot_size)
  plot.heatmap.3(x[ordering,ordering],F,F,"Consensus matrix", ColSideColors=cbPalette[g.Bayes[ordering]], color_gradient = viridis(40))
  dev.off()

  ################## plotting H matrix
  H <- res.Bayes[[2]]
  H <- H[rowSums(H)!=0,]
  K0 <- nrow(H)
  rownames(H) <- paste("EC",seq(1:K0),sep="")
  H.norm <- apply(H,2,function(x) x/sum(x))

  res <- get.sample.association.heatmap(H, g.Bayes, scale=1.0)
  p1 <- res[[1]]
  p2 <- res[[2]]
  p3 <- res[[3]]
  pdf(file=paste(outdir,"H.ordered.pdf",sep=""), width=H_MATRIX_PLOT_WIDTH, height=4)
  plot(p1)
  dev.off()
  png(file=paste(outdir,"H.ordered.png",sep=""), width=H_MATRIX_PLOT_WIDTH, height=4)
  plot(p1)
  dev.off()
  pdf(file=paste(outdir,"H.norm.ordered.pdf",sep=""), width=H_MATRIX_PLOT_WIDTH, height=4)
  plot(p2)
  dev.off()
  png(file=paste(outdir,"H.norm.ordered.png",sep=""), width=H_MATRIX_PLOT_WIDTH, height=4)
  plot(p2)
  dev.off()

  ### Save image
  if(save_image_after_clustering){
    Rimage_bnmf_script_end = paste0(outdir, '/', 'Rsession.image.BNMFclusteringScriptEnd.RData')
    save.image(Rimage_bnmf_script_end)
  }

}


#### Subset samples by sample blacklist or whitelist ###
gene_counts = subset_samples_for_bnmf(gene_counts, SAMPLE_SET, RUN_SUFFIX, subset_set=SUBSET_SET)

#### Run BNMF ####
#Main output dir
if(RUN_MAIN_BNMF){
  outdir.main.bnmf <- paste(RESDIR_SAMPLE_SETS, '/', SAMPLE_SET , '/', BNMF_EXPRESSION_CLUSTERING_SUBDIR, "/",sep="")

  run_gexp_bnmf(gcts=gene_counts, outdir = outdir.main.bnmf,
                n.iter = NUM_BAYESNMF_ITERS_SAMP_BY_SAMP_CLUSTERING,
                num_par = NUM_PARALLEL_BNMF_ITERATIONS)
}


### Subset sample list if applicable for sets
if(RUN_DOWNSAMPLING){
  if(!is.null(DOWNSAMP_SIZE_SINGLE)){
    downsamp_sizes = c(DOWNSAMP_SIZE_SINGLE)
  } else {
    downsamp_sizes = DOWNSAMP_SIZES
  }

  total_ds_runs = length(downsamp_sizes) * DOWNSAMP_RUNS_PER_SIZE * DOWNSAMP_ITERS_PER_RUN
  cat('Running a total iterations for downsampling:', total_ds_runs, '\n')

  for (ds_size in downsamp_sizes){
    for (ds_run in 1:DOWNSAMP_RUNS_PER_SIZE){
      gene_counts_ds = subset_samples_for_bnmf(gene_counts, SAMPLE_SET, RUN_SUFFIX, subset_set=SUBSET_SET, downsamp_n = ds_size, downsamp_seed = ds_run)

      outdir_ds = paste0(RESDIR_SAMPLE_SETS, '/', SAMPLE_SET , '/', BNMF_EXPRESSION_CLUSTERING_SUBDIR, "/downsampling/", 'ds_', ds_size, '_', ds_run, '/')
      print(outdir_ds)
      run_gexp_bnmf(gcts=gene_counts_ds, outdir = outdir_ds,
                    n.iter = DOWNSAMP_ITERS_PER_RUN,
                    num_par = DOWNSAMP_N_PARALLEL_BNMF_ITERS,
                    consensus_clustering_only=F)
    }
  }
}
