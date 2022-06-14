####
# This script uses the H matrix from the expression cluster detection BayesNMF to identify marker genes associated with each cluster
# See Methods in CLLmap (Knisbacher et al, Nature Genetics, 2022) and original method publication by Jaegil Kim (Robertson et al, Cell, 2017)

# Rscript ~/github/Rscripts/rnaseq/get.marker_for_fold.ver3.R -r bc9_v1.5 -d CLL > ~/getmarker.cll.log 2>&1 &
# Rscript ~/github/Rscripts/rnaseq/get.marker_for_fold.ver3.R -r bc9_v1.5_MCLL -d MCLL > ~/getmarker.mcll.log 2>&1 &
# Rscript ~/github/Rscripts/rnaseq/get.marker_for_fold.ver3.R -r bc9_v1.5_UCLL -d UCLL > ~/getmarker.ucll.log 2>&1 &
# Rscript ~/github/Rscripts/rnaseq/get.marker_for_fold.ver3.R -r bc9_v1.5 -d CLL -s clean1 --bnmf_n_runs 100 --n_par 30 --ds_frac 0.8 > ~/get.markers.clean1.par.log 2>&1 &

### Important - for classic runs make sure these are NULL: SUBSET_SET, DOWNSAMP_FRAC

### Defining run - option exists to specify from command line
RUN_SUFFIX = 'bc9_v1.5' #bc9_v1.5
SUBSET_SET = NULL #'clean1' #NULL
DISEASE_LABEL = 'CLL' #UCLL #MCLL   #Define disease label (used for prefix in marker genesets file) #Only "CLL" was used in the CLLmap manuscript
TOTAL_MARKER_BNMF_RUNS = 1 #100 #Total iterations #>1 (e.g. 1000) for parallel
DOWNSAMP_FRAC = NULL #0.8
NUM_PARALLEL_MARKER_BNMF_ITERATIONS = 30 #Total
if(is.null(DISEASE_LABEL))
  DISEASE_LABEL = 'CLL'

suppressPackageStartupMessages(library(optparse))
get_option_list <- function(config=NULL){
  option_list <- list(
    make_option(c("-r", "--run"), action="store", type="character",help="", metavar = 'RUN', default=RUN_SUFFIX),
    make_option(c("-s", "--subset"), action="store", type="character",help="", metavar = 'SUBSET_SET', default=SUBSET_SET),
    make_option(c("-d", "--disease"), action="store", type="character",help="", metavar = 'DISEASE', default=DISEASE_LABEL),
    make_option(c("--bnmf_n_runs"), action="store", type="numeric",help="", metavar = 'TOTAL_MARKER_BNMF_RUNS', default=TOTAL_MARKER_BNMF_RUNS),
    make_option(c("--n_par"), action="store", type="numeric",help="", metavar = 'NUM_PARALLEL_MARKER_BNMF_ITERATIONS', default=NUM_PARALLEL_MARKER_BNMF_ITERATIONS),
    make_option(c("--ds_frac"), action="store", type="numeric",help="", metavar = 'DOWNSAMP_FRAC', default=DOWNSAMP_FRAC)
  )
  return(option_list)
}
option_list = get_option_list()
args <- parse_args(OptionParser(option_list=option_list, usage = "Rscript %prog [options]"), print_help_and_exit=F)
RUN_SUFFIX = args$run
SUBSET_SET = args$subset
DISEASE_LABEL = args$disease
TOTAL_MARKER_BNMF_RUNS = args$bnmf_n_runs
NUM_PARALLEL_MARKER_BNMF_ITERATIONS = args$n_par
DOWNSAMP_FRAC = args$ds_frac
print(args)

####################
#### Functions #####
####################
NMF.W <- function(X,W,tol,K) {
  n.run <- 1
  n.iter <- 1000000
  eps <- 1.e-50
  N <- dim(X)[1]
  M <- dim(X)[2]
  meanX <- mean(X,na.rm=T)
  eps <- 1.e-50
  for (j in 1:n.run) {
    H <- matrix(runif(K * M)*meanX,ncol=M)
    X.ap <- W %*% H
    error.EU <- sum((X-X.ap)^2)
    error.KL <- sum(X*log((X+eps)/(X.ap+eps))+X.ap-X)
    del <- 1
    count <- 1
    while (del >= tol & count < n.iter) {
      H <- H * (t(W) %*% X) / (t(W)%*%(W%*%H) + eps)
      X.ap <- W %*% (H)
      del <- abs(error.EU-sum((X-X.ap)^2))
      error.EU <- sum((X-X.ap)^2)
      error.KL <- sum(X*log((X+eps)/(X.ap+eps))+X.ap-X)
      if (count %% 100 == 0) cat(count,error.EU,error.KL,del,'\n')
      count <- count+1
    }
  }
  return(list(W,H))
}

NMF.H <- function(X,H,tol,K) {
  n.run <- 1
  n.iter <- 1000000
  eps <- 1.e-50
  N <- dim(X)[1]
  M <- dim(X)[2]
  meanX <- mean(X,na.rm=T)
  eps <- 1.e-50
  for (j in 1:n.run) {
    W <- matrix(runif(N * K)*meanX,ncol=K)
    X.ap <- W %*% H
    error.EU <- sum((X-X.ap)^2)
    error.KL <- sum(X*log((X+eps)/(X.ap+eps))+X.ap-X)
    del <- 1
    count <- 1
    while (del >= tol & count < n.iter) {
      W <- W * (X %*% t(H)) / ((W%*%H) %*% t(H)+eps)
      X.ap <- W %*% H
      del <- abs(error.EU-sum((X-X.ap)^2))
      error.EU <- sum((X-X.ap)^2)
      if (count %% 100 == 0) cat(count,error.EU,del,'\n')
      count <- count+1
    }
  }
  return(list(W,H))
}

#Compute BNMF for projection (keeping W of genes constant)
get.single.NMF <- function(X1,W0,tol,K) {
  X1[is.na(X1)] <- 0
  res <- NMF.W(X1,W0,tol,K)
  H1 <- res[[2]]
  g1 <- apply(H1,2,function(x) which.max(x))
  return(list(H1,g1))
}


#Projection based on markers
get.SSEC.fold.short <- function(expr, expr.fold, expr.fold.up, cohort, marker0, W1, group_labels=NULL) {
  n.sample <- ncol(expr)
  comm <- intersect(rownames(expr),marker0)
  W0.tmp <- as.matrix(W1[match(comm,rownames(W1),nomatch=0),])
  W0.tmp.norm <- t(apply(W0.tmp,1,function(x) x/sum(x)))
  K <- ncol(W1)
  H.tmp <- array(0,dim=c(K,n.sample))
  for (i in 1:n.sample) {
    X0.tmp <- as.matrix(expr.fold.up[match(comm,rownames(expr.fold.up),nomatch=0),i])
    x <- get.single.NMF(X0.tmp,W0.tmp,1.e-07,K)
    H.tmp[,i] <- x[[1]]
  }
  rownames(H.tmp) <- colnames(W0.tmp)
  colnames(H.tmp) <- colnames(expr.fold)
  H.tmp.norm <- apply(H.tmp,2,function(x) x/sum(x))
  g.tmp <- apply(H.tmp.norm,2,function(x) which.max(x))
  if(is.null(group_labels)){
    rownames(H.tmp.norm) <- paste("EC",seq(1:ncol(W1)),sep="")
  } else {
    rownames(H.tmp.norm) <- group_labels
  }

  x <- list(g.tmp,H.tmp)
  return(x)
}


get.marker.fold <- function(W1.sig, W1.norm.sig, cut.W, cut.fold, expression) {

  mRNA.fold <- expression
  index.max <- apply(W1.norm.sig,1,function(x) which.max(x))

  summary <- list()

  for(i in 1:ncol(W1.sig)){ # the total number of clusters

    eval(parse(text=paste0("gene",i," <- unique(names(index.max)[index.max==",i,"])")))
    eval(parse(text=paste0("index.sample",i," <- colnames(mRNA.fold)%in%sample",i)))
    eval(parse(text=paste0("index.gene",i," <- rownames(mRNA.fold)%in%gene",i)))
    eval(parse(text=paste0("index.sig",i," <- match(gene",i,",rownames(W1.sig),nomatch=0)")))
    eval(parse(text=paste0("summary.gene",i," <- data.frame(gene",i,",W1.sig[index.sig",i,",",i,"],W1.norm.sig[index.sig",i,",],rowMeans(mRNA.fold[index.gene",i,",index.sample",i,"],na.rm=T),rowMeans(mRNA.fold[index.gene",i,",!index.sample",i,"],na.rm=T))")))

    temp <- ""
    for(j in 1:ncol(W1.sig)){
      temp <- paste0(temp,"\"W.norm.sig",j,"\",")
    }
    eval(parse(text=paste0("colnames(summary.gene",i,") <- c(\"gene\",\"W.sig\",",temp,"\"mean1\",\"mean2\"",")")))

    eval(parse(text=paste0("summary.gene",i,"[,\"fold\"] <- summary.gene",i,"$mean1-summary.gene",i,"$mean2")))
    eval(parse(text=paste0("summary.gene",i," <- summary.gene",i,"[order(summary.gene",i,"$fold,decreasing=T),]")))
    eval(parse(text=paste0("gene",i," <- summary.gene",i,"$gene[(summary.gene",i,"$mean1-summary.gene",i,"$mean2)>cut.fold[",i,"]]")))
    eval(parse(text=paste0("summary[[",i,"]] <- summary.gene",i)))

  }

  temp <- ""
  for(j in 1:ncol(W1.sig)){
    temp <- paste0(temp,"gene",j,",")
  }
  eval(parse(text=paste0("return(list(",temp,"summary))")))

}


#*** Under construction
##To do: convert sample1, sample2 from main() into list too (sample.list)
get.marker.fold.v2 <- function(W1.sig, W1.norm.sig, cut.W, cut.fold, expression) {

  mRNA.fold <- expression
  index.max <- apply(W1.norm.sig,1,function(x) which.max(x))


  gd <- list(genes=list(), index.sample=list(), index.sig=list(), summary.gene=list())
  summary.gene.all <- list()
  for(i in 1:ncol(W1.sig)){ # the total number of clusters

    gd$gene[i] <- unique(names(index.max)[index.max==i])
    gd$index.sample[i] <- colnames(mRNA.fold) %in% sample.list[i] #*** Need to create sample.list (this is a GLOBAL VARIABLE - change this)
    gd$index.gene[i] <-  rownames(mRNA.fold) %in% gd$gene[i]
    gd$index.sig[i] <- match(gd$gene[i], rownames(W1.sig), nomatch=0)

    gd$summary.gene[i] <- data.frame(gd$gene[i],
                                            W1.sig[gd$index.sig[i],i],
                                            W1.norm.sig[gd$index.sig[i],],
                                            rowMeans(mRNA.fold[gd$index.gene[i], gd$index.sample[i]], na.rm=T),
                                            rowMeans(mRNA.fold[gd$index_gene[i], !gd$index.sample[i]],na.rm=T))

    colnames(gd$summary.gene[i]) = c('gene', 'W.sig', paste0('W.norm.sig', 1:col(W1.sig)), 'mean1', 'mean2')

    gd$summary.gene[i][,'fold'] <- gd$summary.gene[i]$mean1 - gd$summary.gene[i]$mean2 #This is fold for log-transformed values
    gd$summary.gene[i] <- gd$summary.gene[i][order(gd$summary.gene[i]$fold, decreasing=T),]
    gd$gene[i] <- gd$summary.gene[i]$gene[(gd$summary.gene[i]$mean1 - gd$summary.gene[i]$mean2) > cut.fold[i]]
    summary.gene.all[i] <- gd$summary.gene[i]

  }

  return(list(gd$gene, summary.gene.all))
}



###############
#### MAIN #####
###############
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(gprofiler2))
suppressPackageStartupMessages(library(msigdbr))
suppressPackageStartupMessages(library(fgsea))
suppressPackageStartupMessages(library(ggpubr))

#### CONSTS ####
GITHUB_DIR = '/Volumes/xchip_cga_home/bknisbac/github'
RSCRIPTS_DIR = paste0(GITHUB_DIR, '/', 'Rscripts')

BNMF_OUTPUT_DIR = '/Volumes/xchip_cga_home/bknisbac/CLL/results_cll1085/sample_sets/core2_tumors'
SUBSETS_ROOT_SUBDIR = 'subsets'
BNMF_EXPRESSION_CLUSTERING_SUBDIR = paste0('bayesnmf_tpm_', RUN_SUFFIX,
                                           ifelse(is.null(SUBSET_SET), '', paste0('/', SUBSETS_ROOT_SUBDIR, '/', SUBSET_SET)))
MARKER_GENES_SUBDIR = paste0('marker_genes_tpms_', RUN_SUFFIX)

## Gene expression matrix loading and transformation
GENE_COUNTS_FILE_FOR_MARKERS = paste0('/Volumes/xchip_cga_home/bknisbac/CLL/results_cll1085/sample_sets/core2_tumors/core2_tumors.rnaseqc_tpm_deseqLog10_', RUN_SUFFIX,'.txt.gz')

SAMPLE_INCLUSION_REGEX = '.*' #sample filtering e.g. to restrict analysis to specific subset

# counts processing and transformation

# GTF for annotation and filtering
#for rnaseqc
GENES_GTF_FILE = '/Volumes/xchip_cga_home/bknisbac/resources/gencode19_noChrPrefix_mitoMT.gtf.gz'
GENE_TYPES_TO_RETAIN = c('protein_coding', 'lincRNA')

CANCER_GENE_CENSUS = '/xchip/cga_home/bknisbac/resources/cancer_genes/cancer_gene_census_all_hg19_cosmic_v90_20200112.tsv'
cgc = fread(CANCER_GENE_CENSUS)

## BNMF biomarker detection consts
FORCE_RERUN_NMF = F

K_FOR_GET_BNMF_RES = NULL #NULL disables

cut.NA = 0.05 #max fraction of columns allowed to be NA
cut.W = 0.75
# cut.fold = c(1, 1, 1, 1, 1) # log2(fold) > cut.fold #one per group (but see function - it uses number 2 for most)
cut.fold.value = log2(1.5) # log2(fold) > cut.fold #one per group (but see function - it uses number 2 for most)
cut.gene = 10 # max number of marker genes per group
EXTRA = FALSE #If to add additional genes to final table (e.g. genes known to be important in this tissue)

##Plotting
#colorblind palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#CC79A7", "#924900", "#490092", "#920000", "#24ff24", "#000000", "#db6d00", "#006ddb")

#### Derived from CONSTS ####
bnmf_clustering_dir = paste(sep='/', BNMF_OUTPUT_DIR, BNMF_EXPRESSION_CLUSTERING_SUBDIR)
bnmf_marker_genes_dir = paste(sep='/', bnmf_clustering_dir, MARKER_GENES_SUBDIR)
if(!is.null(K_FOR_GET_BNMF_RES))
  bnmf_marker_genes_dir = paste0(bnmf_marker_genes_dir, '_', 'k', K_FOR_GET_BNMF_RES)
dir.create(bnmf_marker_genes_dir, showWarnings = F)

rnaseq_source_file = paste0(GITHUB_DIR, '/Rscripts/Sources/RNAseq.src.R') #For get_gene_counts()
suppressPackageStartupMessages(source(rnaseq_source_file))

## Set wd for output files
setwd(bnmf_marker_genes_dir)

#### Get gene counts ####
# Step summary:
# get matrix of numerics + assign Name to rowname
# retain only gene types of interest (e.g. protein_coding)
# deseq normalization
# log2(x+1) transformation
# remove unwanted samples (e.g. Normals)

#Explanation for parameter decisions:
# (a) I'm not using deseq becuase input is already deseq-normalized
# (b) min_value_after_unlog is purposely 1 (this or log_eps need to be >=1; if using deseq norm then specifically min_value_after_unlog must be >=1)
#For standard rnaseqc pipeline
cnts = get_gene_counts(counts_file = GENE_COUNTS_FILE_FOR_MARKERS, genes_gtf_file=GENES_GTF_FILE,
                       return_gene_col='Description', dup_gene_mode='sum', unlog_base = 10, negatives_assign_val=0, deseq_norm=F,
                       log_norm_base = 2, log_increment_val=1, return_class = 'data.frame')

#### Loading BNMF clustering data ####
## Within this code, subsetting of BNMF and counts is done so that same samples (intersect) are in both and dimensions match up ####
#get data for best run
res.Bayes = get_bnmf_results(bnmf_clustering_dir, K=K_FOR_GET_BNMF_RES)
sample_intersect = colnames(cnts)[colnames(cnts) %in% row.names(res.Bayes[[1]])]
# Get W, H matrices from best run
W <- res.Bayes[[1]]
W = W[row.names(W)[row.names(W) %in% sample_intersect],] #subset W for samples with count data
W <- W[,colSums(W)!=0]

H <- res.Bayes[[2]]
H = H[,colnames(H)[colnames(H) %in% sample_intersect]] #subset H for samples with count data
H <- H[rowSums(H)!=0,]
rownames(H) <- paste("EC",seq(1:nrow(H)),sep="")

H.norm <- apply(H,2,function(x) x/sum(x))
g.Bayes <- apply(H.norm,2,function(x) which.max(x))

## Subset only samples present in BNMF clustering
cnts = cnts[,names(g.Bayes)[names(g.Bayes) %in% colnames(cnts)]]

##Compute number of clusters ("groups")
num.groups = length(unique(g.Bayes))
for(i in 1:num.groups){
  eval(parse(text=paste0("sample",i," <- names(g.Bayes)[g.Bayes==",i,"]")))
}


#Split up and down-regulated genes
mRNA.comm = cnts
index.gene <- rowSums(is.na(mRNA.comm))<round(cut.NA*ncol(mRNA.comm)) #Remove genes with fraction NAs > cut.NA
mRNA.comm <- mRNA.comm[index.gene,]
mRNA.fold <- t(apply(mRNA.comm,1,function(x) x-median(x,na.rm=T))) #transform to fold change from median per gene (works for log-transformed counts)
mRNA.fold1 <- mRNA.fold
mRNA.fold2 <- mRNA.fold
mRNA.fold1[is.na(mRNA.fold1)] <- 0
mRNA.fold2[is.na(mRNA.fold2)] <- 0
mRNA.fold1 <- apply(mRNA.fold1,2,function(x) ifelse(x>0,x,0))
mRNA.fold2 <- apply(mRNA.fold2,2,function(x) ifelse(x>0,0,-x))
rownames(mRNA.fold1) <- paste(rownames(mRNA.fold1),"up",sep=".")
rownames(mRNA.fold2) <- paste(rownames(mRNA.fold2),"dn",sep=".")
mRNA.fold0 <- rbind(mRNA.fold1,mRNA.fold2) #contains both up and dn genes as separate rows, all values were tranformed to positive for NMF


############ determine W for tumor.type fold expression and H.norm
H1 <- H
H1.norm <- H.norm

myseed = 1
set.seed(myseed)
if(TOTAL_MARKER_BNMF_RUNS > 1){
  source(paste0(RSCRIPTS_DIR, '/', 'Sources/parallel.src.R'))
  parallel.registerCluster(NUM_PARALLEL_MARKER_BNMF_ITERATIONS) #create and register parallel backend
  parallel_marker_outdir = paste0(bnmf_marker_genes_dir, '/', 'parallel_marker_bnmf_runs', ifelse(!is.null(DOWNSAMP_FRAC), paste0('_ds', DOWNSAMP_FRAC), ''))
  dir.create(parallel_marker_outdir, showWarnings = F)
  foreach(i=1:TOTAL_MARKER_BNMF_RUNS) %dopar% {
    print(paste('Starting iteration', i, 'of', TOTAL_MARKER_BNMF_RUNS))
    set.seed(i)
    run_file = paste0(parallel_marker_outdir,'/',paste("res.NMF2",cut.NA,i,"RData",sep="."))
    if(is.null(DOWNSAMP_FRAC)){
      if(!file.exists(run_file) || FORCE_RERUN_NMF){
        res.NMF2 <- NMF.H(mRNA.fold0, H1.norm, 1.e-05, nrow(H1))
        save(res.NMF2, file=run_file)
      }
    } else {
      samps_selected = sample(colnames(H1.norm), round(ncol(H1.norm) * DOWNSAMP_FRAC))
      mRNA.fold0.subset = mRNA.fold0[,colnames(mRNA.fold0) %in% samps_selected]
      H1.norm.subset = H1.norm[,colnames(H1.norm) %in% samps_selected]

      if(!file.exists(run_file) || FORCE_RERUN_NMF){
        res.NMF2 <- NMF.H(mRNA.fold0.subset, H1.norm.subset, 1.e-05, nrow(H1))
        save(res.NMF2, file=run_file)
      }
    }
  }
  parallel.unregisterClusters() #unregister clusters upon activation
  save.image('r_image_post_marker_detection_parallel.RData')
  stop('Done parallel marker gene BNMF. Planned stop after parallel run - which is run in terminal')
} else if ((!file.exists(paste("res.NMF2",cut.NA,"RData",sep=".")) && !file.exists(paste("res.NMF2",cut.NA,myseed,"RData",sep="."))) || FORCE_RERUN_NMF) {
  res.NMF2 <- NMF.H(mRNA.fold0,H1.norm,1.e-05,nrow(H1))
  save(res.NMF2, file=paste("res.NMF2",cut.NA,"RData",sep="."))
  save.image('r_image_post_marker_detection_single.RData')
} else {
  if (file.exists(paste("res.NMF2",cut.NA,"RData",sep="."))){
    load(paste("res.NMF2",cut.NA,"RData",sep="."))
  } else {
    load(paste("res.NMF2",cut.NA,myseed,"RData",sep="."))
  }
}

# load('r_image_post_marker_detection_single.RData')
set.seed(1)

W1 <- res.NMF2[[1]]
W1 <- W1[index.gene,]
W1.up <- W1[grep("up",rownames(W1)),]
W1.dn <- W1[grep("dn",rownames(W1)),]
W1.norm <- t(apply(W1,1,function(x) x/sum(x)))
W1.norm.up <-  W1.norm[grep("up",rownames(W1.norm)),]
W1.norm.dn <-  W1.norm[grep("dn",rownames(W1.norm)),]

############ compute the maximum association to clusters across genes
############ up-regulated genes
max.up <- apply(W1.norm.up,1,function(x) max(x))
cut.W1 <- quantile(max.up,prob=cut.W,na.rm=T)
gene.sig.up <- rownames(W1.norm.up)[max.up >= cut.W1]
W1.sig.up <- W1.up[match(gene.sig.up,rownames(W1.up),nomatch=0),]
W1.norm.sig.up <- W1.norm.up[match(gene.sig.up,rownames(W1.norm.up),nomatch=0),]
cut.fold = rep(cut.fold.value, num.groups)
x <- get.marker.fold(W1.sig.up,W1.norm.sig.up,cut.W,cut.fold,mRNA.fold1)

for(i in 1:ncol(W1.sig.up)){
  eval(parse(text=paste0("gene",i,".up <- x[[",i,"]]")))
}
summary.gene.up <- x[[ncol(W1.sig.up)+1]]

############ down-regulated genes, we will not use this information
max.dn <- apply(W1.norm.dn,1,function(x) max(x))
cut.W1 <- quantile(max.dn,prob=cut.W,na.rm=T)
gene.sig.dn <- rownames(W1.norm.dn)[max.dn >= cut.W1]
W1.sig.dn <- W1.dn[match(gene.sig.dn,rownames(W1.dn),nomatch=0),]
W1.norm.sig.dn <- W1.norm.dn[match(gene.sig.dn,rownames(W1.norm.dn),nomatch=0),]
cut.fold = rep(cut.fold.value, num.groups)
x <- get.marker.fold(W1.sig.dn,W1.norm.sig.dn,cut.W,cut.fold,mRNA.fold2)

for(i in 1:ncol(W1.sig.up)){
  eval(parse(text=paste0("gene",i,".dn <- x[[",i,"]]")))
}
summary.gene.dn <- x[[ncol(W1.sig.dn)+1]]


#### Actually use only "up" genes ####
for(i in 1:ncol(W1.sig.up)){
  eval(parse(text=paste0("gene",i," <- gsub(\".up\",\"\",gene",i,".up)")))
  eval(parse(text=paste0("if (length(gene",i,")>cut.gene) gene",i," <- gene",i,"[1:cut.gene]")))
  eval(parse(text=paste0("gene",i," <- gene",i,"[gene",i,"%in%rownames(mRNA.comm)]")))
}

#Save markers (UP) as pathways file
markers_up_geneset_filename = "markers_up_as_geneset_list.txt"
if (file.exists(markers_up_geneset_filename))
  file.remove(markers_up_geneset_filename)
markers_up_as_gs = list()
for(i in 1:ncol(W1.sig.up)){#i=1
  eval(parse(text=paste0("print(gene",i,")")))
  markers_up_as_gs[i] = list(eval(parse(text=paste0("gene",i))))
  names(markers_up_as_gs)[i] = paste0(DISEASE_LABEL, '_EC', i, '_MARKERS') #name the gene set
  write(paste(sep='\t', names(markers_up_as_gs)[i], paste(markers_up_as_gs[[i]], collapse='\t')), markers_up_geneset_filename, append=T) #Write markers to file as pathways
}


############################
############################
### Markers up
W1 <- W1.up
rownames(W1) <- gsub(".up","",rownames(W1))

genes_group_order=c()
for(i in 1:ncol(W1.up)){#i=10
  if(eval(parse(text=paste0('length(gene',i,')==1')))){
    eval(parse(text=paste0("W.marker",i," <- as.data.frame(matrix(W1[match(as.character(gene",i,"),rownames(W1),nomatch=0),], nrow=1))")))
    eval(parse(text=paste0("row.names(W.marker",i,") <- as.character(gene",i,")")))
    eval(parse(text=paste0("colnames(W.marker",i,") <- colnames(W1)")))
    eval(parse(text=paste0("W.marker",i,".norm <- t(apply(W.marker",i,",1,function(x) x/sum(x)))")))
    eval(parse(text=paste0("gene",i," <- gene",i,"[order(W.marker",i,".norm[,",i,"],decreasing=T)]")))
    eval(parse(text=paste0("genes_group_order <- c(genes_group_order, gene", i,")")))
  } else if(eval(parse(text=paste0('length(gene',i,')>0')))){
    eval(parse(text=paste0("W.marker",i," <- W1[match(as.character(gene",i,"),rownames(W1),nomatch=0),]")))
    eval(parse(text=paste0("W.marker",i,".norm <- t(apply(W.marker",i,",1,function(x) x/sum(x)))")))
    eval(parse(text=paste0("gene",i," <- gene",i,"[order(W.marker",i,".norm[,",i,"],decreasing=T)]")))
    eval(parse(text=paste0("genes_group_order <- c(genes_group_order, gene", i,")")))
  } else{
    print(paste('Found no sig.up genes for cluster', i))
    eval(parse(text=paste0("gene",i," <- c()")))
  }

}
####

temp <- "marker0 <- unique(c("
for(j in 1:ncol(W1.sig.up)){
  ifelse(j!=ncol(W1.sig.up), temp <- paste0(temp,"as.character(gene",j,"),"),
                             temp <- paste0(temp,"as.character(gene",j,")))")      )
}
eval(parse(text=temp))

W.marker0 <- W1[match(marker0,rownames(W1),nomatch=0),]
write.table(W.marker0, 'W_markers.txt', quote=F)
W.marker0.norm <- t(apply(W.marker0,1,function(x) x/sum(x)))
write.table(W.marker0.norm, 'W_markers_norm.txt', quote=F)

#### Heatmap to visualize biomarkers
pdf('marker_genes_up_W1norm.pdf', height=12, width=10)
heatmap.2(as.matrix(W.marker0.norm[genes_group_order,]), hclustfun=function(c) {hclust(c,method="ward.D")}, distfun=function(c) {dist(c,method="euclidean")},
          na.rm = TRUE, scale="none", dendrogram="none",margins=c(5,5),
          Rowv=F, Colv=F, symbreaks=F, key=TRUE, symkey=F,main="CLL subtype marker genes",
          ColSideColors =cbPalette[1:ncol(W.marker0.norm)],
          density.info="none", trace="none",labCol=colnames(x),labRow=rownames(x),col=greenred(40),cex.lab=0.75,cexRow=1,cexCol=2,
          key.xlab='Gene-Cluster association', key.title = NA,
          keysize=0.75, key.par = list(cex=0.5))
dev.off()


### Heatmap of individual samples
cnts_markers = mRNA.fold[row.names(mRNA.fold) %in% row.names(W.marker0),]
if(all(colnames(cnts_markers)==names(g.Bayes))){
  sample_group_colors = cbPalette[g.Bayes]
  sample_group_colors_ordered = cbPalette[sort(g.Bayes)]
} else {
  print('Check why the columns and g.Bayes are not in same order')
}

max_abs = ceiling(max(abs(min(cnts_markers)),abs(max(cnts_markers))))
num_colors = 40
mybreaks=seq(-max_abs, max_abs, length.out=num_colors+1)
pdf("marker_genes_up_per_sample.log2foldmedian.pdf",38,60)
heatmap.2(as.matrix(cnts_markers[genes_group_order,names(sort(g.Bayes))]), hclustfun=function(c) {hclust(c,method="ward.D")}, distfun=function(c) {dist(c,method="euclidean")},
          na.rm = TRUE, scale="none", dendrogram="none",margins=c(5,5),
          Rowv=F, Colv=F, symbreaks=F, key=TRUE, symkey=F,main="CLL subtype marker genes",
          density.info="none", trace="none",labCol=colnames(x),labRow=rownames(x),col=bluered(num_colors), breaks=mybreaks,
          ColSideColors = sample_group_colors_ordered, cexRow=1, cexCol=0.3, offsetCol=0.5, adjCol=c(1,0.5),
          key.title=NA, key.xlab='log2 fold of median TPM', keysize = 0.5, key.par=list(cex.axis=2, cex.lab=2)) # ,cex.lab=0.2,cexRow=0.2
dev.off()


#####################
#### Here only use "dn" genes ####
for(i in 1:ncol(W1.sig.dn)){
  eval(parse(text=paste0("gene",i," <- gsub(\".dn\",\"\",gene",i,".dn)")))
  eval(parse(text=paste0("if (length(gene",i,")>cut.gene) gene",i," <- gene",i,"[1:cut.gene]")))
  eval(parse(text=paste0("gene",i," <- gene",i,"[gene",i,"%in%rownames(mRNA.comm)]")))
}

for(i in 1:ncol(W1.sig.dn)){
  eval(parse(text=paste0("print(gene",i,")")))
}

### Markers down
W1 <- W1.dn
rownames(W1) <- gsub(".dn","",rownames(W1))

genes_group_order=c()
for(i in 1:ncol(W1.dn)){#i=10
  if(eval(parse(text=paste0('length(gene',i,')==1')))){
    eval(parse(text=paste0("W.marker",i," <- as.data.frame(matrix(W1[match(as.character(gene",i,"),rownames(W1),nomatch=0),], nrow=1))")))
    eval(parse(text=paste0("row.names(W.marker",i,") <- as.character(gene",i,")")))
    eval(parse(text=paste0("colnames(W.marker",i,") <- colnames(W1)")))
    eval(parse(text=paste0("W.marker",i,".norm <- t(apply(W.marker",i,",1,function(x) x/sum(x)))")))
    eval(parse(text=paste0("gene",i," <- gene",i,"[order(W.marker",i,".norm[,",i,"],decreasing=T)]")))
    eval(parse(text=paste0("genes_group_order <- c(genes_group_order, gene", i,")")))
  } else if(eval(parse(text=paste0('length(gene',i,')>0')))){
    eval(parse(text=paste0("W.marker",i," <- W1[match(as.character(gene",i,"),rownames(W1),nomatch=0),]")))
    eval(parse(text=paste0("W.marker",i,".norm <- t(apply(W.marker",i,",1,function(x) x/sum(x)))")))
    eval(parse(text=paste0("gene",i," <- gene",i,"[order(W.marker",i,".norm[,",i,"],decreasing=T)]")))
    eval(parse(text=paste0("genes_group_order <- c(genes_group_order, gene", i,")")))
  } else{
    print(paste('Found no sig.dn genes for cluster', i))
    eval(parse(text=paste0("gene",i," <- c()")))
  }

}
####

temp <- "marker0 <- unique(c("
for(j in 1:ncol(W1.sig.dn)){
  ifelse(j!=ncol(W1.sig.dn), temp <- paste0(temp,"as.character(gene",j,"),"),
         temp <- paste0(temp,"as.character(gene",j,")))")      )
}
eval(parse(text=temp))

W.marker0 <- W1[match(marker0,rownames(W1),nomatch=0),]
write.table(W.marker0, 'W_markers_down.txt', quote=F)
W.marker0.norm <- t(apply(W.marker0,1,function(x) x/sum(x)))
write.table(W.marker0.norm, 'W_markers_norm_down.txt', quote=F)

#Save markers (DOWN) as pathways file
markers_dn_geneset_filename = "markers_dn_as_geneset_list.txt"
if (file.exists(markers_dn_geneset_filename))
  file.remove(markers_dn_geneset_filename)
markers_dn_as_gs = list()
for(i in 1:ncol(W1.sig.dn)){#i=1
  eval(parse(text=paste0("print(gene",i,")")))
  markers_dn_as_gs[i] = list(eval(parse(text=paste0("gene",i))))
  names(markers_dn_as_gs)[i] = paste0(DISEASE_LABEL,'_EC', i, '_MARKERS_DOWN') #name the gene set
  write(paste(sep='\t', names(markers_dn_as_gs)[i], paste(markers_dn_as_gs[[i]], collapse='\t')), markers_dn_geneset_filename, append=T) #Write markers to file as pathways
}


convert_genelist_gconvert <- function(gl){
  convres = tryCatch({
    as.data.table(gconvert(gl, organism = "hsapiens", target = "HGNC", numeric_ns = "", mthreshold = Inf, filter_na = TRUE))
  }, error = function(e) {
    data.table(input=gl, name=gl)
  })
  gl_new = convres$name[match(gl, convres$input)]
  return(gl_new)
}


### Identify markers found in cancer gene census (Note that genes names may have been changed by g:convert!)
#Then write CGC markers to file as pathways
cgc_markers_up_as_gs <- markers_up_as_gs
cgc_markers_dn_as_gs <- markers_dn_as_gs
cgc_upfile = paste0('cgc_', markers_up_geneset_filename)
cgc_dnfile = paste0('cgc_', markers_dn_geneset_filename)
tmp = file.create(cgc_upfile, showWarnings = F)
tmp = file.create(cgc_dnfile, showWarnings = F)
for (i in 1:length(markers_up_as_gs)){ #i=1
    #Get CGC genes for markers-up
    if(length(markers_up_as_gs[[i]])>0){
      cgc_markers_up_as_gs[[i]] = convert_genelist_gconvert(cgc_markers_up_as_gs[[i]])
    } else {
      cgc_markers_up_as_gs[[i]] = c('')
    }
    cgc_markers_up_as_gs[[i]] = cgc_markers_up_as_gs[[i]][cgc_markers_up_as_gs[[i]] %in% cgc$`Gene Symbol`]
    write(paste(sep='\t', names(cgc_markers_up_as_gs)[i], paste(cgc_markers_up_as_gs[[i]], collapse='\t')), cgc_upfile, append=T)
    #Get CGC genes for markers-dn
    if(length(markers_dn_as_gs[[i]])>0){
      cgc_markers_dn_as_gs[[i]] = convert_genelist_gconvert(cgc_markers_dn_as_gs[[i]])
    } else {
      cgc_markers_dn_as_gs[[i]] = c('')
    }
    cgc_markers_dn_as_gs[[i]] = cgc_markers_dn_as_gs[[i]][cgc_markers_dn_as_gs[[i]] %in% cgc$`Gene Symbol`]
    write(paste(sep='\t', names(cgc_markers_dn_as_gs)[i], paste(cgc_markers_dn_as_gs[[i]], collapse='\t')), cgc_dnfile, append=T)
}



################################################
##### Pathway/Gene set enrichment analysis #####
################################################

### GSEA analysis on W matrix
### Set up W matrix for pathway analysis (need to merge results of marker gene detection and set downregulated to negative values)
#Select greatest value from up or down markers
UPDATE_GENES_FOR_GSEA = F
get_larger_abs <- function(x,y){ifelse(abs(x)>abs(y), x, y)}
wpw = get_larger_abs(W1.up, -W1.dn)
row.names(wpw) = sub('.up$', '', row.names(wpw))
w_genes = row.names(wpw)
wpw = as.data.table(wpw)

#Convert gene names to newer version with g:profiler's g:convert
if(UPDATE_GENES_FOR_GSEA){
  convs = as.data.table(gconvert(w_genes, organism = "hsapiens", target = "HGNC",
                                 numeric_ns = "", mthreshold = Inf, filter_na = TRUE))
  wpw$gene = convs$name[match(w_genes, convs$input)]
  wpw = wpw[!is.na(gene),]
  wpw = wpw[!wpw$gene %in% wpw$gene[duplicated(wpw$gene)]] #remove duplicated gene names (could be introduced by ID conversion; shouldn't be in counts data)
} else {
  wpw$gene = w_genes
}

#Get pathways of interest
if(UPDATE_GENES_FOR_GSEA){ #use updated HGNC symbols
  pws = get_msigdb_pathways(msigdb_sets=c('H', 'C5:BP', 'C2:CP:REACTOME'))
} else { #use gencode19-adapted pathways
  PATHWAY_GMT_DIR = '/xchip/cga_home/bknisbac/resources/pathways/msigdb_gencode19'
  pw_file = paste0(PATHWAY_GMT_DIR, '/', 'msigdb_H_C5-BP_C2-CP-REACTOME_gencode19.gmt')
  pws = read_gmt(pw_file)
}

#Run fGSEA
GSEA_FDR_CUTOFF=0.1
gsres_all = NULL
set.seed(1)
pdf('fgsea_plots_per_ec.pdf', 18,18)
for (i in 1:(ncol(wpw)-1)){ #i=1
  setorderv(wpw, colnames(wpw)[i], c(-1))

  ranks = wpw[[i]]
  names(ranks) = wpw$gene
  # ranks = sort(ranks, decreasing=TRUE) #don't need to sort
  fgseaRes <- fgsea(pws, ranks, minSize=12, maxSize=500, nproc=1, eps=0, scoreType='std')
  fgseaRes$group = colnames(wpw)[i]
  fgseaRes_sig = fgseaRes[padj < GSEA_FDR_CUTOFF ,][order(-NES),]
  plot(plotGseaTable(pws[names(pws) %in% fgseaRes_sig$pathway],
                ranks, fgseaRes_sig, gseaParam = 1, colwidths = c(5, 3, 0.8, 1.2, 1.2), render = F), title=colnames(wpw)[i])
  gsres_all = rbindlist(list(gsres_all, fgseaRes))
}
dev.off()
fwrite(gsres_all, 'fgsea_res_per_ec.tsv', sep='\t')
gsres_sig = gsres_all[padj < GSEA_FDR_CUTOFF ,]
fwrite(gsres_sig, 'fgsea_res_per_ec_sig.tsv', sep='\t')

GSEA_COLLAPSE_PATHWAYS = T
if(GSEA_COLLAPSE_PATHWAYS){
  collapsedPathways = collapsePathways(gsres_sig, pathways = pws, stats = ranks, pval.threshold = 0.001)
  gsres_sig_main <- gsres_sig[pathway %in% collapsedPathways$mainPathways]
  fwrite(gsres_sig_main, 'fgsea_res_per_ec_sig_collapsed.tsv', sep='\t')
}


### Table of pathway enrichment by W matrix of ECs
# gsres_sig_for_heatmap = gsres_sig[padj<0.1 & abs(NES)>=1.5,]
if(GSEA_COLLAPSE_PATHWAYS){dt_for_dcast = gsres_sig_main} else {dt_for_dcast = gsres_sig}
gs_mat_df = as.data.frame(dcast.data.table(dt_for_dcast, pathway ~ group, value.var = 'NES', fill=0))
row.names(gs_mat_df) = gs_mat_df$pathway
ec_order = paste0('EC', 1:(ncol(gs_mat_df)-1))
gs_mat_df = gs_mat_df[,ec_order[ec_order %in% colnames(gs_mat_df)] ]

pdf('heatmap_gsea_per_ec_by_W_matrix.pdf', height = dim(gs_mat_df)[1]/5, width=14)
heatmap.2(as.matrix(gs_mat_df), col=colorRampPalette(c('blue', 'yellow'))(256), dendrogram='column', Rowv=T, Colv=F,
          na.color="#E8E8E8", scale="none", density.info="none", trace="none",
          key.title='', key.xlab = 'NES', lhei = c(1,5), margins=c(4,20))
dev.off()

write.table(gs_mat_df, 'fgsea_res_per_ec_sig_NEStable.tsv', sep='\t', quote=F, row.names=T, col.names=T)
fwrite(dcast.data.table(gsres_all, pathway ~ group, value.var = 'NES', fill=0), 'fgsea_res_per_ec_all_NEStable.tsv', sep='\t')


if (any(grepl('^HALLMARK', row.names(gs_mat_df)))){
  gs_mat_df_Hallmarks = gs_mat_df[grepl('^HALLMARK', row.names(gs_mat_df)),]
  pdf('heatmap_gsea_per_ec_by_W_matrix_HallmarksOnly.pdf', height = 8, width=14)
  heatmap.2(as.matrix(gs_mat_df_Hallmarks), col=colorRampPalette(c('blue', 'yellow'))(256), dendrogram='none', Rowv=T, Colv=F,
            na.color="#E8E8E8", scale="none", density.info="none", trace="none",
            key.title='', key.xlab = 'NES', lhei = c(1,5), margins=c(4,22))
  dev.off()
}


#Create table of fgsea results for W matrix, with NES and padj values side by side for all
nes_padj_anysig_dt = dcast.data.table(gsres_all, pathway ~ group, value.var = 'NES', fill=0)[,c('pathway', ec_order),with=F]
setnames(nes_padj_anysig_dt, c('pathway', paste0(ec_order, '_NES')))
padj_dt = dcast.data.table(gsres_all, pathway ~ group, value.var = 'padj', fill=0)[,ec_order,with=F]
setnames(padj_dt, paste0(ec_order, '_padj'))
nes_padj_anysig_dt = cbind(nes_padj_anysig_dt, padj_dt)
rm(padj_dt)
nes_padj_anysig_dt = nes_padj_anysig_dt[apply(nes_padj_anysig_dt[,.SD < 0.1,.SDcols=paste0(ec_order, '_padj')], 1, any)]
nes_padj_anysig_dt[, nes_sdev := apply(.SD, 1, sd), .SDcols=paste0(ec_order, '_NES')]
nes_padj_anysig_dt[, nes_max := apply(.SD, 1, max), .SDcols=paste0(ec_order, '_NES')]
nes_padj_anysig_dt[, nes_min := apply(.SD, 1, min), .SDcols=paste0(ec_order, '_NES')]
nes_padj_anysig_dt[, nes_range := nes_max-nes_min]
fwrite(nes_padj_anysig_dt, 'fgsea_res_per_ec_anysig_NES_padj_table.tsv', sep='\t')

## Save table of top and bottom 10 most extreme NES values per EC
extreme_nes_dt = NULL
n_extreme=25
for (ec in ec_order){
  setorderv(nes_padj_anysig_dt, paste0(ec, '_', 'NES'))
  cols_to_retain = c('pathway', paste0(ec, '_', 'NES'), paste0(ec, '_', 'padj'))
  ec_tmp = rbind(head(nes_padj_anysig_dt[,cols_to_retain,with=F], n=n_extreme), tail(nes_padj_anysig_dt[,cols_to_retain,with=F], n=n_extreme), use.names=F)
  setnames(ec_tmp, c('pathway', 'NES', 'padj'))
  ec_tmp$group = ec
  extreme_nes_dt = rbindlist(list(extreme_nes_dt, ec_tmp))
}
fwrite(extreme_nes_dt, paste0('fgsea_res_per_ec_anysig_topAndBottom', n_extreme, 'extremeNESperEC.tsv'), sep='\t')
fwrite(nes_padj_anysig_dt[pathway %in% extreme_nes_dt$pathway,], 'fgsea_res_per_ec_anysig_topAndBottom_NES_padj_table.tsv', sep='\t')


### Plot limited amount of extreme up or down regulated pathways
extreme_nes_dt_forplot = NULL
n_extreme_forplot=3
for (ec in ec_order){
  setorderv(nes_padj_anysig_dt, paste0(ec, '_', 'NES'))
  cols_to_retain = c('pathway', paste0(ec, '_', 'NES'), paste0(ec, '_', 'padj'))
  ec_tmp = rbind(head(nes_padj_anysig_dt[,cols_to_retain,with=F], n=n_extreme_forplot), tail(nes_padj_anysig_dt[,cols_to_retain,with=F], n=n_extreme_forplot), use.names=F)
  setnames(ec_tmp, c('pathway', 'NES', 'padj'))
  ec_tmp$group = ec
  extreme_nes_dt_forplot = rbindlist(list(extreme_nes_dt_forplot, ec_tmp))
}

gs_mat_df_extreme = gs_mat_df[row.names(gs_mat_df) %in% extreme_nes_dt_forplot$pathway,]
pdf(paste0('heatmap_gsea_per_ec_by_W_matrix_extreme', n_extreme_forplot,'.pdf'), height = 10, width=20)
heatmap.2(as.matrix(gs_mat_df_extreme), col=colorRampPalette(c('blue', 'yellow'))(256), dendrogram='none', Rowv=T, Colv=F, #dendrogram='column',
          na.color="#E8E8E8", scale="none", density.info="none", trace="none",
          key.title='', key.xlab = 'NES', lhei = c(1,5), margins=c(4,50))
dev.off()


#### UMAP based on marker genes for each EC
##prep1: Markers to list of genes
markers_up_genes = sort(unique(c(unlist(markers_up_as_gs))))
markers_dn_genes = sort(unique(c(unlist(markers_dn_as_gs))))
markers_all_genes = sort(unique(c(markers_up_genes, markers_dn_genes)))

##prep2 EC list to datatable
ec_dt = data.table(sample_id=names(g.Bayes), ec_label=paste0('EC', g.Bayes))

library(umap)
get_markers_umap <- function(cnts_in, genes_for_umap=NULL, print_plot=T, pca=NULL, normalize_per_sample=F, standardize_per_gene=F, center_per_gene=T, group_levels=NULL, retmode='umap_dt'){ #, , standardize_per_gene=T
  if(is.null(genes_for_umap)){
    genes_for_umap = row.names(cnts_in)
  }
  cnts_markers = cnts_in[row.names(cnts_in) %in% genes_for_umap,]

  if(normalize_per_sample) #this also transposes
    cnts_markers = apply(cnts_markers,2,function(x){x/sum(x)})
  else
    cnts_markers = cnts_markers

  if(standardize_per_gene)
    cnts_markers_t = t(apply(cnts_markers,2,function(x){(x-mean(x))/sd(x)}))
  else if(center_per_gene)
    cnts_markers_t = t(apply(cnts_markers,2,function(x){x-mean(x)}))
  else
    cnts_markers_t = t(cnts_markers)

  umap_markers = umap(cnts_markers_t, pca=pca)
  umap_markers_dt = as.data.table(umap_markers$layout)
  setnames(umap_markers_dt, c('UMAP1', 'UMAP2'))
  umap_markers_dt$sample_id = row.names(umap_markers$layout)
  umap_markers_dt = merge(umap_markers_dt, ec_dt)
  if(!is.null(group_levels))
    umap_markers_dt$ec_label = factor(umap_markers_dt$ec_label, ordered=T, levels = ec_order)
  else
    umap_markers_dt$ec_label = factor(umap_markers_dt$ec_label, ordered=T, levels = sort(unique(umap_markers_dt$ec_label)))
  if(retmode=='umap_dt')
    return(umap_markers_dt)
  else if(retmode=='markers_dt')
    return(cnts_markers_t)
}


#### Compute and plot umaps
pdf('marker_genes_umap.pdf')
#All markers - normalized samples and centered per genes
umap_all_doubleNorm = get_markers_umap(cnts, markers_all_genes, retmode='umap_dt', normalize_per_sample=T, center_per_gene = T)
norm_p = ggscatter(umap_all_doubleNorm, x = "UMAP1", y = "UMAP2", color = "ec_label", palette = cbPalette, legend.title='Expression cluster',
                   title='UMAP by BNMF markers (all)\nNorm samples & center genes')
norm_p

#All markers - normalized by genes and samples
umap_all_doubleNorm = get_markers_umap(cnts, markers_all_genes, retmode='umap_dt', normalize_per_sample=T, standardize_per_gene = T)
norm_p2 = ggscatter(umap_all_doubleNorm, x = "UMAP1", y = "UMAP2", color = "ec_label", palette = cbPalette, legend.title='Expression cluster',
                   title='UMAP by BNMF markers (all)\nNorm samples & standardize genes')
norm_p2

dev.off()

### Trimmed list of top N marker genes
umap_n_top_markers=20
markers_up_trimmed = sort(unique(c(unlist(lapply(markers_up_as_gs, function(x){return(x[1:umap_n_top_markers])})))))
markers_dn_trimmed = sort(unique(c(unlist(lapply(markers_dn_as_gs, function(x){return(x[1:umap_n_top_markers])})))))
markers_all_trimmed = sort(unique(c(markers_up_trimmed, markers_dn_trimmed)))

pdf(paste0('marker_genes_umap_top_', umap_n_top_markers,'.pdf'))
#All markers
umap_all = get_markers_umap(cnts, markers_all_trimmed)
ggscatter(umap_all, x = "UMAP1", y = "UMAP2", color = "ec_label", palette = cbPalette, legend.title='Expression cluster',
          title='UMAP by BNMF markers (all)')

#Markers (up)
umap_up = get_markers_umap(cnts, markers_up_trimmed)
ggscatter(umap_up, x = "UMAP1", y = "UMAP2", color = "ec_label", palette = cbPalette, legend.title='Expression cluster',
          title='UMAP by BNMF markers (UP)')

#Markers (down)
umap_dn = get_markers_umap(cnts, markers_dn_trimmed)
ggscatter(umap_dn, x = "UMAP1", y = "UMAP2", color = "ec_label", palette = cbPalette, legend.title='Expression cluster',
          title='UMAP by BNMF markers (DOWN)')
dev.off()


### Save image at end of script
if (!file.exists('r_image_end_of_marker_detection_script.RData')){
  save.image('r_image_end_of_marker_detection_script.RData')
}

### Load workspace from last run
# load('r_image_end_of_marker_detection_script.RData')
