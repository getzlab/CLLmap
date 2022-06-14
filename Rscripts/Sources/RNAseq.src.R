###Source file of functions for RNA-seq analysis
#The script contains various utilities for reading/retrieving, filtering and manipulating gene expression matrices, gene set data and BayesNMF results (related to expression clusters)
# Author: Binyamin Knisbacher
library(data.table)
library(DESeq2)
library(msigdbr)
library(stringr)

### Get gene counts table and optionally filter by GTF ###
#dup_gene_mode = 'all', 'sum', 'first'
#Note: the min_value_after_unlog was needed to prevent the final log-transformed result from having negative values (this was after remove_covariates() python method, which leads to negative log-transformed values)
#Expecting counts file/df to have c('Name', 'Description') columns
#negatives_assign_val= NULL, '0' (no reason to choose other value)
#Return value differs if return_gene_col has duplicate entries (it will be assigned to row.names() if return value is data.frame and there are no dups; if there are dups then both gene_id_col and gene_name_col are returned; if return value is data.table then first column will be return_gene_col)
get_gene_counts <- function(counts_file=NULL, gene_counts=NULL, sample_inclusion_regex='.*', sample_inclusion_list=NULL, sample_exclusion_list=NULL,
                            gene_id_col='Name', gene_name_col='Description', dup_gene_mode='sum',
                            genes_gtf_file=NULL, gtf_feature='gene', gene_types_to_retain=c('protein_coding'), unlog_base=NULL, unlog_reduce_res = 1, min_value_after_unlog=NULL, convert_units=NULL,
                            sizefactor_method='ratio', negatives_assign_val=0, negatives_correct_by_shift=F,
                            deseq_norm=FALSE, log_norm_base=NULL, log_increment_val=1,
                            return_class='data.frame', return_gene_col='Name', return_gene_col_merge = 'gene', rename_gene_col=NULL, verbose=T){

  ## Check input
  if(is.null(counts_file) && is.null(gene_counts)){
    warning('Must specify counts_file or gene_counts dataframe - exiting get_gene_counts()')
    return(NULL)
  }

  ## Get genes from file
  if(verbose)
    print('Getting counts')
  if(is.null(gene_counts)){ #counts_file specified
    if(grepl('.gz$', counts_file)){
      gene_counts = fread(cmd = paste('gunzip', '-cq', counts_file), sep='\t')
    } else {
      gene_counts = fread(counts_file, sep='\t')
    }
  } else{ #data.table or data.frame specified in gene_counts
    if(!'data.table' %in% class(gene_counts)){ #data.frame specified
      gene_counts = as.data.table(gene_counts)
    }
  }

  # Retain subset (or retain all using '.*' regex; see above)
  gene_counts = gene_counts[,grepl(sample_inclusion_regex, names(gene_counts)) | !sapply(gene_counts,is.numeric),with=F]

  if(!is.null(sample_inclusion_list))
    gene_counts = gene_counts[, !sapply(gene_counts,is.numeric) | (names(gene_counts) %in% sample_inclusion_list),with=F]
  if(!is.null(sample_exclusion_list))
    gene_counts = gene_counts[, (!names(gene_counts) %in% sample_exclusion_list), with=F]

  #Remove samples that are all NAs or 0 (had this once due to error in running previous algorithm)
  samps_not_all_na_or_zero = apply(gene_counts, 2, function(x){any(!is.na(x)) && any(x>0)})
  if( any(samps_not_all_na_or_zero == F) ){
    print('Removing following columns due to all-NA or all-zero expression:')
    print(colnames(gene_counts)[samps_not_all_na_or_zero==F])
    gene_counts = gene_counts[, samps_not_all_na_or_zero, with=F]
  }

  ### Get GTF data for subsetting
  if(!is.null(genes_gtf_file)){
    if(verbose)
      print('Getting genes from GTF')

    if(grepl('.gz$', genes_gtf_file)){
      genes_gtf = fread(cmd=paste('gunzip -cq', genes_gtf_file), sep='\t', header=F, skip='#')[V3==gtf_feature]
    } else {
      genes_gtf = fread(genes_gtf_file, sep='\t', header=F, skip='#')[V3==gtf_feature]
    }

    gtf_colnames = c('seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute')
    setnames(genes_gtf, gtf_colnames)
    genes_gtf$gene_type = sub('.*gene_type \"([^"]*)\".*', '\\1', genes_gtf$attribute)
    genes_gtf$gene_id = sub('.*gene_id \"([^"]*)\".*', '\\1', genes_gtf$attribute)

    ## Retain only subset of genes of interest (e.g. protein_coding only; see gencode biotypes at https://www.gencodegenes.org/pages/biotypes.html)
    if(!is.null(gene_types_to_retain)){
      ind_to_retain = which(gene_counts[,c(gene_id_col),with=F][[1]] %in% genes_gtf[genes_gtf$gene_type %in% gene_types_to_retain]$gene_id)
      gene_counts <- gene_counts[ind_to_retain]
    }
  }

  #split between numeric and non-numeric columns
  #Save non-numeric columns (important: this must be after gene subsetting, before all numeric calculations)
  char_gct_dt = copy(gene_counts[,!sapply(gene_counts,is.numeric),with=F])
  gene_counts = gene_counts[,sapply(gene_counts,is.numeric),with=F]
  if(ncol(char_gct_dt)>2){
    print(paste('You have', ncol(char_gct_dt), 'non-numeric columns - make sure this is not an error'))
  }

  ## Revert log_N(x+1) transformation applied to input data (N := unlog_base)
  if (!is.null(unlog_base)){ #revert log transformation
    gene_counts = gene_counts[,lapply(.SD, function(x){unlog_base ^ x - unlog_reduce_res})]
    if(!is.null(min_value_after_unlog)){
      min_gct_unlog = min(gene_counts[,lapply(.SD, min)])
      gene_counts = gene_counts - min_gct_unlog + min_value_after_unlog
    }
  }

  ### After numeric actions - add non-numeric columns
  gene_counts = cbind(char_gct_dt, gene_counts)

  ### Merge columns if '|' specified
  if(grepl('\\|', return_gene_col)){ #get two columns and merge into one column
    col1 = strsplit(return_gene_col, '\\|')[[1]][1]
    col2 = strsplit(return_gene_col, '\\|')[[1]][2]
    gene_counts[,return_gene_col_merge] = paste0(gene_counts[,get(col1)], '|', gene_counts[,get(col2)])
    gene_counts = gene_counts[, !c(col1, col2)[!c(col1, col2) %in% c(return_gene_col_merge)], with=F] #drop col1, col2 if not 'gene'
    return_gene_col = return_gene_col_merge
  }

  # Deal with duplicated gene names in GTF and split numeric from non-numeric columns
  if (dup_gene_mode == 'sum'){
    print('Summing counts for duplicate genes')
    gene_counts = gene_counts[,lapply(.SD, sum), by=return_gene_col, .SDcols=sapply(gene_counts, is.numeric)]
    split_cols = c(return_gene_col)
  } else if(dup_gene_mode == 'first'){
    print('Selecting first for duplicate genes')
    gene_counts = unique(gene_counts, by=eval(return_gene_col))
    print(head(colnames(gene_counts)))
    split_cols = c(gene_id_col, gene_name_col)
  } else {
    print('Not taking action or testing for duplicate genes. Either change the option or make sure you are okay with this..')
    split_cols = c(gene_id_col, gene_name_col)
  }

  #Split numeric from non-numeric columns save gene_id and gene_name seperately
  gene_cols_df = gene_counts[, split_cols, with=F]
  gene_counts = gene_counts[, !split_cols, with=F]
  gene_counts = gene_counts[, sapply(gene_counts, is.numeric), with=F] #drop non-numeric

  ## Correct negative values (necessary for DESeq normalization and SignatureAnalyzer ARD)
  neg_val_count = sum(gene_counts < 0)
  if(neg_val_count > 0){
    print(paste0('Total of ', neg_val_count, ' negative counts found! ',
                 round(neg_val_count / (nrow(gene_counts)*ncol(gene_counts)) * 100, 3),'% of values' ))
    negative_sum = sum(sapply(gene_counts,function(x){sum(ifelse(x<0,abs(x),0))}))
    positive_sum = sum(sapply(gene_counts,function(x){sum(ifelse(x>0,x,0))}))
    print(paste0('Total abs value of negative values is ', round(negative_sum, 3), ', Fraction of sum of positive values: ', negative_sum / positive_sum ))

    if(!is.null(negatives_assign_val)){
      print(paste('You have negative counts! assigning ', negatives_assign_val,' to negative values'))
      gene_counts[gene_counts < 0,] = negatives_assign_val
    } else if(negatives_correct_by_shift){ #increment all values so min is zero
      print(paste('Shifting values by', min(gene_counts * -1), 'so minimal value is 0'))
      gene_counts = gene_counts - min(min(gene_counts))
    } else {
      print('You have chosen not to correct negative values')
    }
  }

  ## DESeq2 normalization
  if (deseq_norm){
    if(sum(gene_counts < 0)){
      print('You must remove negative values for deseq2 normalization - returning NULL!')
      return(NULL)
    }

    #Estimate size factors
    if(sizefactor_method=='ratio'){ #regular method using geometric mean
      s <- estimateSizeFactorsForMatrix(gene_counts)
    } else if(sizefactor_method=='poscounts'){ #may be useful for many zeros - option 1 (PREFERRED)
      dds = DESeqDataSetFromMatrix(gene_counts)
      dds <- estimateSizeFactors(dds, type='poscounts')
      s <- sizeFactors(dds)
    } else if(sizefactor_method=='iterate'){ #may be useful for many zeros - option 2 (may be slow and not converge)
      dds = DESeqDataSetFromMatrix(gene_counts)
      dds <- estimateSizeFactors(dds, type='iterate')
      s <- sizeFactors(dds)
    }
    #Use estimated size factors to normalize
    gene_counts <- t(apply(gene_counts, 1, "/", s))
  }

  ## log_N(x+1) transformation (N is specified base)
  if(!is.null(log_norm_base)){
    gene_counts <- log(gene_counts+log_increment_val, base=log_norm_base)
  }

  if(return_class=='data.frame' || return_class=='matrix'){
    #non-unique return_gene_col values exist - can't assign to row names
    if(length(unique(gene_cols_df[,..return_gene_col][[1]])) !=  length(gene_cols_df[,..return_gene_col][[1]]) ){
      if(return_class=='matrix'){ #return matrix
        print('Can not return mat with duplicate row names; returning NULL')
        return(NULL)
      } else { #return data.frame
        gene_counts = as.data.frame(cbind(gene_cols_df, gene_counts))
        if(!is.null(rename_gene_col))
          colnames(gene_counts)[1] = rename_gene_col
      }
    } else { #no duplicates - assing to row.names
      gene_counts = as.data.frame(gene_counts)
      row.names(gene_counts) = gene_cols_df[,..return_gene_col][[1]]
      if(return_class=='matrix') #convert to matrix
        gene_counts = as.matrix(gene_counts)
    }
  } else { #data.table
    gene_counts = cbind(gene_cols_df, gene_counts)
    if(!is.null(rename_gene_col)){
      setnames(gene_counts, return_gene_col, rename_gene_col)
    }
  }
  return(gene_counts)
}

coef.var <- function(x, ignore.na=T){
  if(ignore.na){
    return(100*sd(x[!is.na(x)])/mean(x[!is.na(x)]))
  } else{
    return(100*sd(x)/mean(x))
  }
}

coef.var2 <- function(x, exponent = 1, ignore.na=T){
  if(ignore.na){
    return((sd(x[!is.na(x)])/mean(x[!is.na(x)]))^2)
  } else{
    return((sd(x)/mean(x))^2)
  }
}

subset_gene_expression <- function(gcts, min_quantile_expression_mean=0.1,
                                   min_quantile_expression_sd = 0.75,
                                   max_na_or_zero_samples_frac = 0.1,
                                   min_quantile_expression_cv=0,
                                   min_quantile_expression_mad=0,
                                   top_n_mode=NULL,
                                   top_n=5000,
                                   drop_non_numeric=F){

  original_cols = colnames(gcts)
  if((length(class(gcts))==1 && class(gcts)=='data.frame')){ #data.frame
    gcts = as.data.table(gcts) #make a copy
  } else if('data.table' %in% class(gcts)){ #data.table
    gcts = copy(gcts)
  }

  gcts_sd_cols = colnames(gcts)[sapply(gcts, is.numeric)]
  num_samps = length(gcts_sd_cols)
  gcts[,numZeroOrNA := Reduce(`+`, lapply(.SD,function(x){is.na(x) | x==0})), .SDcols=gcts_sd_cols]
  gcts$meanExp = rowMeans(gcts[, gcts_sd_cols, with=F], na.rm = T)
  gcts$sdExp = rowSds(as.matrix(gcts[, gcts_sd_cols, with=F]), na.rm = T)
  gcts$cvExp = apply(gcts[, gcts_sd_cols, with=F], 1, coef.var)
  gcts$madExp = apply(gcts[, gcts_sd_cols, with=F], 1, mad)

  gcts = gcts[numZeroOrNA < max_na_or_zero_samples_frac * num_samps]
  gcts = gcts[sdExp >= quantile(gcts$sdExp, min_quantile_expression_sd) &
                meanExp >= quantile(gcts$meanExp, min_quantile_expression_mean) &
                cvExp >= quantile(gcts$cvExp, min_quantile_expression_cv) &
                madExp >= quantile(gcts$madExp, min_quantile_expression_mad)]

  if(!is.null(top_n_mode)){
    if(top_n_mode=='cv'){
      setorder(gcts, -cvExp)
    } else if(top_n_mode=='mad'){
      setorder(gcts, -madExp)
    }
    gcts = gcts[1:top_n,]
  }

  #Drop any intermediate compute columns
  gcts = gcts[,original_cols,with=F]

  #subset based on indexes
  if(drop_non_numeric){
    gcts = gcts[,sapply(gcts, is.numeric), with=F]
  }

  return(gcts)
}

### Gene ID conversion
#na_mode=c('original', 'filter') #if to use original or filter out NA (non-matched input genes)
#Usages:
#  convert_genes(my_gene_list, na_mode = 'original') #upgrades hugo symbols with gprofiler's gconvert
#  convert_genes(my_gene_list, gene_table=GENCODE_UPDATE_BED) #upgrades hugo symbols by gene table (main use of the table is for reversion)
#  convert_genes(my_gene_list, gene_table=GENCODE_UPDATE_BED, gt_from = 'gene_name_updated', gt_to='gene_name') #reverts new hugo symbols to gencode19
GENCODE_UPDATE_BED = '/Volumes/xchip_cga_home/bknisbac/resources/gencode19_noChrPrefix_mitoMT.withNewHugo.bed.gz'
update_gene_names <- function(gl, na_mode='original',
                              gene_table=NULL, gt_from='gene_name', gt_to='gene_name_updated',
                              use_default_gene_table=F,
                              gp_target = "HGNC", gp_numeric_ns = "", gp_mthreshold = Inf){
  if(use_default_gene_table){
    gene_table = GENCODE_UPDATE_BED
  }
  if(!is.null(gene_table)){ #filename
    if(length(class(gene_table))==1 && class(gene_table)=='character'){
      gene_dt = fread(gene_table)
    } else { #DF/DT
      gene_dt = copy(gene_table)
    }
    matches = match(gl, gene_dt[[gt_from]])
    gl_new = gene_dt[[gt_to]][matches]

  } else { #use gprofiler
    convres = as.data.table(gconvert(gl, organism = "hsapiens",
                                     target = gp_target,
                                     numeric_ns = gp_numeric_ns,
                                     mthreshold = gp_mthreshold,
                                     filter_na = F)) #NA is filtered later if requested
    matches = match(gl, convres$input)
    gl_new = convres$name[matches]

    gl_new[gl_new=='None'] = NA #convert "None" returned by default for no-matches by gconvert
  }

  if(na_mode=='original'){
    gl_new[is.na(gl_new)] = gl[is.na(gl_new)]
  } else if(na_mode=='filter'){
    gl_new = gl_new[!is.na(gl_new)]
  }

  return(gl_new)
}


#### Pathway enrichment analysis #####
get_msigdb_pathways <- function(msigdb_sets=c('H', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7'), retVal='pathways', downgrade_table=NULL, apply_uniq_afer_downgrade=T){
  msdb_list = NULL
  for (myset in msigdb_sets){
    #parse if contains subcategory specification
    sc = NULL
    if(grepl(':', myset)){
      myset_tmp = sub(':', '#', myset)
      sc = strsplit(myset_tmp, '#')[[1]][2]
      myset = strsplit(myset_tmp, '#')[[1]][1]
    }
    #fetch pathways
    msigdb_df0 = msigdbr(species = "Homo sapiens", category=myset, subcategory = sc)
    msigdb_list0 = msigdb_df0 %>% split(x = .$gene_symbol, f = .$gs_name)
    msdb_list = c(msdb_list, msigdb_list0)
  }

  if(!is.null(downgrade_table)){
    if(length(class(downgrade_table))==1 && class(downgrade_table)=='character'){ #file was provided - avoid reading for each pathway
      gene_dt = fread(downgrade_table)
    }
    for(i in 1:length(msdb_list)){
      msdb_list[[names(msdb_list)[i]]] = update_gene_names(msdb_list[[names(msdb_list)[i]]], gene_table=downgrade_table,
                                                         na_mode = 'original', gt_from = 'gene_name_updated', gt_to='gene_name')
    }
    if(apply_uniq_afer_downgrade)
      msdb_list = sapply(msdb_list, unique)
  }

  if(retVal=='pathways'){
    return(msdb_list)
  } else if(retVal=='gene_list'){
    return(unique(unlist(msdb_list)))
  }
}

find_genes_missing_in_pathways <- function(gene_list, msigdb_sets=c('H', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7'),
                                           pathways=NULL, retVal='', verbose=T){
  if(!is.null(pathways)){ #test missing in list of pathways
    missing_genes = sort(unique(unlist(pathways))[!unique(unlist(pathways)) %in% gene_list])
    missing_genes_pc = round(length(unique(missing_genes)) / length(unique(unlist(pathways))) * 100, 2)
    if(verbose)
      print(paste0('missing percent for custom pathways: ', missing_genes_pc, '%'))
    if(retVal=='genesAndPercent'){
      return(list(missing_genes, missing_genes_pc))
    } else if(retVal=='genes'){
      return(missing_genes)
    }
  } else { #test missing in each msigdb specified
    for (myset in msigdb_sets){
      msigdb_df0 = msigdbr(species = "Homo sapiens", category=myset)
      pwl = unique(unlist(msigdb_df0 %>% split(x = .$gene_symbol, f = .$gs_name)))
      missing_genes = sort(unique(unlist(pwl))[!unique(unlist(pwl)) %in% gene_list])
      missing_genes_pc = round(length(unique(missing_genes)) / length(unique(unlist(pwl))) * 100, 2)
      print(paste0('missing percent for ', myset, ': ', missing_genes_pc, '%'))
      rm(msigdb_df0, pwl)
    }
  }
}


# pathway_list_to_gmt(pathways, 'test.gmt')
pathway_list_to_gmt <- function(pws, outfile=NULL, col2='name', file_format='gmt'){
  pws_char = rep('', length(pws))

  if(file_format=='flat' || is.null(col2)){ #no second column before actual genelists
    for (i in 1:length(pws)){
      pwchar=paste(names(pws)[i], paste(pws[[i]], collapse='\t'), sep='\t')
      pws_char[i] = pwchar
    }
  } else if(col2=='name'){ #print name in first two columns
    for (i in 1:length(pws)){
      pwchar=paste(names(pws)[i], names(pws)[i], paste(pws[[i]], collapse='\t'), sep='\t')
      pws_char[i] = pwchar
    }
  } else if(col2==''){ #second column empty
    for (i in 1:length(pws)){
      pwchar=paste(names(pws)[i], '', paste(pws[[i]], collapse='\t'), sep='\t')
      pws_char[i] = pwchar
    }
  } #future option: add other options here

  if(is.null(outfile)){ #return as gmt char (for testing)
    return(pws_char)
  } else { #to file
    fh<-file(outfile)
    writeLines(pws_char, fh)
    close(fh)
  }
}

#Returns data.table
pathway_list_to_dt <- function(pw_list, join_str=NULL){
  pwdt_list = list()
  if(!is.null(join_str)){
    res = lapply(pw_list, function(pw){paste0(pw, collapse=join_str)})
    pwdt = data.table(pathway=names(res), genes=unlist(res))
  } else {
    for(i in 1:length(pw_list)){
      pwdt0 = data.table(pathway=names(pw_list)[i],
                         gene=pw_list[[i]])
      pwdt_list[[i]] = pwdt0
    }
    pwdt = rbindlist(pwdt_list)
  }
  return(pwdt)
}

#default mode truncates column 2 of GMT
#file_format='flat'
read_gmt <- function(gmt_file, retmode='gslist', file_format='gmt'){
  pathway_list = readLines(gmt_file)
  mynames = unlist(str_extract(pattern = '([^\t]+)', pathway_list))
  if(file_format=='gmt'){
    pathway_list = lapply(pathway_list, function(x){return(sub(pattern='^[^\t]+\t[^\t]*\t', '', x))})
  } else if (file_format=='flat'){
    pathway_list = lapply(pathway_list, function(x){return(sub(pattern='^[^\t]*\t', '', x))})
  } else {
    print('Bad file format specified')
    return(NULL)
  }
  if(retmode=='gslist'){
    pathway_list = str_split(pathway_list, '\t')
    names(pathway_list) = mynames
    return(pathway_list)
  } else if(retmode=='genelist'){
    glist = unique(unlist(str_split(pathway_list, '\t')))
    return(glist)
  }
}

dt_to_geneset_list <- function(gs_dt, gs_col='geneset', gene_col='gene'){
  gs_list = list()
  for (gs in unique(gs_dt[,get(gs_col)])){
    gs_list[[gs]] = gs_dt[gs_dt[[gs_col]]==gs, get(gene_col)]
  }
  return(gs_list)
}

#get genes intersecting a list of genesets
get_geneset_intersect <- function(genelist, gslist=NULL, gmt_file=NULL, gs_names=NULL, join_str='|', ret_type='list'){
  if(!is.null(gmt_file)){
    gslist = read_gmt(gmt_file)
  }
  if(!is.null(gs_names)){
    gslist = gslist[names(gslist) %in% gs_names]
  }
  res = lapply(gslist, function(pw, genes=genelist){paste0(genes[genes %in% pw], collapse=join_str)})

  if (ret_type=='data.table'){
    res = data.table(pathway=names(res), genes=unlist(res))
  }
  return(res)
}

get_clusters_from_H_matrix <- function(H) {
  H <- H[rowSums(H)!=0,]
  H.norm <- apply(H,2,function(x) x/sum(x))
  g.Bayes <- apply(H.norm,2,function(x) which.max(x)) # select most probable cluster per sample
  return(g.Bayes)
}

get_bnmf_results <- function(bnmf_dir=CLUSTERING_BNMF_DIR ,verbose=F, retval='best_res', test_n=Inf, K=NULL,
                             all_file_suffix='res.L1EU.Bayes.all.txt', extract_regex = "res.L1EU.Bayes.(\\d+).RData"){

  bnmf_all_file = paste0(bnmf_dir ,'/', all_file_suffix)

  n.iter = max(as.numeric(sub(extract_regex, '\\1', list.files(bnmf_dir, pattern=extract_regex))))
  n.iter = min(n.iter, test_n)
  if(retval=='clusters_mat'){
    clusters_mat = NULL
  }

  if(!file.exists(bnmf_all_file) || retval=='clusters_mat'){
    # Retrieve df information
    tmpK <- rep(0,n.iter)
    tmpE <- rep(0,n.iter)
    for (i in 1:n.iter) {
      bnmf_run_res_file = paste(bnmf_dir, paste("res.L1EU.Bayes",i,"RData",sep="."),sep="/")
      if (file.exists(bnmf_run_res_file)){ #Files that aren't the best run may have been deleted
        if(verbose)
          print(paste('Loading file:', bnmf_run_res_file))
        #Get run scores
        load(file=bnmf_run_res_file)
        lambda <- res.Bayes[[5]]
        lambda <- unlist(lambda[length(lambda)])
        lambda <- lambda-min(lambda)
        if(verbose)
          cat(lambda,sum(lambda!=0),'\n')
        tmpK[i] <- sum(lambda > 0)
        tmpE[i] <- res.Bayes[[4]][length(res.Bayes[[4]])]

        #get cluster assignments
        if(retval=='clusters_mat'){
          if(is.null(clusters_mat)){
            n_samples = ncol(res.Bayes[[2]]) #cols of H == num samples
            clusters_mat = matrix(rep(0, n.iter * n_samples), ncol=n.iter, nrow=n_samples)
            row.names(clusters_mat) = colnames(res.Bayes[[2]])
          }
          clusters_mat[,i] = get_clusters_from_H_matrix(res.Bayes[[2]])
        }
      } else {
        if(verbose)
          print(paste('Skipping non-existing file:', bnmf_run_res_file))
      }
    }

    #Construct and name DF from run results
    df <- data.frame(seq(1:n.iter),tmpK,unlist(tmpE))
    colnames(df) <- c("run","K","evid")
  } else {
    df <- read.delim(bnmf_all_file, header=T, sep='\t')
  }
  if(retval=='run_df'){
    return(df)
  }

  #Select best
  if(retval=='best_res' || retval=='best_res_run'){
    df <- df[df$K>0,] #only runs that had Rdata files (I previously deleted these intermediate files of non-best solutions)
    df <- df[order(df$evid,decreasing=T),]

    if(!is.null(K)){
      k_selected = K
      out_print_pref = paste0('Loading data for specified K=', K, ', best run is: ')
    } else{
      best_ks = as.numeric(names(table(df$K))[table(df$K)==max(table(df$K))])
      if(length(best_ks)>1){
        cat("Multiple K values tied for most frequent - selecting minimal among:", best_ks, '\n')
      }
      k_selected = min(best_ks) #min selected to reduce complexity if several K were tied for most frequent
      out_print_pref = paste0('Loading data for run most frequent K (',k_selected,'), best run is: ')
    }
    if(!k_selected %in% df$K){
      cat('No run found for', K, '\n')
      return(NULL)
    }

    run.K = min(subset(df,evid==min(subset(df,K==k_selected)$evid))$run) #select by lowest evid; if tied for lowest, take first
    if(retval=='best_res_run')
      return(run.K)
    cat(out_print_pref, run.K, '\n')
    load(file=paste(bnmf_dir, paste("res.L1EU.Bayes",run.K,"RData",sep="."),sep="/"))
    return(res.Bayes)
  } else if(retval=='clusters_mat'){
    return(clusters_mat)
  }
}
