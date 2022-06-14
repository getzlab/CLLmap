# Utilities for assembling, manipulating and analyzing comut tables (sample by mutation tables)
# Author: Binyamin Knisbacher

library(data.table)
manual_gconvert <- function(glist, gl_from=NULL, gl_to=NULL, revert=F){
  if(is.null(gl_from)){
    gl_from = manual_gconvert_from
    gl_to = manual_gconvert_to
  }
  if(!revert)
    glist[glist %in% gl_from] = gl_to[match(glist, gl_from, nomatch=F)]
  else
    glist[glist %in% gl_to] = gl_from[match(glist, gl_to, nomatch=F)]
  return(glist)
}

get_drivers <- function(files, modes, mutsig_q_cutoff=0.1, clumps_q_cutoff=0.1, clumps_q_restricted_cutoff=0.1, gistic_q_cutoff=0.1,
                        retmode='uniq', verbose=T){
  drivers_list = list()
  for(i in 1:length(files)){

    driver_dt = fread(files[i])
    if(modes[i]=='mutsig'){
      drivers_list[[i]] = driver_dt[q < mutsig_q_cutoff,]$gene
    } else if (modes[i]=='clumps'){
      drivers_list[[i]] = driver_dt[CLUMPS_Q_FULL < clumps_q_cutoff,]$GENE_NAMES
    } else if (modes[i]=='clumps_or_cgc'){
      drivers_list[[i]] = driver_dt[CLUMPS_Q_FULL < clumps_q_cutoff | (CLUMPS_Q_RESTRICTED < clumps_q_restricted_cutoff & CLUMPS_Q_RESTRICTED > 0) ,]$GENE_NAMES
    }
    if(verbose){
      print(paste('Found', length(drivers_list[[i]]), 'drivers in', modes[i], '-', files[i]))
    }
  }

  if(retmode=='uniq'){
    uniq_drivers = unique(do.call(c, drivers_list))
    if(verbose){
      print(paste('Returning unique list of', length(uniq_drivers)))
    }
    return(uniq_drivers)
  } else if(retmode=='multilist') {
    return(drivers_list)
  } else if(retmode=='nonuniq') { #don't remove redundancy
    return(do.call(c, drivers_list))
  }
}

get_clumps_sig_per_sample <- function(clumps_file, clumps_mut_file, sig_type='clumps',
                                      clumps_q_cutoff=0.1, clumps_q_restricted_cutoff=0.1, sample_label='wes_tumor_sample_id',
                                      id_map=NULL, id_map_return_cols=c('participant_id', 'wes_tumor_sample_id') ,
                                      drop_dups_by_sample_gene_pair=F){
  #Get drivers
  clumps = fread(clumps_file)
  if (sig_type=='clumps'){
    clumps = clumps[CLUMPS_Q_FULL < clumps_q_cutoff,]
  } else if (sig_type=='clumps_cgc'){
    clumps = clumps[CLUMPS_Q_RESTRICTED < clumps_q_restricted_cutoff & CLUMPS_Q_RESTRICTED > 0 ,]
  } else if (sig_type=='clumps_or_cgc'){
    clumps = clumps[CLUMPS_Q_FULL < clumps_q_cutoff | (CLUMPS_Q_RESTRICTED < clumps_q_restricted_cutoff & CLUMPS_Q_RESTRICTED > 0) ,]
  }

  #Get mutations (has sample information)
  muts = fread(clumps_mut_file)
  muts$uniprot_id = sub('([^:]+):.*', '\\1', muts$uniprot_change)
  # muts$pos = sub('[^:]+:[A-Z](\\d+).*', '\\1', muts$uniprot_change) #shouldn't need this, just by uniprot_id
  muts = muts[muts$uniprot_id %in% clumps$UNIPROT_ID]
  muts = merge(muts[,c('patient','uniprot_id', 'uniprot_change')], clumps[,c('UNIPROT_ID', 'GENE_NAMES')], by.x='uniprot_id', by.y='UNIPROT_ID')
  muts = unique(muts[,c('patient', 'GENE_NAMES', 'uniprot_id', 'uniprot_change')]) #this could still return multiple events per single sample
  setnames(muts, c('patient', 'GENE_NAMES'), c(ifelse(!is.null(sample_label), sample_label, 'patient'), 'gene'))

  if(drop_dups_by_sample_gene_pair){
    muts = unique(muts, by=c(sample_label, 'gene'))
  }

  if(!is.null(id_map)){
    muts = merge(id_map[,unique(c(sample_label, id_map_return_cols)),with=F], muts, all.y=T, by=sample_label)
    if(! sample_label %in% id_map_return_cols)
      muts[,(sample_label):= NULL]
    muts = muts[order(gene, get(id_map_return_cols[1])),]
    muts = cbind(muts[,id_map_return_cols[1], with=F], muts[,!id_map_return_cols[1], with=F])
  } else { #different sort
    muts = muts[order(gene, get(sample_label)),]
  }

  return(muts)
}



preprocess_comut <- function(comut_file, id_map,
                             mutsig_sig_file=NULL, arm_sig_file=NULL,
                             q_cutoff=0.1,
                             id_map_return_cols = c('participant_id'),
                             drop_pair_col=T,
                             wes_pairs_to_participant=NULL,
                             sample_label='wes_tumor_sample_id',
                             id_map_linker=NULL,
                             drop_nontarget_id_cols=T,
                             annot_cols=1:9,
                             focal_input_type='default',
                             focal_name_mode='default',
                             arm_input_type='default',
                             comut_type='genes',
                             convert_to_ones=T){ #comut_type='focal_cnv'
  ### Allow for file or dt input for wes_pairs_to_participant
  if(length(class(wes_pairs_to_participant))==1 && class(wes_pairs_to_participant)=='character'){
    if(file.exists(wes_pairs_to_participant))
      wes_pairs_to_participant = fread(wes_pairs_to_participant)
    else
      stop(paste('pairs_to_participant file not found:', wes_pairs_to_participant))
  }

  if(comut_type=='genes'){
    comut0 = fread(comut_file)
    if(!is.null(mutsig_sig_file)){
      mutsig_sig_genes = fread(mutsig_sig_file)[q<q_cutoff,]$gene
      comut0 = comut0[comut0$gene %in% mutsig_sig_genes,]
    }
    cm_names=comut0[['gene']]
    setnames(comut0, 'gene', 'event')
    del_cols = c('event')
    comut_col_backup = comut0[,del_cols,with=F]
    id_map_linker = sample_label
  } else if (comut_type=='focal_cnv'){
    if(is.null(focal_input_type) || focal_input_type=='default'){
      comut0 = fread(comut_file)
      setnames(comut0, names(comut0)[annot_cols], gsub(' ', '.', names(comut0)[annot_cols]))
      comut0 = comut0[comut0$`q.values` < q_cutoff,]
      comut0 = comut0[Amplitude.Threshold != 'Actual Copy Change Given',] #remove duplicate entries holding actual values from bottom of file
      del_cols = colnames(comut0)[annot_cols]
      comut_col_backup = comut0[,del_cols,with=F]
      if(!is.null(focal_name_mode) && focal_name_mode=='event_col'){ #use Event col as is
        cm_names = comut_col_backup$Event
        cm_names = sub('^[Gg]ain|[Aa]mp', 'amp', cm_names)
        cm_names = sub('^[Dd]el|^[Ll]oss', 'del', cm_names)
      } else { #build from uniqe.name and descriptor cols
        cm_names = gsub('Amplification','Amp', gsub('Deletion', 'Del', gsub(' ', '', paste(sep="_", comut_col_backup$Unique.Name, gsub('_', '.', comut_col_backup$Descriptor) ))))
      }
    } else { #Ziao's fixed format that does not have the gistic event description columns
      #No changes needed
      comut0 = fread(comut_file)
      cm_names = colnames(comut0) #effectively no change
      del_cols = c() #nothing to del (this format has no description columns)
    }

    if(is.null(id_map_linker))
      id_map_linker = 'wes_participant' #sample_label is pair_id (used below)
  } else if (comut_type=='arm_cnv' || comut_type=='chrom_cnv'){ #do not change these mode names, as they're used in preprocess function
    if(arm_input_type != 'pair_by_event'){ #raw gistic input
      comut0 = preprocess_gistic_arm_level(comut_file, arm_sig_file=arm_sig_file, q_cutoff=q_cutoff, retval=comut_type)
      del_cols = colnames(comut0)[annot_cols] #*** Check if should be 1:3
    } else { #non-raw input (processed by Ziao)
      comut0 = preprocess_gistic_arm_level2(comut_file, arm_sig_file=arm_sig_file, q_cutoff=q_cutoff, retval=comut_type)
      del_cols = colnames(comut0)[1:3] #so length is 1
    }
    cm_names=comut0[['event']]
    if(is.null(id_map_linker))
      id_map_linker = 'wes_participant' #sample_label is pair_id (used below)
  }

  #transpose (from sample columns to sample rows)
  if(focal_input_type != 'fixed_no_event_cols'){
    if(length(del_cols)>0){
      comut = transpose(comut0)[-1:-length(del_cols),]
    } else {
      comut = transpose(comut0)[-1,]
    }
    setnames(comut, cm_names)

    comut[, names(comut) := lapply(.SD, as.numeric)] #make all cols numeric

    #link to sample data
    if(sample_label=='pair_id'){ #current format for arm_cnv and chrom_cnv (need special file to link pair_id to master table)
      comut[,'pair_id'] = colnames(comut0)[-1:-length(del_cols)]
      comut = merge(comut, wes_pairs_to_participant[,c('pair_id', 'participant'),with=F], by='pair_id')
      if(drop_nontarget_id_cols)
        comut[,pair_id:=NULL]
      setnames(comut, 'participant', id_map_linker)
    } else { #use master table
      comut[,id_map_linker] = colnames(comut0)[-1:-length(del_cols)]
      #do nothing, use sample_label below (which must be valid column of master table)
    }

  } else {
    comut = comut0
    id_col_dt = comut[,'pair_id']
    comut = comut[,-'pair_id',with=F]
    comut[, names(comut) := lapply(.SD, as.numeric)] #make all cols numeric
    comut = cbind(id_col_dt, comut)
    setnames(comut, 'pair_id', id_map_linker)
  }

  if(convert_to_ones){
    cols = colnames(comut)[sapply(comut, is.numeric)]
    comut[ , (cols) := lapply(.SD, function(x){ifelse(x>0,1,0)}), .SDcols = cols]
  }

  #Add id_map columns of interest
  comut = merge(comut, id_map[,unique(c(id_map_linker, id_map_return_cols)),with=F], by=id_map_linker, all.x=T)

  if(drop_nontarget_id_cols){
    cols_to_drop = colnames(comut)[colnames(comut) %in% c(id_map_linker, 'wgs_participant', 'wes_participant', 'wes_participant_id', 'wes_tumor_sample_id') &
                                     (!colnames(comut) %in% id_map_return_cols)]
    if(length(cols_to_drop)>0){
      if(length(cols_to_drop)==1)
        comut[, (cols_to_drop):=NULL]
      else
        comut[, cols_to_drop := NULL]
    }
  }

  #place target label in first column
  comut = cbind(comut[,id_map_return_cols,with=F], comut[,-id_map_return_cols,with=F])

  return(comut)
}


parse_maf <- function(maf_file, maf_dt=NULL, drop_silent=F,
                      id_map_file=NULL, id_linker_col='wes_tumor_sample_id', id_target_col='participant_id',
                      class_replace_pairs=NULL, no_class_underscores=T,
                      classes_to_drop=c('Targeted_Region', 'Unknown'),
                      classes_to_drop_renamed=c('Synonymous', 'Non-coding'),
                      classes_to_retain_renamed=NULL,
                      variant_classes_ordered=NULL,
                      maf_table_as_int=F, maf_table_fill_val=0, sample_by_gene=T,
                      maf_cols=c('Hugo_Symbol', 'Variant_Classification'), gene_col='Hugo_Symbol', gene_col_return_name='gene',
                      silent_regex="Silent|^Synonymous|3'\\-?UTR|5'\\-?UTR|3'\\-?Flank|5'\\-?Flank|IGR|Intron|RNA",
                      gene_retain_list=NULL, maf_sample_col='Tumor_Sample_Barcode',
                      retval='maf', verbose=T){

  if(is.null(maf_dt)){ #read maf_file only if maf_dt is not given
    if(!is.null(maf_cols)) #drop duplicates by maf_cols
      maf = fread(maf_file)[,unique(c(maf_sample_col, maf_cols)),with=F]
    else #no duplicate filtering upon reading
      maf = fread(maf_file)
  } else {
    maf = maf_dt
  }

  #Dealt with differently
  if(drop_silent){
    maf = maf[!grepl(silent_regex, maf$Variant_Classification),]
  }

  if(!is.null(id_map_file)){
    id_map = fread(id_map_file)
    maf = merge(id_map[,unique(c(id_target_col, id_linker_col)),with=F], maf, by.x=id_linker_col, by.y=maf_sample_col, all.y=T)
    if(id_linker_col != id_target_col)
      maf[,(id_linker_col) := NULL]
  }

  if(!is.null(classes_to_drop)){
    for (to_drop in classes_to_drop){
      n_pre = nrow(maf)
      maf = maf[maf$Variant_Classification != to_drop,]
      n_post = nrow(maf)
      if(verbose) cat('Dropped from MAF', n_pre-n_post , 'instances of the following (before rename):', to_drop, '\n')
    }
  }

  if(!is.null(class_replace_pairs)){
    for(sub_pair in class_replace_pairs){
      maf$Variant_Classification = sub(sub_pair[1], sub_pair[2], maf$Variant_Classification, ignore.case = T)
    }
  }

  if(!is.null(classes_to_drop_renamed)){
    for (to_drop in classes_to_drop_renamed){
      n_pre = nrow(maf)
      maf = maf[maf$Variant_Classification != to_drop,]
      n_post = nrow(maf)
      if(verbose) cat('Dropped from MAF', n_pre-n_post , 'instances of the following (after rename):', to_drop, '\n')
    }
  }

  if(!is.null(classes_to_retain_renamed)){
      n_pre = nrow(maf)
      maf = maf[maf$Variant_Classification %in% classes_to_retain_renamed,]
      n_post = nrow(maf)
      if(verbose) cat('Dropped from MAF', n_pre-n_post , 'instances not in retain category list:', paste(classes_to_retain_renamed, collapse=', '), '\n')
  }

  if(gene_col %in% colnames(maf)){
    setnames(maf, gene_col, gene_col_return_name) #rename gene column
    if(!is.null(gene_retain_list)) #retain gene names of interest
      maf = maf[maf[[gene_col_return_name]] %in% gene_retain_list]
  } else {
    print(paste('No ',gene_col,'column found in maf - returning NULL'))
    return(NULL)
  }

  if(no_class_underscores){
    maf$Variant_Classification = gsub('_', ' ', maf$Variant_Classification)
  }

  if(retval=='maf_table'){ #comut but with string describing mutation type
    if(!is.null(variant_classes_ordered)){ #choose by precedence list (if not found in classification - chooses randomly)
      if(no_class_underscores)
        variant_classes_ordered = gsub('_', ' ', variant_classes_ordered)
      maf_table = dcast.data.table(maf, formula=as.formula(ifelse(sample_by_gene, 'participant_id ~ gene', 'gene ~ participant_id')),
                        value.var='Variant_Classification', fill=maf_table_fill_val,
                        fun.aggregate = function(x){ifelse(maf_table_as_int,
                                                           min(match(x, variant_classes_ordered)),
                                                           variant_classes_ordered[min(match(x, variant_classes_ordered))])})
      return(maf_table)
    } else { #choose 1 randomly
      maf_table = dcast.data.table(maf, formula=as.formula(ifelse(sample_by_gene, 'participant_id ~ gene', 'gene ~ participant_id')),
                        value.var='Variant_Classification',
                        fun.aggregate = function(x){sample(x,1)}, fill=0)
      return(maf_table)
    }
  } else if (retval=='comut') {
    maf$Variant_Classification = 1 #return 1 instead of label
    if(sample_by_gene){
      comut = dcast.data.table(maf, formula=as.formula('gene ~ participant_id'), value.var='Variant_Classification',
                    fun.aggregate = function(x){x[1]}, fill=0)
    } else {
      comut = dcast.data.table(maf, formula=as.formula('participant_id ~ gene'), value.var='Variant_Classification',
                    fun.aggregate = function(x){x[1]}, fill=0)
    }
    return(comut)
  } else if (retval=='maf_syn_counts'){
    maf$metaclass = ifelse(grepl(silent_regex , maf$Variant_Classification), 'Syn.', 'Non Syn.')
    maf_syn_counts = dcast.data.table(maf[,.N,by=c('participant_id', 'metaclass')], 'participant_id ~ metaclass', fill = 0, value.var = 'N')
    return(maf_syn_counts)
  } else if (retval=='maf'){
    return(maf)
  } else {
    return(maf)
  }
}

#q_cutoff used only if arm_sig_file is provided
preprocess_gistic_arm_level <- function(raw_arm_level_file, amp_cutoff=0, del_cutoff=0, arm_sig_file=NULL, q_cutoff=0.1, retval='both'){
  if(!is.null(arm_sig_file)){
    arm_sig = fread(arm_sig_file)
    arm_sig_amp = arm_sig[arm_sig$`Amp q-value` < q_cutoff, ]$Arm
    arm_sig_del = arm_sig[arm_sig$`Del q-value` < q_cutoff, ]$Arm
  }

  #convert all vals to -1, 0, 1
  cm = fread(raw_arm_level_file)
  numeric_cols = names(cm[,sapply(cm,is.numeric),with=F])

  # cm = cm[, lapply(.SD, function(x){ifelse(x==0, 0, x/x*sign(x))}), .SDcols=numeric_cols]
  #retain all amps as 1 in amp dt, other to 0
  cm_amp = copy(cm)
  if(!is.null(arm_sig_file)){
    cm_amp = cm_amp[cm_amp$`Chromosome Arm` %in% arm_sig_amp,]
  }
  cm_amp$event = paste0('amp(',cm_amp[['Chromosome Arm']], ')')
  # cm_amp = cm_amp[, lapply(.SD, function(x){ifelse(x==1, 1, 0)}), .SDcols=numeric_cols]
  for (j in numeric_cols)
    data.table::set(cm_amp,which(cm_amp[[j]] <= amp_cutoff),j,0)

  #convert all dels to 1 in del dt, other to 0
  cm_del = copy(cm)
  if(!is.null(arm_sig_file)){
    cm_del = cm_del[cm_del$`Chromosome Arm` %in% arm_sig_del,]
  } else {
    cat('No significance filter in Arm/Chrom level CNV!\n')
  }
  cm_del$event = paste0('del(',cm_del[['Chromosome Arm']], ')')
  # cm_del = cm_del[, lapply(.SD, function(x){ifelse(x==-1, 1, 0)}), .SDcols=numeric_cols]
  for (j in numeric_cols)
    data.table::set(cm_del,which(cm_del[[j]] >= del_cutoff),j, 0)

  cm_del$event = paste0('del(',cm_del[['Chromosome Arm']], ')')

  #join arm-level amp and del
  cm = rbind(cm_amp, cm_del)
  for (j in numeric_cols)
    cm[[j]] = ifelse(cm[[j]]==0, 0, cm[[j]]/cm[[j]])

  cm$chrom = sub('[^0-9]*(\\d+).*', '\\1', cm$event)
  cm$cnv_type = ifelse(grepl('^[Dd]el|^[Ll]oss', cm$event), 'del', 'amp')
  cm_chrom = cm[,lapply(.SD, sum), by=.(chrom, cnv_type), .SDcols=which(sapply(cm, is.numeric))]

  cm_chrom$event = paste0(ifelse(cm_chrom$cnv_type=='amp', 'amp', 'del'), '(', cm_chrom$chrom, ')')
  for (j in numeric_cols)
    cm_chrom[[j]] = ifelse(cm_chrom[[j]]==2, 1, 0)

  #place non-numeric columns first
  cm = cm[,c(names(cm[,!sapply(cm, is.numeric), with=F]), names(cm[,sapply(cm, is.numeric), with=F])),with=F]
  cm_chrom = cm_chrom[,c(names(cm_chrom[,!sapply(cm_chrom, is.numeric), with=F]), names(cm_chrom[,sapply(cm_chrom, is.numeric), with=F])),with=F]

  if(retval=='both'){
    return(list(cm_arm=cm, cm_chrom=cm_chrom))
  } else if(retval=='arm_cnv'){
    return(cm)
  } else if(retval=='chrom_cnv'){
    return(cm_chrom)
  }
}


#The format for preprocess_gistic_arm_level2 is pair-id by event in format as in "del_1p" "gain_1p" (all values are 0 or 1)
preprocess_gistic_arm_level2 <- function(modified_arm_level_file, arm_sig_file=NULL, q_cutoff=0.1, retval='both'){

  if(!is.null(arm_sig_file)){
    arm_sig = fread(arm_sig_file)
    arm_sig_amp = paste0('amp_', arm_sig[arm_sig$`Amp q-value` < q_cutoff, ]$Arm)
    arm_sig_del = paste0('del_', arm_sig[arm_sig$`Del q-value` < q_cutoff, ]$Arm)
    arm_sig_all = c(arm_sig_amp, arm_sig_del)
  }

  cm0 = fread(modified_arm_level_file)
  event_names = colnames(cm0)[-1]
  cm = transpose(cm0) #there is only one annotation column to drop in this format
  setnames(cm, as.character(cm[1,]))
  cm = cm[-1,]
  cm[, names(cm) := lapply(.SD, as.numeric)] #make all cols numeric
  cm = cbind(event=event_names, cm)
  cm$event = sub('gain_', 'amp_', cm$event)
  # print('Printing events')
  # print(cm$event)
  numeric_cols = names(cm[,sapply(cm,is.numeric),with=F]) #for chrom events

  if(!is.null(arm_sig_file)){
    cm = cm[cm$event %in% arm_sig_all,]
  } else {
    cat('No significance filter in Arm/Chrom level CNV!\n')
  }

  cm$chrom = sub('[^0-9]*(\\d+).*', '\\1', cm$event)
  cm$cnv_type = ifelse(grepl('^[Dd]el|^[Ll]oss', cm$event), 'del', 'amp')
  cm_chrom = cm[,lapply(.SD, sum), by=.(chrom, cnv_type), .SDcols=which(sapply(cm, is.numeric))]
  # cm_chrom$event = paste0(ifelse(cm_chrom$cnv_type=='amp', 'gain', 'loss'), '(', cm_chrom$chrom, ')')
  for (j in numeric_cols)
    cm_chrom[[j]] = ifelse(cm_chrom[[j]]==2, 1, 0)
  cm_chrom$event = paste0(ifelse(cm_chrom$cnv_type=='amp', 'amp', 'del'), '(', cm_chrom$chrom, ')')

  #place non-numeric columns first
  cm = cm[,c(names(cm[,!sapply(cm, is.numeric), with=F]), names(cm[,sapply(cm, is.numeric), with=F])),with=F]
  cm_chrom = cm_chrom[,c(names(cm_chrom[,!sapply(cm_chrom, is.numeric), with=F]), names(cm_chrom[,sapply(cm_chrom, is.numeric), with=F])),with=F]
  if(retval=='both'){
    return(list(cm_arm=cm, cm_chrom=cm_chrom))
  } else if(retval=='arm_cnv'){
    return(cm)
  } else if(retval=='chrom_cnv'){
    return(cm_chrom)
  }
}

#This works only if the columns are unique in the two comuts other than the 'by' column
merge_comuts <- function(comut_list, participant_col='participant_id', filter_nonnumeric=T){
  for(i in 1:length(comut_list)){print(paste('comut', i, 'in merge has dim:', dim(comut_list[[i]])));  print(comut_list[[i]][1:5,1:min(5,ncol(comut_list[[i]]))])} #***
  mymerge <- function(x,y){merge(x,y,all=T, by=participant_col)}
  comut_merged = Reduce(merge,comut_list)
  if(any(is.na(comut_merged))){
    print('Check why you have NAs in your merge - incomplete overlap in participant_id per comut!')
    print('Not merging - returning list of comuts!')
    return(comut_list)
  }

  if(filter_nonnumeric){
    cols_to_retain = c(participant_col, colnames(comut_merged[,sapply(comut_merged, is.numeric),with=F])[colSums(comut_merged[,sapply(comut_merged, is.numeric),with=F])!=0])
    comut_merged = comut_merged[,..cols_to_retain]
  }
  return(comut_merged)
}


#Note hard coding for 'participant_id' and 'wes_participant' columns in id_map
#Currently requires to have 'participant_id' column in id_map, if it is used and manual participant name is conducted
get_cnv_comut <- function(focal_file=NULL, arm_file=NULL, arm_sig_file=NULL, q_cutoff = 0.1, arm_level=T, chrom_level=T,
                          chrom_level_list=NULL,
              id_map=NULL, id_map_return_cols=c('participant_id', 'wes_participant'), id_map_linker=NULL,
              sample_label='pair_id', wes_pairs_to_participant=NULL,
              focal_annot_col_indexes=1:9, arm_annot_col_indexes=1:4, chrom_annot_col_indexes=1:3,
              focal_input_type='default',
              focal_name_mode='default', arm_input_type='default',

              manual_patid_convert_from=NULL, manual_patid_convert_to=NULL,
              comut_event_drop=NULL,
              verbose=T){
  if(length(class(wes_pairs_to_participant))==1 && class(wes_pairs_to_participant)=='character'){
    if(file.exists(wes_pairs_to_participant))
      wes_pairs_to_participant = fread(wes_pairs_to_participant)
    else
      stop(paste('wes_pairs_to_participant file not found:', wes_pairs_to_participant))
  }

  if(length(class(id_map))==1 && class(id_map)=='character'){
    if(file.exists(id_map))
      id_map = fread(id_map)
    else
      stop(paste('id_map file not found:', id_map))
  }

  if(is.null(focal_file) && is.null(arm_file)){
    stop('In get_cnv_comut: focal_file or arm_file must not be NULL')
  }
  comut_list = list()
  #focal
  if(!is.null(focal_file)){
    if(verbose)
      print('Parsing focal CNV')
    comut_focal = preprocess_comut(focal_file, q_cutoff=q_cutoff, id_map=id_map, id_map_return_cols = id_map_return_cols,
                                   sample_label=sample_label, annot_cols=focal_annot_col_indexes, comut_type='focal_cnv',
                                   focal_input_type=focal_input_type,
                                   focal_name_mode=focal_name_mode,
                                   id_map_linker=id_map_linker,
                                   wes_pairs_to_participant=wes_pairs_to_participant)
    if(!is.null(manual_patid_convert_from) && !is.null(manual_patid_convert_to)){
      mymatch = match(comut_focal$wes_participant[is.na(comut_focal$participant_id)], manual_patid_convert_from)
      if(any(!is.na(mymatch))){
        print(paste('Converting participant ids!', paste(comut_focal$wes_participant[is.na(comut_focal$participant_id)], collapse=" ")))
        comut_focal$participant_id[is.na(comut_focal$participant_id)] = manual_patid_convert_to[mymatch]
      } else {
      print(paste('Not converting participant ids! - WES samples with NA participants are: ', paste(comut_focal$wes_participant[is.na(comut_focal$participant_id)], collapse=" ")))
      }
    }
    if('wes_participant' %in% colnames(comut_focal)){
      comut_focal[,c('wes_participant'):=NULL]
    }
    if('wgs_participant' %in% colnames(comut_focal)){
      comut_focal[,c('wgs_participant'):=NULL]
    }
    comut_list[[length(comut_list)+1]] = comut_focal
  }

  if(!is.null(arm_file)){
    if(verbose)
      print('Parsing arm-level CNV')
    #arm
    if(arm_level){
      comut_arm = preprocess_comut(arm_file, arm_sig_file=arm_sig_file, q_cutoff=q_cutoff, id_map=id_map, id_map_return_cols = id_map_return_cols,
                                   sample_label=sample_label, annot_cols=arm_annot_col_indexes, comut_type='arm_cnv', arm_input_type=arm_input_type,
                                   id_map_linker=id_map_linker,
                                   wes_pairs_to_participant=wes_pairs_to_participant)
      if(!is.null(manual_patid_convert_from) && !is.null(manual_patid_convert_to)){
        mymatch = match(comut_arm$wes_participant[is.na(comut_arm$participant_id)], manual_patid_convert_from)
        if(any(!is.na(mymatch)))
          comut_arm$participant_id[is.na(comut_arm$participant_id)] = manual_patid_convert_to[mymatch]
      }
      if('wes_participant' %in% colnames(comut_arm)){
        comut_arm[,c('wes_participant'):=NULL]
      }
      if('wgs_participant' %in% colnames(comut_arm)){
        comut_arm[,c('wgs_participant'):=NULL]
      }
      comut_list[[length(comut_list)+1]] = comut_arm
    }

    #chrom
    if(chrom_level){
      if(verbose)
        print('Generalizing arm-level CNV to chrom CNV')
      comut_chrom = preprocess_comut(arm_file, arm_sig_file=arm_sig_file, q_cutoff=q_cutoff, id_map=id_map, id_map_return_cols = id_map_return_cols,
                                     sample_label=sample_label, annot_cols=chrom_annot_col_indexes, comut_type='chrom_cnv', arm_input_type=arm_input_type,
                                     id_map_linker=id_map_linker,
                                     wes_pairs_to_participant=wes_pairs_to_participant)

      if(!is.null(chrom_level_list)){
        chrom_to_del = grep('^([Aa]mp|[Dd]el|[gG]ain|[Ll]oss)', colnames(comut_chrom), value=T)
        chrom_to_del = chrom_to_del[!chrom_to_del %in% chrom_level_list]
        comut_chrom = comut_chrom[,-chrom_to_del,with=F]
      }

      if(!is.null(manual_patid_convert_from) && !is.null(manual_patid_convert_to)){
        mymatch = match(comut_chrom$wes_participant[is.na(comut_chrom$participant_id)], manual_patid_convert_from)
        if(any(!is.na(mymatch)))
          comut_chrom$participant_id[is.na(comut_chrom$participant_id)] = manual_patid_convert_to[mymatch]
      }
      if('wes_participant' %in% colnames(comut_chrom)){
        comut_chrom[,c('wes_participant'):=NULL]
      }
      if('wgs_participant' %in% colnames(comut_chrom)){
        comut_chrom[,c('wgs_participant'):=NULL]
      }
      comut_list[[length(comut_list)+1]] = comut_chrom
    }
  }

  # Merge CNV comuts
  if(length(comut_list)>1){
    merged_comut = merge_comuts(comut_list)
    if(!is.null(comut_event_drop))
      merged_comut = merged_comut[, colnames(merged_comut)[!colnames(merged_comut) %in% comut_event_drop], with=F]
    return(merged_comut)
  } else {
    if(!is.null(comut_event_drop)) #typically to filter redundant events (e.g. 12q,12p amp when there is 12 amp)
      comut_list[[1]] = comut_list[[1]][, colnames(comut_list[[1]])[!colnames(comut_list[[1]]) %in% comut_event_drop], with=F]
    return(comut_list[[1]])
  }
}


#Needs both to have
merge_and_parse_mafs <- function(wes_maf_file, wgs_maf_file, targeted_maf_file = NULL,
                                            id_map_file=NULL, gene_retain_list=NULL,
                                            class_replace_pairs=NULL, variant_classes_ordered=NULL,
                                            prefer_dna_type='wes', retval='comut', sample_by_gene=F, verbose=T){

  if(verbose)
    print('Parsing WES maf')
  wes_maf_dt = parse_maf(wes_maf_file, id_map_file = id_map_file,
                        class_replace_pairs=class_replace_pairs, no_class_underscores = T,
                        variant_classes_ordered=variant_classes_ordered, retval='maf',
                        gene_retain_list=gene_retain_list, id_linker_col='wes_tumor_sample_id')
  print("dim(wes_maf_dt)")
  print(dim(wes_maf_dt))

  if(verbose)
    print('Parsing WGS maf')
  wgs_maf_dt = parse_maf(wgs_maf_file, id_map_file = id_map_file,
                         class_replace_pairs=class_replace_pairs, no_class_underscores = T,
                         variant_classes_ordered=variant_classes_ordered, retval='maf',
                         # classes_to_drop_renamed=c('Synonymous'), #Do not drop Non-coding
                         gene_retain_list=gene_retain_list, id_linker_col='wgs_tumor_sample_id')

  print("dim(wgs_maf_dt)")
  print(dim(wgs_maf_dt))

  if(!is.null(targeted_maf_file)){
    if(verbose)
      print('Parsing targeted maf')
    targeted_dt = parse_maf(targeted_maf_file,
                           class_replace_pairs=class_replace_pairs, no_class_underscores = T,
                           variant_classes_ordered=variant_classes_ordered, retval='maf',
                           classes_to_drop_renamed=c('Synonymous'), #Do not drop Non-coding
                           maf_sample_col='participant_id',
                           gene_retain_list=gene_retain_list) #no need to convert input
  } else {
    targeted_dt = NULL
  }

  # Check that ID map conversion worked for all samples
  if(any(is.na(wgs_maf_dt$participant_id)) || any(is.na(wes_maf_dt$participant_id))){
    print('Conversion map for MAF is incomplete - fix this! Returning NULL')
    return(NULL)
  }

  ## What data to use ("prefer") for participant with WES and WGS
  if(prefer_dna_type=='wes'){ #Prefer WES - subset WGS results
    wgs_maf_dt = wgs_maf_dt[wgs_maf_dt$participant_id %in% wes_maf_dt$participant_id,]
  } else if(prefer_dna_type=='wgs'){ #Prefer WGS - subset WES results
    wes_maf_dt = wes_maf_dt[wes_maf_dt$participant_id %in% wgs_maf_dt$participant_id,]
  } else if(prefer_dna_type=='both') { #No subsetting
  } else { #No subsetting
  }

  merged_maf = rbind(wes_maf_dt, wgs_maf_dt, targeted_dt) #handles NULL properly too

  if(verbose)
    print('Parsing merged WES + WGS maf')
  merged_maf_parsed = parse_maf(maf_file=NULL, maf_dt=merged_maf, drop_silent=F,
            id_map_file=NULL, #already converted
            class_replace_pairs=class_replace_pairs, no_class_underscores = T,
            variant_classes_ordered=variant_classes_ordered,
            classes_to_drop_renamed=c('Synonymous'), #Do not drop Non-coding
            sample_by_gene=sample_by_gene, #To enable merge by participant_id column
            maf_cols=c('gene', 'Variant_Classification'), #already modified
            gene_col='gene',
            gene_col_return_name='gene',
            maf_sample_col='participant_id',
            retval=retval, #retval='comut'
            gene_retain_list = NULL, #already applied above
            verbose=T)

  return(merged_maf_parsed)
}


### Function to merge focal_del and focal_amp from "fixed" output from Ziao's algorithm that corrects potential ZD bug
merge_focal_fixed_output <- function(focal_amp_file, focal_del_file, negative_to_zero_cutoff=0,
                                     id_col = 'pair_id',
                                     amp_retain_regex='pair_id|gain',
                                     del_retain_regex='pair_id|del',
                                     assign_zero_negative_cutoff=0){
  ### Get focal amps
  amps = fread(focal_amp_file)
  amps = amps[, grepl(amp_retain_regex, colnames(amps)), with=F]
  for (nc in colnames(amps)[sapply(amps,is.numeric)])
    data.table::set(amps,which(amps[[nc]] < assign_zero_negative_cutoff),nc,0)
  colnames(amps) = sub('^[Gg]ain', 'amp', colnames(amps))

  ### Get focal dels
  dels = fread(focal_del_file)
  dels = dels[, grepl(del_retain_regex, colnames(dels)), with=F]
  for (nc in colnames(dels)[sapply(dels,is.numeric)])
    data.table::set(dels,which(dels[[nc]] < assign_zero_negative_cutoff),nc,0)


  #Merge and Write to file
  merged_cnv_comut = merge(amps, dels, by='pair_id')
  outfile = sub('_amp_', '_amp_and_del_', focal_amp_file)
  fwrite(merged_cnv_comut, outfile, sep='\t')

  return(outfile)
}

#Will work if supply only WES or WGS too
get_patient_counts_and_rates <- function(wes_counts_file=NULL, wgs_counts_file=NULL, prefer='wes', id_map_file=NULL,
                                         wes_linker='wes_tumor_sample_id', wgs_linker='wgs_tumor_sample_id',
                                         wes_blacklist=NULL, wgs_blacklist=NULL,
                                         id_map_target='participant_id', with_dnaseq_type=T, dnaseq_col='dnaseq_type', name_col_rename='participant_id',
                                         cols_original_keep=c('rate_sil', 'rate_non', 'rate_tot', 'N_tot'),
                                         add_cols_for_comut=T, log_base=10,
                                         outfile=NULL, retval='dt'){
  #Get WES rates
  if(!is.null(wes_counts_file)){
    wes_counts = fread(wes_counts_file)
    if(!is.null(wes_blacklist)){
      wes_counts = wes_counts[!wes_counts$name %in% wes_blacklist]
    }
    if(!is.null(id_map_file)){
      id_map = fread(id_map_file)[,c(wes_linker, id_map_target),with=F]
      wes_counts[,'name'] = id_map[match(wes_counts[['name']], id_map[[wes_linker]])][[id_map_target]]
    }
    if(with_dnaseq_type)
      wes_counts[, dnaseq_col] = 'WES'
  } else {
    wes_counts = NULL
  }

  #Get WGS rates
  if(!is.null(wgs_counts_file)){
    wgs_counts = fread(wgs_counts_file)
    if(!is.null(wgs_blacklist)){
      wgs_counts = wgs_counts[!wgs_counts$name %in% wgs_blacklist]
    }
    if(!is.null(id_map_file)){
      id_map = fread(id_map_file)[,c(wgs_linker, id_map_target),with=F]
      wgs_counts[,'name'] = id_map[match(wgs_counts[['name']], id_map[[wgs_linker]])][[id_map_target]]
    }
    if(with_dnaseq_type)
      wgs_counts[, dnaseq_col] = 'WGS'
  } else {
    wgs_counts = NULL
  }

  if(prefer=='wes' && !is.null(wes_counts) && !is.null(wgs_counts)){
    all_counts = rbind(wes_counts, wgs_counts[!wgs_counts$name %in% wes_counts$name])
  } else if (prefer=='wgs' && !is.null(wes_counts)  && !is.null(wgs_counts)) { #prefer wgs
    all_counts = rbind(wgs_counts, wes_counts[!wes_counts$name %in% wgs_counts$name])
  } else { #concatenate with redundancy (unless only one of WES/WGS is specified)
    all_counts = rbind(wes_counts, wgs_counts)
  }
  if(!is.null(name_col_rename))
    setnames(all_counts, 'name', id_map_target)

  ## Keep cols of interest
  if(with_dnaseq_type)
    cols_keep = c(id_map_target, cols_original_keep, dnaseq_col)
  else
    cols_keep = c(id_map_target, cols_original_keep, dnaseq_col)
  all_counts = all_counts[, cols_keep, with=F]

  ## Add cols for comut (calculations for viz)
  if(add_cols_for_comut){
    all_counts[['rate_sil_frac']] = all_counts$rate_sil / all_counts$rate_tot
    all_counts[['rate_non_frac']] = all_counts$rate_non / all_counts$rate_tot
    raise_power = max(ceiling(-log(all_counts$rate_tot, log_base)))+1
    all_counts$incremented_by = log_base^raise_power
    all_counts[['raise_power']] = raise_power
    all_counts[['rate_tot_inc_log']] = log(all_counts$rate_tot * log_base^raise_power, base=log_base)
    #used for stacked bar plot in log scale (need to subtract raise_power in axis ticks)
    all_counts[['rate_tot_inc_log_sil_frac']] = all_counts[['rate_tot_inc_log']] * all_counts[['rate_sil_frac']]
    all_counts[['rate_tot_inc_log_non_frac']] = all_counts[['rate_tot_inc_log']] * all_counts[['rate_non_frac']]
    #used for point plots (corrected: removed the raised power)
    all_counts[['rate_tot_inc_log_per_mb']] = log(all_counts$rate_tot * log_base^raise_power, base=log_base) - raise_power + log(1e+6, log_base)
    all_counts[['rate_sil_inc_log_per_mb']] = log(all_counts$rate_sil * log_base^raise_power, base=log_base) - raise_power + log(1e+6, log_base)
    all_counts[['rate_non_inc_log_per_mb']] = log(all_counts$rate_non * log_base^raise_power, base=log_base) - raise_power + log(1e+6, log_base)
    per_mb_cols = c('rate_tot_inc_log_per_mb', 'rate_sil_inc_log_per_mb', 'rate_non_inc_log_per_mb')
    all_counts[, (per_mb_cols) := replace(.SD, .SD == -Inf, 0), .SDcols = per_mb_cols] #correct any log(0)

    #For overall mutation rate per mb (log10)
    all_counts[['rate_tot_per_mb_log']] = log(all_counts[['rate_tot']], base=log_base) - log(1e+6, log_base)
  }

  if(!is.null(outfile))
    fwrite(all_counts, outfile, sep='\t')

  if(!is.null(retval) && retval=='dt')
    return(all_counts)
}


### Filter comuts with following optional filters:
#event blacklist
#participant blacklist
#participant inclusion file (e.g., if want to include only U-CLL/M-CLL patients)
filter_comut <- function(comut, comut_label='comut', blacklist=NULL,
                         pat_blacklist=NULL, pat_inclusion_file=NULL,
                         pat_col='participant_id',
                         pat_extend_missing = NULL,
                         drop_zero_sum_events=F,
                         cna_regex = '^[Aa]mp|^[Tt]ri|^[Gg]ain|^[Dd]el|^[Ll]oss',
                         min_cna_pc = NULL,
                         min_driver_pc = NULL,
                         min_pc_force_retain_cols = NULL,
                         to_binary=F,
                         outfile=NULL, retval='comut'){
  comut_filt = copy(comut)

  if(!is.null(blacklist)){
    print(paste(comut_label, '-', 'Filtering following events found in blacklist:'))
    print(colnames(comut_filt)[colnames(comut_filt) %in% blacklist])
    comut_filt = comut_filt[, !colnames(comut_filt) %in% blacklist, with=F]
  }

  if(!is.null(pat_blacklist)){
    print(paste(comut_label, '-', 'Filtering following participants found in blacklist:'))
    print(comut_filt[[pat_col]][comut_filt[[pat_col]] %in% pat_blacklist])
    comut_filt = comut_filt[!comut_filt[[pat_col]] %in% pat_blacklist]
  }

  if(!is.null(pat_inclusion_file)){
    print(paste(comut_label, '-', 'Filtering following participants NOT found in participant inclusion file:'))
    print(comut_filt[[pat_col]][!comut_filt[[pat_col]] %in% pat_inclusion_file])
    comut_filt = comut_filt[comut_filt[[pat_col]] %in% pat_inclusion_file]
  }

  if(!is.null(pat_extend_missing)){
    comut_filt = merge(comut_filt, data.table(participant_id=pat_extend_missing), all=T, by='participant_id')
    comut_filt[is.na(comut_filt)] = 0
  }

  if(drop_zero_sum_events){
    participant_col = which(colnames(comut_filt)=='participant_id')
    comut_filt_numeric = comut_filt[,sapply(comut_filt, is.numeric), with=F]
    sum_zero_cols = colnames(comut_filt_numeric)[colSums(comut_filt_numeric)==0]
    if(length(sum_zero_cols)>0){
      print(paste(comut_label, '-', 'Dropping following zero-occurrence events from comut:'))
      print(sum_zero_cols)
      comut_filt = comut_filt[, -sum_zero_cols, with=F]
    }
  }

  if(is.null(min_pc_force_retain_cols)) #for following code to work
      min_pc_force_retain_cols = c()

  if(!is.null(min_cna_pc)){
    nonzero_freq = sapply(comut_filt, function(x){sum(!(x==0 | x=='0' | is.na(x)))}) / nrow(comut_filt) * 100
    cna_cols = grepl(cna_regex, colnames(comut_filt))
    comut_filt = comut_filt[,(nonzero_freq >= min_cna_pc | !cna_cols | colnames(comut_filt)%in%min_pc_force_retain_cols),with=F]
  }

  if(!is.null(min_driver_pc)){
    #note: nonzero is used to retain participant_col too
    nonzero_freq = sapply(comut_filt, function(x){sum(!(x==0 | x=='0' | is.na(x)))}) / nrow(comut_filt) * 100
    cna_cols = grepl(cna_regex, colnames(comut_filt))
    comut_filt = comut_filt[, (nonzero_freq >= min_driver_pc | cna_cols | colnames(comut_filt)%in%min_pc_force_retain_cols) ,with=F]
  }

  if(to_binary){
    pat_backup = comut_filt[[pat_col]]
    comut_filt = comut_filt[,-c(pat_col),with=F]
    comut_filt[!is.na(comut_filt) & comut_filt!=0] = 1
    comut_filt = cbind(pat_backup, comut_filt)
    setnames(comut_filt, 'pat_backup', pat_col)
  }

  if(!is.null(outfile)){
    print(paste('Writing', comut_label, 'to:', outfile))
    fwrite(comut_filt, outfile, sep='\t')
  }

  if(!is.null(retval) &&  retval=='comut')
    return(comut_filt)
}


#Reformat events (used this initially to change chromosome events from e.g. amp(12) to amp_12, like the focal format I had)
reformat_events <- function(events, regex_from=c('\\(', '\\)$'), regex_to=c('_', ''),
                            sub_file=NULL, from_col=1, to_col=2){
  ### Reformat with list of regexes
  if(!is.null(regex_from)){
    for(i in 1:length(regex_from)){
      events = sub(regex_from[i], regex_to[i], events)
    }
  }

  ### Replace values by renaming file (default: sub names in 1st col with vals from 2nd col)
  if(!is.null(sub_file)){
    if(file.exists(sub_file)){
      subber = fread(sub_file)
      match_bool = !is.na(match(events, subber[[from_col]]))
      match_ind = match(events, subber[[from_col]])
      match_ind = match_ind[!is.na(match_ind)]

      events[match_bool] = subber[[to_col]][match_ind] #replace matches

    }
  }

  return(events)
}

### Get CNV event information from raw file, generate unique events in Ziao format and create "pretty" names
get_focal_event_annotation <- function(lesions_file, n_annot_cols = 9, blacklist=NULL,
                                      rename_from=c('Residual.q.values.after.removing.segments.shared.with.higher.peaks'),
                                      rename_to=c('residual.q'),
                                      drop_cols=c('Broad.or.Focal', 'Amplitude.Threshold'),
                                      verbose=T){
  #Get events
  lesions = fread(lesions_file) #lesions_file = args$gistic.focal.raw.file
  colnames(lesions) = gsub(' ', '.', colnames(lesions))

  #For sanity check
  if(verbose){
    cat(paste0('Annotation cols in lesion file:\n', paste(colnames(lesions)[1:n_annot_cols], collapse='\n') ))
    print(paste('First column in lesions file dropped:', colnames(lesions)[n_annot_cols + 1]))
  }

  #drop value rows and extra columns
  lesions = lesions[!grepl('CN values', lesions$Unique.Name), 1:n_annot_cols, with=F]

  #make unique as in Ziao's format
  #e.g. del_1p36_11 (add _1 for second event in same cytoband)
  #Necessary to split amp and del if same cytogenetic band has amp and del
  #rename amp
  amp_bool = grepl('^Amp', lesions$Unique.Name)
  event_id_amp = paste0('amp_', sub('\\.', '_', make.unique(as.character(lesions$Descriptor[amp_bool]), sep='_')))
  #rename del
  del_bool = grepl('^Del', lesions$Unique.Name)
  event_id_del = paste0('del_', sub('\\.', '_', make.unique(as.character(lesions$Descriptor[del_bool]), sep='_')))
  #merge
  lesions$event_id = c(event_id_amp, event_id_del)

  #Filter blacklist to preventa adding a/b suffixes in event_pretty where one is filtered out
  if(!is.null(blacklist)){
    lesions = lesions[!lesions$event_id %in% blacklist]
  }

  ## Create pretty names e.g. for plotting
  lesions$event_pretty = paste0(sub('(Del|Amp).*', '\\1', lesions$Unique.Name) , '(', lesions$Descriptor, ')')
  # Deal with non-unique (add a, b, c suffix)
  nonuniq = table(lesions$event_pretty)[table(lesions$event_pretty)>1]
  for(i in 1:length(nonuniq)){
    event = names(nonuniq)[i]
    lesions$event_pretty[lesions$event_pretty==event] = paste0(event, letters[1:nonuniq[i]])
  }

  if(!is.null(rename_from)){
    setnames(lesions, rename_from, rename_to)
  }

  if(!is.null(drop_cols)){
    lesions = lesions[,-drop_cols,with=F]
  }

  make_first = c('event_id', 'event_pretty')
  lesions = cbind(lesions[,make_first,with=F], lesions[,-make_first,with=F])

  #remove . from colnames
  colnames(lesions) = gsub('\\.', '_', colnames(lesions))

  return(lesions)
}


preprocess_counts_and_rates <- function(infile, log_base=10, outfile=NULL, retval=NULL){
  rates_dt = fread(infile)

  rates_dt[['rate_sil_frac']] = rates_dt$rate_sil / rates_dt$rate_tot
  rates_dt[['rate_non_frac']] = rates_dt$rate_non / rates_dt$rate_tot
  raise_power = max(ceiling(-log(rates_dt$rate_tot, log_base)))+1
  rate_raised_power_label = paste0('rate_tot_raised',log_base,'to', raise_power)
  rates_dt[[rate_raised_power_label]] = log(rates_dt$rate_tot * log_base^raise_power, base=log_base)
  rates_dt[[paste0(rate_raised_power_label,'_sil_frac')]] = rates_dt[[rate_raised_power_label]] * rates_dt[['rate_sil_frac']]
  rates_dt[[paste0(rate_raised_power_label,'_non_frac')]] = rates_dt[[rate_raised_power_label]] * rates_dt[['rate_non_frac']]

  if(!is.null(outfile)){
    fwrite(rates_dt, outfile, sep='\t')
  }

  if(!is.null(retval))
    return(rates_dt)
}



cna_binary_to_text_comut <- function(cm, del_label='CN-del', amp_label='CN-amp',
                                     del_regex='^[Dd]el|^[Ll]oss', amp_regex='^[Aa]mp|^[Gg]ain|^[Tt]ri'){
  cm[, names(cm) := lapply(.SD, as.character)]
  for(col in grep(del_regex, colnames(cm))){
    data.table::set(cm, which(cm[[col]]=='1'), col, del_label)
  }
  for(col in grep(amp_regex, colnames(cm))){
    data.table::set(cm, which(cm[[col]]=='1'), col, amp_label)
  }

  return(cm)
}

serialize_comut <- function(cm, txt_labs, serial_order=NULL, zero_base=F, cols_skip=c('participant_id'), inplace=F){
  if(!inplace)
    cm = copy(cm)

  if(is.null(serial_order)){ #default serials (no specific order)
    serial_order = 1:length(txt_labs)
    if(zero_base){
      serial_order = serial_order - 1
    }
  }

  cm[, names(cm) := lapply(.SD, as.character)]
  cols_for_replace = colnames(cm)[!colnames(cm) %in% cols_skip]

  for(col in cols_for_replace){
    for(i in 1:length(txt_labs)){
      data.table::set(cm, which(cm[[col]]==txt_labs[i]), col, serial_order[i])
    }
  }

  return(cm)
}


comut_to_freq <- function(cm_file=NULL, cm=NULL, discard_cols=c('MCLL', 'UCLL'), ret='pc_round', dig=3){
  if(!is.null(cm_file)){
    cm = fread(cm_file)
  }
  if(!is.null(discard_cols))
    cm = cm[,-discard_cols,with=F]

  if(ret=='pc_round' || ret=='pc'){
    cm = cm[,lapply(.SD, sum), .SDcols=sapply(cm, is.numeric)] / nrow(cm) * 100
  } else if (ret=='counts'){
    cm = cm[,lapply(.SD, sum), .SDcols=sapply(cm, is.numeric)]
  }

  cm = as.data.table(t(cm), keep.rownames = T)
  setnames(cm, c('event_id', 'freq'))

  if(ret=='pc_round')
    cm$freq = round(cm$freq, dig)

  return(cm)
}

comut_fisher_per_group <- function(comut=NULL, comut_file=NULL, group_col='ec_label', min_mut_events_all=10, min_mut_events_group=0, blacklist=NULL, group_subset=NULL){
  if(!is.null(comut_file)){
    comut = fread(comut_file)
  }

  if(!is.null(blacklist))
    comut = comut[, colnames(comut)[!colnames(comut) %in% blacklist], with=F]

  if(any(duplicated(comut$participant_id))){
    print('Dropping duplicates')

  }
  comut = comut[!duplicated(participant_id)][,-c('participant_id')]

  #Reformat for analyses
  cmlong = melt(comut, id.vars = c(group_col))
  setnames(cmlong, c(group_col,'event', 'is_mut'))
  cmlong$is_mut[cmlong$is_mut>0] = 1 #set all non-zero to 1

  cmlong_count = cmlong[,list(freq=.N), by=c(group_col, 'event', 'is_mut')]
  cmstat = dcast.data.table(cmlong_count, paste(group_col, '+ event ~ is_mut'), value.var = 'freq', fill=0)
  setnames(cmstat, c('0','1'), c('wt', 'mut'))

  cmstat[,'n_all'] = nrow(comut)
  mut_totals = cmlong_count[is_mut==1,][,.(mut_all=sum(freq)),by='event']

  cmstat = merge(cmstat, mut_totals, by='event')
  cmstat[,mut_pc_of_all_mut := round(mut/mut_all*100,2),] #fraction of mutated in this EC of all mutated
  cmstat[,`:=` ('wt_other' = n_all - wt - mut_all, 'mut_other' = mut_all - mut)]
  cmstat[,mut_pc_in_group := round((mut/(wt+mut))*100,2),] #fraction in this EC that are mutated

  cmstat[, fisher_pval := mapply(function(a,a2,b,b2){runChisq(a,a2,b,b2,increment_odds = 0.5,fisher = T, retVal=2)[1]}, mut,wt,mut_other,wt_other)]
  cmstat[, odds_ratio := mapply(function(a,a2,b,b2){runChisq(a,a2,b,b2,increment_odds = 0.5,fisher = T, retVal=2)[2]}, mut,wt,mut_other,wt_other)]

  if(!is.null(group_subset)){
    cmstat = cmstat[cmstat[[group_col]] %in% group_subset]
  }

  cmstat = cmstat[mut_all >= min_mut_events_all & mut >= min_mut_events_group]

  cmstat[,'fisher_minusLogP'] = -log10(cmstat$fisher_pval)
  cmstat[,'fisher_minusLogP_signed'] = cmstat[,'fisher_minusLogP'] * ifelse(cmstat$odds>=1, 1, -1)
  cmstat[,'fisher_padj_BY'] = p.adjust(cmstat[['fisher_pval']], method='BY')
  cmstat[,'fisher_padj_BH'] = p.adjust(cmstat[['fisher_pval']], method='BH')

  cmstat_w_signp = as.data.frame(dcast(cmstat[cmstat$mut_all>=min_mut_events_all], paste('event ~', group_col), value.var='fisher_minusLogP_signed'))
  row.names(cmstat_w_signp) = cmstat_w_signp$event
  cmstat_w_freq = as.data.frame(dcast(cmstat[cmstat$mut_all>=min_mut_events_all], paste('event ~', group_col), value.var='mut'))
  row.names(cmstat_w_freq) = cmstat_w_signp$event

  return(list(cmstat=cmstat, cmstat_w_signp=cmstat_w_signp, cmstat_w_freq=cmstat_w_freq))
}
