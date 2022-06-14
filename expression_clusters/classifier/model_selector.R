### Model selection for expression cluster classifier
### Author: Binyamin Knisbacher
### Provided for review and educational purposes only.

<<<<<<< HEAD
source('~/github/Rscripts/Sources/CLL1085/CLL1100_current.src.R', proj<-new.env())
=======
source('~/github/Rscripts/Sources/CLLmap_project.src.R', proj<-new.env())
>>>>>>> 5892dc2afe20ea72e74026403ee0dddf16754e3a
source('~/github/Rscripts/Sources/graphics.src.R')

suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(pROC))
suppressPackageStartupMessages(library(PRROC))

EC_PRED_RUN = 'bc_tpm_genes' #'tpm_genes'
ec_preds = proj$get_ec_preds_args(opt=EC_PRED_RUN)
outdir = ec_preds$pred_outdir

OPEN_TEST_SET = T
DEFAULT_GS_SIZES = c(5, 10, 20, 50)
PER_EC_TH_METHOD = "pr_diff_vs_best_other" #'pr_within_ec'

if(is.null(ec_preds$geneset_subdir_suffix)){
  ec_preds$geneset_subdir_suffix = 'with_genesets'
}

run_all_models_per_test_set = ifelse(EC_PRED_RUN %in% c('tpm_genes', 'bc_tpm_genes'), T, F)


### Predictor functions
get_f1 <- function(pr, rec, digits=2){
  f1 = 2*pr*rec / (pr+rec)
  if(!is.null(digits)){
    f1 = round(f1, digits)
  }
  return(f1)
}

get_pred_metrics <- function(real, pred, no_per_class=F){
  #some code from: https://blog.revolutionanalytics.com/2016/03/com_class_eval_metrics_r.html

  #Create matrix, but in case some level isn't predicted, need to pad (adding class_levels) with all classes and then reomve the padding (subtract diag of 1)
  class_levels = sort(unique(real))
  cm = as.matrix(table(Actual = c(real, class_levels), Predicted = c(pred, class_levels) )) # create the confusion matrix
  cm = cm - diag(length(class_levels))

  n = sum(cm) # number of instances
  nc = nrow(cm) # number of classes
  diag = diag(cm) # number of correctly classified instances per class
  rowsums = apply(cm, 1, sum) # number of instances per class
  colsums = apply(cm, 2, sum) # number of predictions per class
  p = rowsums / n # distribution of instances over the actual classes
  q = colsums / n # distribution of instances over the predicted classes

  accuracy = round(sum(diag) / n, 3)

  precision = diag / colsums
  precision[is.na(precision)] = 0
  recall = diag / rowsums
  recall[is.na(recall)] = 0
  f1 = (1 + 1^2) * precision * recall / (1^2 * precision + recall)
  f1[is.na(f1)] = 0
  f2 = (1 + 2^2) * precision * recall / (2^2 * precision + recall)
  f2[is.na(f2)] = 0
  f05 = (1 + 0.5^2) * precision * recall / (0.5^2 * precision + recall)
  f05[is.na(f05)] = 0

  per_class_df = round(data.frame(precision, recall, f1, f2, f05), 3)
  per_class_df[is.na(per_class_df)] = 0

  macroPrecision = mean(precision)
  macroRecall = mean(recall)
  macroF1 = mean(f1)
  macroF2 = mean(f2)
  macroF05 = mean(f05)

  macros = round(c(macroPrecision, macroRecall, macroF1, macroF2, macroF05), 3)
  names(macros) = c('macroPrecision', 'macroRecall', 'macroF1', 'macroF2', 'macroF05')

  if(no_per_class){
    return( c(accuracy=accuracy, macros))
  } else {
    return(list(accuracy=accuracy, per_class=per_class_df, macro=macros))
  }
}

get_best_f1_from_roc <- function(roc_obj, select_method='last', ret='best'){
  pr_th_dt = as.data.table(coords(roc_obj, ret=c('threshold', 'recall', 'precision')))
  pr_th_dt$f1 = 2*pr_th_dt$precision * pr_th_dt$recall / (pr_th_dt$precision + pr_th_dt$recall)
  best_pr_th_dt = pr_th_dt[f1==max(pr_th_dt$f1, na.rm = T)]

  #selects last among multiple results with max F1
  if(select_method=='first'){
    best_pr_th_dt = best_pr_th_dt[1,]
  } else if(select_method=='last') {
    best_pr_th_dt = best_pr_th_dt[nrow(best_pr_th_dt),]
  }

  if(ret=='best')
    return(best_pr_th_dt)
}

get_pr_curve_and_roc_metrics <- function(corrects, incorrects, th=NULL){
  #Note: corrects/incorrects are the metric values to which threshold (th) is applied

  #Get basic results for series
  roc.res = roc(case=corrects, control=incorrects, direction="<")

  # Precision-recall AUC
  pr_curve_res <- pr.curve(scores.class0 = corrects, scores.class1 = incorrects, curve = T)
  pr_curve_auc = pr_curve_res$auc.integral

  #ROC AUC and youden-based metrics
  roc_youden = coords(roc.res, "best", best.method = 'youden')
  roc_auc_ec = roc.res$auc

  # Precision-recall max F1
  proc_max_f1 = round(get_best_f1_from_roc(roc.res)$f1,2)
  # Precision, Recall and F1 per threshold
  if(!is.null(th)){
    proc_res = coords(roc.res, th, ret=c('threshold', 'recall', 'precision'))
    proc_res$f1 = get_f1(proc_res$precision, proc_res$recall)
    proc_res$precision = round(proc_res$precision,2)
    proc_res$recall = round(proc_res$recall,2)

    pr_and_roc_mets = as.list(c(round(c(prcurve_auc=pr_curve_auc, prcurve_max_f1=proc_max_f1, prcurve_precomputed_th=th, prcurve_th_precision=proc_res$precision[1], prcurve_th_recall=proc_res$recall[1], prcurve_th_f1=proc_res$f1[1]),2),
                                round(c(roc_auc=roc_auc_ec, roc_youden_th=roc_youden$threshold, roc_youden_sensitivity=roc_youden$sensitivity, roc_youden_specificity=roc_youden$specificity),2)))
  } else {
    pr_and_roc_mets = as.list(c(round(c(prcurve_auc=pr_curve_auc, prcurve_max_f1=proc_max_f1),2),
                                round(c(roc_auc=roc_auc_ec, roc_youden_th=roc_youden$threshold, roc_youden_sensitivity=roc_youden$sensitivity, roc_youden_specificity=roc_youden$specificity),2)))
  }

  return(pr_and_roc_mets)
}

#input dt must have a column ec_label and continuous value columns per EC
get_ec_diffs_best_other <- function(mydt, ec_list=proj$EC_LIST){
  for (ec in ec_list){ #ec='EC1'
    get_best_f = function(x){x[names(x)==ec] - max(x[names(x)!=ec])}
    mydt[[paste0(ec, '_diff_best_other')]] = apply(mydt[,ec_list,with=F], 1, get_best_f)
  }

  ## Get one column for the ec_label vs other best (for all-EC stats)
  get_diff_best_other_col <- function(x){return(x[paste0(x['ec_label'], '_diff_best_other')])}
  mydt[['this_diff_best_other']] = apply(mydt[,c('ec_label', paste0(paste0(ec_list, '_diff_best_other'))), with=F], 1, get_diff_best_other_col)

  return(mydt)
}

#add_roc_and_youden: don't use this option (it over-estimates in multi-class)
get_binarized_pr_and_roc <- function(mydt, ec_thresholds, no_per_ec=F, no_ec_means=F, ec_list=proj$EC_LIST, add_roc_and_youden=F, all_ret_fmt='list'){
  #AUC and 'optimal' prediction per EC using PR curve thresholds
  per_ec_dt = data.table()
  for (ec in ec_list){
    roc.one_ec = roc(case=mydt[ec_label==ec][[paste0(ec, '_diff_best_other')]], control=mydt[ec_label!=ec][[paste0(ec, '_diff_best_other')]], direction = "<")

    proc_max_f1_ec = round(get_best_f1_from_roc(roc.one_ec)$f1, 2)
    proc_res_one_ec = coords(roc.one_ec, ec_thresholds[[ec]], ret=c('threshold', 'recall', 'precision'))
    proc_res_one_ec$f1 = get_f1(proc_res_one_ec$precision, proc_res_one_ec$recall)
    proc_res_one_ec$precision = round(proc_res_one_ec$precision,2)
    proc_res_one_ec$recall = round(proc_res_one_ec$recall,2)

    pr.curve.ec <- pr.curve(scores.class0 = mydt[ec_label==ec][[paste0(ec, '_diff_best_other')]],
                            scores.class1 = mydt[ec_label!=ec][[paste0(ec, '_diff_best_other')]], curve = T)
    pr.curve.ec.auc = pr.curve.ec$auc.integral

    per_ec_dt0 = as.list(c(ec_label=ec,
                           round(c(prcurve_auc=pr.curve.ec.auc, prcurve_max_f1=proc_max_f1_ec,
                                   prcurve_th=ec_thresholds[[ec]], prcurve_th_precision=proc_res_one_ec$precision,
                                   prcurve_th_recall=proc_res_one_ec$recall, prcurve_th_f1=proc_res_one_ec$f1),2)))

    #ROC curve metrics
    if(add_roc_and_youden){
      youden_ec = coords(roc.one_ec, "best", best.method = 'youden')
      roc_auc_ec = roc.one_ec$auc
      youden_ec_metric_list = as.list(round(c(roc_auc=roc_auc_ec, roc_youden_th=youden_ec$threshold, roc_youden_sensitivity=youden_ec$sensitivity, roc_youden_specificity=youden_ec$specificity),2))
      per_ec_dt0 = append(per_ec_dt0, youden_ec_metric_list)
    }

    per_ec_dt = rbindlist(list(per_ec_dt, per_ec_dt0))
  }

  #Macro-Metrics per full dataset (mean of per-EC)
  per_ec_dt = merge(per_ec_dt, mydt[,.N, by='ec_label'])[order(ec_label)]
  mean_col_drop = c('ec_label', 'N', 'prcurve_th', 'roc_youden_th')
  mean_col_drop = mean_col_drop[mean_col_drop %in% colnames(per_ec_dt)]
  mets_ec_mean = sapply(per_ec_dt[,-mean_col_drop,with=F], function(x){mean(as.numeric(x))})
  my_weights = per_ec_dt$N
  mets_ec_mean_w = sapply(per_ec_dt[,-mean_col_drop,with=F], function(x){weighted.mean(as.numeric(x), as.numeric(my_weights))})

  if(all_ret_fmt=='dt'){
    ec_means = data.table(metric=names(mets_ec_mean), ec_mean=round(mets_ec_mean, 2), weighted_ec_mean=round(mets_ec_mean_w, 2))
  } else { #as list
    names(mets_ec_mean_w) = paste0('ec_w_mean_', names(mets_ec_mean_w))
    names(mets_ec_mean) = paste0('ec_mean_', names(mets_ec_mean))
    ec_means = as.list(c(round(mets_ec_mean, 2), round(mets_ec_mean_w, 2)))
  }

  if(no_ec_means){
    return(per_ec_dt)
  } else if(no_per_ec){
    return(ec_means)
  } else {
    return(list(ec_means=all_mean_ret, per_ec_dt=per_ec_dt))
  }
}

get_pr_curves_per_ec <- function(mydt, ec_list=proj$EC_LIST, col_suffix='_diff_best_other', direction = "<", plotfile=NULL, prc_palette=rev(heat.colors(100)), ret='list'){
  pr_curves_per_ec = list()
  for (ec in ec_list){
    roc.one_ec = roc(case=mydt[ec_label==ec][[paste0(ec, col_suffix)]], control=mydt[ec_label!=ec][[paste0(ec, col_suffix)]], direction = direction)

    pr.curve.ec <- pr.curve(scores.class0 = mydt[ec_label==ec][[paste0(ec, '_diff_best_other')]],
                            scores.class1 = mydt[ec_label!=ec][[paste0(ec, '_diff_best_other')]], curve = T)

    pr_curves_per_ec[[ec]] = pr.curve.ec
  }

  if(!is.null(plotfile)){
    pdf(plotfile)
    for(prc in pr_curves_per_ec){
      plot(prc, scale.color=prc_palette)
    }
    dev.off()
  }

  if(!is.null(ret)){
    return(pr_curves_per_ec)
  }
}


### SCREEN ALL MODELS TO IDENTIFY BEST ####
if(is.null(ec_preds$gs_sizes)){
  gs_sizes = DEFAULT_GS_SIZES
} else {
  gs_sizes = ec_preds$gs_sizes
}


rf_opts = c('rf', 'rf2')
scaling_opts = c('', '_scaled')
fold_inds = 0:4
oversamp_opts = c('', 'ADASYN', 'BorderlineSMOTE', 'SMOTE', 'SVMSMOTE')

modres = data.table() #model results
for (gs_size in gs_sizes){#gs_size=gs_sizes[1]
  for (rf_opt in rf_opts){ #rf_opt=rf_opts[1]
    for (scaling_opt in scaling_opts){ #scaling_opt='_scaled'
      for (fold_ind in fold_inds){ #fold_ind=0
        for (oversamp in oversamp_opts){ #oversamp='ADASYN'
          is_scaled = (scaling_opt=='_scaled')
          file_re = paste0(ec_preds$pred_rootdir, '/results/', 'topBot_', gs_size,'_a_', ec_preds$geneset_subdir_suffix, '/soft_clust/', rf_opt, '_', oversamp, '_split_', fold_ind, scaling_opt, '_', oversamp,ifelse(oversamp!='', '_', ''),  '[0-9.]*.csv')
          f = Sys.glob(file_re)
          if(length(f)==0){
            print(paste('No result for:', file_re))
          } else {
            # print(paste('Found file:', f))
            mod_dt = proj$parse_ec_pred_files(f, parse_opt = 'train_test')
            modres0 = get_pred_metrics(mod_dt[test_sample==1]$ec_real, mod_dt[test_sample==1]$ec_pred, no_per_class = T)
            modres0 = as.list(c(gs_size=gs_size, rf_opt=rf_opt, oversamp=oversamp, is_scaled=is_scaled, fold=fold_ind+1, modres0, file=f))

            modres = rbindlist(list(modres, modres0))
          }
        }
      }
    }
  }
}

## convert to numeric
models_dt = data.table(modres)
metric_cols = c('accuracy', 'macroPrecision', 'macroRecall', 'macroF1', 'macroF2', 'macroF05')
for (col in metric_cols){
  models_dt[[col]] = as.numeric(models_dt[[col]])
}

## Identify mean scores and select the best
mod_means = models_dt[, lapply(.SD, mean), by=list(gs_size, rf_opt, oversamp, is_scaled), .SDcols=metric_cols]
mod_means$hm_accuracy_macroF1 = (2 * mod_means$accuracy * mod_means$macroF1) / (mod_means$accuracy + mod_means$macroF1)
setorder(mod_means, -hm_accuracy_macroF1)
mod_means$best = c(1, rep(0, nrow(mod_means)-1))

model_mean_metrics_outfile = paste0(ec_preds$pred_outdir, '/', 'model_metrics_mean_of_fold.',ec_preds$output_timestamp,'.tsv')
fwrite(mod_means, model_mean_metrics_outfile, sep='\t')

### show correlation of major selection metrics
model_metrics_scatter_pdf = paste0(outdir, '/', 'model_metrics_scatter',ec_preds$output_timestamp,'.pdf')
p_mod_scatter = ggscatter(mod_means[is_scaled==FALSE], x='accuracy', y='macroF1',
                          add='reg.line', conf.int = TRUE, cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))
ggsave(model_metrics_scatter_pdf, p_mod_scatter, height = 4, width = 4)

best_folds = merge(models_dt, mod_means[,c('gs_size', 'rf_opt', 'oversamp', 'is_scaled', 'best')])[best==1]
best_fold = best_folds[order(-accuracy)][1,]

best_per_fold_outfile = paste0(ec_preds$pred_outdir, '/', 'best_model_metrics_per_fold.',ec_preds$output_timestamp,'.tsv')
fwrite(best_folds, best_per_fold_outfile, sep='\t')

mod_best_dt = proj$parse_ec_pred_files(best_fold$file, parse_opt = 'train_test')
modres_best = get_pred_metrics(mod_best_dt[test_sample==1]$ec_real, mod_best_dt[test_sample==1]$ec_pred, no_per_class = T)

### Get folds per best and show boxplot of correct vs. incorrect scores
pred_dt = proj$parse_ec_pred_files(file_list=best_folds$file, parse_opt='train_test')
pred_dt$correct = as.integer(pred_dt$ec_real==pred_dt$ec_pred)
pred_dt$ec_label = paste0('EC', pred_dt$ec_real)

correct_vs_incorrect_boxplot_pdf = paste0(outdir, '/', 'correct_vs_incorrect_boxplot_per_fold.',ec_preds$output_timestamp,'.pdf')
pdf(correct_vs_incorrect_boxplot_pdf, useDingbats = F, height=4, width=5)
p_box_correct = ggboxplot(pred_dt[test_sample==1], x='correct', y='diff_top2', facet.by='cv_fold', nrow=1) +
  geom_jitter(aes(color=ec_label), size=1) +
  scale_color_manual(values=proj$get_colors('ecs'), name='Expression cluster') + labs(x='Correct', y='Prediction margin (Best - 2nd)') +
  stat_compare_means(label='p.format', size=3)
p_box_correct
dev.off()

### Get per-class stats (for best model)
per_class_dt = data.table()
for (i in fold_inds+1){ #i=1
  f = best_folds$file[i]
  mod_folds_dt = proj$parse_ec_pred_files(f, parse_opt = 'train_test')
  per_class_dt0 = get_pred_metrics(mod_folds_dt[test_sample==1]$ec_real, mod_folds_dt[test_sample==1]$ec_pred, no_per_class = F)$per_class
  per_class_dt0$ec_label = paste0('EC', row.names(per_class_dt0))
  per_class_dt0$cv_fold = paste0('fold', i)

  #these should be the same for all
  per_class_dt0$gs_size = best_folds$gs_size[i]
  per_class_dt0$rf_opt = best_folds$rf_opt[i]
  per_class_dt0$oversamp = best_folds$oversamp[i]
  per_class_dt0$is_scaled = best_folds$is_scaled[i]
  per_class_dt0$file = best_folds$file[i]

  per_class_dt0 = as.data.table(per_class_dt0)

  per_class_dt = rbindlist(list(per_class_dt, per_class_dt0))
}

per_class_metrics_outfile = paste0(ec_preds$pred_outdir, '/', 'per_class_metrics.',ec_preds$output_timestamp,'.tsv')
fwrite(per_class_dt, per_class_metrics_outfile, sep='\t')

### Get best thresholds per all and per-EC based on folds
folds_test_dt = copy(pred_dt[test_sample==1]) #only test samples for all folds
folds_test_dt$correct = ifelse(folds_test_dt$correct==1, 'yes', 'no')

fold_test_thresholds = list()
# Get best threshold per EC
if(PER_EC_TH_METHOD=='pr_within_ec'){
  for (ec in proj$EC_LIST){
    folds.test.roc.per_ec = roc(case=folds_test_dt[folds_test_dt$correct=='yes' & ec_label==ec]$diff_top2, control=folds_test_dt[folds_test_dt$correct=='no' & ec_label==ec]$diff_top2, direction = "<")
    folds.test.per_ec.best_f1 = get_best_f1_from_roc(folds.test.roc.per_ec, ret='best')
    print(ec)
    print(folds.test.per_ec.best_f1)
    fold_test_thresholds[[ec]] = folds.test.per_ec.best_f1$threshold
  }
} else { #based on per-EC scores, comparing the EC to all others (PER_EC_TH_METHOD=="pr_diff_vs_best_other")
  folds_test_dt = get_ec_diffs_best_other(folds_test_dt)
  for (ec in proj$EC_LIST){
    folds.test.roc.per_ec = roc(case=folds_test_dt[ec_label==ec][[paste0(ec, '_diff_best_other')]],
                                control=folds_test_dt[ec_label!=ec][[paste0(ec, '_diff_best_other')]], direction = "<")
    folds.test.per_ec.best_f1 = get_best_f1_from_roc(folds.test.roc.per_ec, ret='best')
    print(ec)
    print(folds.test.per_ec.best_f1)
    fold_test_thresholds[[ec]] = folds.test.per_ec.best_f1$threshold
  }
}

#### Analysis showing performance across range of N markers
#Get best model per gs_size
best_per_gs_size_dt = mod_means[order(gs_size, -hm_accuracy_macroF1)][!duplicated(gs_size)]
best_per_gs_size_dt$gs_size = as.integer(best_per_gs_size_dt$gs_size)
setorder(best_per_gs_size_dt, gs_size)


### Get performance on test set - code to be used only when finalized
if(OPEN_TEST_SET){
  test_dt = proj$get_ec_preds(opt=EC_PRED_RUN, filter_opt='master_tumors_mainCohorts', ret_mode='fold_sum_max_dt')
  test_dt$ec_pred_numeric = as.numeric(sub('EC', '', test_dt$ec_pred))
  ec_test_dt = test_dt[sample_id %in% proj$get_sample_set('ec_test')]
  ec_test_dt = merge(ec_test_dt, proj$get_ecs()[,c('sample_id', 'expression_cluster')], by='sample_id')
  ec_test_dt$ec_label = paste0('EC', ec_test_dt$expression_cluster)
  test_mets = get_pred_metrics(ec_test_dt$expression_cluster, ec_test_dt$ec_pred_numeric)
  test_mets_per_ec_dt = data.table(test_mets$per_class)
  test_mets_per_ec_dt$ec_label = proj$EC_LIST

  #create confusion matrix
  ec_confusion = matrix(table(ec_test_dt$ec_label, ec_test_dt$ec_pred), nrow=length(unique(ec_test_dt$ec_label)))
  colnames(ec_confusion) = proj$EC_LIST
  row.names(ec_confusion) = proj$EC_LIST

  #compute dominance score
  ec_confusion2 = ec_confusion
  for (i in 1:nrow(ec_confusion2)){
    for(j in 1:ncol(ec_confusion2)){
      ec_confusion2[i,j] = ec_confusion[i,j] / (sum(ec_confusion[i,]) + sum(ec_confusion[,j]) - ec_confusion[i,j])
    }
  }

  p_confusion = Heatmap(ec_confusion2, cluster_rows=F, cluster_columns=F,
                        row_title='Real', column_title='Predicted', col=c('#FFFF66', '#FF8B01', '#A40001'),
                        heatmap_legend_param = list(title = "Dominance"))
  test_set_confusion_mat_pdf = paste0(outdir, '/', 'test_set_confusion_matrix_heatmap.',ec_preds$output_timestamp,'.pdf')
  pdf(test_set_confusion_mat_pdf)
  draw(p_confusion)
  dev.off()

  test_all_dt = proj$get_ec_preds(opt=EC_PRED_RUN, filter_opt='master_tumors_mainCohorts', ret_mode='fold_sum_dt')
  test_all_dt = dcast.data.table(test_all_dt[,-c('participant_id'),with=F], 'sample_id ~ ec_pred', value.var = 'ec_pred_p')
  test_all_dt$diff_top2 =  apply(test_all_dt[,-c('sample_id'),with=F], 1, function(x){return(sort(x, decreasing=T)[1] - sort(x, decreasing=T)[2])})
  test_all_dt = merge(test_all_dt, ec_test_dt)
  test_all_dt$correct = factor(ifelse(test_all_dt$ec_pred_numeric==test_all_dt$expression_cluster, 'yes', 'no'), levels=c('yes', 'no'), ordered=T)
  test_all_dt = get_ec_diffs_best_other(test_all_dt)


  ### Diff2-based plots
  # 1. Correct vs. Incorrect diff2 boxplots
  test_set_correct_vs_incorrect_boxplot_pdf = paste0(outdir, '/', 'correct_vs_incorrect_boxplot_test_set.',ec_preds$output_timestamp,'.pdf')
  p_box_correct_test_set = ggboxplot(test_all_dt, x='correct', y='diff_top2', nrow=1) +
    geom_jitter(aes(color=ec_label), size=1) +
    scale_color_manual(values=proj$get_colors('ecs'), name='Expression cluster') + labs(x='Correct', y='Prediction margin (Best - 2nd)') +
    stat_compare_means(label='p.format', label.x.npc='center') #+ geom_hline(yintercept = proc_res_ft$threshold, linetype=2)
  ggsave(test_set_correct_vs_incorrect_boxplot_pdf, p_box_correct_test_set, height=4, width=3)

  # 2. Correct vs. Incorrect ROC curve
  test_roc.c_vs_ic = roc(case=test_all_dt[test_all_dt$correct=='yes']$diff_top2, control=test_all_dt[test_all_dt$correct=='no']$diff_top2, direction = "<")
  test_roc.c_vs_ic.auc = round(as.numeric(sub('.* ', '', test_roc.c_vs_ic$auc)), 2)
  roc_testset_pdf = paste0(outdir, '/', 'roc_curve_test_set.correct_vs_incorrect.',ec_preds$output_timestamp,'.pdf')
  pdf(roc_testset_pdf, height=4, width=4)
  plot(test_roc.c_vs_ic, main=paste0('AUC=', test_roc.c_vs_ic.auc)) #print.thres="best", print.thres.best.method="youden",
  dev.off()

  ### PR curve per EC
  #AUC and 'optimal' prediction per EC using PR curve thresholds
  selected_model_per_ec_dt = data.table()
  roc.ts_ec.list = list() #for ROC AUC plots
  for (ec in proj$EC_LIST){
    roc.ts_ec = roc(case=test_all_dt[ec_label==ec][[paste0(ec, '_diff_best_other')]], control=test_all_dt[ec_label!=ec][[paste0(ec, '_diff_best_other')]], direction = "<")
    roc.ts_ec.list[[ec]] = roc.ts_ec

    proc_max_f1_ec = round(get_best_f1_from_roc(roc.ts_ec)$f1,2)
    proc_res_ft_ec = coords(roc.ts_ec, fold_test_thresholds[[ec]], ret=c('threshold', 'recall', 'precision'))
    proc_res_ft_ec$f1 = get_f1(proc_res_ft_ec$precision, proc_res_ft_ec$recall)
    proc_res_ft_ec$precision = round(proc_res_ft_ec$precision,2)
    proc_res_ft_ec$recall = round(proc_res_ft_ec$recall,2)

    test.pr.curve.ec <- pr.curve(scores.class0 = test_all_dt[ec_label==ec][[paste0(ec, '_diff_best_other')]],
                              scores.class1 = test_all_dt[ec_label!=ec][[paste0(ec, '_diff_best_other')]],
                              curve = T)
    test.pr.curve.auc.ec = test.pr.curve.ec$auc.integral

    #ROC curve metrics
    youden_ec_test = coords(roc.ts_ec, "best", best.method = 'youden')

    selected_model_per_ec_dt0 = as.list(c(ec_label=ec,
                                          round(c(prcurve_auc=test.pr.curve.auc.ec, prcurve_max_f1=proc_max_f1_ec,
                                                  prcurve_cv_th=fold_test_thresholds[[ec]], prcurve_cv_th_precision=proc_res_ft_ec$precision,
                                                  prcurve_cv_th_recall=proc_res_ft_ec$recall, prcurve_cv_th_f1=proc_res_ft_ec$f1),2),
                                          round(c(roc_auc=roc.ts_ec$auc, roc_youden_th=youden_ec_test$threshold, roc_youden_sensitivity=youden_ec_test$sensitivity, roc_youden_specificity=youden_ec_test$specificity),2)))
    selected_model_per_ec_dt = rbindlist(list(selected_model_per_ec_dt, selected_model_per_ec_dt0))
  }

  ### Collect/merge and write per-EC metrics
  #get N, N_correct, N_incorrect per best model per test set
  n_per_ec_test = merge(test_all_dt[,.N,by='ec_label'], dcast.data.table(test_all_dt[,.N,by=c('ec_label', 'correct')], 'ec_label ~ correct', value.var = 'N'), by='ec_label')
  setnames(n_per_ec_test, c('N', 'yes', 'no'), c('N', 'correct', 'incorrect'))

  selected_model_per_ec_merged = merge(merge(n_per_ec_test, test_mets_per_ec_dt), selected_model_per_ec_dt)
  selected_model_per_ec_file = paste0(outdir, '/', 'performance_best_model_on_test_set_per_ec.',ec_preds$output_timestamp,'.tsv')
  fwrite(selected_model_per_ec_merged, selected_model_per_ec_file, sep='\t')

  #Macro-Metrics per full dataset (mean of per-EC)
  selected_model_per_ec_dt = merge(selected_model_per_ec_dt, test_all_dt[,.N, by='ec_label'])[order(ec_label)]
  my_weights = selected_model_per_ec_dt$N
  metrics_all_weighted_mean_of_per_ec = sapply(selected_model_per_ec_dt[,-c('ec_label'),with=F], function(x){weighted.mean(as.numeric(x), as.numeric(my_weights))})
  metrics_all_mean_of_per_ec = sapply(selected_model_per_ec_dt[,-c('ec_label'),with=F], function(x){mean(as.numeric(x))})
  metrics_all_mean_of_per_ec_dt = data.table(metric=names(metrics_all_mean_of_per_ec), ec_mean=round(metrics_all_mean_of_per_ec, 2), weighted_ec_mean=round(metrics_all_weighted_mean_of_per_ec, 2))
  selected_model_mean_of_per_ec_file = paste0(outdir, '/', 'performance_best_model_on_test_set_mean_of_per_ec.',ec_preds$output_timestamp,'.tsv')
  fwrite(metrics_all_mean_of_per_ec_dt, selected_model_mean_of_per_ec_file, sep='\t')

  #Plot pr curves per EC
  pr_curve_per_ec_pdf = paste0(outdir, '/', 'precision_recall_curves_per_ec.',ec_preds$output_timestamp,'.pdf')
  get_pr_curves_per_ec(test_all_dt, plotfile=pr_curve_per_ec_pdf, ret=NULL)

  #Plot correct vs incorrect ROC curves per EC
  roc_testset_ec_pdf = paste0(outdir, '/', 'roc_curve_test_set.correct_vs_incorrect.per_ec.',ec_preds$output_timestamp,'.pdf')
  pdf(roc_testset_ec_pdf, height=4, width=4)
  for(ec in proj$EC_LIST){
    plot(roc.ts_ec.list[[ec]], main=paste0('AUC=', roc.ts_ec.list[[ec]]$auc))
  }
  dev.off()


  ####### Test each model on the test set
  if(run_all_models_per_test_set){
    tt_fold_sum_res = data.table() #model results
    for (gs_size in gs_sizes){
      for (rf_opt in rf_opts){
        for (scaling_opt in scaling_opts){ #scaling_opt='_scaled'; scaling_opt=""
          for (oversamp in oversamp_opts){
            is_scaled = (scaling_opt=='_scaled')
            #Note: this regex won't work for >=10-fold cross validation
            file_re = paste0(ec_preds$pred_rootdir, '/results/', 'topBot_', gs_size,'_a_', ec_preds$geneset_subdir_suffix, '/', rf_opt, '_', oversamp, '_split_', '[0-9]', scaling_opt, "_", oversamp, ifelse(oversamp!='', '_', ''), 'full_predict_fold_[0-9].csv')
            f = Sys.glob(file_re)

            if(length(f)==0){
              print(paste('No result for:', file_re))
            } else {
              tt_dt = proj$parse_ec_pred_files(files_regex = f, parse_opt='full')[,-c('ec_pred'),with=F]
              tt_dt = tt_dt[sample_id %in% proj$get_sample_set('ec_test')]
              tt_dt = tt_dt[, lapply(.SD, function(x){sum(x)/length(unique(tt_dt$cv_fold))}), by=c('sample_id'), .SDcols=proj$EC_LIST]
              tt_dt$diff_top2 = apply(tt_dt[,proj$EC_LIST,with=F], 1, function(x){return(sort(x, decreasing=T)[1] - sort(x, decreasing=T)[2])})
              tt_dt$ec_pred_p = apply(tt_dt[,proj$EC_LIST,with=F], 1, max)
              tt_dt$ec_pred = apply(tt_dt[,proj$EC_LIST,with=F], 1, which.max)
              tt_dt = merge(tt_dt, proj$get_ecs()[,c('sample_id', 'expression_cluster')])
              setnames(tt_dt, 'expression_cluster', 'ec_real')
              setcolorder(tt_dt, c('sample_id', proj$EC_LIST, 'diff_top2','ec_pred', 'ec_real', 'ec_pred_p'))

              tt_dt$ec_label = paste0('EC', tt_dt$ec_real)
              tt_dt = get_ec_diffs_best_other(tt_dt)

              tt_binarized_mets = get_binarized_pr_and_roc(tt_dt, fold_test_thresholds, no_per_ec = T, add_roc_and_youden = F)

              tt_fold_sum_res0 = get_pred_metrics(tt_dt$ec_real, tt_dt$ec_pred, no_per_class = T)
              tt_fold_sum_res0 = as.list(c(gs_size=gs_size, rf_opt=rf_opt, oversamp=oversamp, is_scaled=is_scaled, tt_fold_sum_res0, tt_binarized_mets))

              tt_fold_sum_res = rbindlist(list(tt_fold_sum_res, tt_fold_sum_res0), fill=T)
            }
          }
        }
      }
    }

    tt_fold_sum_res$gs_size = as.integer(tt_fold_sum_res$gs_size)

    #Results on test set per all models
    metrics_per_model_on_test_set_file = paste0(outdir, '/', 'performance_per_model_fold_sum_on_test_set.',ec_preds$output_timestamp,'.tsv')
    fwrite(tt_fold_sum_res, metrics_per_model_on_test_set_file, sep='\t')

    #Get stats (including AUC on test) per model based on best per gs_size on TRAINING set
    model_cols = c('gs_size', 'rf_opt', 'oversamp', 'is_scaled')
    mod_means_best_per_size = mod_means[order(gs_size, -hm_accuracy_macroF1)][!duplicated(gs_size)][,model_cols,with=F]
    tt_fold_sum_res$gs_size = as.character(tt_fold_sum_res$gs_size)
    tt_fold_sum_res$is_scaled = as.character(tt_fold_sum_res$is_scaled)
    tt_fold_sum_res_best = merge(tt_fold_sum_res, mod_means_best_per_size, by=model_cols)

    metrics_per_model_on_test_set_bestModelPerSize_file = paste0(outdir, '/', 'performance_per_model_fold_sum_on_test_set.best_model_per_size_by_training.',ec_preds$output_timestamp,'.tsv')
    fwrite(tt_fold_sum_res_best, metrics_per_model_on_test_set_bestModelPerSize_file, sep='\t')

    METRICS_FOR_PER_SIZE_LINEPLOT = c('accuracy', 'macroF1', 'ec_w_mean_prcurve_auc', 'ec_w_mean_prcurve_th_f1')
    METRICS_FOR_PER_SIZE_LINEPLOT_NAMES = c("Accuracy", "F1 (EC mean)", "pr-curve AUC\n(EC weighted average)", "pr-curve F1\n(EC weighted average)\n(at threshold per EC)")

    tt_fold_sum_res_best_long = melt(tt_fold_sum_res_best[,c('gs_size', METRICS_FOR_PER_SIZE_LINEPLOT),with=F], id.vars=c('gs_size'), value.vars=METRICS_FOR_PER_SIZE_LINEPLOT)
    tt_fold_sum_res_best_long$gs_size = factor(tt_fold_sum_res_best_long$gs_size, levels=gs_sizes, ordered=T)
    tt_fold_sum_res_best_long$value = as.numeric(tt_fold_sum_res_best_long$value)
    #plot with facet
    gs_size_comp = ggline(tt_fold_sum_res_best_long, 'gs_size', "value", facet.by="variable",
                          xlab='N genes per EC (Up+Down)', ylab='Value',
                          panel.labs = list(variable = METRICS_FOR_PER_SIZE_LINEPLOT_NAMES))
    gs_size_comp = ggpar(gs_size_comp, ylim=c(0,1)) + font("xy.text", size = 10)
    gs_size_comp_pdf = paste0(outdir, '/', 'performance_on_test_set_per_best_per_N_markers_used_faceted.',ec_preds$output_timestamp,'.pdf')
    ggsave(gs_size_comp_pdf, gs_size_comp, height=3, width=7, useDingbats=FALSE)

    #no facet
    gs_size_comp2 = ggline(tt_fold_sum_res_best_long, 'gs_size', "value", color="variable",
                           xlab='N genes per EC (Up+Down)', ylab='', legend='right') +
      scale_colour_colorblind(labels=METRICS_FOR_PER_SIZE_LINEPLOT_NAMES)
    gs_size_comp2 = ggpar(gs_size_comp2, ylim=c(0,1))
    gs_size_comp_pdf2 = paste0(outdir, '/', 'performance_on_test_set_per_best_per_N_markers_used.',ec_preds$output_timestamp,'.pdf')
    ggsave(gs_size_comp_pdf2, gs_size_comp2, height=4, width=6, useDingbats=FALSE)
  }
}
