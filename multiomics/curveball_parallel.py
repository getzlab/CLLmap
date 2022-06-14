#!/usr/bin/env python
# coding: utf-8

## Code for performing Curveball method for associations between mutations and molecular features
#Author: Ziao Lin

#from curveball import*
import pandas as pd
import numpy as np
import scipy.stats as stats
import statsmodels.stats.multitest
import copy
import random
import matplotlib.pyplot as plt
from random import sample, randint, shuffle
from multiprocessing import Pool


comut = pd.read_csv("/home/zlin/mcll/comut_for_curveball_wes_wgs_noContam_withMBL_withIGHV_20210106.tsv", sep="\t", index_col=False)

comut=comut[comut['MCLL'] ==1]
#comut_np = np.array(comut[comut.columns[1:204]])

#get contingency dict from the original matrix.
comut_np = np.array(comut[comut.columns[1:204]])


def whole_process(input_):
    input_matrix = input_[0]
    iter_ = int(input_[1])
    #merge
    contigency_dict_init = {}
    for i in range(input_matrix.shape[1]):
        print(i)
        contigency_dict_init[i] = {}
        for j in range(i+1, input_matrix.shape[1]):
            contigency_dict_init[i][j] = [sum((input_matrix[:, i] == 1) & (input_matrix[:, j] == 1))]
            contigency_dict_init[i][j].append(sum((input_matrix[:, i] == 0) & (input_matrix[:, j] == 1)))
            contigency_dict_init[i][j].append(sum((input_matrix[:, i] == 1) & (input_matrix[:, j] == 0)))
            contigency_dict_init[i][j].append(sum((input_matrix[:, i] == 0) & (input_matrix[:, j] == 0)))
    p_value_diff_mat_positive = np.array([]).reshape(0,20503)
    p_value_diff_mat_negative = np.array([]).reshape(0,20503)
    odds_ratio_mat_permute = np.array([]).reshape(0, 20503)
    r_hp = find_presences(input_matrix)
    comut_ = random_permute(input_matrix, r_hp, input_matrix.shape[1] * (iter_+1))
    contigency_dict_starting_point = {}
    for i in range(comut_.shape[1]):
        print(i)
        contigency_dict_starting_point[i] = {}
        for j in range(i+1, comut_.shape[1]):
            contigency_dict_starting_point[i][j] = [sum((comut_[:, i] == 1) & (comut_[:, j] == 1))]
            contigency_dict_starting_point[i][j].append(sum((comut_[:, i] == 0) & (comut_[:, j] == 1)))
            contigency_dict_starting_point[i][j].append(sum((comut_[:, i] == 1) & (comut_[:, j] == 0)))
            contigency_dict_starting_point[i][j].append(sum((comut_[:, i] == 0) & (comut_[:, j] == 0)))
    contigency_dict_ = copy.deepcopy(contigency_dict_starting_point)

    #get fisher p values for each pair.
    fisher_p_value_negative_asso_list = []
    fisher_p_value_positive_asso_list = []
    odds_ratio_list_permute = []
    for j in range(input_matrix.shape[1]):
        #print(i)
        for k in range(j+1, input_matrix.shape[1]):
            oddsratio, pvalue_ = stats.fisher_exact([[contigency_dict_[j][k][0], contigency_dict_[j][k][1]], [contigency_dict_[j][k][2], contigency_dict_[j][k][3]]], alternative='less')
            fisher_p_value_negative_asso_list.append(pvalue_)
            odds_ratio_list_permute.append(oddsratio)
            oddsratio, pvalue_ = stats.fisher_exact([[contigency_dict_[j][k][0], contigency_dict_[j][k][1]], [contigency_dict_[j][k][2], contigency_dict_[j][k][3]]], alternative='greater')
            fisher_p_value_positive_asso_list.append(pvalue_)
    p_value_diff_mat_positive = np.concatenate((p_value_diff_mat_positive, np.array(fisher_p_value_positive_asso_list).reshape(1, len(fisher_p_value_positive_asso_list))), axis = 0)
    p_value_diff_mat_negative = np.concatenate((p_value_diff_mat_negative, np.array(fisher_p_value_negative_asso_list).reshape(1, len(fisher_p_value_negative_asso_list))), axis = 0)
    odds_ratio_mat_permute = np.concatenate((odds_ratio_mat_permute, np.array(odds_ratio_list_permute).reshape(1, len(odds_ratio_list_permute))), axis = 0)
    #comut_ = copy.deepcopy(input_matrix)
    for i in range(5000):
        print (i)
        #ucll_contigency_dict_ = copy.deepcopy(ucll_contigency_dict_init)
        #r_hp=find_presences(ucll_comut)
        RM, contigency_dict=curve_ball(comut_, contigency_dict_)
        #p_value_diff = fisher_p_comparison_m_vs_u(RM, ucll_column_index, ucll_column_index)
        #p_value_diff_mat = np.concatenate((p_value_diff_mat, p_value_diff.reshape(1, len(p_value_diff))), axis = 0)
        comut_ = copy.deepcopy(RM)
        contigency_dict_ = copy.deepcopy(contigency_dict)
        #get fisher p values for each pair.
        fisher_p_value_negative_asso_list = []
        fisher_p_value_positive_asso_list = []
        odds_ratio_list_permute = []
        for j in range(input_matrix.shape[1]):
            #print(i)
            for k in range(j+1, input_matrix.shape[1]):
                oddsratio, pvalue_ = stats.fisher_exact([[contigency_dict_[j][k][0], contigency_dict_[j][k][1]], [contigency_dict_[j][k][2], contigency_dict_[j][k][3]]], alternative='less')
                fisher_p_value_negative_asso_list.append(pvalue_)
                odds_ratio_list_permute.append(oddsratio)
                oddsratio, pvalue_ = stats.fisher_exact([[contigency_dict_[j][k][0], contigency_dict_[j][k][1]], [contigency_dict_[j][k][2], contigency_dict_[j][k][3]]], alternative='greater')
                fisher_p_value_positive_asso_list.append(pvalue_)

        p_value_diff_mat_positive = np.concatenate((p_value_diff_mat_positive, np.array(fisher_p_value_positive_asso_list).reshape(1, len(fisher_p_value_positive_asso_list))), axis = 0)
        p_value_diff_mat_negative = np.concatenate((p_value_diff_mat_negative, np.array(fisher_p_value_negative_asso_list).reshape(1, len(fisher_p_value_negative_asso_list))), axis = 0)
        odds_ratio_mat_permute = np.concatenate((odds_ratio_mat_permute, np.array(odds_ratio_list_permute).reshape(1, len(odds_ratio_list_permute))), axis = 0)
    #save matrix.
    np.save('merge_positive_permute_mat_iter_' + str(iter_), p_value_diff_mat_positive)
    np.save('merge_negative_permute_mat_iter_' + str(iter_), p_value_diff_mat_negative)
    np.save('merge_odds_ratio_mat_iter_' + str(iter_), odds_ratio_mat_permute)

    #input matrix's fisher_exact p value.
    merge_p_value_negative_init = []
    merge_p_value_positive_init = []
    for j in range(input_matrix.shape[1]):
        print(j)
        for k in range(j+1, input_matrix.shape[1]):
            oddsratio, pvalue_ = stats.fisher_exact([[contigency_dict_init[j][k][0], contigency_dict_init[j][k][1]], [contigency_dict_init[j][k][2], contigency_dict_init[j][k][3]]], alternative='less')
            merge_p_value_negative_init.append(pvalue_)
            oddsratio, pvalue_ = stats.fisher_exact([[contigency_dict_init[j][k][0], contigency_dict_init[j][k][1]], [contigency_dict_init[j][k][2], contigency_dict_init[j][k][3]]], alternative='greater')
            merge_p_value_positive_init.append(pvalue_)

    #get empirical p values from curveball.
    empirical_p_value_negative_list = []
    count = -1
    for j in range(input_matrix.shape[1]):
        for k in range(j+1, input_matrix.shape[1]):
            count = count +1
            empirical_p_value_negative_list.append(sum(merge_p_value_negative_init[count] >= p_value_diff_mat_negative[:, count])/len(p_value_diff_mat_negative[:, count]))
    empirical_p_value_positive_list = []
    count = -1
    for j in range(input_matrix.shape[1]):
        for k in range(j+1, input_matrix.shape[1]):
            count = count +1
            empirical_p_value_positive_list.append(sum(merge_p_value_positive_init[count] >= p_value_diff_mat_positive[:, count])/len(p_value_diff_mat_positive[:, count]))

    #plot.
    plt.figure(figsize = (8,8))
    plt.hist(empirical_p_value_negative_list, bins = 50)
    plt.title("distribution of empirical pvalues for negative association (mutual exclusivity)")
    plt.xlabel("empirical p values of negative association")
    plt.ylabel("freq")
    plt.savefig('histogram_of_empirical_p_values_negative_association_' + str(iter_) + '.png', bbox_inches='tight')

    plt.figure(figsize = (8,8))
    plt.hist(empirical_p_value_positive_list, bins = 50)
    plt.title("distribution of empirical pvalues for positive association (co occurence)")
    plt.xlabel("empirical p values of positive association")
    plt.ylabel("freq")
    plt.savefig('histogram_of_empirical_p_values_positive_association_' + str(iter_) + '.png', bbox_inches='tight')

    #append pair name.
    pair_name_list = []
    for j in range(input_matrix.shape[1]):
        for k in range(j+1, input_matrix.shape[1]):
            pair_name_list.append(str(comut.columns[j+1]) + "_" + str(comut.columns[k+1]))

    #mutual exclusivity

    negative_significant_pair_index = np.where(statsmodels.stats.multitest.multipletests(np.array(empirical_p_value_negative_list), alpha = 0.10, method = "fdr_bh")[0])[0]
    mutual_exclusivity_pair = [pair_name_list[i] for i in negative_significant_pair_index]

    with open('merge_mutual_exclusivity_pairs_iter1000_starting_' + str(iter_) + '.txt', 'w') as filehandle:
        for listitem in mutual_exclusivity_pair:
            filehandle.write('%s\n' % listitem)

    #co-occurence

    positive_significant_pair_index = np.where(statsmodels.stats.multitest.multipletests(np.array(empirical_p_value_positive_list), alpha = 0.10, method = "fdr_bh")[0])[0]
    co_occurence_pair = [pair_name_list[i] for i in positive_significant_pair_index]

    with open('merge_co_occurence_pairs_iter1000_starting_' + str(iter_) + '.txt', 'w') as filehandle:
        for listitem in co_occurence_pair:
            filehandle.write('%s\n' % listitem)
    return (1)


def random_permute(input_matrix_, r_hp_, num_iterations=-1):
    num_rows, num_cols = input_matrix_.shape
    l = range(len(r_hp_))
    num_iters = 5*min(num_rows, num_cols) if num_iterations == -1 else num_iterations
    for rep in range(num_iters):
        AB = sample(l, 2)
        a = AB[0]
        b = AB[1]
        ab = set(r_hp_[a])&set(r_hp_[b]) # common elements
        l_ab=len(ab)
        l_a=len(r_hp_[a])
        l_b=len(r_hp_[b])
        if l_ab not in [l_a,l_b]:
            tot=list(set(r_hp_[a]+r_hp_[b])-ab)
            ab=list(ab)
            shuffle(tot)
            L=l_a-l_ab
            r_hp_[a] = ab+tot[:L]
            r_hp_[b] = ab+tot[L:]
    out_mat = np.zeros(input_matrix_.shape, dtype='int8') if num_cols >= num_rows else np.zeros(input_matrix_.T.shape, dtype='int8')
    for r in range(min(num_rows, num_cols)):
        out_mat[r, r_hp_[r]] = 1
    result = out_mat if num_cols >= num_rows else out_mat.T
    return result


def find_presences(input_matrix):
    num_rows, num_cols = input_matrix.shape
    hp = []
    iters = num_rows if num_cols >= num_rows else num_cols
    input_matrix_b = input_matrix if num_cols >= num_rows else np.transpose(input_matrix)
    for r in range(iters):
        hp.append(list(np.where(input_matrix_b[r] == 1)[0]))
    return hp


def curve_ball(input_matrix, contigency_dict_, num_iterations=-1):
    num_rows, num_cols = input_matrix.shape
    l = range(input_matrix.shape[1])
    num_iters = 5*min(num_rows, num_cols) if num_iterations == -1 else num_iterations
    out_mat = copy.deepcopy(input_matrix)
    r_hp=find_presences(out_mat)
    contigency_dict = copy.deepcopy(contigency_dict_)
    for rep in range(num_iters):
        AB = sample(l, 2)
        a = AB[0]
        b = AB[1]
        #min_ = min(a, b)
        #max_ = max(a, b)
        ab = set(r_hp[a])&set(r_hp[b]) # common elements
        l_ab=len(ab)
        l_a=len(r_hp[a])
        l_b=len(r_hp[b])
        if l_ab not in [l_a,l_b]:
            #tot=list(set(r_hp[a]+r_hp[b])-ab)
            ab=list(ab)
            #shuffle(tot)
            #L=l_a-l_ab
            #randomly pick one unique event from a.
            unique_a = set(r_hp[a]) - set(r_hp[b])
            sample_a = random.sample(unique_a, 1)
            #print(sample_a)
            #randomly pick one unique event from b.
            unique_b = set(r_hp[b]) - set(r_hp[a])
            sample_b = random.sample(unique_b, 1)
            #new_a = ab+tot[:L]
            #new_b = ab+tot[L:]
            #check which index are different from before, update the contigency table for those indices.
            #set_diff_a = new_a - r_hp[a]
            #set_diff_b = new_b - r_hp[b]
            #update the contigency table (4 values: upleft (min =1, max = 1), upright (min = 0, max = 1), bottomleft (min = 1, max = 0), bottomright (min = 0, max = 0)).
            for i in range(out_mat.shape[1]):
                if i != a and i != b:
                    #update contingency table for pairs involving a.
                    min_ = min(i, a)
                    max_ = max(i, a)
                    if out_mat[sample_a[0], i] ==1 and out_mat[sample_b[0], i] ==0:
                        contigency_dict[min_][max_][0] = contigency_dict[min_][max_][0] - 1
                        if contigency_dict[min_][max_][0] == -1:
                            print ([out_mat[sample_a[0], i], out_mat[sample_b[0], i], out_mat[sample_a[0], a], out_mat[sample_b[0], a]])
                            #print (out_mat[sample_b[0], i])
                            #print (out_mat[sample_a[0], a])
                            #print (out_mat[sample_b[0], a])
                        contigency_dict[min_][max_][1] = contigency_dict[min_][max_][1] + 1
                        contigency_dict[min_][max_][2] = contigency_dict[min_][max_][2] + 1
                        contigency_dict[min_][max_][3] = contigency_dict[min_][max_][3] - 1
                    elif out_mat[sample_a[0], i] ==0 and out_mat[sample_b[0], i] ==1:
                        contigency_dict[min_][max_][0] = contigency_dict[min_][max_][0] + 1
                        contigency_dict[min_][max_][1] = contigency_dict[min_][max_][1] - 1
                        contigency_dict[min_][max_][2] = contigency_dict[min_][max_][2] - 1
                        contigency_dict[min_][max_][3] = contigency_dict[min_][max_][3] + 1
                    #update contingency table for pairs involving b.
                    min_ = min(i, b)
                    max_ = max(i, b)
                    if out_mat[sample_a[0], i] ==1 and out_mat[sample_b[0], i] ==0:
                        contigency_dict[min_][max_][0] = contigency_dict[min_][max_][0] + 1
                        contigency_dict[min_][max_][1] = contigency_dict[min_][max_][1] - 1
                        contigency_dict[min_][max_][2] = contigency_dict[min_][max_][2] - 1
                        contigency_dict[min_][max_][3] = contigency_dict[min_][max_][3] + 1
                    elif out_mat[sample_a[0], i] ==0 and out_mat[sample_b[0], i] ==1:
                        contigency_dict[min_][max_][0] = contigency_dict[min_][max_][0] - 1
                        if contigency_dict[min_][max_][0] == -1:
                            print ([out_mat[sample_a[0], i], out_mat[sample_b[0], i], out_mat[sample_a[0], b], out_mat[sample_b[0], b]])
                            #print (out_mat[sample_b[0], i])
                            #print (out_mat[sample_a[0], a])
                            #print (out_mat[sample_b[0], a])
                        contigency_dict[min_][max_][1] = contigency_dict[min_][max_][1] + 1
                        contigency_dict[min_][max_][2] = contigency_dict[min_][max_][2] + 1
                        contigency_dict[min_][max_][3] = contigency_dict[min_][max_][3] - 1
            #update the new out_mat.
            out_mat[sample_a[0],a]=0
            out_mat[sample_b[0],a]=1
            out_mat[sample_a[0],b]=1
            out_mat[sample_b[0],b]=0
            r_hp=find_presences(out_mat)
            #r_hp[a] = ab+tot[:L]
            #r_hp[b] = ab+tot[L:]



    #out_mat = np.zeros(input_matrix.shape, dtype='int8') if num_cols >= num_rows else np.zeros(input_matrix.T.shape, dtype='int8')
    #out_mat[]
    #for r in range(min(num_rows, num_cols)):
    #out_mat[r, r_hp[r]] = 1
    #result = out_mat if num_cols >= num_rows else out_mat.T
    #result = out_mat
    return out_mat, contigency_dict


### MAIN ###
input_pairs = [[comut_np, i] for i in range(7)]

#start running
with Pool(13) as p:
    output = p.map(whole_process, input_pairs)
