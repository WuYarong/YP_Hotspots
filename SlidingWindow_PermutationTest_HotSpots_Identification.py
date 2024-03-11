#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import print_function

import argparse
import collections
import math
from multiprocessing import Pool, Queue
import numpy as np
import os
import random
import re
from scipy import stats
from statsmodels.stats.multitest import multipletests
import sys
import time


__author__ = 'wyr'


def read_snp_matrix(file, region_len):
    snp_count = 0
    snp_loc = []
    with open(file) as f:
        for line in f.readlines():
            line = line.replace('\n', '')
            if line.startswith('#') or re.match('\t', line):
                continue
            else:
                tmp = line.split('\t')
                snp_count += 1
                snp_loc.append(int(tmp[0].split('_')[-1]))
    mut_rate = float(snp_count) / region_len

    return mut_rate, snp_loc


def calcu_exp_snp_num(mut_rate, window_len, p_value=0.05):
    binorm_exp_snp= stats.binom.ppf(1-p_value, window_len, mut_rate)
    return binorm_exp_snp


def sliding_window_hotregions(snp_loc, mut_rate, p_value, window=500):
    hot_regions = []
    exp_snp_num = calcu_exp_snp_num(mut_rate, window, p_value)

    for i in range(0, len(snp_loc)):
        real_snp_num = 1
        real_snp_set = [snp_loc[i]]
        for j in range(i+1, len(snp_loc)):
            if snp_loc[j] - snp_loc[i] + 1 <= window:
                real_snp_num += 1
                real_snp_set.append(snp_loc[j])
            else:
                break

        if real_snp_num > exp_snp_num and real_snp_num >= 2:
            hot_regions.append([min(real_snp_set), max(real_snp_set), max(real_snp_set)- min(real_snp_set) + 1, exp_snp_num, real_snp_num, real_snp_set])

    return hot_regions


def merge_hotregions(hot_regions, mut_rate, snp_loc, p_value, merge_dist=1000):
    merge_hot_region = collections.OrderedDict()

    if hot_regions == '' or hot_regions == []:
        print('No hot regions were found. Quit.')
        sys.exit(-1)

    if len(merge_hot_region.keys()) == 0:
        merge_hot_region[str(hot_regions[0][0])] = hot_regions[0]

    for item in hot_regions:
        flag = 0
        for key in merge_hot_region.keys():
            # if there is an intersection
            if set(item[5]) & set(merge_hot_region[key][5]):
                flag = 1
                snp_set = list(set(item[5]) | set(merge_hot_region[key][5]))
                merge_hot_region[key][5] = sorted(snp_set)
                merge_hot_region[key][0] = str(min(int(item[0]), int(merge_hot_region[key][0])))
                merge_hot_region[key][1] = str(max(int(item[1]), int(merge_hot_region[key][1])))
                merge_hot_region[key][2] = int(merge_hot_region[key][1]) - int(merge_hot_region[key][0]) + 1
                merge_hot_region[key][3] = calcu_exp_snp_num(mut_rate, merge_hot_region[key][2], p_value)
                merge_hot_region[key][4] = len(snp_set)
            # If the distance between two hotspots is less than merge_dist and the merged SNP count remains higher than the expected value of the binomial distribution, merge; otherwise, don't merge
            elif abs(int(item[0]) - int(merge_hot_region[key][1])) <= merge_dist:
                snp_set = []
                for snp in snp_loc:
                    if snp >= min(int(item[0]), int(merge_hot_region[key][0])) and snp <= max(int(item[1]), int(merge_hot_region[key][1])):
                        snp_set.append(snp)
                exp_snp_num = calcu_exp_snp_num(mut_rate, max(int(item[1]), int(merge_hot_region[key][1])) - min(int(item[0]), int(merge_hot_region[key][0])) + 1, p_value)

                if len(snp_set) > exp_snp_num:
                    flag = 2
                    merge_hot_region[key][0] = str(min(int(item[0]), int(merge_hot_region[key][0])))
                    merge_hot_region[key][1] = str(max(int(item[1]), int(merge_hot_region[key][1])))
                    merge_hot_region[key][5] = sorted(snp_set)
                    merge_hot_region[key][2] = int(merge_hot_region[key][1]) - int(merge_hot_region[key][0]) + 1
                    #merge_hot_region[key][3] = exp_snp_num
                    merge_hot_region[key][3] = calcu_exp_snp_num(mut_rate, merge_hot_region[key][2], p_value)
                    merge_hot_region[key][4] = len(snp_set)

        if flag == 0:
            merge_hot_region[str(item[0])] = item

    out = open('tmp_merged_hot_region.txt', 'w')
    out.write('#Index\tStart\tEnd\tHR_len\tExpSNPNum\tRealSNPNum\tSNPSet\n')
    count = 0   
    for key in merge_hot_region.keys():
        count += 1
        out.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(count, merge_hot_region[key][0], merge_hot_region[key][1], merge_hot_region[key][2], calcu_exp_snp_num(mut_rate, merge_hot_region[key][2], p_value), merge_hot_region[key][4], merge_hot_region[key][5]))
    out.close()

    return merge_hot_region


def localtime():
    local_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
    return local_time


def fake_genome(snp_loc, chr_len, permutation_snp_set):
    tmp_snp_set = []
    while len(tmp_snp_set) < len(snp_loc):
        loc = random.randint(1, chr_len)
        if loc not in tmp_snp_set:
            tmp_snp_set.append(loc)
    permutation_snp_set.append(sorted(tmp_snp_set))

    return permutation_snp_set

def fake_genome_permutation_test(snp_loc, chr_len, repeat_num=1000, cpu_num=10):
    permutation_snp_set = []
    result = []

    p = Pool(cpu_num)
    for i in range(0, repeat_num):
        # if (i+1) % 1000 == 0:
        print('Run{}'.format(i+1))
        result.append(p.apply_async(fake_genome, args=(snp_loc, chr_len, permutation_snp_set)))
    p.close()
    p.join()

    for item in result:
        permutation_snp_set.append(item.get()[0])

    out = open('tmp_fake_genome_permutation_test.txt', 'w') 
    for item in permutation_snp_set:
        out.write('{}\n'.format('\t'.join(map(lambda x: str(x), item))))
    out.close()

    return permutation_snp_set


def get_freq_for_hot_region(merge_hot_region, permutation_snp_set_item, i, hr_permutation_freq):
    i = str(i)
    hr_permutation_freq[i] = dict()
    for k in merge_hot_region.keys():
        hr_len = merge_hot_region[k][2]
        hr_permutation_freq[i][k] = 0

        for j in range(0, len(permutation_snp_set_item)):
            count = 1
            for m in range(j+1, len(permutation_snp_set_item)):
                if int(permutation_snp_set_item[m]) - int(permutation_snp_set_item[j]) + 1 <= hr_len:
                    count += 1
                else:
                    break
            if count >= merge_hot_region[k][4]:
                hr_permutation_freq[i][k] += 1
    return hr_permutation_freq


def hot_region_probability(merge_hot_region, permutation_snp_set, outfile, p_value, mt_method='fdr_bh', cpu_num=10):
    permu_hot_region = []
    orign_pvalue = []
    hr_permutation_freq = dict()
    result = []

    p = Pool(cpu_num)
    for i in range(0, len(permutation_snp_set)):
        # if (i+1) % 1000 == 0:
        print('Rerun{}'.format(i+1))
        result.append(p.apply_async(get_freq_for_hot_region, args=(merge_hot_region, permutation_snp_set[i], i, hr_permutation_freq)))
    p.close()
    p.join()

    for item in result:
        hr_permutation_freq.update(item.get().items())

    tmp_feq = collections.defaultdict(lambda: 0)

    for k in hr_permutation_freq.keys():
        for sk in hr_permutation_freq[k].keys():
            tmp_feq[sk] += hr_permutation_freq[k][sk]

    for k in merge_hot_region.keys():
        if len(permutation_snp_set) == 0:
            continue
        prob = tmp_feq[k] / float(len(permutation_snp_set))
        orign_pvalue.append(prob)
        item = merge_hot_region[k]
        item.append(tmp_feq[k])
        item.append(prob)
        permu_hot_region.append(item)

    out = multipletests(pvals=orign_pvalue, alpha=p_value, method=mt_method)

    for item in range(0, len(permu_hot_region)):
        permu_hot_region[item].append(out[1][item])

    tmp_hr_out = open('tmp_hr_region_add_qvalue.txt', 'w')
    for item in permu_hot_region:
        tmp_hr_out.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(item[0], item[1], item[2], item[3], item[4], ','.join(map(lambda x: str(x), item[5])), item[6], item[7], item[8]))
    tmp_hr_out.close()

    hr_out = open(outfile, 'w')
    hr_out.write('\tStart\tEnd\tHR_len\tExpSNPNum\tRealSNPNum\tSNPSet\tfreq\tp-value\tq-value\n')

    count = 0
    for item in permu_hot_region:
        if item[-1] < p_value:
            count += 1
    tmp_d = int(math.log(count))

    count = 0
    for item in permu_hot_region:
        if item[-1] < p_value:
            count += 1
            hr_out.write('HR{:0{}d}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(count, tmp_d, item[0], item[1], item[2], item[3], item[4], ','.join(map(lambda x: str(x), item[5])), item[6], item[7], item[8]))
    hr_out.close()


def parameters():
    parser = argparse.ArgumentParser()

    parser.add_argument("-s", "--snp_matrix", required=True, help="file of SNP-INDEL matrix (tab-delimited)")
    parser.add_argument("-l", "--chr_len", required=True, type=int, help="chromosome length")
    parser.add_argument("-o", "--outfile", required=True, help="chromosome length")

    parser.add_argument("-w", "--window_len", default=1000, type=int, help="sliding window size")
    parser.add_argument("-d", "--merge_distance", default=10000, type=int, help="merged nearby hot regions with in this distance")
    parser.add_argument("-r", "--repeat_num", default=10000, type=int, help="repeat number for permutation test")
    parser.add_argument("-t", "--threads", default=10, type=int, help="threads")
    parser.add_argument("-p", "--pvalue", default=0.05, type=float, help="pvalue")

    args = parser.parse_args()
    help = parser.format_help()

    return args, help


def main():
    args, help = parameters()

    mut_rate, snp_loc = read_snp_matrix(args.snp_matrix, args.chr_len)
    print('Step1: sliding window to identify hot regions with a window size of {} bp, Time: {}'.format(args.window_len, localtime()))
    hot_regions = sliding_window_hotregions(snp_loc, mut_rate, p_value=args.pvalue, window=args.window_len)
    print('Step2: merging nearby hot regions within {} bp, Time: {}'.format(args.merge_distance, localtime()))
    merge_hot_regions = merge_hotregions(hot_regions, mut_rate, snp_loc, p_value=args.pvalue, merge_dist=args.merge_distance)
    print('Step3: starting permutation test, Time: {}'.format(localtime()))
    permutation_snp_set = fake_genome_permutation_test(snp_loc, args.chr_len, repeat_num=args.repeat_num, cpu_num=args.threads)
    print('Step4: obtaining the probability for each merged hot regions according to permutation test, Time: {}'.format(localtime()))
    hot_region_probability(merge_hot_regions, permutation_snp_set, args.outfile, p_value=args.pvalue, mt_method='fdr_bh', cpu_num=args.threads)
    print('Successly finished all steps, Time: {}'.format(localtime()))


if __name__ == '__main__':
    main()

