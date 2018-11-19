# Copyright (C) 2017 University of California, Los Angeles (UCLA)
# Emad Bahrami-Samani, Yi Xing
#
# Authors: Emad Bahrami-Samani, Yi Xing
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.

# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along with
# this program. If not, see http://www.gnu.org/licenses/.


import re,os,sys,warnings,numpy,scipy,math,itertools;

from scipy import stats;
from numpy import *;
from multiprocessing import Pool;
from scipy.optimize import fmin_cobyla
from scipy.optimize import fmin_l_bfgs_b
from math import log;
import fisher
import itertools as it
from collections import defaultdict

import mne
from mne import io
from mne.datasets import sample
from mne.stats import bonferroni_correction, fdr_correction

import progressbar

class Model :
  def __init__(self, minimum_coverage, rho, cutoff, snp_counter) :
    self.minimum_coverage = minimum_coverage
    self.rho = rho
    self.cutoff = cutoff
    self.snp_counter = snp_counter

  def perform_test(self, genotype_info, clip_reads, clip_coverage, rna_reads, rna_coverage) :
    asprin_test = defaultdict(lambda: defaultdict(list))
    asprin_pvalues = []
    asprin_odds_ratio = []

    bar = progressbar.ProgressBar(maxval=self.snp_counter, \
         widgets=[progressbar.Bar('=','[',']'), ' ', progressbar.Percentage()])
    progress_counter = 0
    bar.start()

    for chrom in genotype_info:
      for pos in genotype_info[chrom]:
        if (genotype_info[chrom][pos][2] != "none" and \
            clip_coverage[chrom][pos] >= self.minimum_coverage and \
            rna_coverage[chrom][pos] >= self.minimum_coverage) :
          asprin_test[chrom][pos] = \
            stats.fisher_exact([\
              [clip_reads[chrom][pos][genotype_info[chrom][pos][0]],\
               clip_reads[chrom][pos][genotype_info[chrom][pos][1]]],\
              [rna_reads[chrom][pos][genotype_info[chrom][pos][0]],\
               rna_reads[chrom][pos][genotype_info[chrom][pos][1]]]],\
              'two-sided')
          asprin_pvalues.append(asprin_test[chrom][pos][1])
          asprin_odds_ratio.append(asprin_test[chrom][pos][0])

          progress_counter += 1
          bar.update(progress_counter)

    bar.finish()

    alpha = 0.1
    reject_fdr, asprin_qvalues = fdr_correction(asprin_pvalues, \
                                              alpha=alpha, method='indep')

    return asprin_qvalues, asprin_odds_ratio

  def logit(self, x):
    if x<0.01:
      x=0.01;
    if x>0.99:
      x=0.99;
    return(log(x/(1-x)));
  
  def mats_clip(self, x,\
                clip_ref_count,\
                clip_alt_count,\
                rna_ref_count,\
                rna_alt_count, var):
  
    binomial_sum = -1 * (clip_ref_count*log(x + self.cutoff) +\
                         clip_alt_count*log(1 - (x + self.cutoff)) +\
                         rna_ref_count*log(x) +\
                         rna_alt_count*log(1 - x));
  
    multivar_sum = 0.1 * 0.5 * (pow(self.rho,2))/(1-pow(self.rho,2)) *\
                               (pow(stats.norm.ppf(x+self.cutoff),2) +\
                                pow(stats.norm.ppf(x),2) -\
                                2 * self.rho * stats.norm.ppf(x+self.cutoff) * \
                                stats.norm.ppf(x))
  
    return(binomial_sum + multivar_sum);
  
  def mats_clip_derivitive(self, x,\
                           clip_ref_count,\
                           clip_alt_count,\
                           rna_ref_count,\
                           rna_alt_count, var):
  
    res1 = -1 * (clip_ref_count/(x + self.cutoff) - clip_alt_count/(1 - (x + self.cutoff)));
    res1 += 0.1 * 0.5 * (pow(self.rho,2)) / (1-pow(self.rho,2)) *\
                        (2 * stats.norm.ppf(x + self.cutoff) - 2 * self.rho * stats.norm.ppf(x)) /\
                                               (stats.norm.pdf(stats.norm.ppf(x + self.cutoff)))
  
    res2 = -1 * (rna_ref_count / x - rna_alt_count / (1 - x));
    res2 += 0.1 * 0.5 * (pow(self.rho,2)) / (1-pow(self.rho,2)) *\
                        (2 * stats.norm.ppf(x) - 2 * self.rho * stats.norm.ppf(x + self.cutoff)) /\
                                               (stats.norm.pdf(stats.norm.ppf(x)));
    return(numpy.array(res1 + res2));
  
  def mats_rna(self, x,\
               clip_ref_count,\
               clip_alt_count,\
               rna_ref_count,\
               rna_alt_count, var):
  
    binomial_sum = -1 * (clip_ref_count * log(x) + \
                         clip_alt_count * log(1-x) + \
                         rna_ref_count * log((x + self.cutoff)) + \
                         rna_alt_count * log(1 - (x + self.cutoff)));
  
    multivar_sum = 0.1 * 0.5 * (pow(self.rho,2)) / (1-pow(self.rho,2)) * \
                   (pow(stats.norm.ppf(x), 2) + \
                    pow(stats.norm.ppf(x+self.cutoff), 2) - 2 * \
                    self.rho * stats.norm.ppf(x) * \
                    stats.norm.ppf(x + self.cutoff))
  
    return(binomial_sum+multivar_sum);
  
  def mats_rna_derivitive(self, x,\
                          clip_ref_count,\
                          clip_alt_count,\
                          rna_ref_count,\
                          rna_alt_count, var):
  
    res1 = -1 * (clip_ref_count / x - clip_alt_count / \
                                         (1 - x));
    res1 += 0.1 * 0.5 * (pow(self.rho,2)) / (1-pow(self.rho,2)) * \
                        (2 * stats.norm.ppf(x) - 2 * self.rho * \
                         stats.norm.ppf((x+self.cutoff))) / \
                              stats.norm.pdf(stats.norm.ppf(x))
  
    res2 = -1 * (rna_ref_count / (x + self.cutoff) - rna_alt_count / \
                                             (1-(x + self.cutoff)));
  
    res2 += 0.1 * 0.5 * (pow(self.rho,2)) / (1-pow(self.rho,2)) * \
                        (2 * stats.norm.ppf((x+self.cutoff)) - 2 * self.rho * \
                         stats.norm.ppf(x)) / \
                               stats.norm.pdf(stats.norm.ppf((x+self.cutoff)));
  
    return(numpy.array(res1+res2));
  
  def mats_individual(self, x,\
                      clip_ref_count,\
                      clip_alt_count,\
                      rna_ref_count,\
                      rna_alt_count, var):
  
    binomial_sum = -1 * (clip_ref_count * log(x[0]) + \
                         clip_alt_count * log(1 - x[0]) + \
                         rna_ref_count * log(x[1]) + \
                         rna_alt_count * log(1 - x[1]));
  
    multivar_sum = 0.1 * 0.5 * (pow(self.rho,2)) / (1 - pow(self.rho,2)) * \
                               (pow(stats.norm.ppf(x[0]),2) + \
                                pow(stats.norm.ppf(x[1]),2) - \
                                2 * self.rho * stats.norm.ppf(x[0]) * \
                                stats.norm.ppf(x[1]))
  
    return(binomial_sum+multivar_sum);
  
  def mats_individual_der(self, x,\
                          clip_ref_count,\
                          clip_alt_count,\
                          rna_ref_count,\
                          rna_alt_count, var):
  
    res1 = -1 * (clip_ref_count / x[0] -\
                 clip_alt_count / (1 - x[0]));
  
    res1 += 0.1 * 0.5 * (pow(self.rho,2)) / (1 - pow(self.rho,2)) *\
                        (2 * stats.norm.ppf(x[0]) -\
                         2 * self.rho * stats.norm.ppf(x[1])) /\
                        stats.norm.pdf(stats.norm.ppf(x[0]))
  
    res2 = -1 * (rna_ref_count / x[1] - \
                 rna_alt_count / (1 - x[1]));
  
    res2 += 0.1 * 0.5 * (pow(self.rho,2)) / (1-pow(self.rho,2)) * \
                        (2 * stats.norm.ppf(x[1]) - \
                         2 * self.rho * stats.norm.ppf(x[0])) / \
                        stats.norm.pdf(stats.norm.ppf(x[1]));
  
    return(numpy.array([res1, res2]));
  
  def mats_likelihood(self, x,\
                      clip_ref_count,\
                      clip_alt_count,\
                      rna_ref_count,\
                      rna_alt_count, var):
  
    sum = 0;
    N1 = clip_ref_count + clip_alt_count;
    N2 = rna_ref_count + rna_alt_count;
    sum += -0.5 * ((clip_ref_count - N1 * x[0]) *\
                   (clip_ref_count - N1 * x[0]) /\
                   (N1 * x[0]) +\
                   (clip_alt_count - N1 * (1 - x[0])) *\
                   (clip_alt_count - N1 * (1 - x[0])) /\
                   (N1 * (1 - x[0])));
  
    sum += -0.5 * ((rna_ref_count - N2 * x[1]) *\
                   (rna_ref_count - N2 * x[1]) /\
                   (N2 * x[1]) +\
                   (rna_alt_count - N2 * (1 - x[1])) *\
                   (rna_alt_count - N2 * (1 - x[1])) /\
                   (N2 * (1 - x[1])));
  
    sum += pow(self.logit(self.raf(clip_ref_count,clip_alt_count)) -\
               self.logit(self.raf(rna_ref_count,rna_alt_count)) -\
               self.logit(x[0]) + self.logit(x[1]), 2);
  
    return(sum);
  
  def MLE_iteration_constrain(self, clip_ref_count,\
                              clip_alt_count,\
                              rna_ref_count,\
                              rna_alt_count):
  
    clip_raf = self.raf(clip_ref_count,clip_alt_count);
    rna_raf = self.raf(rna_ref_count,rna_alt_count);
    iter_cutoff=1;
    iter_maxrun=100;
    count=0;
    previous_sum=0;
  
    while((iter_cutoff>0.01) and (count<=iter_maxrun)):
      count+=1;
      var1=0;
      var2=0;
      current_sum=0;
      likelihood_sum=0;
      if (clip_raf)>(rna_raf):
        xopt = fmin_l_bfgs_b(self.mats_clip,[rna_raf],\
                             self.mats_clip_derivitive,args=[clip_ref_count,\
                                                             clip_alt_count,\
                                                             rna_ref_count,\
                                                             rna_alt_count, var1],\
                             bounds=[[0.001,0.999-self.cutoff]], iprint=-1)
        theta2 = max(min(float(xopt[0]),1-self.cutoff),0);
        theta1 = theta2 + self.cutoff;
      else:
        xopt = fmin_l_bfgs_b(self.mats_rna,[clip_raf],\
                             self.mats_rna_derivitive,args=[clip_ref_count,\
                                                            clip_alt_count,\
                                                            rna_ref_count,\
                                                            rna_alt_count, var1],\
                             bounds=[[0.001,0.999-self.cutoff]], iprint=-1)
        theta1 = max(min(float(xopt[0]),1-self.cutoff),0);
        theta2 = theta1+self.cutoff;
      current_sum += float(xopt[1]);
      clip_raf = theta1;
      rna_raf = theta2;
      if count>1:
        iter_cutoff = abs(previous_sum - current_sum) / abs(previous_sum);
      previous_sum = current_sum;
    return([current_sum,[clip_raf,rna_raf,clip_raf,rna_raf,var1,var2]]);
  
  def MLE_iteration(self, clip_ref_count,\
                    clip_alt_count,\
                    rna_ref_count,\
                    rna_alt_count):
  
    clip_raf = self.raf(clip_ref_count,clip_alt_count);
    rna_raf = self.raf(rna_ref_count,rna_alt_count);
    iter_cutoff=1;
    iter_maxrun=100;
    count=0;
    previous_sum=0;
    while((iter_cutoff>0.01) and (count<=iter_maxrun)):
      count+=1;
      var1=0;
      var2=0;
      current_sum=0;
      likelihood_sum=0;
      xopt = fmin_l_bfgs_b(self.mats_individual,[clip_raf,rna_raf],\
                           self.mats_individual_der,args = [clip_ref_count,\
                                                       clip_alt_count,\
                                                       rna_ref_count,\
                                                       rna_alt_count, var1],\
                           bounds=[[0.01,0.99],[0.01,0.99]], iprint=-1);
  
      new_clip_raf = float(xopt[0][0]);
      current_sum += float(xopt[1]);
      new_rna_raf = float(xopt[0][1]);
      likelihood_sum += self.mats_likelihood([new_clip_raf,new_rna_raf],\
                                              clip_ref_count,\
                                              clip_alt_count,\
                                              rna_ref_count,\
                                              rna_alt_count, var1);
      clip_raf = new_clip_raf;
      rna_raf = new_rna_raf;
      if count > 1:
        iter_cutoff = abs(previous_sum-current_sum)/abs(previous_sum);
      previous_sum = current_sum;
    if count > iter_maxrun:
      return([current_sum,[clip_raf,rna_raf,0,0,var1,var2]]);
    return([current_sum,[clip_raf,rna_raf,clip_raf,rna_raf,var1,var2]]);
  
  def likelihood_test(self, clip_ref_count,\
                      clip_alt_count,\
                      rna_ref_count,\
                      rna_alt_count):
  
    res = self.MLE_iteration(clip_ref_count,\
                             clip_alt_count,\
                             rna_ref_count,\
                             rna_alt_count);
    if abs(res[1][2]-res[1][3]) <= self.cutoff:
      return(1);
    else:
      res_constrain = self.MLE_iteration_constrain(clip_ref_count,\
                                                   clip_alt_count,\
                                                   rna_ref_count,\
                                                   rna_alt_count);
      return(1-scipy.stats.chi2.cdf(2*(abs(res_constrain[0]-res[0])),1));
  
  def mats_model(self, n_original_diff):
    clip_ref_count = n_original_diff[0];
    clip_alt_count = n_original_diff[1];
    rna_ref_count = n_original_diff[2];
    rna_alt_count = n_original_diff[3];
    P = self.likelihood_test(clip_ref_count,\
                             clip_alt_count,\
                             rna_ref_count,\
                             rna_alt_count);
    return(P);
  
  def raf(self, reference_allele_count, alternative_allele_count):
    return(float(reference_allele_count /\
                (reference_allele_count + alternative_allele_count)))
  
  def perform_test_2(self, genotype_info, clip_reads, clip_coverage, rna_reads, rna_coverage) :
    asprin_pvalues = []
    bar = progressbar.ProgressBar(maxval=self.snp_counter, \
         widgets=[progressbar.Bar('=','[',']'), ' ', progressbar.Percentage()])
    progress_counter = 0
    bar.start()
  
    for chrom in genotype_info:
      for pos in genotype_info[chrom]:
        if (genotype_info[chrom][pos][2] != "none" and \
            clip_coverage[chrom][pos] >= self.minimum_coverage and \
            rna_coverage[chrom][pos] >= self.minimum_coverage) :
          asprin_pvalues.append(self.mats_model([\
              clip_reads[chrom][pos][genotype_info[chrom][pos][0]],\
              clip_reads[chrom][pos][genotype_info[chrom][pos][1]],\
              rna_reads[chrom][pos][genotype_info[chrom][pos][0]],\
              rna_reads[chrom][pos][genotype_info[chrom][pos][1]]]))
          progress_counter += 1
          bar.update(progress_counter)

    bar.finish()

    alpha = 0.1
    reject_fdr, asprin_qvalues = fdr_correction(asprin_pvalues, \
                                                alpha=alpha, method='indep')

    return asprin_qvalues


