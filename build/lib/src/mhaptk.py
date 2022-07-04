import math
import tabix
import pandas as pd
import matplotlib
import os

import matplotlib.pyplot as plt
import argparse, sys, os, time
import numpy as np
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
import scipy.stats  as ss
import seaborn as sns
import matplotlib.gridspec as gs
import matplotlib as mpl
import matplotlib.ticker as ticker
from pathlib import Path
import glob,re
from tqdm import tqdm
from scipy.stats import binom,chi2
import gzip
import matplotlib.gridspec as gridspec
import copy
import matplotlib.colors as mcolors
import random

matplotlib.use('Agg')
plt.switch_backend('agg')
version = '1.0'

class cpgAnnotation:
    def __init__(self, iChr, posArray):
        self.iChr = iChr
        self.posArray = posArray

class HT:
    def __init__(self, hChr, hStart, hEnd, HapMet, count, strand):
        self.hChr = hChr
        self.hStart = int(hStart)
        self.hEnd = int(hEnd)
        self.HapMet = HapMet
        self.count = int(count)
        self.WC = strand

class Stat:
    def __init__(self,
                 mhap_path: str,
                 cpg_path: str,
                 maxK: int,
                 minK: int,
                 stranded: str,
                 K: int):
        self.mhap_path = mhap_path
        self.cpg_path = cpg_path
        self.stranded = stranded
        self.K = K
        self.startK = minK
        self.endK = maxK
        self.data = tabix.open(mhap_path)
        self.statslist = ["MM", 'CHALM', 'PDR', 'MHL', 'MBS', 'MCR', 'Entropy', 'R2']
        self.statsbool = {"MM": False, 'CHALM': False, 'PDR': False, 'MHL': False, 'MBS': False, 'MCR': False,
                          'Entropy': False}
    def getBed(self,bed_path, line=0):
        self.bed = pd.read_csv(bed_path, sep='\t', header=None).iloc[line, :]
        self.start = int(self.bed[1])
        self.end = int(self.bed[2])
        self.Chr = self.bed[0]
        self.len = self.end - self.start
        self.name = str(self.bed[0]) + ':' + str(self.bed[1]) + '-' + str(self.bed[2])
    def Region(self,string):
        self.Chr = string.split(':')[0]
        self.start = int(string.split(':')[1].split('-')[0])-1
        self.end = int(string.split(':')[1].split('-')[1])
        self.len= self.end - self.start
        self.name = self.Chr + ":" + str(self.start) + "-" + str(self.end)

    def getrecord(self, strand: str = 'both'):

        self.record = []
        self.strand = strand

        try:
            records = self.data.query(self.Chr, self.start, self.end)
        except:
        #     print('record query error')
            print(self.Chr, self.start, self.end, 'the region can not be loaded')
            return self.record
        for rd in records:

            if strand == 'plus' and rd[5] != '+':
                continue
            if strand == 'minus' and rd[5] != '-':
                continue

            self.record.append(HT(rd[0], rd[1], rd[2], rd[3], rd[4], rd[5]))


        return self.record

    def aimming(self, stats_):

        for k in self.statslist:
            if k in stats_:
                self.statsbool[k] = True
            else:
                self.statsbool[k] = False
        #print('config reset as: ', self.statsbool)
        return self.statsbool
    def Kmer_gen(self, s, k=4):  # r:str
        '''
        To cut the kmers from a long sequense.
        - input:
         - s: input sequence, type : str
         - k: kmer length , type:int , default=k=4
        - output
         - ret_ : return as a bag of kmer(string type) in a list.
        '''
        ret_ = []
        for i in range(len(s) - k + 1):
            ret_.append(s[i: i + k])
        return ret_
    def calculating(self):
        self.kmer = self.K  # kmer片段长度
        self.result = {}  # 结果dictionay

        for a in self.statsbool.keys():
            if self.statsbool[a] == True:
                self.result[a] = ''

        if self.record ==[]:
            return self.result
        tBase = 0
        mBase = 0
        K4plus = 0
        nDR = 0
        nMR = 0
        cBase = 0
        nReads = 0

        # if self.record != []:
        #     for val in self.record:
        #         nReads += val.count
        #         tBase += val.count * len(val.HapMet)
        #         mBase += val.HapMet.count('1') * val.count
        #         if '1' in val.HapMet:
        #             cBase += val.HapMet.count('0') * val.count
        #         if len(val.HapMet) >= self.kmer:
        #             K4plus += val.count
        #             if '1' in val.HapMet:
        #                 nMR += val.count
        #
        #             if len(set(val.HapMet)) > 1:
        #                 nDR += val.count

        self.Total_dict = {}
        self.epiallell_dic = {}

        if self.record != []:
            for val in self.record:
                nReads += val.count
                tBase += val.count * len(val.HapMet)
                mBase += val.HapMet.count('1') * val.count
                if '1' in val.HapMet:
                    cBase += val.HapMet.count('0') * val.count
                if len(val.HapMet) >= self.kmer:
                    K4plus += val.count
                    if '1' in val.HapMet:
                        nMR += val.count

                    if len(set(val.HapMet)) > 1:
                        nDR += val.count
                if self.statsbool['Entropy']:
                    if len(val.HapMet) >= self.kmer:
                        rk_ = self.Kmer_gen(val.HapMet, self.kmer)
                        for rk in rk_:
                            try:
                                self.epiallell_dic[rk] += val.count
                            except:
                                self.epiallell_dic[rk] = val.count
            self.nReads, self.mBase, self.tBase, self.K4plus, self.nDR, self.nMR, self.cBase = nReads, mBase, tBase, K4plus, nDR, nMR, cBase

        else:
            self.nReads, self.mBase, self.tBase, self.K4plus, self.nDR, self.nMR, self.cBase=0, 0, 0, 0, 0, 0,0
        '''
        ——————————按照statsbool的情况计算统计量，要分情况计算————————————
        '''
        if self.statsbool['PDR']:
            if self.K4plus == 0:
                self.result['PDR'] = ''
            else:
                self.result['PDR'] = self.nDR /self.K4plus #if self.K4plus != 0 else 0
        if self.statsbool['MHL']:
            cpgAnn = self.tabixCPG(self.cpg_path)
            self.buildBinaryMatrix()
            self.result['MHL'] = self.getMHL()

        if self.statsbool['MM']:  # 简单
            if self.tBase ==0:
                self.result['MM'] = ''
            else:
                self.result['MM'] = self.mBase / self.tBase #if self.tBase != 0 else 0 # 甲基化位点数/所有位点数
        if self.statsbool['CHALM']:  # 简单
            if self.K4plus == 0:
                self.result['CHALM'] = ''
            else:
                self.result['CHALM'] = self.nMR / self.K4plus #if self.K4plus != 0 else 0 # 含甲基化的片段比例
        if self.statsbool['MBS']:
            self.result['MBS'] = self.get_mark_MBS(self.record)
        if self.statsbool['MCR']:
            if self.tBase ==0:
                self.result['MCR'] = ''
            else:
                self.result['MCR'] = self.cBase/self.tBase #if self.tBase != 0 else 0
        if self.statsbool['Entropy']:
            temp = 0
            all = sum(self.epiallell_dic.values())
            for i in self.epiallell_dic.values():
                temp += i/all*np.log2(i/all)

            self.result['Entropy'] = abs(-1/4 * temp)


        return self.result

    def Kmer(self, s):
        ret_ = []
        for i in range(len(s) - self.K + 1):
            ret_.append(s[i: i + self.K])
        return ret_

    def tabixCPG(self, htGZ, shift=500):
        self.cpg_path = htGZ
        tb = tabix.open(htGZ)
        records = tb.query(self.Chr, self.start - shift, self.end + shift)
        posArray = []
        for rd in records:
            if len(rd) < 3:
                continue
            posArray.append(int(rd[1]))
        self.cpgAnn = cpgAnnotation(self.Chr, posArray)

        return self.cpgAnn

    def get_r_MBS(self, r):
        Li = len(r)
        r_split = r.split("0")
        res = 0
        for lij in r_split:
            res = res + len(lij) ** 2
        res = res / Li ** 2
        return res

    def get_mark_MBS(self, r, cut=4):
        res = 0
        n_r = 0
        for i in r:
            if len(i.HapMet) < cut:
                continue
            res = res + self.get_r_MBS(i.HapMet) * int(i.count)
            n_r = n_r + int(i.count)
        try:
            res = res / n_r
            return res
        except:
            return  -1

    def getMHL(self):
        Nc = np.shape(self.MC)[1]
        obsK = max(np.sum(self.MC, 1))
        maxK = min([obsK, self.endK])
        if (self.startK > maxK):
            print("Error: startK is too large.\n")
            sys.exit(-1)

        uCount = np.zeros(maxK, dtype='int')
        mCount = np.zeros(maxK, dtype='int')
        tCount = np.zeros(maxK, dtype='int')

        for k in range(maxK):
            for j in range(Nc - k):

                xM0 = self.M0[..., j:(j + k + 1)]
                xM1 = self.M1[..., j:(j + k + 1)]
                xMC = self.MC[..., j:(j + k + 1)]
                uCount[k] += sum(((np.sum(xM0, 1) == (k + 1)) + 0) * self.count)
                mCount[k] += sum(((np.sum(xM1, 1) == (k + 1)) + 0) * self.count)
                tCount[k] += sum(((np.sum(xMC, 1) == (k + 1)) + 0) * self.count)

        mHL = sumFactor = 0.0
        for k in range(self.startK - 1, maxK):
            mHL += (k + 1.0) * mCount[k] / tCount[k]
            sumFactor += k + 1.0
        mHL = round(mHL / sumFactor * 100000) / 100000.0

        return mHL

    def info_to_file(self, k=4):
        tBase = 0
        mBase = 0
        K4plus = 0
        nDR = 0
        nMR = 0
        cBase = 0
        nReads = 0

        if self.record != []:
            for val in self.record:
                nReads += val.count
                tBase += val.count * len(val.HapMet)
                mBase += val.HapMet.count('1') * val.count
                if '1' in val.HapMet:
                    cBase += val.HapMet.count('0') * val.count
                if len(val.HapMet) >= k:
                    K4plus += val.count
                    if '1' in val.HapMet:
                        nMR += val.count

                    if len(set(val.HapMet)) > 1:
                        nDR += val.count

            self.nReads, self.mBase, self.tBase, self.K4plus, self.nDR, self.nMR, self.cBase = nReads, mBase, tBase, K4plus, nDR, nMR, cBase
            df = pd.DataFrame([self.Chr, self.start, self.end, nReads, mBase, tBase, K4plus, nDR, nMR]).T
            df.columns = ['chr', 'start', 'end', 'nReads', 'mBase', 'tBase', 'K4plus', 'nDR', 'nMR']
        else:
            self.nReads, self.mBase, self.tBase, self.K4plus, self.nDR, self.nMR = 0, 0, 0, 0, 0, 0
            df = pd.DataFrame([self.Chr, self.start, self.end, 0, 0, 0, 0, 0, 0]).T
            df.columns = ['chr', 'start', 'end', 'nReads', 'mBase', 'tBase', 'K4plus', 'nDR', 'nMR']

        return df

    def buildBinaryMatrix(self):
        posArray = self.cpgAnn.posArray
        posDict = dict()
        Nc = len(posArray)
        Nr = len(self.record)

        for i in range(Nc):
            posDict[self.cpgAnn.iChr + ":" + str(posArray[i])] = i

        self.MC = np.zeros((Nr, Nc), dtype='int')
        self.M0 = np.zeros((Nr, Nc), dtype='int')
        self.M1 = np.zeros((Nr, Nc), dtype='int')

        self.count = np.zeros(Nr, dtype='int')
        strand = np.zeros(Nr, dtype='int')

        for i in range(Nr):
            ht = self.record[i]
            pos = ht.hChr + ":" + str(ht.hStart)

            self.count[i] = ht.count
            if ht.WC == "+":
                strand[i] = 1
            # print(pos, posDict)
            if pos not in posDict:
                print("Haplotype positions are out of range.")
                sys.exit(-1)
            Idx = posDict[pos]
            HapMet = ht.HapMet
            for j in range(len(HapMet)):
                self.MC[i, Idx + j] = 1
                if HapMet[j] == "0":
                    self.M0[i, Idx + j] = 1
                elif HapMet[j] == "1":
                    self.M1[i, Idx + j] = 1
                else:
                    print("Haplotypes must be 0/1 only.")
                    sys.exit(-1)
        return [self.MC, self.M0, self.M1, self.count, strand]

    def count2pos(self, MC, M0, M1, count, pi, pj):
        if pi >= np.shape(MC)[1] or pj >= np.shape(MC)[1]:
            print("pi or pj are out of range.")
            sys.exit(-1)

        iAll = MC[..., pi] + MC[..., pj] == 2
        iN00 = M0[..., pi] + M0[..., pj] == 2
        iN01 = M0[..., pi] + M1[..., pj] == 2
        iN10 = M1[..., pi] + M0[..., pj] == 2
        iN11 = M1[..., pi] + M1[..., pj] == 2

        N00 = sum(count[np.logical_and(iAll, iN00)])
        N01 = sum(count[np.logical_and(iAll, iN01)])
        N10 = sum(count[np.logical_and(iAll, iN10)])
        N11 = sum(count[np.logical_and(iAll, iN11)])

        iAlli = MC[..., pi] == 1
        iAllj = MC[..., pj] == 1
        iNi0 = M0[..., pi] == 1
        iNi1 = M1[..., pi] == 1
        iNj0 = M0[..., pj] == 1
        iNj1 = M1[..., pj] == 1

        Ni0 = sum(count[np.logical_and(iAlli, iNi0)])
        Ni1 = sum(count[np.logical_and(iAlli, iNi1)])
        Nj0 = sum(count[np.logical_and(iAllj, iNj0)])
        Nj1 = sum(count[np.logical_and(iAllj, iNj1)])
        return [N00, N01, N10, N11], [Ni0, Ni1, Nj0, Nj1]

    def get_r2(self, N00, N01, N10, N11, Ni0, Ni1, Nj0, Nj1, cov=10):

        n = N00 + N01 + N10 + N11

        if n < cov:
            try:
                Pi = Ni1 / (Ni1 + Ni0)
                Pj = Nj1 / (Nj1 + Nj0)
                return (np.nan, np.nan)
            except:
                return (np.nan, np.nan)

        Pa = (N10 + N11) / n
        Pb = (N01 + N11) / n
        Pi = Ni1 / (Ni1 + Ni0)
        Pj = Nj1 / (Nj1 + Nj0)

        if (Pa * (1 - Pa) * Pb * (1 - Pb) == 0):
            return (0, 0.5)

        r = (N11 / n - Pa * Pb) / np.sqrt(Pa * (1 - Pa) * Pb * (1 - Pb))
        if r >= 1:
            return (r ** 2, 1)
        if r <= -1:
            return (r ** 2, 0)
        t = r * np.sqrt((n - 1) / (1 - r ** 2))

        pval = ss.t.cdf(t, n - 1)
        return (r ** 2, pval)


class GenomeWide:
    def __init__(self,
                 mhap_path:str,
                 cpg_path:str,
                 maxK:int,
                 minK:int,
                 stranded:str,
                 K:int):
        self.mhap_path = mhap_path
        self.cpg_path = cpg_path
        self.K = K
        self.startK = minK
        self.endK = maxK
        self.statslist = ["MM", 'CHALM', 'PDR', 'MHL', 'MBS', 'MCR', 'Entropy','R2']

        if stranded == 'both':
            self.stranded = 'both'
        else:
            self.stranded = '+' if stranded == 'plus' else '-'

    def Kmer(self, s):
        ret_ = []
        for i in range(len(s) - self.K + 1):
            ret_.append(s[i: i + self.K])
        return ret_

    def Kmer_gen(self, s, k=4):  # r:str
        '''
        To cut the kmers from a long sequense.
        - input:
         - s: input sequence, type : str
         - k: kmer length , type:int , default=k=4
        - output
         - ret_ : return as a bag of kmer(string type) in a list.
        '''
        ret_ = []
        for i in range(len(s) - k + 1):
            ret_.append(s[i: i + k])
        return ret_
    def get_pdr_reads(self, read):
        dict_i = {}
        read_ = self.Kmer(read)
        for r in read_:
            try:
                dict_i[r] += 1
            except:
                dict_i[r] = 1
        if dict_i == {}:
            return -1
        try:
            p0 = dict_i['0' * self.K]
        except:
            p0 = 0
        try:
            p1 = dict_i['1' * self.K]
        except:
            p1 = 0
        pdr = 1 - (p0 + p1) / sum(dict_i.values())
        return 1 - pdr

    def E4(self, r):
        result = {}
        for i in range(len(r)):
            result[i] = {}

        for i in range(len(r) - 3):
            frac = r[i:i + 4]
            for j in range(i, i + 4):
                try:
                    result[j][frac] += 1
                except:
                    result[j][frac] = 1
        return result

    def M4(self, r, count):
        l = len(r)
        fenzi_ =  np.zeros(self.endK)
        fenmu_ = np.zeros(self.endK)
        for i in range(min(l, self.endK)):
            tmp = self.Kmer_gen(r, i+1)
            fenmu_[i] += len(tmp) * count
            fenzi_[i] += tmp.count('1'*(i+1)) * count * (i+1)
        return fenzi_,fenmu_

    def MM(self):
        cpg_matrix = pd.read_csv(self.cpg_path, sep='\t', header=None).iloc[:, :2]
        cpg_matrix['stats'] = 0
        cpg_matrix['count'] = 0

        d_result = {}
        chr_i = 'chr1'
        cpg_index = 0
        mat_chr_i = cpg_matrix.loc[cpg_matrix.iloc[:, 0] == chr_i, :]

        for file in tqdm(gzip.open(self.mhap_path, 'rb')):
            file_info = file.decode().split('\t')
            chain = file_info[5]
            if self.stranded != 'both' :
                if chain != self.stranded:
                    continue
            chr_i_ = file_info[0]
            if chr_i_ != chr_i:
                chr_i = chr_i_
                mat_chr_i = cpg_matrix.loc[cpg_matrix.iloc[:, 0] == chr_i, :]
                cpg_index = mat_chr_i.index[0]
            cpg_i_s, cpg_i_e = file_info[1:3]

            while int(mat_chr_i.loc[cpg_index, 1]) < int(cpg_i_s):
                cpg_index += 1

            for l in range(len(file_info[3])):
                try:
                    d_result[cpg_index + l][0] += int(file_info[4]) * (1 if file_info[3][l] == '1' else 0)
                    d_result[cpg_index + l][1] += int(file_info[4])
                except:
                    d_result[cpg_index + l] = [int(file_info[4]) * (1 if file_info[3][l] == '1' else 0),
                                               int(file_info[4])]

        cpg_matrix.loc[d_result.keys(), ['stats', 'count']] = np.array(pd.DataFrame(d_result).T)
        cpg_matrix['MM'] = cpg_matrix['stats'] / cpg_matrix['count']
        return cpg_matrix

    def CHALM(self):
        cpg_matrix = pd.read_csv(self.cpg_path, sep='\t', header=None).iloc[:, :2]
        cpg_matrix['stats'] = 0
        cpg_matrix['count'] = 0

        d_result = {}
        chr_i = 'chr1'
        cpg_index = 0
        mat_chr_i = cpg_matrix.loc[cpg_matrix.iloc[:, 0] == chr_i, :]

        for file in tqdm(gzip.open(self.mhap_path, 'rb')):
            file_info = file.decode().split('\t')
            chain = file_info[5]
            if self.stranded == '+' and chain != '+':
                continue
            if self.stranded == '-' and chain != '-':
                continue
            if len(file_info[3]) < 4:
                continue

            chr_i_ = file_info[0]
            if chr_i_ != chr_i:
                chr_i = chr_i_
                mat_chr_i = cpg_matrix.loc[cpg_matrix.iloc[:, 0] == chr_i, :]
                cpg_index = mat_chr_i.index[0]
            cpg_i_s, cpg_i_e = file_info[1:3]

            while int(mat_chr_i.loc[cpg_index, 1]) < int(cpg_i_s):
                cpg_index += 1
            nMR = int(file_info[4]) if '1' in file_info[3] else 0
            K4plus = int(file_info[4])
            for l in range(len(file_info[3])):
                try:
                    d_result[cpg_index + l][0] += nMR
                    d_result[cpg_index + l][1] += K4plus
                except:
                    d_result[cpg_index + l] = [nMR, K4plus]


        cpg_matrix.loc[d_result.keys(), ['stats', 'count']] = np.array(pd.DataFrame(d_result).T)
        cpg_matrix['CHALM'] = cpg_matrix['stats'] / cpg_matrix['count']
        return cpg_matrix

    def PDR(self):
        cpg_matrix = pd.read_csv(self.cpg_path, sep='\t', header=None).iloc[:, :2]
        cpg_matrix['stats'] = 0
        cpg_matrix['count'] = 0

        d_result = {}
        chr_i = 'chr1'
        cpg_index = 0
        mat_chr_i = cpg_matrix.loc[cpg_matrix.iloc[:, 0] == chr_i, :]

        for file in tqdm(gzip.open(self.mhap_path, 'rb')):
            file_info = file.decode().split('\t')
            chain = file_info[5]
            if self.stranded != 'both':
                if chain != self.stranded:
                    continue
            if len(file_info[3]) < self.K:
                continue

            chr_i_ = file_info[0]
            if chr_i_ != chr_i:
                chr_i = chr_i_
                mat_chr_i = cpg_matrix.loc[cpg_matrix.iloc[:, 0] == chr_i, :]  # 所有染色体chr_i
                cpg_index = mat_chr_i.index[0]  # mat_chr_i开始的索引 mat包含了所有染色体chr

            cpg_i_s, cpg_i_e = file_info[1:3]

            while int(mat_chr_i.loc[cpg_index, 1]) < int(cpg_i_s):
                cpg_index += 1

            if '1' in file_info[3] and '0' in file_info[3]:
                nDR = int(file_info[4])
                K4plus = int(file_info[4])
            else:
                nDR = 0
                K4plus = int(file_info[4])
            for l in range(len(file_info[3])):
                try:
                    d_result[cpg_index + l][0] += nDR
                    d_result[cpg_index + l][1] += K4plus
                except:
                    d_result[cpg_index + l] = [nDR, K4plus]

        cpg_matrix.loc[d_result.keys(), ['stats', 'count']] = np.array(pd.DataFrame(d_result).T)
        cpg_matrix['meanPDR'] = cpg_matrix['stats'] / cpg_matrix['count']
        return cpg_matrix

    def MHL(self):

        cpg_matrix = pd.read_csv(self.cpg_path, sep='\t', header=None).iloc[:, :2]
        cpg_matrix['MHL'] = 0

        d_result = {}
        chr_i = 'chr1'
        cpg_index = 0
        mat_chr_i = cpg_matrix.loc[cpg_matrix.iloc[:, 0] == chr_i, :]

        for file in tqdm(gzip.open(self.mhap_path, 'rb')):
            file_info = file.decode().split('\t')
            chain = file_info[5]
            if self.stranded != 'both':
                if chain != self.stranded:
                    continue
            chr_i_ = file_info[0]
            if chr_i_ != chr_i:
                chr_i = chr_i_
                mat_chr_i = cpg_matrix.loc[cpg_matrix.iloc[:, 0] == chr_i, :]
                cpg_index = mat_chr_i.index[0]
            cpg_i_s, cpg_i_e = file_info[1:3]
            # print(mat_chr_i.loc[cpg_index, 1],cpg_i_s, cpg_i_e,file_info[3])
            while int(mat_chr_i.loc[cpg_index, 1]) < int(cpg_i_s):
                cpg_index += 1

            fenzi, fenmu = self.M4(file_info[3], int(file_info[4]))

            for l in range(len(file_info[3])):
                try:
                    d_result[cpg_index + l][0] = d_result[cpg_index + l][0] + fenzi
                    d_result[cpg_index + l][1] = d_result[cpg_index + l][1] + fenmu
                except:

                    d_result[cpg_index + l] = [fenzi, fenmu]


        result = {}
        for key in d_result.keys():
            result[key] = d_result[key][0] / d_result[key][1]
            sumfactor = self.endK - sum(np.isnan(result[key]))

            result[key] = sum(np.nan_to_num(result[key])) / sum(range(self.startK, sumfactor + 1))

        cpg_matrix.loc[d_result.keys(), 'MHL'] = np.array(pd.DataFrame(result, index=[0]).T)
        cpg_matrix = cpg_matrix.loc[d_result.keys()]
        return cpg_matrix

    def MBS(self):
        cpg_matrix = pd.read_csv(self.cpg_path, sep='\t', header=None).iloc[:, :2]
        cpg_matrix['stats'] = 0
        cpg_matrix['total'] = 0
        d_result = {}
        chr_i = 'chr1'
        cpg_index = 0
        mat_chr_i = cpg_matrix.loc[cpg_matrix.iloc[:, 0] == chr_i, :]

        for file in tqdm(gzip.open(self.mhap_path, 'rb')):
            file_info = file.decode().split('\t')
            chain = file_info[5]
            if self.stranded != 'both':
                if chain != self.stranded:
                    continue
            if len(file_info[3]) < self.K:
                continue

            chr_i_ = file_info[0]
            if chr_i_ != chr_i:
                chr_i = chr_i_
                mat_chr_i = cpg_matrix.loc[cpg_matrix.iloc[:, 0] == chr_i, :]
                cpg_index = mat_chr_i.index[0]
            cpg_i_s, cpg_i_e = file_info[1:3]

            while int(mat_chr_i.loc[cpg_index, 1]) < int(cpg_i_s):
                cpg_index += 1

            def get_r_MBS(r):
                Li = len(r)
                r_split = r.split("0")
                res = 0
                for lij in r_split:
                    res = res + len(lij) ** 2
                res = res / Li ** 2
                return res

            info = get_r_MBS(file_info[3])

            for l in range(len(file_info[3])):
                try:
                    d_result[cpg_index + l][0] += int(file_info[4]) * get_r_MBS(file_info[3])
                    d_result[cpg_index + l][1] += int(file_info[4])
                except:
                    d_result[cpg_index + l] = [int(file_info[4]) * get_r_MBS(file_info[3]), int(file_info[4])]


        cpg_matrix.loc[d_result.keys(), ['stats', 'count']] = np.array(pd.DataFrame(d_result).T)
        cpg_matrix['MBS'] = cpg_matrix['stats'] / cpg_matrix['count']

        return cpg_matrix

    def MCR(self):
        cpg_matrix = pd.read_csv(self.cpg_path, sep='\t', header=None).iloc[:, :2]
        cpg_matrix['stats'] = 0
        cpg_matrix['count'] = 0

        d_result = {}
        chr_i = 'chr1'
        cpg_index = 0
        mat_chr_i = cpg_matrix.loc[cpg_matrix.iloc[:, 0] == chr_i, :]

        for file in tqdm(gzip.open(self.mhap_path, 'rb')):
            file_info = file.decode().split('\t')
            chain = file_info[5]
            if self.stranded != 'both':
                if chain != self.stranded:
                    continue
            chr_i_ = file_info[0]
            if chr_i_ != chr_i:
                chr_i = chr_i_
                mat_chr_i = cpg_matrix.loc[cpg_matrix.iloc[:, 0] == chr_i, :]
                cpg_index = mat_chr_i.index[0]
            cpg_i_s, cpg_i_e = file_info[1:3]

            while int(mat_chr_i.loc[cpg_index, 1]) < int(cpg_i_s):
                cpg_index += 1

            if '1' in file_info[3]:
                cBase = int(file_info[4]) * file_info[3].count('0')
                tBase = int(file_info[4]) * len(file_info[3])
            else:
                cBase = 0
                tBase = int(file_info[4]) * len(file_info[3])
            for l in range(len(file_info[3])):
                try:
                    d_result[cpg_index + l][0] += cBase
                    d_result[cpg_index + l][1] += tBase
                except:
                    d_result[cpg_index + l] = [cBase, tBase]

        cpg_matrix.loc[d_result.keys(), ['stats', 'count']] = np.array(pd.DataFrame(d_result).T)
        cpg_matrix['MCR'] = cpg_matrix['stats'] / cpg_matrix['count']
        return cpg_matrix

        cpg_matrix.loc[d_result.keys(), ['stats', 'count']] = np.array(pd.DataFrame(d_result).T)
        cpg_matrix['MCR'] = cpg_matrix['stats'] / cpg_matrix['count']
        return cpg_matrix

    def Entropy(self, k=4):
        cpg_matrix = pd.read_csv(self.cpg_path, sep='\t', header=None).iloc[:, :2]
        d_result = {}
        cpg_matrix['Entropy'] = 0
        chr_i = 'chr1'
        cpg_index = 0

        mat_chr_i = cpg_matrix.loc[cpg_matrix.iloc[:, 0] == chr_i, :]
        for file in tqdm(gzip.open(self.mhap_path, 'rb')):
            file_info = file.decode().split('\t')
            chain = file_info[5]
            if self.stranded != 'both':
                if chain != self.stranded:
                    continue
            if len(file_info[3]) < self.K:
                continue

            chr_i_ = file_info[0]
            if chr_i_ != chr_i:
                chr_i = chr_i_
                mat_chr_i = cpg_matrix.loc[cpg_matrix.iloc[:, 0] == chr_i, :]
                cpg_index = mat_chr_i.index[0]

            cpg_i_s, cpg_i_e = file_info[1:3]
            while int(mat_chr_i.loc[cpg_index, 1]) < int(cpg_i_s):
                cpg_index += 1
            frac = self.Kmer(file_info[3])
            for val in frac:
                for l in range(len(file_info[3])):
                    try:
                        d_result[cpg_index + l]
                    except:
                        d_result[cpg_index + l] = {}
                    try:
                        d_result[cpg_index + l][val] += int(file_info[4])
                    except:
                        d_result[cpg_index + l][val] = int(file_info[4])


        result = {}

        def cal(total, dic):
            res = 0
            for val in dic.values():
                res += val / total * np.log2(val / total)
            res = -1 / 4 * res
            return abs(res)

        for key in d_result.keys():
            total = sum(c for c in d_result[key].values())
            result[key] = cal(total, d_result[key])

        cpg_matrix.loc[result.keys(), 'Entropy'] = np.array(pd.DataFrame(result, index=[0]).T)
        cpg_matrix = cpg_matrix.loc[result.keys()]
        return cpg_matrix

    def R2(self):
        cpg_matrix = pd.read_csv(self.cpg_path, sep='\t', header=None).iloc[:, :2]
        cpg_matrix['R2'] = None
        d_result = {}
        result = {}
        check_chr = 'cls'
        chr_i = 'cls'
        cpg_index = 0
        mat_chr_i = cpg_matrix.loc[cpg_matrix.iloc[:, 0] == chr_i, :]

        def cal(lst):

            n = sum(lst)
            if n == 0:
                return None
            n11, n00, n01, n10 = lst

            pa = (n11 + n10) / n
            pb = (n11 + n01) / n
            den = pa * pb * (1 - pa) * (1 - pb)
            if den == 0:
                return None

            pab = n11 / n - pa * pb
            return pab ** 2 / den if pab >= 0 else -1 * pab ** 2 / den

        def cal2():

            for key in d_result.keys():
                d_result[key]['right1'] = cal(d_result[key]['right1'])
                d_result[key]['right2'] = cal(d_result[key]['right2'])
                try:
                    result[key]
                    if d_result[key]['right1'] is not None:
                        result[key].append(d_result[key]['right1'])
                    if d_result[key]['right2'] is not None:
                        result[key].append(d_result[key]['right2'])
                except:
                    # print('chuangjian')
                    result[key] = []
                    if d_result[key]['right1'] is not None:
                        result[key].append(d_result[key]['right1'])
                    if d_result[key]['right2'] is not None:
                        result[key].append(d_result[key]['right2'])

                try:
                    d_result[key - 1]
                    if d_result[key - 1]['right1'] is not None:
                        result[key].append(d_result[key - 1]['right1'])
                except:
                    pass
                try:
                    d_result[key - 2]
                    if d_result[key - 2]['right2'] is not None:
                        result[key].append(d_result[key - 2]['right2'])
                except:
                    pass
                # print(result)
                try:
                    if result[key] == []:
                        result[key] = None
                    else:
                        result[key] = sum(result[key]) / len(result[key])
                except:
                    result[key] = None
            # print(d_result)
            cpg_matrix.loc[result.keys(), 'R2'] = np.array(pd.DataFrame(result, index=[0]).T)

        for file in tqdm(gzip.open(self.mhap_path, 'rb')):
            file_info = file.decode().split('\t')

            chr_i_ = file_info[0]
            if chr_i_ != chr_i:
                if check_chr == 'cls':
                    check_chr = chr_i_

                if check_chr != chr_i_:
                    cal2()
                    result = {}
                    d_result = {}

                chr_i = chr_i_
                mat_chr_i = cpg_matrix.loc[cpg_matrix.iloc[:, 0] == chr_i, :]
                cpg_index = mat_chr_i.index[0]

            cpg_i_s, cpg_i_e = file_info[1:3]
            # print(cpg_index, chr_i)
            # cpg_location = int(mat_chr_i.loc[cpg_index, 1])
            while int(mat_chr_i.loc[cpg_index, 1]) < int(cpg_i_s):
                cpg_index += 1

            c = int(file_info[4])

            for l in range(len(file_info[3])):
                try:
                    d_result[cpg_index + l]
                except:
                    d_result[cpg_index + l] = {'right1': [0, 0, 0, 0],
                                               'right2': [0, 0, 0, 0]}  # n11, n00, n01, n11

                try:
                    file_info[3][l + 1]  # 看有没有点

                    if file_info[3][l:l + 2] == '11':
                        d_result[cpg_index + l]['right1'][0] += c
                    elif file_info[3][l:l + 2] == '00':
                        d_result[cpg_index + l]['right1'][1] += c
                    elif file_info[3][l:l + 2] == '01':
                        d_result[cpg_index + l]['right1'][2] += c
                    else:
                        d_result[cpg_index + l]['right1'][3] += c
                except:
                    pass

                try:
                    file_info[3][l + 2]  # 看有没有点

                    if file_info[3][l] == '1' and file_info[3][l + 2] == '1':
                        d_result[cpg_index + l]['right2'][0] += c
                    elif file_info[3][l] == '0' and file_info[3][l + 2] == '0':
                        d_result[cpg_index + l]['right2'][1] += c
                    elif file_info[3][l] == '0' and file_info[3][l + 2] == '1':
                        d_result[cpg_index + l]['right2'][2] += c
                    else:
                        d_result[cpg_index + l]['right2'][3] += c
                except:
                    pass

        cal2()

        return cpg_matrix

class Tanghulu:
    def __init__(self,
                 mhap_path: str,
                 cpg_path: str):
        self.mhap_path = mhap_path
        self.cpg_path = cpg_path
        self.data = tabix.open(mhap_path)
        # print(self.data)
        # self.data.query('chr1',14434,14653)

    def Region(self, string):
        self.Chr = string.split(':')[0]
        self.start = int(string.split(':')[1].split('-')[0]) - 1
        self.end = int(string.split(':')[1].split('-')[1])
        self.len = self.end - self.start
        self.name = self.Chr + ":" + str(self.start) + "-" + str(self.end)

    def getrecord(self, strand: str = 'both'):

        self.record = []
        self.strand = strand
        # print(self.Chr, self.start, self.end)
        try:
            records = self.data.query(self.Chr, self.start, self.end)
        except:
            print(self.Chr, self.start, self.end, 'the region can not be loaded')
            return self.record
        for rd in records:

            if strand == 'plus' and rd[5] != '+':
                continue
            if strand == 'minus' and rd[5] != '-':
                continue

            self.record.append(HT(rd[0], rd[1], rd[2], rd[3], rd[4], rd[5]))

        return self.record

    def paint_tanghulu_plot(self, args, MC, M0, M1, count, strand, samplename, outpdf):
        cpgAnn = self.cpgAnn
        def scatter(df, ypos, color="black"):
            zero_position = list(df[df["zero"] == 1]["posArray"])
            one_position = list(df[df["one"] == 1]["posArray"])
            ax.scatter(zero_position, [ypos] * len(zero_position), c='white', edgecolors=color, s=102, zorder=2)
            ax.scatter(one_position, [ypos] * len(one_position), c=color, s=100, zorder=2)

        def line(array, ypos, color='black'):
            mc_string = ''.join(list(map(str, array)))
            over_start = cpgAnn.posArray[mc_string.find('1')]
            over_end = cpgAnn.posArray[mc_string.rfind('1')]
            ax.plot([over_start, over_end], [ypos, ypos], c=color, zorder=1)

        def ref_text(ref_posArray, xinches, yinches, a, b, xchange):
            text_y = (-1) / (xchange * a)
            text_minus_x = (-1) / (2 * b)
            for pos in ref_posArray:
                ax.text(pos + text_minus_x, text_y, '%1.0f' % (pos), rotation=45, color="grey", zorder=3)

        def count_text2(array, ypos, count_read, ref_posArray, a, b, color='black'):
            mc_string = ''.join(list(map(str, array)))
            over_end = cpgAnn.posArray[mc_string.rfind('1')]
            text_y = -0.2 / (4 * a)
            text_plus_x = (ref_posArray[-1] - ref_posArray[0]) / 50
            ax.text(over_end + text_plus_x, ypos + text_y, int(count_read), color=color)

        def legend_x(plt, ax):
            plt.scatter([], [], c='black', s=100, label='Methylated,+ strand')
            plt.scatter([], [], c='b', s=100, label='Methylated,- strand')
            plt.scatter([], [], c='white', edgecolors='black', s=102, label='Unmethylated,+ strand')
            plt.scatter([], [], c='white', edgecolors='b', s=102, label='Unmethylated,- strand')

            plt.legend(loc=3, bbox_to_anchor=(1.005, 0.5))
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['bottom'].set_visible(False)

        plt.clf()
        fig = plt.gcf()
        xinches = 20
        if args.merge:
            yinches = 6 / 24 * len(MC)
            yinches = 6 if yinches < 6 else yinches
        else:
            yinches = 6 / 24 * sum(count)
            # yinches = (6/24 * sum(count)) if sum(count) > 50 else (6/18 * sum(count))
            yinches = 8 if yinches < 8 else yinches

        fig.set_size_inches(xinches, yinches)
        # ax = plt.subplot(111, xlabel='pos')
        ax = plt.axes([0.1, 0.1, .7, .8], xlabel='Genomic position')
        ref_posArray = [cpgAnn.posArray[int(i)] for i in np.where(sum(MC) != 0)[0]]
        ax.scatter(ref_posArray, [0] * len(ref_posArray), c="#c0c0c0", s=100, zorder=2)
        ax.plot(ref_posArray, [0] * len(ref_posArray), c="#c0c0c0", zorder=1)
        try:
            if args.merge:
                a = yinches / len(MC)
            else:
                a = yinches / sum(count)
        except ZeroDivisionError as e:
            print("there's no records in this region from input mhap file, please choose another one\n" + e.message)

        b = xinches / (ref_posArray[-1] - ref_posArray[0])


        if args.merge:
            ref_text(ref_posArray, xinches, yinches, a, b, 1.3)
        else:
            ref_text(ref_posArray, xinches, yinches, a, b, 1.5)
        plt.xticks([])

        plt.yticks(range(sum(count)+1))
        name = self.Chr + ':' + str(self.start+1) + '-' +str(self.end)
        ax.set_title(name+ ' '+ '(' + samplename + ')')
        ypos = 1
        for i in range(0, len(MC)):
            df = pd.DataFrame(zip(cpgAnn.posArray, M0[i], M1[i], ), columns=["posArray", "zero", 'one'])
            if args.merge:
                if strand[i] > 0:
                    line(MC[i], ypos)
                    scatter(df, ypos)
                    count_text2(MC[i], ypos, str(count[i]), ref_posArray, a, b)
                else:
                    line(MC[i], ypos, color='b')
                    scatter(df, ypos, color='b')
                    count_text2(MC[i], ypos, str(count[i]), ref_posArray, a, b, color='b')
                ypos = ypos + 1
            else:
                for j in range(0, count[i]):
                    if strand[i] > 0:
                        line(MC[i], ypos)
                        scatter(df, ypos)
                    else:
                        line(MC[i], ypos, color='b')
                        scatter(df, ypos, color='b')
                    ypos = ypos + 1
        if args.merge:
            ax.axis(ymin=((-1) / (1.1 * a)))
        else:
            ax.axis(ymin=((-1) / (1 * a)))
        legend_x(plt, ax)

        fig.savefig(outpdf)
        plt.close()

    def MM(self):
        tBase = 0
        mBase = 0
        for val in self.record:
            tBase += val.count * len(val.HapMet)
            mBase += val.HapMet.count('1') * val.count
        #         print(tBase)
        if tBase == 0:
            mm = ''
        else:
            mm = mBase / tBase
        return round(mm, 8)

    def buildBinaryMatrix(self):
        posArray = self.cpgAnn.posArray
        posDict = dict()
        Nc = len(posArray)
        Nr = len(self.record)

        for i in range(Nc):
            posDict[self.cpgAnn.iChr + ":" + str(posArray[i])] = i

        self.MC = np.zeros((Nr, Nc), dtype='int')
        self.M0 = np.zeros((Nr, Nc), dtype='int')
        self.M1 = np.zeros((Nr, Nc), dtype='int')

        count = np.zeros(Nr, dtype='int')
        strand = np.zeros(Nr, dtype='int')

        for i in range(Nr):
            ht = self.record[i]
            pos = ht.hChr + ":" + str(ht.hStart)

            count[i] = ht.count
            if ht.WC == "+":
                strand[i] = 1
            # print(pos, posDict)
            if pos not in posDict:

                print("Haplotype positions are out of range.")
                sys.exit(-1)
            Idx = posDict[pos]
            HapMet = ht.HapMet
            for j in range(len(HapMet)):
                self.MC[i, Idx + j] = 1
                if HapMet[j] == "0":
                    self.M0[i, Idx + j] = 1
                elif HapMet[j] == "1":
                    self.M1[i, Idx + j] = 1
                else:
                    print("Haplotypes must be 0/1 only.")
                    sys.exit(-1)
        return [self.MC, self.M0, self.M1, count, strand]

    def tabixCPG(self, htGZ, shift=500):
        self.cpg_path = htGZ
        tb = tabix.open(htGZ)
        records = tb.query(self.Chr, self.start - shift, self.end + shift)
        posArray = []
        for rd in records:
            if len(rd) < 3:
                continue
            posArray.append(int(rd[1]))
        self.cpgAnn = cpgAnnotation(self.Chr, posArray)

        return self.cpgAnn

    def simulate(self):

        ref_posArray = [x for x in self.cpgAnn.posArray if (x >= self.start) and (x <= self.end)]
        if len(ref_posArray) > 20:
            print('nums of cpg are larger than 20 ',print(len(ref_posArray)))
            return

        CpgHpMat = np.zeros((len(self.record), len(self.cpgAnn.posArray)))
        cpgIdDic = {}
        for i, idx in enumerate(self.cpgAnn.posArray):
            cpgIdDic[idx] = i
        for i, rd in enumerate(self.record):
            idx = cpgIdDic[int(rd.hStart)]
            if int(rd.hStart) not in cpgIdDic:
                print('cpg not in dic')
            else:
                for j, cpg in enumerate(rd.HapMet):
                    if cpg == '0':
                        CpgHpMat[i, idx + j] = -1
                    if cpg == '1':
                        CpgHpMat[i, idx + j] = 1
        s_idx = self.cpgAnn.posArray.index(ref_posArray[0])
        e_idx = self.cpgAnn.posArray.index(ref_posArray[-1])
        CpgHpMat = CpgHpMat[:, s_idx:e_idx + 1]

        MM = []
        for i in range(CpgHpMat.shape[1]):
            MM.append(CpgHpMat[CpgHpMat[:, i] == 1].shape[0] / CpgHpMat[CpgHpMat[:, i] != 0].shape[0])
        MM = np.array(MM)

        dic = {}
        for i in range(CpgHpMat.shape[1]):
            if i == 0:
                to_left = [[0, 0], [0, 0]]  # w2w,w2b,b2w,b2b
            else:
                w2w = CpgHpMat[(CpgHpMat[:, i] == -1) & (CpgHpMat[:, i - 1] == -1), :].shape[0] \
                      / CpgHpMat[(CpgHpMat[:, i] == -1) & (CpgHpMat[:, i - 1] != 0), :].shape[0] \
                    if CpgHpMat[(CpgHpMat[:, i] == -1) & (CpgHpMat[:, i - 1] != 0), :].shape[0] != 0 else 0
                w2b = CpgHpMat[(CpgHpMat[:, i] == -1) & (CpgHpMat[:, i - 1] == 1), :].shape[0] \
                      / CpgHpMat[(CpgHpMat[:, i] == -1) & (CpgHpMat[:, i - 1] != 0), :].shape[0] \
                    if CpgHpMat[(CpgHpMat[:, i] == -1) & (CpgHpMat[:, i - 1] != 0), :].shape[0] != 0 else 0
                b2w = CpgHpMat[(CpgHpMat[:, i] == 1) & (CpgHpMat[:, i - 1] == -1), :].shape[0] \
                      / CpgHpMat[(CpgHpMat[:, i] == 1) & (CpgHpMat[:, i - 1] != 0), :].shape[0] \
                    if CpgHpMat[(CpgHpMat[:, i] == 1) & (CpgHpMat[:, i - 1] != 0), :].shape[0] != 0 else 0
                b2b = CpgHpMat[(CpgHpMat[:, i] == 1) & (CpgHpMat[:, i - 1] == 1), :].shape[0] \
                      / CpgHpMat[(CpgHpMat[:, i] == 1) & (CpgHpMat[:, i - 1] != 0), :].shape[0] \
                    if CpgHpMat[(CpgHpMat[:, i] == 1) & (CpgHpMat[:, i - 1] != 0), :].shape[0] != 0 else 0
                to_left = [[w2w, w2b], [b2w, b2b]]

            if i == CpgHpMat.shape[1] - 1:
                to_right = [[0, 0], [0, 0]]
            else:
                w2w = CpgHpMat[(CpgHpMat[:, i] == -1) & (CpgHpMat[:, i + 1] == -1), :].shape[0] \
                      / CpgHpMat[(CpgHpMat[:, i] == -1) & (CpgHpMat[:, i + 1] != 0), :].shape[0] \
                    if CpgHpMat[(CpgHpMat[:, i] == -1) & (CpgHpMat[:, i + 1] != 0), :].shape[0] != 0 else 0
                w2b = CpgHpMat[(CpgHpMat[:, i] == -1) & (CpgHpMat[:, i + 1] == 1), :].shape[0] \
                      / CpgHpMat[(CpgHpMat[:, i] == -1) & (CpgHpMat[:, i + 1] != 0), :].shape[0] \
                    if CpgHpMat[(CpgHpMat[:, i] == -1) & (CpgHpMat[:, i + 1] != 0), :].shape[0] != 0 else 0
                b2w = CpgHpMat[(CpgHpMat[:, i] == 1) & (CpgHpMat[:, i + 1] == -1), :].shape[0] \
                      / CpgHpMat[(CpgHpMat[:, i] == 1) & (CpgHpMat[:, i + 1] != 0), :].shape[0] \
                    if CpgHpMat[(CpgHpMat[:, i] == 1) & (CpgHpMat[:, i + 1] != 0), :].shape[0] != 0 else 0
                b2b = CpgHpMat[(CpgHpMat[:, i] == 1) & (CpgHpMat[:, i + 1] == 1), :].shape[0] \
                      / CpgHpMat[(CpgHpMat[:, i] == 1) & (CpgHpMat[:, i + 1] != 0), :].shape[0] \
                    if CpgHpMat[(CpgHpMat[:, i] == 1) & (CpgHpMat[:, i + 1] != 0), :].shape[0] != 0 else 0
                to_right = [[w2w, w2b], [b2w, b2b]]

            dic[i] = dict(to_right=to_right, to_left=to_left)

        #             print(dic)

        def dfs(x, y):
            if y < 0 or y >= CpgHpMat.shape[1]:
                return
            if CpgHpMat[x][y] != 0:
                return

            if y - 1 >= 0 and CpgHpMat[x][y - 1] != 0:
                check = 0 if CpgHpMat[x][y - 1] == -1 else 1
                CpgHpMat[x][y] = -1 if dic[y - 1]['to_right'][check][0] > random.uniform(0, 1) else 1
                dfs(x, y + 1)

            if y + 1 < CpgHpMat.shape[1] and CpgHpMat[x][y + 1] != 0:
                check = 0 if CpgHpMat[x][y - 1] == -1 else 1
                CpgHpMat[x][y] = -1 if dic[y + 1]['to_left'][check][0] > random.uniform(0, 1) else 1
                dfs(x, y - 1)

        for i in range(CpgHpMat.shape[0]):
            for j in range(CpgHpMat.shape[1]):
                dfs(i, j)

        CpgHpMat = (CpgHpMat == 1).astype(int).astype(str)

        Cpgframe = pd.DataFrame(CpgHpMat.astype(int))

        frame = Cpgframe.sample(n=10)

        # Simulation Annealing
        MM_i = []

        for i in range(frame.shape[1]):
            MM_i.append((sum(frame.iloc[:, i]) / frame.shape[0]))
        MM_i = np.array(MM_i)
        loss = sum((MM - MM_i) ** 2)

        while frame.shape[0] < 20:
            new_frame = Cpgframe.sample(n=1)
            new_frame = pd.concat([frame, new_frame])

            MM_i = []

            for i in range(new_frame.shape[1]):
                MM_i.append((sum(new_frame.iloc[:, i]) / new_frame.shape[0]))
            MM_i = np.array(MM_i)
            new_loss = sum((MM - MM_i) ** 2)
            if new_loss <= loss:
                frame = new_frame
                loss = new_loss
            else:
                if random.uniform(0, 1) < 0.1:
                    frame = new_frame
        MM_i = []
        for i in range(frame.shape[1]):
            MM_i.append((sum(frame.iloc[:, i]) / frame.shape[0]))

        frame_count = new_frame
        def mm_line(x):
            x = list(x)
            a = x.count(0)
            b = x.count(1)
            return a / b if b != 0 else 0

        frame_count['mm'] = frame_count.apply(mm_line, axis=1)

        frame_count.sort_values(by='mm', ascending=True, inplace=True)


        return np.matrix(frame_count.iloc[:, :-1]), loss, MM, MM_i, ref_posArray

class R2_c:
    def __init__(self,
                 mhap_path:str,
                 cpg_path:str,
                 stranded:str,):
        self.mhap_path = mhap_path
        self.cpg_path = cpg_path
        self.data = tabix.open(mhap_path)

        if stranded == 'both':
            self.stranded = 'both'
        else:
            self.stranded = '+' if stranded == 'plus' else '-'

    def tabixCPG(self, htGZ, shift=500):
        self.cpg_path = htGZ
        tb = tabix.open(htGZ)
        records = tb.query(self.Chr, self.start - shift, self.end + shift)
        posArray = []
        for rd in records:
            if len(rd) < 3:
                continue
            posArray.append(int(rd[1]))
        self.cpgAnn = cpgAnnotation(self.Chr, posArray)

        return self.cpgAnn

    def Region(self, string):
        self.Chr = string.split(':')[0]
        self.start = int(string.split(':')[1].split('-')[0]) - 1
        self.end = int(string.split(':')[1].split('-')[1])
        self.len = self.end - self.start
        self.name = self.Chr + ":" + str(self.start) + "-" + str(self.end)

    def getrecord(self, strand: str = 'both'):

        self.record = []
        self.strand = strand

        try:
            records = self.data.query(self.Chr, self.start, self.end)
        except:
        #     print('record query error')
            print(self.Chr, self.start, self.end, 'the region can not be loaded')
            return self.record
        for rd in records:

            if strand == 'plus' and rd[5] != '+':
                continue
            if strand == 'minus' and rd[5] != '-':
                continue

            self.record.append(HT(rd[0], rd[1], rd[2], rd[3], rd[4], rd[5]))

        # if self.record ==[] or len(self.record) <= 5:
        #     print('Get Nothing in this line of bed or get no more than 5 reads')
        #     return []
        # #print('Get record successfully')
        return self.record

    def hp_CPG(self, resultPath):
        # print(self.cpgAnn.posArray)
        ref_posArray = [x for x in self.cpgAnn.posArray if (x>=self.start) and (x<=self.end)]
        # print(ref_posArray)
        CpgHpMat = np.zeros((len(self.record), len(self.cpgAnn.posArray)))
        cpgIdDic = {}
        for i,idx in enumerate(self.cpgAnn.posArray):
            cpgIdDic[idx] = i
        for i,rd in enumerate(self.record):
            idx = cpgIdDic[int(rd.hStart)]
            if int(rd.hStart) not in cpgIdDic:
                print('cpg not in dic')
            else:
                for j,cpg in enumerate(rd.HapMet):
                    if cpg == '0':
                        CpgHpMat[i, idx+j] = -1
                    if cpg == '1':
                        CpgHpMat[i, idx+j] = 1
        s_idx = self.cpgAnn.posArray.index(ref_posArray[0])
        e_idx = self.cpgAnn.posArray.index(ref_posArray[-1])
        CpgHpMat = CpgHpMat[:,s_idx:e_idx+1]
        # print(ref_posArray)
        fig = plt.figure()
        ax = fig.gca()
        im = sns.heatmap(CpgHpMat, cmap='Greys',cbar=False)
        # mpl.colors.to_rgb('#D3D3D3')
        cmap, norm = mcolors.from_levels_and_colors([-1,0,1], ['white','#D3D3D3'])
        plt.pcolor(CpgHpMat, cmap=cmap, norm=norm)
        cb = im.figure.colorbar(im.collections[0])
        cb.set_ticks([-1,0,1])
        # im = ax.imshow(CpgHpMat, cmap='Greys')

        # cb1 = plt.colorbar(im, fraction=0.03, pad=0.05)
        # tick_locator = ticker.MaxNLocator(nbins=5)  # colorbar上的刻度值个数
        # cb1.locator = tick_locator
        # cb1.set_ticks([-1, 0, 1])
        # cb1.update_ticks()
        plt.title('hp of cpg')
        plt.xticks([])
        plt.yticks([])

        plt.savefig(resultPath + 'hp_cpg.pdf', dpi=1200)
        return  CpgHpMat

    def get_r2(self, N00, N01, N10, N11, Ni0, Ni1, Nj0, Nj1, cov=10):

        n = N00 + N01 + N10 + N11

        if n < cov:
            try:
                Pi = Ni1 / (Ni1 + Ni0)
                Pj = Nj1 / (Nj1 + Nj0)
                return (np.nan, np.nan)
            except:
                return (np.nan, np.nan)

        Pa = (N10 + N11) / n
        Pb = (N01 + N11) / n
        Pi = Ni1 / (Ni1 + Ni0)
        Pj = Nj1 / (Nj1 + Nj0)

        if (Pa * (1 - Pa) * Pb * (1 - Pb) == 0):
            return (0, 0.5)

        r = (N11 / n - Pa * Pb) / np.sqrt(Pa * (1 - Pa) * Pb * (1 - Pb))
        if r >= 1:
            return (r ** 2, 1)
        if r <= -1:
            return (r ** 2, 0)
        t = r * np.sqrt((n - 1) / (1 - r ** 2))

        pval = ss.t.cdf(t, n - 1)
        return (r ** 2, pval)

    def D_r2(self, n_11, n_10, n_01, n_00):

        n = n_11 + n_10 + n_01 + n_00
        if n == 0:
            return [np.nan, np.nan]

        P_A = (n_10 + n_11) / n
        P_B = (n_01 + n_11) / n

        # calculate pval
        D = n_11 / n - P_A * P_B
        Num = D * D

        Den = P_A * (1 - P_A) * P_B * (1 - P_B)

        if Den == 0:
            r2 = np.nan
        else:
            r2 = Num / Den
            if D < 0:
                r2 = -1 * r2  # add sign to r2

        # p-value
        pval = 1-binom.cdf(k=n_11-1, n=n, p=P_A * P_B)
        return [r2, pval]

    def count2pos(self,  MC, M0, M1, count, pi, pj):
        if pi >= np.shape(MC)[1] or pj >= np.shape(MC)[1]:
            print("pi or pj are out of range.")
            sys.exit(-1)

        iAll = MC[..., pi] + MC[..., pj] == 2
        iN00 = M0[..., pi] + M0[..., pj] == 2
        iN01 = M0[..., pi] + M1[..., pj] == 2
        iN10 = M1[..., pi] + M0[..., pj] == 2
        iN11 = M1[..., pi] + M1[..., pj] == 2

        N00 = sum(count[np.logical_and(iAll, iN00)])
        N01 = sum(count[np.logical_and(iAll, iN01)])
        N10 = sum(count[np.logical_and(iAll, iN10)])
        N11 = sum(count[np.logical_and(iAll, iN11)])

        iAlli = MC[..., pi] == 1
        iAllj = MC[..., pj] == 1
        iNi0 = M0[..., pi] == 1
        iNi1 = M1[..., pi] == 1
        iNj0 = M0[..., pj] == 1
        iNj1 = M1[..., pj] == 1

        Ni0 = sum(count[np.logical_and(iAlli, iNi0)])
        Ni1 = sum(count[np.logical_and(iAlli, iNi1)])
        Nj0 = sum(count[np.logical_and(iAllj, iNj0)])
        Nj1 = sum(count[np.logical_and(iAllj, iNj1)])
        return [N00, N01, N10, N11], [Ni0,Ni1,Nj0,Nj1]

    def buildBinaryMatrix(self):
        posArray = self.cpgAnn.posArray
        posDict = dict()
        Nc = len(posArray)
        Nr = len(self.record)

        for i in range(Nc):
            posDict[self.cpgAnn.iChr + ":" + str(posArray[i])] = i

        self.MC = np.zeros((Nr, Nc), dtype='int')
        self.M0 = np.zeros((Nr, Nc), dtype='int')
        self.M1 = np.zeros((Nr, Nc), dtype='int')

        count = np.zeros(Nr, dtype='int')
        strand = np.zeros(Nr, dtype='int')

        for i in range(Nr):
            ht = self.record[i]
            pos = ht.hChr + ":" + str(ht.hStart)

            count[i] = ht.count
            if ht.WC == "+":
                strand[i] = 1
            # print(pos, posDict)
            if pos not in posDict:

                print("Haplotype positions are out of range.")
                sys.exit(-1)
            Idx = posDict[pos]
            HapMet = ht.HapMet
            for j in range(len(HapMet)):
                self.MC[i, Idx + j] = 1
                if HapMet[j] == "0":
                    self.M0[i, Idx + j] = 1
                elif HapMet[j] == "1":
                    self.M1[i, Idx + j] = 1
                else:
                    print("Haplotypes must be 0/1 only.")
                    sys.exit(-1)
        return [self.MC, self.M0, self.M1, count, strand]

    def paint_rsquare_plot(self, samplename, R_matrix, xx_matrix, MC, outpdf):
        cpgAnn = self.cpgAnn
        # print(R_matrix)
        def ref_text(ax, ref_posArray, xinches, yinches, a, b):
            text_y = (-7) / (4 * a)
            text_minus_x = (-1) / (3 * b)
            for pos in ref_posArray:
                ax.text(pos + text_minus_x, text_y, '%1.0f' % (pos), rotation=45, color="grey", zorder=3)

        def shape(ax, nn, m, l, ref_posArray, R_matrix, xx_matrix):
            grid = np.mgrid[ref_posArray[0]:ref_posArray[0] + nn * m:m, 0:nn * m:m].reshape(2, -1).T
            patches = []
            colors = []
            for i in range(0, len(grid)):
                xx = int(i / nn)
                yy = i % nn
                if grid[i][1] - m * xx >= 0:
                    down_kuang = mpatches.Rectangle(grid[i] + [m * yy + m, -m * xx], l, l, angle=45, ec="white")
                    patches.append(down_kuang)
                    # print(R_matrix[xx][yy + 1])
                    colors.append(R_matrix[xx][yy + 1] / 2)
                    center_word_x = grid[i][0] + m * yy + 0.7 * m
                    center_word_y = grid[i][1] + m * (1 - xx)
                    if xx_matrix[xx][yy + 1] == 1:
                        ax.text(center_word_x, center_word_y, '%1.4f' % (R_matrix[xx][yy + 1]), color="black")

            # collection = PatchCollection(patches, cmap=plt.cm.Beds, ec="black", alpha = 0.75)
            collection = PatchCollection(patches, cmap=plt.cm.seismic, ec="black")
            # colors.append(1)
            # colors.append(-1)
            colors = np.array(colors)
            # colors = np.linspace(0, 1, len(patches))
            collection.set_array(colors)
            collection.set_clim(-1, 1)
            ax.add_collection(collection)

        ref_posArray_tmp = [cpgAnn.posArray[int(i)] for i in np.where(sum(MC) != 0)[0]]
        ref_posArray = [x for x in ref_posArray_tmp if (x > self.start) and (x <= self.end)]
        n = len(ref_posArray)
        plt.clf()
        fig = plt.gcf()
        yinches = n if n >= 8 else 8
        xinches = 1.5 * yinches
        fig.set_size_inches(xinches, yinches)
        ax = plt.subplot(111, xlabel='pos')

        nn = n - 1  # n shows how much sites, while nn shows Rsquare rows, l shows the square line length, m shows xiebian/2
        l = (ref_posArray[-1] - ref_posArray[0]) / math.sqrt(2) / nn
        m = math.sqrt(2) / 2 * l
        a = yinches / (ref_posArray[-1] - ref_posArray[0])
        b = xinches / (ref_posArray[-1] - ref_posArray[0])

        ave_posArray = np.arange(ref_posArray[0], ref_posArray[-1] + 1, (ref_posArray[-1] - ref_posArray[0]) / nn)

        ax.scatter(ave_posArray, [0] * len(ave_posArray), c="#c0c0c0", s=100)
        ax.scatter(ref_posArray, [-(1 / a)] * len(ref_posArray), c="#c0c0c0", s=100)
        ax.plot(ref_posArray, [-(1 / a)] * len(ref_posArray), c="#c0c0c0", zorder=1)
        for i in range(0, n):
            ax.plot([ave_posArray[i], ref_posArray[i]], [0, -(1 / a)], c="#c0c0c0", zorder=1)
            pass
        ref_text(ax, ref_posArray, xinches, yinches, a, b)

        plt.xticks([])
        plt.yticks([])
        ax.set_title('Rsquare plot(' + samplename + ') ' + self.name)
        shape(ax, nn, m, l, ref_posArray, R_matrix, xx_matrix)

        ax.axis(ymax=yinches / xinches * (ax.axis()[1] - ax.axis()[0]))
        ax.axis(ymin=(n - 2) * m - yinches / xinches * (ax.axis()[1] - ax.axis()[0]))

        fig.savefig(outpdf)

    def paint_rsquare_heatmap(self, R_matrix, MC, outpdf):
        cpgAnn = self.cpgAnn
        ref_posArray_tmp = [cpgAnn.posArray[int(i)] for i in np.where(sum(MC) != 0)[0]]
        ref_posArray = [x for x in ref_posArray_tmp if (x > self.start) and (x <= self.end)]

        n = len(ref_posArray)

        plt.clf()
        fig = plt.gcf()
        yinches = n if n >= 8 else 8
        xinches = 1.5 * yinches
        fig.set_size_inches(xinches, yinches)
        ax = plt.subplot(111, xlabel='pos')
        nn = n - 1  # n shows how much sites, while nn shows Rsquare rows, l shows the square line length, m shows xiebian/2
        mask = np.zeros_like(R_matrix)
        for i in range(1, len(mask)):
            for j in range(0, i + 1):
                mask[i][j] = True
        p1 = sns.heatmap(R_matrix, cmap='Reds',square=True, mask=mask)
        ax.set_xticklabels(ref_posArray)
        ax.set_yticklabels(ref_posArray, rotation='horizontal')


        s1 = p1.get_figure()
        s1.savefig(outpdf)

    def join_pic(self,R_matrix, xx_matrix, MC, CpgHpMat,outpdf):
        cpgAnn = self.cpgAnn

        # print(R_matrix)
        def ref_text(ax, ref_posArray, xinches, yinches, a, b):
            text_y = (-7) / (4 * a)
            text_minus_x = (-1) / (3 * b)
            # for pos in ref_posArray:
            #     ax.text(pos + text_minus_x, text_y, '%1.0f' % (pos), rotation=45, color="grey", zorder=3)

        def shape(ax, nn, m, l, ref_posArray, R_matrix, xx_matrix):
            grid = np.mgrid[ref_posArray[0]:ref_posArray[0] + nn * m:m, 0:nn * m:m].reshape(2, -1).T
            patches = []
            colors = []
            for i in range(0, len(grid)):
                xx = int(i / nn)
                yy = i % nn
                if grid[i][1] - m * xx >= 0:
                    down_kuang = mpatches.Rectangle(grid[i] + [m * yy + m, -m * xx], l, l, angle=45, ec="white")
                    patches.append(down_kuang)
                    # print(R_matrix[xx][yy + 1])
                    colors.append(R_matrix[xx][yy + 1])
                    center_word_x = grid[i][0] + m * yy + 0.7 * m
                    center_word_y = grid[i][1] + m * (1 - xx)
                    # if xx_matrix[xx][yy + 1] == 1:
                    #     ax.text(center_word_x, center_word_y, '%1.4f' % (R_matrix[xx][yy + 1]), color="black")

            # collection = PatchCollection(patches, cmap=plt.cm.Beds, ec="black", alpha = 0.75)
            collection = PatchCollection(patches, cmap='RdBu_r', ec="white")
            # colors.append(1)
            # colors.append(-1)
            colors = np.array(colors)
            # print(colors)
            # colors = np.linspace(0, 1, len(patches))
            collection.set_array(colors)
            collection.set_clim(-1,1)
            ax.add_collection(collection)

        ref_posArray_tmp = [cpgAnn.posArray[int(i)] for i in np.where(sum(MC) != 0)[0]]
        ref_posArray = [x for x in ref_posArray_tmp if (x > self.start) and (x <= self.end)]
        n = len(ref_posArray)
        plt.clf()
        yinches = n if n >= 8 else 8
        xinches = 2 * yinches
        gs = gridspec.GridSpec(yinches * 3, xinches + 3)
        fig = plt.figure()
        ax = plt.subplot(gs[yinches:, :-3])
        ax2 = plt.subplot(gs[:yinches, :-3])
        plt.title(self.Chr+":"+str(ref_posArray[0])+"-"+str(ref_posArray[-1]),fontdict={"size":80})
        sns.heatmap(CpgHpMat, cmap='Greys', cbar=None, square=False, linewidths=0)
        cmap, norm = mcolors.from_levels_and_colors([-1,0,1], ['white','#D3D3D3'])
        plt.pcolor(CpgHpMat, cmap=cmap, norm=norm)

        # plt.figure(figsize=(xinches,yinches))
        # cb1 = plt.colorbar(im, fraction=0.03, pad=0.05)
        # tick_locator = ticker.MaxNLocator(nbins=5)  # colorbar上的刻度值个数
        # cb1.locator = tick_locator
        # cb1.set_ticks([-1, 0, 1])
        # cb1.update_ticks()
        # plt.title('hp of cpg')
        plt.xticks([])
        plt.yticks([])
        plt.subplots_adjust(wspace=0)

        fig.set_size_inches(xinches, yinches*2)


        nn = n - 1  # n shows how much sites, while nn shows Rsquare rows, l shows the square line length, m shows xiebian/2
        l = (ref_posArray[-1] - ref_posArray[0]) / math.sqrt(2) / nn
        m = math.sqrt(2) / 2 * l
        a = yinches / (ref_posArray[-1] - ref_posArray[0])
        b = xinches / (ref_posArray[-1] - ref_posArray[0])

        ave_posArray = np.arange(ref_posArray[0], ref_posArray[-1] + 1, (ref_posArray[-1] - ref_posArray[0]) / nn)

        ax.scatter(ave_posArray, [0] * len(ave_posArray), c="white", s=0)
        # ax.scatter(ref_posArray, [-(1 / a)] * len(ref_posArray), c="#c0c0c0", s=100)
        # ax.plot(ref_posArray, [-(1 / a)] * len(ref_posArray), c="blue", zorder=1)
        # for i in range(0, n):
        #     ax.plot([ave_posArray[i]], [0], c="#c0c0c0", zorder=1)
        #     pass
        # ref_text(ax, ref_posArray, xinches, yinches, a, b)

        # ax.set_title('Rsquare plot(' + samplename + ') ' + self.name)
        shape(ax, nn, m, l, ref_posArray, R_matrix, xx_matrix)
        # print(ax.axis())
        # print(ax2.axis())
        ax.axis(ymax=yinches / xinches * (ax.axis()[1] - ax.axis()[0])+yinches)
        ax.axis(ymin=0)
        ax.axis(xmin=ref_posArray[0],xmax=ref_posArray[-1])
        # print(ax.axis())
        ax.axis('off')
        # [ax.spines[loc_axis].set_visible(False) for loc_axis in ['top', 'right', 'bottom', 'left']]
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        plt.xticks([])
        plt.yticks([])
        plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=1)
        ax.invert_yaxis()

        #——————————————————————————————————————————————————————————————

        ax3 = plt.subplot(gs[yinches+int(0.3*yinches):-int(0.3*yinches),-1:])
        norm = mpl.colors.Normalize(vmin=-1, vmax=1)

        plt.subplots_adjust(bottom=0.5)
        cb1 = plt.colorbar(
            mpl.cm.ScalarMappable(norm=norm, cmap="RdBu_r"),
            cax=ax3
        )
        # tick_locator = ticker.MaxNLocator(nbins=3)
        # cb1.locator = tick_locator
        cb1.set_ticks([-1, 0, 1])
        cb1.update_ticks()
        cb1.ax.tick_params(labelsize=20)
        fig.savefig(outpdf)

class IR:

    def __init__(self, iChr, iStart, iEnd, fullText):
        self.iChr = iChr
        self.iStart = int(iStart)
        self.iEnd = int(iEnd)
        self.fullText = fullText

    def toString(self):
        self.name = self.iChr + ":" + str(self.iStart) + "-" + str(self.iEnd)
        return self.name

class HT_mbh:

    def __init__(self, hChr, hStart, hEnd, HapMet, count, strand):
        self.hChr = hChr
        self.hStart = int(hStart)
        self.hEnd = int(hEnd)
        self.HapMet = HapMet
        self.count = int(count)
        self.WC = strand

    def toString(self):
        return self.hChr + '\t' + str(self.hStart) + '\t' + str(self.hEnd) + '\t' + self.HapMet + '\t' + str(self.count) + '\t' + str(self.WC)

class MHB:
    def __init__(self):
        pass

    def loadIR(self, iFile, BED=False):

        irList = []
        f = open(iFile, "r")
        temp = f.read().splitlines()
        f.close()

        for text in temp:
            chunks = text.split("\t")
            if len(chunks) < 3:
                print("Interval not in correct format.")
                sys.exit(-1)
            if BED:
                queryIR = IR(chunks[0], int(chunks[1]) + 1, chunks[2], text)
            else:
                queryIR = IR(chunks[0], chunks[1], chunks[2], text)
            irList.append(queryIR)
        return irList

    def tabixHT(self, htGZ, queryIR, strand='both'):
        tb = tabix.open(htGZ)
        htList = []
        try:
            records = tb.query(queryIR.iChr, queryIR.iStart, queryIR.iEnd)
        except:
            return htList
        for rd in records:
            if len(rd) < 5:
                continue
            if strand == 'plus' and rd[5] != '+':
                continue
            if strand == 'minus' and rd[5] != '-':
                continue
            htList.append(HT_mbh(rd[0], rd[1], rd[2], rd[3], rd[4], rd[5]))
        return htList

    def tabixCPG(self, htGZ, queryIR, shift=500):

        tb = tabix.open(htGZ)
        records = tb.query(queryIR.iChr, queryIR.iStart - shift, queryIR.iEnd + shift)
        posArray = []

        for rd in records:
            if len(rd) < 3:
                continue
            posArray.append(int(rd[1]))

        return cpgAnnotation(queryIR.iChr, posArray)

    def buildBinaryMatrix(self, htList, cpgAnn):
        posArray = cpgAnn.posArray
        posDict = dict()
        Nc = len(posArray)
        Nr = len(htList)

        for i in range(Nc):
            posDict[cpgAnn.iChr + ":" + str(posArray[i])] = i

        MC = np.zeros((len(htList), Nc), dtype='int')
        M0 = np.zeros((len(htList), Nc), dtype='int')
        M1 = np.zeros((len(htList), Nc), dtype='int')

        count = np.zeros(Nr, dtype='int')
        strand = np.zeros(Nr, dtype='int')

        for i in range(Nr):
            ht = htList[i]
            pos = ht.hChr + ":" + str(ht.hStart)
            count[i] = ht.count
            if ht.WC == "+":
                strand[i] = 1
            if pos not in posDict:
                print("Haplotype positions are out of range.")
                sys.exit(-1)
            Idx = posDict[pos]
            HapMet = ht.HapMet
            for j in range(len(HapMet)):
                MC[i, Idx + j] = 1
                if HapMet[j] == "0":
                    M0[i, Idx + j] = 1
                elif HapMet[j] == "1":
                    M1[i, Idx + j] = 1
                else:
                    print("Haplotypes must be 0/1 only.")
                    sys.exit(-1)

        return [MC, M0, M1, count, strand]

    def count2pos(self,MC, M0, M1, count, pi, pj):
        if pi >= np.shape(MC)[1] or pj >= np.shape(MC)[1]:
            print("pi or pj are out of range.")
            sys.exit(-1)

        iAll = MC[..., pi] + MC[..., pj] == 2
        iN00 = M0[..., pi] + M0[..., pj] == 2
        iN01 = M0[..., pi] + M1[..., pj] == 2
        iN10 = M1[..., pi] + M0[..., pj] == 2
        iN11 = M1[..., pi] + M1[..., pj] == 2

        N00 = sum(count[np.logical_and(iAll, iN00)])
        N01 = sum(count[np.logical_and(iAll, iN01)])
        N10 = sum(count[np.logical_and(iAll, iN10)])
        N11 = sum(count[np.logical_and(iAll, iN11)])

        return [N00, N01, N10, N11]

    def get_r2(self, n_11, n_10, n_01, n_00):
        n = n_11 + n_10 + n_01 + n_00

        if n == 0:
            return (None, None)

        P_A = (n_10 + n_11) / n
        P_B = (n_01 + n_11) / n

        # calculate pval
        D = n_11 / n - P_A * P_B
        Num = D * D

        Den = P_A * (1 - P_A) * P_B * (1 - P_B)

        if Den == 0:
            r2 = None
        else:
            r2 = Num / Den

        if D < 0:
            r2 = -1 * r2  # add sign to r2

        # calculate pval
        pval = 1 - binom.cdf(k=n_11 - 1, n=n, p=P_A * P_B)
        return (r2, pval)

    # Get r square from array [n_11, n_10, n_01, n_00]
    def LD(self, array):
        return self.get_r2(array[0], array[1], array[2], array[3])

    def getMHB(self, MC, M0, M1, count, cpgAnn, window=3, r_square=0.5, p_value=0.05):
        Nc = np.shape(MC)[1]
        start = 0
        end = 0
        mhb = []
        while end < (Nc - 1):
            end += 1
            extend = True
            # print("running", end)
            for i in range(1, window):
                xPOS = end - i
                if xPOS < 0:
                    continue
                [N00, N01, N10, N11] = self.count2pos(MC, M0, M1, count, xPOS, end)
                [r2, pval] = self.LD([N11, N10, N01, N00])
                if r2 is None or r2 < r_square or pval > p_value:
                    extend = False
                    break
            if not extend:
                candidate = [start, end - 1]
                size = (end - start)
                start = end
                if size >= window:
                    mhb.append(candidate)
        if end - start >= window:
            mhb.append([start, end - 1])
        irMHB = []
        for blocks in mhb:
            irMHB.append(IR(cpgAnn.iChr, cpgAnn.posArray[blocks[0]], cpgAnn.posArray[blocks[1]], ''))
        return irMHB



def main():

    parser = argparse.ArgumentParser()

    subparsers1 =  parser.add_subparsers()
    gw = subparsers1.add_parser('genomeWide', help='calculate methylation metrics for mHaps that cover each CpG site across the genome')
    gw.add_argument('--tag', type=str, required=True, help='prefix of the output file(s)')
    gw.add_argument('--mhapPath', type=str,required=True, help='input file, mhap.gz format, generated by mHapTools and indexed')
    gw.add_argument('--cpgPath', type=str, required=True, help='genomic CpG file, gz format and indexed')
    gw.add_argument('--metrics', required=True, nargs='*',help='mHap-level metrics, including MM, PDR, CHALM, MHL, MCR, MBS, Entropy, and R2')
    gw.add_argument('--outputDir', type=str, required=True, help='output directory, created in advance')
    gw.add_argument('--minK',  default=np.inf,type=int, help='minimum k-mer length for MHL [1]')
    gw.add_argument('--maxK',default=np.inf,type=int, help='maximum k-mer length for MHL [10]')
    gw.add_argument('--K',  default=np.inf,type=int, help='k-mer length for entropy, PDR, and CHALM, can be 3, 4, or 5 [4]')
    gw.add_argument('--strand', type=str, default='both', help='strand information, one of plus, minus and both [both]')
    gw.set_defaults(func='genomeWide')



    stat = subparsers1.add_parser('stat', help='calculate methylation metrics for mHaps that cover predefined regions')
    stat.add_argument('--metrics', nargs='*',default=None, help='mHap-level metrics, including MM, PDR, CHALM, MHL, MCR, MBS, and Entropy [None]')
    stat.add_argument('--mhapPath', type=str, required=True, help='input file, mhap.gz format, generated by mHapTools and indexed')
    stat.add_argument('--cpgPath', type=str, required=True, help='genomic CpG file, gz format and Indexed')
    stat.add_argument('--region', type=str, default=None, help='one region, in the format of chr:start-end')
    stat.add_argument('--bedPath', type=str, default=None, help='input BED file')
    stat.add_argument('--outputFile', type=str, required=True, help='output file name')
    stat.add_argument('--minK',  default=np.inf,type=int, help='minimum k-mer length for MHL [1]')
    stat.add_argument('--maxK', default=np.inf,type=int, help='maximum k-mer length for MHL [10]')
    stat.add_argument('--K',  default=np.inf,type=int,
                    help='k-mer length for entropy, PDR, and CHALM, can be 3, 4, or 5 [4]')
    stat.add_argument('--strand', type=str, default='both', help='strand information, one of plus, minus and both [both]')
    stat.set_defaults(func='stat')


    tanghulu = subparsers1.add_parser('tanghulu', help='plot the DNA methylation status for mHaps in a region')
    tanghulu.add_argument('--mhapPath', type=str,required=True ,help='input file, mhap.gz format, generated by mHapTools and indexed')
    tanghulu.add_argument('--simulation',action='store_true',help='indicates whether mHaps should be simulated')
    tanghulu.add_argument('--cpgPath', type=str,  required=True, help='genomic CpG file, gz format and Indexed')
    tanghulu.add_argument('--region', type=str, required=True, help='one region, in the format of chr:start-end')
    tanghulu.add_argument('--merge', action='store_true', help='indicates whether identical mHaps should be merged')
    tanghulu.add_argument( "--outcut", type=int, default=2000, help="the max length of region to plot, default is 2000")
    tanghulu.add_argument('--outputFile', type=str,required=True, help='output file name')
    tanghulu.set_defaults(func='tanghulu')

    R2 = subparsers1.add_parser('R2', help='calculate linkage disequilibrium between CpG sites within predefined regions')
    R2.add_argument('--tag', type=str, required=True, help='prefix of the output file(s)')
    R2.add_argument('--mhapPath', type=str, required=True, help='input file, mhap.gz format, sorted by samtools')
    R2.add_argument('--cpgPath', type=str, required=True, help='genomic CpG file, gz format and Indexed')
    R2.add_argument('--region', type=str, required=True, help='one region, in the format of chr:start-end')
    R2.add_argument('--outputDir', type=str, required=True, help='output directory name')
    R2.add_argument('--mHapView',action='store_true',  help='plot linkage disequilibrium patterns of pair-wise CpGs')
    R2.add_argument('--longrange', action='store_true', help='indicates whether generate a file in longrange format')
    R2.add_argument('--strand', type=str, default='both', help='strand information, one of plus, minus and both [both]')
    R2.set_defaults(func='R2')


    MHBDiscovery = subparsers1.add_parser('MHBDiscovery', help='identification of methylation haplotype blocks within a region or genome-wide')
    MHBDiscovery.add_argument('--mhapPath', type=str,required=True, help='input file, mhap.gz format, generated by mHapTools and indexed')
    MHBDiscovery.add_argument('--cpgPath', type=str, required=True, help='genomic CpG file, gz format and Indexed')
    MHBDiscovery.add_argument('--region', type=str, help='one region, in the format of chr:start-end')
    MHBDiscovery.add_argument('--bedPath', type=str,default=None, help='input BED file')
    MHBDiscovery.add_argument('--outputFile', type=str, required=True, help='output file name')
    MHBDiscovery.add_argument("--window", type=int, default=5, required=False, help="size of core window [5]")
    MHBDiscovery.add_argument("--r_square", type=float, default=0.5, required=False, help="R-square cutoff [0.5]")
    MHBDiscovery.add_argument("--p_value", type=float, default=0.05, required=False, help="P-value cutoff [0.05]")
    MHBDiscovery.set_defaults(func='MHBDiscovery')

    args = parser.parse_args()
    try:
        args.func
    except:
        print('mHapTk:A comprehensive tool kit for analysis of DNA methylation haplotypes')
        print('version:', version)
        sys.exit()

    if args.func == 'genomeWide':
        assert args.strand == 'both' or args.strand == 'plus' or args.strand == 'minus', '--stranded should be both plus minus'

        if 'MHL' not in args.metrics and (args.maxK != np.inf or args.minK != np.inf):
            print('Warning: --maxK and --minK is only for mhl')
        if not('PDR' in args.metrics or 'CHALM' in args.metrics or 'Entropy' in args.metrics) and args.K != np.inf:
            print('Warning: --K is only for PDR CHALM Entropy')
        if args.maxK == np.inf:
            args.maxK = 10
        if args.minK == np.inf:
            args.minK = 1
        if args.K == np.inf:
            args.K = 4
        assert isinstance(args.maxK, int), 'maxK should be int'
        assert isinstance(args.K, int), 'K should be int'
        assert isinstance(args.minK, int), 'minK should be int'
        assert Path(args.outputDir).exists(), 'outDir does not exist'
        assert 3<= args.K <=5, 'K：the default is 4 and values must be between 3 and 5'
        assert  args.maxK > args.minK, 'maxK should be larger than minK'
        assert 1 <= args.maxK <= 10, 'maxK should be in 1 to 10'
        assert 1 <= args.minK <= 10, 'minK should be in 1 to 10'
        resultPath = args.outputDir + '/' + args.tag + '_'
        print('the stat u chose is ', args.metrics)
        gw = GenomeWide(args.mhapPath,
                        args.cpgPath,
                        args.maxK,
                        args.minK,
                        args.strand,
                        args.K)
        for stat in args.metrics:
            print(stat)
            if stat not in ['MM', 'CHALM', 'PDR', 'MHL', 'MBS', 'MCR', 'Entropy', 'R2']:
                print('you input a wrong stat')
                print('the right input like', gw.statslist)
            else:
                if stat == 'MM':
                    Time = time.time()
                    MM = gw.MM()
                    print('MM was done')
                    # MM.to_csv(resultPath + 'MM GW.csv', sep='\t', index=False, header=None)
                    MM_new = pd.concat([MM.iloc[:, 0], MM.iloc[:, 1] - 1, MM.iloc[:, 1], MM.iloc[:, -1]], axis=1).round(8)
                    MM_new[(1 - np.isnan(MM_new.iloc[:, 3])).astype(np.bool_)].to_csv(resultPath + 'MM.bedGraph',
                                                                                      index=False, header=None,
                                                                                      sep='\t')
                    print('MM time span:', time.time() - Time)
                if stat == 'CHALM':
                    Time = time.time()
                    CHALM = gw.CHALM()
                    print('CHALM was done')
                    CHALM_new = pd.concat(
                        [CHALM.iloc[:, 0], CHALM.iloc[:, 1] - 1, CHALM.iloc[:, 1], CHALM.iloc[:, -1]], axis=1).round(8)
                    CHALM_new[(1 - np.isnan(CHALM_new.iloc[:, 3])).astype(np.bool)].to_csv(
                        resultPath + 'CHALM.bedGraph', index=False, header=None, sep='\t')
                    print('CHALM time span:', time.time() - Time)
                if stat == 'PDR':
                    Time = time.time()
                    PDR = gw.PDR()
                    print('PDR was done ')
                    # PDR.to_csv(resultPath + 'PDR GW.csv', sep='\t', index=False, header=None)
                    PDR_new = pd.concat([PDR.iloc[:, 0], PDR.iloc[:, 1] - 1, PDR.iloc[:, 1], PDR.iloc[:, -1]],
                                        axis=1).round(8)
                    PDR_new[(1 - np.isnan(PDR_new.iloc[:, 3])).astype(np.bool_)].to_csv(
                        resultPath + 'PDR.bedGraph', index=False, header=None, sep='\t')
                    print('PDR time span:', time.time() - Time)
                if stat == 'MHL':
                    Time = time.time()
                    MHL = gw.MHL()
                    print('MHL was done')
                    # MHL.to_csv(resultPath + 'MHL GW.csv', sep='\t', index=False, header=None)
                    MHL_new = pd.concat([MHL.iloc[:, 0], MHL.iloc[:, 1] - 1, MHL.iloc[:, 1], MHL.iloc[:, -1]],
                                        axis=1).round(8)
                    MHL_new[(1 - np.isnan(MHL_new.iloc[:, 3])).astype(np.bool_)].to_csv(
                        resultPath + 'MHL.bedGraph', index=False, header=None, sep='\t')
                    print('MHL time span:', time.time() - Time)
                if stat == 'MBS':
                    Time = time.time()
                    MBS = gw.MBS()
                    print('MBS was done ')
                    # MBS.to_csv(resultPath + 'MBS GW.csv', sep='\t', index=False, header=None)
                    MBS_new = pd.concat([MBS.iloc[:, 0], MBS.iloc[:, 1] - 1, MBS.iloc[:, 1], MBS.iloc[:, -1]],
                                        axis=1).round(8)
                    MBS_new[(1 - np.isnan(MBS_new.iloc[:, 3])).astype(np.bool_)].to_csv(
                        resultPath + 'MBS.bedGraph', index=False, header=None, sep='\t')
                    print('MBS time span:', time.time() - Time)
                if stat == 'MCR':
                    Time = time.time()
                    MCR = gw.MCR()
                    print('MCR was done')
                    # MCR.to_csv(resultPath + 'MCR GW.csv', sep='\t', index=False, header=None)
                    MCR_new = pd.concat([MCR.iloc[:, 0], MCR.iloc[:, 1] - 1, MCR.iloc[:, 1], MCR.iloc[:, -1]],
                                        axis=1).round(8)
                    MCR_new[(1 - np.isnan(MCR_new.iloc[:, 3])).astype(np.bool_)].to_csv(
                        resultPath + 'MCR.bedGraph', index=False, header=None, sep='\t')
                    print('MCR time span:', time.time() - Time)
                if stat == 'Entropy':
                    Time = time.time()
                    Entropy = gw.Entropy()
                    print('Entropy was done')
                    # Entropy.to_csv(resultPath + 'Entropy GW.csv', sep='\t', index=False, header=None)
                    Entropy_new = pd.concat(
                        [Entropy.iloc[:, 0], Entropy.iloc[:, 1] - 1, Entropy.iloc[:, 1], Entropy.iloc[:, -1]],
                        axis=1).round(8)
                    Entropy_new[(1 - np.isnan(Entropy_new.iloc[:, 3])).astype(np.bool_)].to_csv(
                        resultPath + 'Entropy.bedGraph', index=False, header=None, sep='\t')
                    print('Entropy time span:', time.time() - Time)
                if stat == 'R2':
                    Time = time.time()
                    R2 = gw.R2()
                    print('R2 was done')
                    # R2.to_csv(resultPath + 'R2 GW.csv', sep='\t', index=False, header=None)
                    # R2_new = pd.concat([R2.iloc[:, 0], R2.iloc[:, 1] - 1, R2.iloc[:, 1], R2.iloc[:, -1]], axis=1)
                    # print(R2_new)
                    # print(np.isnan(R2.iloc[:, 3]))
                    R2_new = pd.concat([R2.iloc[:, 0], R2.iloc[:, 1] - 1, R2.iloc[:, 1], R2.iloc[:, -1]], axis=1)
                    R2_new = R2_new.dropna().round(8)
                    R2_new.to_csv(resultPath + 'R2.bedGraph', index=False, header=None, sep='\t')
                    # R2_new[(1 - np.isnan(R2.iloc[:, 3])).astype(np.bool_)].to_csv(
                    #     resultPath + 'R2.bedGraph', index=False, header=None, sep='\t')
                    print('R2 time span:', time.time() - Time)

    if args.func == 'R2':
        assert Path(args.outputDir).exists(), 'outDir does not exist'
        assert args.strand == 'both' or args.strand == 'plus' or args.strand == 'minus', '--stranded should be both plus minus'

        resultPath = args.outputDir + '/' + args.tag + '_'
        M = R2_c(args.mhapPath,
                     args.cpgPath,
                     args.strand)
        M.Region(args.region)
        M.tabixCPG(args.cpgPath, shift=500)

        M.getrecord(strand=args.strand)

        position = M.Chr + "_" + str(M.start + 1) + "_" + str(M.end)
        cpgAnn = M.tabixCPG(args.cpgPath)
        [MC, M0, M1, count, strand_] = M.buildBinaryMatrix()

        samplename = os.path.basename(args.mhapPath)
        outtxt = resultPath + position + ".cpg_sites_rsquare.txt"
        outpdf1 = resultPath + position + ".cpg_sites_rsquare.pdf"
        outpdf2 = resultPath + position + ".cpg_sites_rsquare_hp.pdf"
        if args.longrange:
            outlongrange = resultPath  + M.Chr + '_' + str( M.start) + '_' + str(M.end) + '_' +'longrange'


        sites_pos = np.where(sum(MC) != 0)[0]

        xx_matrix = np.zeros([len(sites_pos), len(sites_pos)], dtype=int)
        pval_matrix = np.zeros([len(sites_pos), len(sites_pos)], dtype=float)
        ref_posArray = [cpgAnn.posArray[int(i)] for i in np.where(sum(MC) != 0)[0]]
        R_matrix = np.zeros([len(sites_pos), len(sites_pos)], dtype=float)
        with open(outtxt, 'w+') as f:
            if args.longrange:
                f2 = open(outlongrange, 'w+')
            f.write("\t".join(['chr', 'posi', 'posj','N00', 'N01', 'N10', 'N11', 'Chr', 'r2', 'pvalue']))
            f.write("\n")
            for i in range(0, len(sites_pos)):
                pi = sites_pos[i]

                if not M.start < ref_posArray[i] <= M.end:
                    continue
                for j in [j for j in range(0, len(sites_pos)) if sites_pos[j] >= sites_pos[i]]:
                    if not M.start < ref_posArray[j] <= M.end:
                        continue
                    pj = sites_pos[j]


                    (N00, N01, N10, N11), (Ni0, Ni1, Nj0, Nj1) = M.count2pos(MC, M0, M1, count, pi, pj)
                    r2, pval = M.D_r2(N11, N10, N01, N00)
                    R_matrix[i][j] = r2
                    pval_matrix[i][j] = pval
                    xx_matrix[i][j] = 1
                    if i != j:
                        f.write(
                            "\t".join([cpgAnn.iChr, str(ref_posArray[i]), str(ref_posArray[j]), str(N00), str(N01), str(N10), str(N11),
                                        str(round(r2,8)), str(pval)]))
                        if args.longrange:
                            f2.write(f'{cpgAnn.iChr}\t{ref_posArray[i]}\t{ref_posArray[i]+1}\t{cpgAnn.iChr}:{ref_posArray[j]}-{ref_posArray[j]+1},{round(r2,8)}')
                            f2.write('\n')
                        f.write("\n")

        rd1 = []
        rd2 = []
        for i in range(R_matrix.shape[0]):
            if sum(R_matrix[i]) != 0:
                rd1.append(i)
            if sum(R_matrix[:, i]) != 0:
                rd2.append(i)

        R_matrix = R_matrix[rd1][:, rd2]
        xx_matrix = xx_matrix[rd1][:, rd2]

        M.paint_rsquare_plot(samplename, R_matrix, xx_matrix, MC, outpdf1)
        M.paint_rsquare_heatmap( R_matrix, MC, outpdf2)

        if args.mHapView:
            if args.cpgPath is None:
                print('lack of cpg path')
            else:
                M.tabixCPG(args.cpgPath)
                hp_cpg = M.hp_CPG(resultPath)
            outpdf_join = resultPath + position + "_join.pdf"
            M.join_pic( R_matrix, xx_matrix, MC, hp_cpg, outpdf_join)

    if args.func == 'stat':
        if (args.metrics is not None) and ('MHL' not in args.metrics)  and (args.maxK != np.inf or args.minK != np.inf):
            print('Warning: --maxK and --minK is only for mhl')
        if (args.metrics is not None) and not('PDR' in args.metrics or 'CHALM' in args.metrics or 'Entropy' in args.metrics) and args.K != np.inf:
            print('Warning: --K is only for PDR CHALM Entropy')
        if args.maxK == np.inf:
            args.maxK = 10
        if args.minK == np.inf:
            args.minK = 1
        if args.K == np.inf:
            args.K = 4
        assert isinstance(args.maxK, int), 'maxK should be int'
        assert isinstance(args.K, int), 'K should be int'
        assert isinstance(args.minK, int), 'minK should be int'
        assert Path(os.path.dirname(args.outputFile)).exists(), 'outDir does not exist'
        assert  args.maxK > args.minK, 'maxK should be larger than minK'
        assert 3 <= args.K <= 5, 'K：the default is 4 and values must be between 3 and 5'
        assert  args.region is  None or args.bedPath is  None, 'U should only inuput bedPath or region'
        assert (args.region is not None) or (args.bedPath is not None), 'U should input bedPath or region'
        assert 1<= args.maxK <=10, 'maxK should be in 1 to 10'
        assert 1 <= args.minK <= 10, 'minK should be in 1 to 10'
        assert args.strand == 'both' or args.strand == 'plus' or args.strand == 'minus', '--strand should be both plus minus'

        resultPath = args.outputFile
        M = Stat(args.mhapPath,
                        args.cpgPath,
                        args.maxK,
                        args.minK,
                        args.strand,
                        args.K)
        if args.bedPath is not None:

            lines = pd.read_csv(args.bedPath, sep='\t', header=None).shape[0]

            f_info = open(resultPath, 'w+')
            f_info.write('chr\tstart\tend\tnReads\tmBase\tcBase\ttBase\tK4plus\tnDR\tnMR')

            if args.metrics:
                if 'MHl' in args.metrics:
                    cpgAnn = M.tabixCPG(args.cpgPath)
                    M.buildBinaryMatrix()
                for i in args.metrics:
                    f_info.write('\t')
                    f_info.write(i)
            f_info.write('\n')

            for i in tqdm(range(lines), desc='Read bed lines'):
                M.getBed(args.bedPath, i)
                M.getrecord(strand=args.strand)
                # --------------------把mhap的数据计算储存---------------------------------------------
                # columns = ['chr', 'start', 'end', 'nReads', 'mBase', 'tBase', 'K4plus', 'nDR', 'nMR']
                M.info_to_file()
                f_info.write(f'{M.Chr}\t{M.start}\t{M.end}\t{M.nReads}\t{M.mBase}\t{M.cBase}\t{M.tBase}\t{M.K4plus}\t{M.nDR}\t{M.nMR}')


                # ---------------------选择计算的统计量--------------------------------------------

                if args.metrics is not None:
                    stats_list = []
                    for stats in args.metrics:
                        stats_list.append(stats)
                    M.aimming(stats_list)
                    if 'MHL' in args.metrics:
                        dic_stat = M.calculating()
                    else:
                        dic_stat = M.calculating()

                    for i, val in enumerate(args.metrics):

                            f_info.write(f'\t{round(dic_stat[val],8)}')

                f_info.write('\n')
            f_info.close()

        elif args.region is not None:
            M.Region(args.region)
            M.getrecord(args.strand)
            M.info_to_file()

            f = open(resultPath, 'w+')
            f.write('chr\tstart\tend\tnReads\tmBase\tcBase\ttBase\tK4plus\tnDR\tnMR')
            if args.metrics:
                stats_list = []
                for stats in args.metrics:
                    stats_list.append(stats)
                M.aimming(stats_list)
                dic_stat = M.calculating()
                print(dic_stat)
                for key in dic_stat:
                    f.write('\t' + key)

            f.write('\n')
            f.write(f'{M.Chr}\t{M.start}\t{M.end}\t{M.nReads}\t{M.mBase}\t{M.cBase}\t{M.tBase}\t{M.K4plus}\t{M.nDR}\t{M.nMR}')
            if args.metrics:
                for i,key in enumerate(dic_stat):
                        f.write(f'\t{round(dic_stat[key],8)}')
            f.write('\n')

            f.close()

        else:
            print('you should input the region you need if you dont know what to input.Please use help.')

    if args.func == 'tanghulu':
        assert Path(os.path.dirname(args.outputFile)).exists(), 'outDir does not exist'
        resultPath = args.outputFile
        M = Tanghulu(args.mhapPath,args.cpgPath)
        M.Region(args.region)
        M.getrecord()
        M.tabixCPG(args.cpgPath, shift=500)

        position = M.Chr + "_" + str(M.start + 1) + "_" + str(M.end)
        if not args.simulation:
            if M.len > args.outcut:
                print("The region is larger than " + str(
                    args.outcut) + ", it's not recommanded to do tanghulu plotting and system will exit right now...")
                time.sleep(0.1)
                sys.exit()
            [MC, M0, M1, count, strand_] = M.buildBinaryMatrix()
            samplename = os.path.basename(args.mhapPath)

            outpdf = resultPath


            M.paint_tanghulu_plot(args, MC, M0, M1, count, strand_, samplename, outpdf)
        else:
            a, _, _, _,ref_posArray = M.simulate()
            mm = M.MM()
            plot_x = []
            plot_y = []
            meth = []

            for i in range(a.shape[0]):
                for j in range(a.shape[1]):
                    plot_x.append(j)
                    plot_y.append(i)
                    meth.append(a[i, j])

            plt.clf()
            fig = plt.gcf()
            plt.figure(dpi=200, figsize=(8, 15))
            xinches = 10
            yinches = 6 / 24 * 40
            ax = plt.axes([0.1, 0.1, .7, .8])
            plt.title(f'Average methylation:{mm}', fontdict={'size': 22})
            for x, y, m in zip(plot_x, plot_y, meth):
                inner = 'w' if m == 0 else 'k'
                ax.scatter(x, y, c=inner, linewidths=1, edgecolors='k', s=800, zorder=2)
            for i in range(a.shape[0]):
                for j in range(a.shape[1] - 1):
                    ax.plot([j, j + 1], [i, i], zorder=1, c='k')
            for i in range(a.shape[1]):
                ax.text(i , -2, ref_posArray[i], rotation=60, fontdict={'size': 13},horizontalalignment='center')
            ax.set_xticks([])
            ax.set_yticks([])
            plt.savefig(resultPath, dpi=200)

    if args.func == "MHBDiscovery":
        assert Path(os.path.dirname(args.outputFile)).exists(), 'outDir does not exist'
        Time = time.time()
        resultPath = args.outputFile
        M = MHB()
        irList = []
        if args.region is not None:
            chunks = re.split(':|-', args.region)
            queryIR = IR(chunks[0], chunks[1], chunks[2], '')
            irList.append(queryIR)
        else:
            irList = M.loadIR(args.bedPath, BED=True)

        OUT = open(resultPath , "w")
        for queryIR in irList:
            htList = M.tabixHT(args.mhapPath, queryIR)
            if len(htList) == 0:
                continue
            cpgAnn = M.tabixCPG(args.cpgPath, queryIR, shift=500)
            [MC, M0, M1, count, strand] = M.buildBinaryMatrix(htList, cpgAnn)
            irHB = M.getMHB(MC, M0, M1, count, cpgAnn, window=args.window, r_square=args.r_square, p_value=args.p_value)
            for ir in irHB:
                outString = '\t'.join([ir.iChr, str(ir.iStart), str(ir.iEnd)])
                # print(outString)
                OUT.write(outString + "\n")
        OUT.close()
        print('MHB time span:', time.time() - Time)







if __name__ == '__main__':

    main()