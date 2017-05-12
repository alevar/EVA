import os
import argparse
import sys
import subprocess
import random
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA
import matplotlib.pyplot as plt
import multiprocessing
import glob
import signal
import pandas as pd
import seaborn as sbn
from ast import literal_eval
from sklearn import preprocessing
from scipy.stats import boxcox
from pandas.tools.plotting import scatter_matrix
from matplotlib import animation
from scipy.stats import t

# lets also calculate pearson r correlation coefficient
def pearson(dataT):
    dataTC = pd.DataFrame(dataT.groupby(by=["ID","sf"])["cov"].mean()).reset_index()
    dataTC["std"] = pd.DataFrame(dataT.groupby(by=["ID","sf"])["cov"].std()).reset_index()["cov"]
    dataTC["count"] = pd.DataFrame(dataT.groupby(by=["ID","sf"])["cov"].count()).reset_index()["cov"]
    dataTC["tpm"] = pd.DataFrame(dataT.groupby(by=["ID","sf"])["tpm"].mean()).reset_index()["tpm"]
    dataTC["tpm"] = dataTC.apply(lambda row: [row["tpm"]],axis=1)
    dataG = pd.DataFrame(dataTC.groupby(by=["cov"])["tpm"].sum()).reset_index()
    dataG["tpmM"] = dataG.apply(lambda row: sum(row[0])/float(len(row[0])),axis=1)
    del dataTC
    return np.corrcoef(dataG["cov"],dataG["tpmM"])[0][1]

# Here we shall calculate whiskers
def calcWhisk(row,data,param):
    if row["sf"] == 1.0:
        return [100,100,[]]
    iqr = row[param+"_q75"] - row[param+"_q25"]
    lowWhisker = float(row[param+"_q25"])-1.5*float(iqr)
    highWhisker = float(row[param+"_q75"])+1.5*float(iqr)

    rDF = data[data["sf"] == row['sf']]

    wiskhi = np.max(rDF[rDF[param+"Q50"]<=highWhisker][param+"Q50"])
    wisklo = np.min(rDF[rDF[param+"Q50"]>=lowWhisker][param+"Q50"])
    extremesHigh = rDF[rDF[param+"Q50"]>wiskhi][param+"Q50"].tolist()
    extremesLow = rDF[rDF[param+"Q50"]<wisklo][param+"Q50"].tolist()

    del rDF
    return [wisklo,wiskhi,extremesLow+extremesHigh]

def calcWhiskSTD(row,data,param):
    if row["sf"] == 1.0:
        return [0,0,[0]]
    iqr = row[param+"_q75"] - row[param+"_q25"]
    lowWhisker = float(row[param+"_q25"])-1.5*float(iqr)
    highWhisker = float(row[param+"_q75"])+1.5*float(iqr)

    rDF = data[data["sf"] == row['sf']]

    wiskhi = np.max(rDF[rDF["tpmSTD"]<=highWhisker]["tpmSTD"])
    wisklo = np.min(rDF[rDF["tpmSTD"]>=lowWhisker]["tpmSTD"])
    extremesHigh = rDF[rDF["tpmSTD"]>wiskhi]["tpmSTD"].tolist()
    extremesLow = rDF[rDF["tpmSTD"]<wisklo]["tpmSTD"].tolist()

    del rDF
    return [wisklo,wiskhi,extremesLow+extremesHigh]

def calcWhiskCV(row,data,param):
    if row["sf"] == 1.0:
        return [0,0,[0]]
    iqr = row[param+"_q75"] - row[param+"_q25"]
    lowWhisker = float(row[param+"_q25"])-1.5*float(iqr)
    highWhisker = float(row[param+"_q75"])+1.5*float(iqr)

    rDF = data[data["sf"] == row['sf']]

    wiskhi = np.max(rDF[rDF["tpmCV"]<=highWhisker]["tpmCV"])
    wisklo = np.min(rDF[rDF["tpmCV"]>=lowWhisker]["tpmCV"])
    extremesHigh = rDF[rDF["tpmCV"]>wiskhi]["tpmCV"].tolist()
    extremesLow = rDF[rDF["tpmCV"]<wisklo]["tpmCV"].tolist()

    del rDF
    return [wisklo,wiskhi,extremesLow+extremesHigh]

def KendalTau(df,dfBASE,topF=1.0,orderTop=True):
    dfBASE = dfBASE.sort_values(by=["tpm","RankSampleID"],ascending=False)
    dfBASE = dfBASE.reset_index().drop("index",1)
    uniqueCOMBINATION = pd.unique(dfBASE["RankSampleID"])
    outDF = pd.DataFrame(uniqueCOMBINATION.reshape(uniqueCOMBINATION.shape[0],1),columns=["RankSampleID"])
    outDF[str(1.0)] = pd.DataFrame(np.array(range(1,uniqueCOMBINATION.shape[0]+1)).reshape(uniqueCOMBINATION.shape[0],1))

    df = df.sort_values(by=["tpm","RankSampleID"],ascending=False)
    df = df.reset_index().drop("index",1)
    uniqueCOMBINATION = pd.unique(df["RankSampleID"])
    df["rank"] = pd.DataFrame(np.array(range(1,uniqueCOMBINATION.shape[0]+1)).reshape(uniqueCOMBINATION.shape[0],1))

    outDF = pd.merge(outDF,df[["RankSampleID","rank"]],on='RankSampleID',how='outer')

    if orderTop:
        top = int(len(outDF)*topF)
        tau = outDF[:top].corr(method="kendall")
    else:
        bottom = len(outDF)-int(len(outDF)*topF)
        tau = outDF[bottom:].corr(method="kendall")

    del dfBASE
    del df
    del uniqueCOMBINATION
    return tau["rank"]["1.0"]

def readStatsSFRange(data1,data2,outDir,gene=False):
    data=data2
    data3=data1
    cols = ['ID','covBase','tpmBase','falseNegative','tpmMEAN','paMEAN','sf']

    dataOff = pd.unique(data3[(data3["sf"] == 1.0) & (data3["tpm"] == 0.0)]["ID"])
    uniqueID1Maintain = pd.unique(data3[~data3["ID"].isin(dataOff)]["ID"])
    del dataOff
    dataN = data3[~data3["ID"].isin(uniqueID1Maintain)]
    data3 = data3[data3["ID"].isin(uniqueID1Maintain)]
    data3["lost"] = data3.apply(lambda row: row["tpm"] == 0.0,axis=1)

    setTrueNeg = set(dataN[(dataN["sf"]==1.0)&(dataN["tpm"]==0.0)]["ID"].unique()) # Lets try counting the number of false Positives
    setTruePos = set(dataN[(dataN["sf"]==1.0)&(dataN["tpm"]!=0.0)]["ID"].unique()) # Lets try counting the number of false Positives
    setFalsePos = list(setTrueNeg.difference(setTruePos)) # Lets try counting the number of falsePositives

    dataSF = pd.DataFrame(data[cols].groupby(['sf']).mean()).reset_index()
    dataSF["falsePositives"] = np.nan
    dataSF["recall"] = np.nan

    dataSF["falseNegatives"] = dataSF.apply(lambda row: len(data3[(data3["sf"] == row["sf"])&(data3["lost"])]),axis=1)
    dataLostAll = pd.DataFrame(data3.groupby(["ID","sf"]).mean()).reset_index()
    dataSF["falseNegativesFull"] = dataSF.apply(lambda row: len(dataLostAll[(dataLostAll["sf"] == row["sf"])&(dataLostAll["lost"] == 1.0)]),axis=1)
    dataSF["NumTranscripts"] = pd.DataFrame(data3.groupby(["sf"],as_index=False)["tpm"].count())["tpm"]
    dataSF["precision"] = (dataSF["NumTranscripts"]-dataSF["falseNegatives"])/dataSF["NumTranscripts"]
    dictBase = pd.Series(data3[data3["sf"]==1.0].tpm.values,index=data3[data3["sf"]==1.0].ID).to_dict()

    dataSF["pa_q25"] = pd.DataFrame(data.groupby(["sf"])["paQ50"].quantile(0.25)).reset_index()["paQ50"]
    dataSF["pa_median"] = pd.DataFrame(data.groupby(["sf"])["paQ50"].quantile(0.50)).reset_index()["paQ50"]
    dataSF["pa_q75"] = pd.DataFrame(data.groupby(["sf"])["paQ50"].quantile(0.75)).reset_index()["paQ50"]
    dataSF["pa_mean"] = pd.DataFrame(data.groupby(["sf"])["paQ50"].mean()).reset_index()["paQ50"]
    dataSF[['pa_whiskLow','pa_whiskHigh','pa_extremes']] = pd.DataFrame([x for x in dataSF.apply(lambda row: calcWhisk(row,data,"pa"),axis=1)])
    dataSF["pa_numExtremes"] = dataSF.apply(lambda row: len(row["pa_extremes"]),axis=1)
    dataSF["pa_weightedNumExtremes"] = dataSF.apply(lambda row: 0 if row["sf"] == 1.0 else abs(np.array(row["pa_extremes"])-row["pa_median"]).mean()*row["pa_numExtremes"],axis=1)
    dataSF["pa_weightedNormalizedNumExtremes"] = dataSF.apply(lambda row: row["pa_weightedNumExtremes"]/row["NumTranscripts"],axis=1)
    dataSF["pa_std"] = pd.DataFrame(data.groupby(["sf"])["paQ50"].std()).reset_index()["paQ50"]
    dataSF["pa_cv"] = dataSF.apply(lambda row: (row["pa_std"]/row['pa_mean'])*100,axis=1)

    # dataSF["tpm_q25"] = pd.DataFrame(data.groupby(["sf"])["tpmMEAN"].quantile(0.25)).reset_index()["tpmMEAN"]
    # dataSF["tpm_median"] = pd.DataFrame(data.groupby(["sf"])["tpmMEAN"].quantile(0.50)).reset_index()["tpmMEAN"]
    # dataSF["tpm_q75"] = pd.DataFrame(data.groupby(["sf"])["tpmMEAN"].quantile(0.75)).reset_index()["tpmMEAN"]
    # dataSF["tpm_mean"] = pd.DataFrame(data.groupby(["sf"])["tpmMEAN"].mean()).reset_index()["tpmMEAN"]
    # dataSF[['tpm_whiskLow','tpm_whiskHigh','tpm_extremes']] = pd.DataFrame([x for x in dataSF.apply(lambda row: calcWhisk(row,data,"tpm"),axis=1)])
    # dataSF["tpm_numExtremes"] = dataSF.apply(lambda row: len(row["tpm_extremes"]),axis=1)
    # dataSF["tpm_weightedNumExtremes"] = dataSF.apply(lambda row: 0 if row["sf"] == 1.0 else abs(np.array(row["tpm_extremes"])-row["tpm_median"]).mean()*row["tpm_numExtremes"],axis=1)
    # dataSF["tpm_weightedNormalizedNumExtremes"] = dataSF.apply(lambda row: row["tpm_weightedNumExtremes"]/row["NumTranscripts"],axis=1)
    # dataSF["tpm_std"] = pd.DataFrame(data.groupby(["sf"])["tpmMEAN"].std()).reset_index()["tpmMEAN"]
    # dataSF["tpm_cv"] = dataSF.apply(lambda row: (row["tpm_std"]/row['tpm_mean'])*100,axis=1)

    dataSF["std_q25"] = pd.DataFrame(data.groupby(["sf"])["tpmSTD"].quantile(0.25)).reset_index()["tpmSTD"]
    dataSF["std_median"] = pd.DataFrame(data.groupby(["sf"])["tpmSTD"].quantile(0.50)).reset_index()["tpmSTD"]
    dataSF["std_q75"] = pd.DataFrame(data.groupby(["sf"])["tpmSTD"].quantile(0.75)).reset_index()["tpmSTD"]
    dataSF["std_mean"] = pd.DataFrame(data.groupby(["sf"])["tpmSTD"].mean()).reset_index()["tpmSTD"]
    dataSF[['std_whiskLow','std_whiskHigh','std_extremes']] = pd.DataFrame([x for x in dataSF.apply(lambda row: calcWhiskSTD(row,data,"std"),axis=1)])

    dataSF["cv_q25"] = pd.DataFrame(data.groupby(["sf"])["tpmCV"].quantile(0.25)).reset_index()["tpmCV"]
    dataSF["cv_median"] = pd.DataFrame(data.groupby(["sf"])["tpmCV"].quantile(0.50)).reset_index()["tpmCV"]
    dataSF["cv_q75"] = pd.DataFrame(data.groupby(["sf"])["tpmCV"].quantile(0.75)).reset_index()["tpmCV"]
    dataSF["cv_mean"] = pd.DataFrame(data.groupby(["sf"])["tpmCV"].mean()).reset_index()["tpmCV"]
    dataSF[['cv_whiskLow','cv_whiskHigh','cv_extremes']] = pd.DataFrame([x for x in dataSF.apply(lambda row: calcWhiskCV(row,data,"cv"),axis=1)])

    data3["RankSampleID"] = data3["ID"]+":"+data3["sample"].astype(str)

    dataSF["tauFull"] = dataSF.apply(lambda row: KendalTau(data3[data3["sf"] == row["sf"]][["RankSampleID","tpm"]],data3[data3["sf"] == 1.0][["RankSampleID","tpm"]],1.0,True),axis=1)
    dataSF["tauTop10"] = dataSF.apply(lambda row: KendalTau(data3[data3["sf"] == row["sf"]][["RankSampleID","tpm"]],data3[data3["sf"] == 1.0][["RankSampleID","tpm"]],0.1,True),axis=1)
    dataSF["tauTop20"] = dataSF.apply(lambda row: KendalTau(data3[data3["sf"] == row["sf"]][["RankSampleID","tpm"]],data3[data3["sf"] == 1.0][["RankSampleID","tpm"]],0.2,True),axis=1)
    dataSF["tauTop50"] = dataSF.apply(lambda row: KendalTau(data3[data3["sf"] == row["sf"]][["RankSampleID","tpm"]],data3[data3["sf"] == 1.0][["RankSampleID","tpm"]],0.5,True),axis=1)
    dataSF["tauBottom10"] = dataSF.apply(lambda row: KendalTau(data3[data3["sf"] == row["sf"]][["RankSampleID","tpm"]],data3[data3["sf"] == 1.0][["RankSampleID","tpm"]],0.1,False),axis=1)
    dataSF["tauBottom20"] = dataSF.apply(lambda row: KendalTau(data3[data3["sf"] == row["sf"]][["RankSampleID","tpm"]],data3[data3["sf"] == 1.0][["RankSampleID","tpm"]],0.2,False),axis=1)
    dataSF["tauBottom50"] = dataSF.apply(lambda row: KendalTau(data3[data3["sf"] == row["sf"]][["RankSampleID","tpm"]],data3[data3["sf"] == 1.0][["RankSampleID","tpm"]],0.5,False),axis=1)
    print(" << Done Ranking and Tau coeeficient calculation")
    if gene==True:
        dataSF[["sf",
                "falsePositives",
                "falseNegatives",
                "falseNegativesFull",
                "NumTranscripts",
                "pa_q25",
                "pa_median",
                "pa_q75",
                "pa_mean",
                "pa_whiskLow",
                "pa_whiskHigh",
                "pa_weightedNumExtremes",
                "pa_weightedNormalizedNumExtremes",
                "pa_std",
                "pa_cv",
                "std_q25",
                "std_median",
                "std_q75",
                "std_mean",
                "std_whiskLow",
                "std_whiskHigh",
                "cv_q25",
                "cv_median",
                "cv_q75",
                "cv_mean",
                "cv_whiskLow",
                "cv_whiskHigh",
                "tauFull",
                "tauTop10",
                "tauTop20",
                "tauTop50",
                "tauBottom10",
                "tauBottom20",
                "tauBottom50",
                "recall",
                "precision"]].to_csv(outDir+"/csv/groupedGeneBySF.csv")
    else:
        dataSF[["sf",
                "falsePositives",
                "falseNegatives",
                "falseNegativesFull",
                "NumTranscripts",
                "pa_q25",
                "pa_median",
                "pa_q75",
                "pa_mean",
                "pa_whiskLow",
                "pa_whiskHigh",
                "pa_weightedNumExtremes",
                "pa_weightedNormalizedNumExtremes",
                "pa_std",
                "pa_cv",
                "std_q25",
                "std_median",
                "std_q75",
                "std_mean",
                "std_whiskLow",
                "std_whiskHigh",
                "tauFull",
                "tauTop10",
                "tauTop20",
                "tauTop50",
                "tauBottom10",
                "tauBottom20",
                "tauBottom50",
                "recall",
                "precision"]].to_csv(outDir+"/csv/groupedTranscriptBySF.csv")

    del dataSF
    # del data
    print("< Done grouping transcripts by downsampling factor")

def readStatsSFFull(data1,data2,outDir,gene=False):
    data=data2
    data3=data1
    cols = ['ID','covBase','tpmBase','falseNegative','tpmMEAN','paMEAN','sf']

    dataOff = pd.unique(data3[(data3["sf"] == 1.0) & (data3["tpm"] == 0.0)]["ID"])
    uniqueID1Maintain = pd.unique(data3[~data3["ID"].isin(dataOff)]["ID"])
    del dataOff
    dataN = data3[~data3["ID"].isin(uniqueID1Maintain)]
    data3 = data3[data3["ID"].isin(uniqueID1Maintain)]
    data3["lost"] = data3.apply(lambda row: row["tpm"] == 0.0,axis=1)

    setTrueNeg = set(dataN[(dataN["sf"]==1.0)&(dataN["tpm"]==0.0)]["ID"].unique()) # Lets try counting the number of false Positives
    setTruePos = set(dataN[(dataN["sf"]==1.0)&(dataN["tpm"]!=0.0)]["ID"].unique()) # Lets try counting the number of false Positives
    setFalsePos = list(setTrueNeg.difference(setTruePos)) # Lets try counting the number of falsePositives
    dataT = dataN[dataN['ID'].isin(setFalsePos)]
    del dataN

    dataSF = pd.DataFrame(data[cols].groupby(['sf']).mean()).reset_index()
    dataSF["falsePositives"] = dataSF.apply(lambda row: len(dataT[(dataT['ID'].isin(setFalsePos))&(dataT["tpm"]!=0.0)&(dataT["sf"]==row["sf"])]),axis=1)
    setTruePos2 = set(data3[(data3["sf"]==1.0)]["ID"].unique()) # Lets try counting the number of true Positives
    dataSF["truePositives"] = dataSF.apply(lambda row: len(data3[(data3['ID'].isin(setTruePos2))&(data3["tpm"]!=0.0)&(data3["sf"]==row["sf"])]),axis=1)
    dataSF["precision"] = dataSF['truePositives']/(dataSF["truePositives"]+dataSF["falsePositives"])

    dataSF["falseNegatives"] = dataSF.apply(lambda row: len(data3[(data3["sf"] == row["sf"])&(data3["lost"])]),axis=1)
    dataSF["recall"] = dataSF['truePositives']/(dataSF["truePositives"]+dataSF["falseNegatives"])

    dataLostAll = pd.DataFrame(data3.groupby(["ID","sf"]).mean()).reset_index()
    dataSF["falseNegativesFull"] = dataSF.apply(lambda row: len(dataLostAll[(dataLostAll["sf"] == row["sf"])&(dataLostAll["lost"] == 1.0)]),axis=1)
    # 3 Calculating the total number of transcripts at a particular coverage point
    dataSF["NumTranscripts"] = pd.DataFrame(data3.groupby(["sf"],as_index=False)["tpm"].count())["tpm"]
    # First we need to express tpm as percentage deviation from the baseTPM
    dictBase = pd.Series(data3[data3["sf"]==1.0].tpm.values,index=data3[data3["sf"]==1.0].ID).to_dict()

    dataSF["pa_q25"] = pd.DataFrame(data.groupby(["sf"])["paQ50"].quantile(0.25)).reset_index()["paQ50"]
    dataSF["pa_median"] = pd.DataFrame(data.groupby(["sf"])["paQ50"].quantile(0.50)).reset_index()["paQ50"]
    dataSF["pa_q75"] = pd.DataFrame(data.groupby(["sf"])["paQ50"].quantile(0.75)).reset_index()["paQ50"]
    dataSF["pa_mean"] = pd.DataFrame(data.groupby(["sf"])["paQ50"].mean()).reset_index()["paQ50"]
    dataSF[['pa_whiskLow','pa_whiskHigh','pa_extremes']] = pd.DataFrame([x for x in dataSF.apply(lambda row: calcWhisk(row,data,"pa"),axis=1)])
    dataSF["pa_numExtremes"] = dataSF.apply(lambda row: len(row["pa_extremes"]),axis=1)
    dataSF["pa_weightedNumExtremes"] = dataSF.apply(lambda row: 0 if row["sf"] == 1.0 else abs(np.array(row["pa_extremes"])-row["pa_median"]).mean()*row["pa_numExtremes"],axis=1)
    dataSF["pa_weightedNormalizedNumExtremes"] = dataSF.apply(lambda row: row["pa_weightedNumExtremes"]/row["NumTranscripts"],axis=1)
    dataSF["pa_std"] = pd.DataFrame(data.groupby(["sf"])["paQ50"].std()).reset_index()["paQ50"]
    dataSF["pa_cv"] = dataSF.apply(lambda row: (row["pa_std"]/row['pa_mean'])*100,axis=1)

    # dataSF["tpm_q25"] = pd.DataFrame(data.groupby(["sf"])["tpmMEAN"].quantile(0.25)).reset_index()["tpmMEAN"]
    # dataSF["tpm_median"] = pd.DataFrame(data.groupby(["sf"])["tpmMEAN"].quantile(0.50)).reset_index()["tpmMEAN"]
    # dataSF["tpm_q75"] = pd.DataFrame(data.groupby(["sf"])["tpmMEAN"].quantile(0.75)).reset_index()["tpmMEAN"]
    # dataSF["tpm_mean"] = pd.DataFrame(data.groupby(["sf"])["tpmMEAN"].mean()).reset_index()["tpmMEAN"]
    # dataSF[['tpm_whiskLow','tpm_whiskHigh','tpm_extremes']] = pd.DataFrame([x for x in dataSF.apply(lambda row: calcWhisk(row,data,"tpm"),axis=1)])
    # dataSF["tpm_numExtremes"] = dataSF.apply(lambda row: len(row["tpm_extremes"]),axis=1)
    # dataSF["tpm_weightedNumExtremes"] = dataSF.apply(lambda row: 0 if row["sf"] == 1.0 else abs(np.array(row["tpm_extremes"])-row["tpm_median"]).mean()*row["tpm_numExtremes"],axis=1)
    # dataSF["tpm_weightedNormalizedNumExtremes"] = dataSF.apply(lambda row: row["tpm_weightedNumExtremes"]/row["NumTranscripts"],axis=1)
    # dataSF["tpm_std"] = pd.DataFrame(data.groupby(["sf"])["tpmMEAN"].std()).reset_index()["tpmMEAN"]
    # dataSF["tpm_cv"] = dataSF.apply(lambda row: (row["tpm_std"]/row['tpm_mean'])*100,axis=1)

    dataSF["std_q25"] = pd.DataFrame(data.groupby(["sf"])["tpmSTD"].quantile(0.25)).reset_index()["tpmSTD"]
    dataSF["std_median"] = pd.DataFrame(data.groupby(["sf"])["tpmSTD"].quantile(0.50)).reset_index()["tpmSTD"]
    dataSF["std_q75"] = pd.DataFrame(data.groupby(["sf"])["tpmSTD"].quantile(0.75)).reset_index()["tpmSTD"]
    dataSF["std_mean"] = pd.DataFrame(data.groupby(["sf"])["tpmSTD"].mean()).reset_index()["tpmSTD"]
    dataSF[['std_whiskLow','std_whiskHigh','std_extremes']] = pd.DataFrame([x for x in dataSF.apply(lambda row: calcWhiskSTD(row,data,"std"),axis=1)])

    dataSF["cv_q25"] = pd.DataFrame(data.groupby(["sf"])["tpmCV"].quantile(0.25)).reset_index()["tpmCV"]
    dataSF["cv_median"] = pd.DataFrame(data.groupby(["sf"])["tpmCV"].quantile(0.50)).reset_index()["tpmCV"]
    dataSF["cv_q75"] = pd.DataFrame(data.groupby(["sf"])["tpmCV"].quantile(0.75)).reset_index()["tpmCV"]
    dataSF["cv_mean"] = pd.DataFrame(data.groupby(["sf"])["tpmCV"].mean()).reset_index()["tpmCV"]
    dataSF[['cv_whiskLow','cv_whiskHigh','cv_extremes']] = pd.DataFrame([x for x in dataSF.apply(lambda row: calcWhiskCV(row,data,"cv"),axis=1)])

    data3["RankSampleID"] = data3["ID"]+":"+data3["sample"].astype(str)

    dataSF["tauFull"] = dataSF.apply(lambda row: KendalTau(data3[data3["sf"] == row["sf"]][["RankSampleID","tpm"]],data3[data3["sf"] == 1.0][["RankSampleID","tpm"]],1.0,True),axis=1)
    dataSF["tauTop10"] = dataSF.apply(lambda row: KendalTau(data3[data3["sf"] == row["sf"]][["RankSampleID","tpm"]],data3[data3["sf"] == 1.0][["RankSampleID","tpm"]],0.1,True),axis=1)
    dataSF["tauTop20"] = dataSF.apply(lambda row: KendalTau(data3[data3["sf"] == row["sf"]][["RankSampleID","tpm"]],data3[data3["sf"] == 1.0][["RankSampleID","tpm"]],0.2,True),axis=1)
    dataSF["tauTop50"] = dataSF.apply(lambda row: KendalTau(data3[data3["sf"] == row["sf"]][["RankSampleID","tpm"]],data3[data3["sf"] == 1.0][["RankSampleID","tpm"]],0.5,True),axis=1)
    dataSF["tauBottom10"] = dataSF.apply(lambda row: KendalTau(data3[data3["sf"] == row["sf"]][["RankSampleID","tpm"]],data3[data3["sf"] == 1.0][["RankSampleID","tpm"]],0.1,False),axis=1)
    dataSF["tauBottom20"] = dataSF.apply(lambda row: KendalTau(data3[data3["sf"] == row["sf"]][["RankSampleID","tpm"]],data3[data3["sf"] == 1.0][["RankSampleID","tpm"]],0.2,False),axis=1)
    dataSF["tauBottom50"] = dataSF.apply(lambda row: KendalTau(data3[data3["sf"] == row["sf"]][["RankSampleID","tpm"]],data3[data3["sf"] == 1.0][["RankSampleID","tpm"]],0.5,False),axis=1)
    print(" << Done Ranking and Tau coeeficient calculation")
    if gene==True:
        dataSF[["sf",
                "falsePositives",
                "falseNegatives",
                "falseNegativesFull",
                "NumTranscripts",
                "pa_q25",
                "pa_median",
                "pa_q75",
                "pa_mean",
                "pa_whiskLow",
                "pa_whiskHigh",
                "pa_weightedNumExtremes",
                "pa_weightedNormalizedNumExtremes",
                "pa_std",
                "pa_cv",
                "std_q25",
                "std_median",
                "std_q75",
                "std_mean",
                "std_whiskLow",
                "std_whiskHigh",
                "cv_q25",
                "cv_median",
                "cv_q75",
                "cv_mean",
                "cv_whiskLow",
                "cv_whiskHigh",
                "tauFull",
                "tauTop10",
                "tauTop20",
                "tauTop50",
                "tauBottom10",
                "tauBottom20",
                "tauBottom50",
                "recall",
                "precision"]].to_csv(outDir+"/csv/groupedGeneBySF.csv")
    else:
        dataSF[["sf",
                "falsePositives",
                "falseNegatives",
                "falseNegativesFull",
                "NumTranscripts",
                "pa_q25",
                "pa_median",
                "pa_q75",
                "pa_mean",
                "pa_whiskLow",
                "pa_whiskHigh",
                "pa_weightedNumExtremes",
                "pa_weightedNormalizedNumExtremes",
                "pa_std",
                "pa_cv",
                "std_q25",
                "std_median",
                "std_q75",
                "std_mean",
                "std_whiskLow",
                "std_whiskHigh",
                "tauFull",
                "tauTop10",
                "tauTop20",
                "tauTop50",
                "tauBottom10",
                "tauBottom20",
                "tauBottom50",
                "recall",
                "precision"]].to_csv(outDir+"/csv/groupedTranscriptBySF.csv")

    del dataSF
    print("< Done grouping transcripts by downsampling factor")

def readStatsID(dataPrime,outDir,gene=False):
    #group by ID
    unique=dataPrime["sf"].unique().tolist()
    frames = []
    for sf in unique:
        df = pd.DataFrame([])
        dfTPM = pd.DataFrame([])
        df[["ID","covBase","tpmBase","sample"]] = dataPrime[(dataPrime["sf"]==1.0)&~(dataPrime["cov"]==0.0)].reset_index().drop("index",1)[["ID","cov","tpm",'sample']]
        df[["ID","covBase","tpmBase","sample","tpmSF"]] = pd.merge(df,dataPrime[dataPrime["sf"]==sf].reset_index().drop("index",1)[["ID","tpm","sample"]],on=["ID","sample"],how="left")
        df["percentAway"] = df["tpmSF"]/df["tpmBase"]*100
        dfG = pd.DataFrame(df.groupby(by=["ID"])[["covBase","tpmBase"]].max()).reset_index()[["ID","covBase","tpmBase"]]
        df['tpmSF'].replace(0,np.nan,inplace=True)
        dfG["hits"] = pd.DataFrame(df.groupby(by='ID')["tpmSF"].count()).reset_index()['tpmSF']
        df['tpmSF'].replace(np.nan,0,inplace=True)
        dfG["count"] = pd.DataFrame(df.groupby(by='ID')["tpmSF"].count()).reset_index()['tpmSF']
        dfG["falseNegative"] = dfG["count"]-dfG["hits"]
        
        dfG["tpmMEAN"] = pd.DataFrame(df.groupby(by=["ID"])["tpmSF"].mean()).reset_index()["tpmSF"]
        dfTPM[["ID","tpmSTD"]] = pd.DataFrame(df.groupby(by=["ID"])["tpmSF"].std()).reset_index()[["ID","tpmSF"]]
        dfG = pd.merge(dfG,dfTPM[["ID","tpmSTD"]],on=["ID"],how='left')
        dfTPM[["ID","tpmQ25"]] = pd.DataFrame(df.groupby(by=["ID"])["tpmSF"].quantile(0.25)).reset_index()[["ID","tpmSF"]]
        dfG = pd.merge(dfG,dfTPM[["ID","tpmQ25"]],on=["ID"],how='left')
        dfTPM[["ID","tpmQ50"]] = pd.DataFrame(df.groupby(by=["ID"])["tpmSF"].quantile(0.50)).reset_index()[["ID","tpmSF"]]
        dfG = pd.merge(dfG,dfTPM[["ID","tpmQ50"]],on=["ID"],how='left')
        dfTPM[["ID","tpmQ75"]] = pd.DataFrame(df.groupby(by=["ID"])["tpmSF"].quantile(0.75)).reset_index()[["ID","tpmSF"]]
        dfG = pd.merge(dfG,dfTPM[["ID","tpmQ75"]],on=["ID"],how='left')
        dfG["tpmCV"] = dfG["tpmSTD"]/dfG["tpmMEAN"]*100
        dfG["tpmNORM"] = (3*(dfG["tpmMEAN"]-dfG["tpmQ50"]))/dfG["tpmSTD"]
        dfG.replace([np.inf, -np.inf], 0,inplace=True)
        dfG["tpmIQR"] = dfG["tpmQ75"]-dfG["tpmQ25"]
        
        dfG["paMEAN"] = pd.DataFrame(df.groupby(by=["ID"])["percentAway"].mean()).reset_index()["percentAway"]
        dfTPM[["ID","paSTD"]] = pd.DataFrame(df.groupby(by=["ID"])["percentAway"].std()).reset_index()[["ID","percentAway"]]
        dfG = pd.merge(dfG,dfTPM[["ID","paSTD"]],on=["ID"],how='left')
        dfTPM[["ID","paQ25"]] = pd.DataFrame(df.groupby(by=["ID"])["percentAway"].quantile(0.25)).reset_index()[["ID","percentAway"]]
        dfG = pd.merge(dfG,dfTPM[["ID","paQ25"]],on=["ID"],how='left')
        dfTPM[["ID","paQ50"]] = pd.DataFrame(df.groupby(by=["ID"])["percentAway"].quantile(0.50)).reset_index()[["ID","percentAway"]]
        dfG = pd.merge(dfG,dfTPM[["ID","paQ50"]],on=["ID"],how='left')
        dfTPM[["ID","paQ75"]] = pd.DataFrame(df.groupby(by=["ID"])["percentAway"].quantile(0.75)).reset_index()[["ID","percentAway"]]
        dfG = pd.merge(dfG,dfTPM[["ID","paQ75"]],on=["ID"],how='left')
        dfG["paCV"] = dfG["paSTD"]/dfG["paMEAN"]*100
        dfG["paNORM"] = (3*(dfG["paMEAN"]-dfG["paQ50"]))/dfG["paSTD"]
        dfG.replace([np.inf, -np.inf], 0,inplace=True)
        dfG["paIQR"] = dfG["paQ75"]-dfG["paQ25"]
        
        dfG = dfG.sort_values(by=["covBase"]).reset_index().drop("index",1)
        dfG.replace("inf",0,inplace=True)
        dfG.replace(np.nan,0,inplace=True)
        dfG["sf"]=sf
        frames.append(dfG[['ID',
                            'covBase',
                            'tpmBase',
                            'falseNegative',
                            'tpmMEAN',
                            'tpmSTD',
                            'tpmQ25',
                            'tpmQ50',
                            'tpmQ75',
                            'tpmCV',
                            'tpmIQR',
                            'tpmNORM',
                            'paMEAN',
                            'paSTD',
                            'paQ25',
                            'paQ50',
                            'paQ75',
                            'paCV',
                            'paIQR',
                            'paNORM',
                            'sf']])
    data=pd.concat(frames)
    data=data.reset_index(drop=True)
    if gene==True:
        data.to_csv(outDir+"/csv/groupedGeneByID.csv")
    else:
        data.to_csv(outDir+"/csv/groupedTranscriptByID.csv")

    return data

def readStatsStudentTest(data,outDir,gene=False):
    data0 = data[~(data['tpm']==0.0)]
    numSamples = len(data['sample'].unique().tolist())
    numSF = len(data['sf'].unique().tolist())
    data05 = pd.DataFrame([])
    data05[["ID","count"]] = pd.DataFrame(data0.groupby("ID")["tpm"].count()).reset_index()[["ID","tpm"]]
    del data0
    data5 = data05[data05["count"]==numSamples*numSF]
    del data05
    dataT = data[data["ID"].isin(data5["ID"].unique().tolist())]
    del data5

    dataT.sort_values(by="tpm",ascending=False,inplace=True)
    dataIDs = dataT[(dataT["sf"]==1.0)&(dataT["sample"]==0)]["ID"].unique().tolist()
    dataT = dataT[dataT["ID"].isin(dataIDs)]
    dataT.reset_index(inplace=True)
    dataT.drop("index",axis=1,inplace=True)
    data1 = pd.DataFrame(dataT[dataT["sf"]==1.0].groupby(["ID"]).mean()).reset_index()
    dataT = dataT[~(dataT["sf"]==1.0)]
    dataT.reset_index(inplace=True)
    dataT.drop("index",axis=1,inplace=True)

    df = pd.DataFrame([])
    df[["ID","sf","tpmMean"]] = pd.DataFrame(dataT.groupby(["ID","sf"]).mean()).reset_index()[["ID","sf","tpm"]]
    df[["ID","sf","tpmMean","covBase","tpmBase"]] = pd.merge(df,data1[["ID","cov","tpm"]],on="ID",how="outer")
    del data1
    df["std"] = pd.DataFrame(dataT.groupby(["ID","sf"])["tpm"].std()).reset_index()["tpm"]
    df["n"] = pd.DataFrame(dataT.groupby(["ID","sf"])["tpm"].count()).reset_index()["tpm"]
    del dataT
    df["df"] = df["n"]-1
    df["isf"] = t.isf(0.01,df["df"])
    df["denominator"] = df["std"]/np.sqrt(df["n"])
    df["score"] = ((df["isf"]*df["denominator"])-df["tpmMean"])*(-1)
    df["scoreRev"] = df['tpmMean']+(df["tpmMean"]-df["score"])
    df["fold"] = df["scoreRev"]/df["tpmMean"]
    # df = df[(df["tpmBase"]>500)&(df["tpmBase"]<1000)]
    df.sort_values(by="tpmBase",inplace=True)
    df.reset_index(inplace=True)
    df.drop("index",axis=1,inplace=True)
    if gene==True:
        df.to_csv(outDir+"/csv/groupedGeneByIDStudent.csv")
    else:
        df.to_csv(outDir+"/csv/groupedTranscriptByIDStudent.csv")
    del df

def main(args):
    if not os.path.exists(os.path.abspath(args.out)):
        os.makedirs(os.path.abspath(args.out))
        os.makedirs(os.path.abspath(args.out)+'/csv')

    if not args.input == None:
        dtypeC={'refID':'object',
                'tissue':'object',
                'sample':'int',
                'sf':'float',
                'cov':'float',
                'tpm':'float',
                'chr':'object',
                'start':'int',
                'end':'int',
                'tID':'object'}
        names =['tissue',
                'chr',
                'refID',
                'tID',
                'start',
                'end',
                'sample',
                'sf',
                'cov',
                'tpm']

        data = pd.read_csv(os.path.abspath(args.input),sep="\t",skiprows=1,names=names,dtype=dtypeC)
        data['tissue']=data['tissue'].str.strip()
        data["ID"] = data["tissue"]+":"+data["chr"]+":"+data['refID']+":"+data["tID"]+":"+data["start"].astype(str)+"-"+data["end"].astype(str)
        full = True
        try:
            iqrC = float(args.coverage)
            # lets try identifying upper outliers in covBase
            q25,q50,q75 = data['cov'].quantile([0.25,0.5,0.75])
            iqr = q75-q25
            thw = q75+iqrC*iqr
            tlw = q25-iqrC*iqr
            ahw = data[data["cov"]<thw]["cov"].max()
            alw = data[data["cov"]>tlw]["cov"].min()
            transcs = data[(data['cov']<ahw)&(data['cov']>alw)&(data["sf"]==1.0)]["ID"].unique()
            data = data[data["ID"].isin(transcs)]
            data.reset_index(inplace=True)
            data = data.drop("index",axis=1)
            full = False
        except:
            if args.coverage == "full":
                full = True
            else:
                try:
                    bounds = args.coverage.split(":")
                    transcripts = data[(data['cov']<float(bounds[1]))&(data['cov']>float(bounds[0]))&(data["sf"]==1.0)]["ID"].unique()
                    data = data[data["ID"].isin(transcripts)]
                    data.reset_index(inplace=True)
                    data = data.drop("index",axis=1)
                    full = False
                except:
                    print("Seems the coverage parameter is specified incorectly: ", sys.exc_info())

        if not args.top == None:
            data.sort_values(by="tpm",ascending=False,inplace=True)
            dataIDs = data[(data["sf"]==1.0)&(data["sample"]==0)].iloc[:args.top]["ID"].unique().tolist()
            data = data[data["ID"].isin(dataIDs)]
            data.sort_values(by="sf",ascending=True,inplace=True)
            data.reset_index(inplace=True)
            data.drop("index",axis=1,inplace=True)

        dataIDBASE = readStatsID(data,os.path.abspath(args.out))

        if full and args.top == None:
            readStatsSFFull(data,dataIDBASE,os.path.abspath(args.out))
        else:
            readStatsSFRange(data,dataIDBASE,os.path.abspath(args.out))

        if args.de == True:
            readStatsStudentTest(data,os.path.abspath(args.out))

    if not args.gene == None: # analyze gene-level expression
        data = pd.read_csv(os.path.abspath(args.gene)).drop(["Unnamed: 0","Strand","Start","End","FPKM"],axis=1)
        columns = ["geneID","geneName","chr","cov","tpm","tissue","sf","sample"]
        data.columns = columns
        data["ID"] = data["tissue"]+":"+data["chr"]+":"+data['geneID']+":"+data["geneName"]
        full = True
        try:
            iqrC = float(args.coverage)
            # lets try identifying upper outliers in covBase
            q25,q50,q75 = data['cov'].quantile([0.25,0.5,0.75])
            iqr = q75-q25
            thw = q75+iqrC*iqr
            tlw = q25-iqrC*iqr
            ahw = data[data["cov"]<thw]["cov"].max()
            alw = data[data["cov"]>tlw]["cov"].min()
            transcs = data[(data['cov']<ahw)&(data['cov']>alw)&(data["sf"]==1.0)]["ID"].unique()
            data = data[data["ID"].isin(transcs)]
            data.reset_index(inplace=True)
            data = data.drop("index",axis=1)
            full = False
        except:
            if args.coverage == "full":
                full = True
            else:
                try:
                    bounds = args.coverage.split(":")
                    transcripts = data[(data['cov']<float(bounds[1]))&(data['cov']>float(bounds[0]))&(data["sf"]==1.0)]["ID"].unique()
                    data = data[data["ID"].isin(transcripts)]
                    data.reset_index(inplace=True)
                    data = data.drop("index",axis=1)
                    full = False
                except:
                    print("Seems the coverage parameter is specified incorectly: ", sys.exc_info())

        if not args.top == None:
            data.sort_values(by="tpm",ascending=False,inplace=True)
            dataIDs = data[(data["sf"]==1.0)&(data["sample"]==0)].iloc[:args.top]["ID"].unique().tolist()
            data = data[data["ID"].isin(dataIDs)]
            data.sort_values(by="sf",ascending=True,inplace=True)
            data.reset_index(inplace=True)
            data.drop("index",axis=1,inplace=True)
        
        dataIDBASE = readStatsID(data,os.path.abspath(args.out),True)

        if full and args.top == None:
            readStatsSFFull(data,dataIDBASE,os.path.abspath(args.out),True)
        else:
            readStatsSFRange(data,dataIDBASE,os.path.abspath(args.out),True)

        if args.de == True:
            readStatsStudentTest(data,os.path.abspath(args.out),True)