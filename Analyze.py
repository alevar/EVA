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
def calcWhisk(row,data):
    if row["sf"] == 1.0:
        return [100,100,[]]
    iqr = row["q75"] - row["q25"]
    lowWhisker = float(row["q25"])-1.5*float(iqr)
    highWhisker = float(row["q75"])+1.5*float(iqr)

    rDF = data[data["sf"] == row['sf']]

    wiskhi = np.max(rDF[rDF["percentAway"]<=highWhisker]["percentAway"])
    wisklo = np.min(rDF[rDF["percentAway"]>=lowWhisker]["percentAway"])
    extremesHigh = rDF[rDF["percentAway"]>wiskhi]["percentAway"].tolist()
    extremesLow = rDF[rDF["percentAway"]<wisklo]["percentAway"].tolist()

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

def readStatsSF(path,outDir):
    print("> Begin grouping transcripts by downsampling factor")
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
    data = pd.read_csv(path,sep="\t",skiprows=1,names=['tissue',
                                                        'chr',
                                                        'refID',
                                                        'tID',
                                                        'start',
                                                        'end',
                                                        'sample',
                                                        'sf',
                                                        'cov',
                                                        'tpm'],dtype=dtypeC)

    # Aggregate a column for uniqueID
    # ID column will be used for quantification of results and building the results DF
    data["ID"] = data["tissue"]+":"+data["chr"]+":"+data['refID']+":"+data["tID"]+":"+data["start"].astype(str)+"-"+data["end"].astype(str)

    # First reject all transcripts for which tpm==0.0 at sf==1.0
    dataOff = pd.unique(data[(data["sf"] == 1.0) & (data["tpm"] == 0.0)]["ID"])
    uniqueID1Maintain = pd.unique(data[~data["ID"].isin(dataOff)]["ID"])
    del dataOff
    dataN = data[~data["ID"].isin(uniqueID1Maintain)]
    data = data[data["ID"].isin(uniqueID1Maintain)]
    data["lost"] = data.apply(lambda row: row["tpm"] == 0.0,axis=1)

    # 7 Computes the number of novel mistakes at each coverage point
    setNo = set(dataN[(dataN["sf"]==1.0)&(dataN["tpm"]==0.0)]["ID"].unique()) # Lets try counting the number of novel mistake
    setYes = set(dataN[(dataN["sf"]==1.0)&(dataN["tpm"]!=0.0)]["ID"].unique()) # Lets try counting the number of novel mistake
    setFin = list(setNo.difference(setYes)) # Lets try counting the number of novel mistake
    dataT = dataN[dataN['ID'].isin(setFin)]
    del dataN
    dataSF = pd.DataFrame(dataT.groupby(['sf']).mean()).reset_index()
    print("> Begin counting false positives")
    dataSF["falsePositives"] = dataSF.apply(lambda row: len(dataT[(dataT['ID'].isin(setFin))&(dataT["tpm"]!=0.0)&(dataT["sf"]==row["sf"])]),axis=1)
    print("< Done counting false positives")
    # 2 Here we shall also count the number of losses
    print("> Begin counting false negatives")
    dataSF["falseNegatives"] = dataSF.apply(lambda row: len(data[(data["sf"] == row["sf"])&(data["lost"])]),axis=1)
    print("< Done counting false negatives")
    # 2 here we shall add information about total losses
    dataLostAll = pd.DataFrame(data.groupby(["ID","sf"]).mean()).reset_index()
    dataSF["falseNegativesFull"] = dataSF.apply(lambda row: len(dataLostAll[(dataLostAll["sf"] == row["sf"])&(dataLostAll["lost"] == 1.0)]),axis=1)
    # 3 Calculating the total number of transcripts at a particular coverage point
    dataSF["NumTranscripts"] = pd.DataFrame(data.groupby(["sf"],as_index=False)["tpm"].count())["tpm"]
    # First we need to express tpm as percentage deviation from the baseTPM
    dictBase = pd.Series(data[data["sf"]==1.0].tpm.values,index=data[data["sf"]==1.0].ID).to_dict()
    print("> Begin calculating percent away from the original tpm estimation")
    # dataBASE = pd.DataFrame([])
    # dataBASE[["ID","baseTPM"]] = data[data["sf"]==1.0][["ID","tpm"]]
    # data = pd.merge(data, dataBASE,on='ID',how='outer')
    # data = data.drop_duplicates()
    # del dataBASE
    # data = data.drop("tpm_y",axis=1)
    # data.columns = ["tissue","chr","refID","tID","start","end","sample","sf","cov","tpm","ID","lost","baseTPM","percentAway"]
    # data["percentAway"] = (data["tpm"]/data["baseTPM"])*100
    # data.to_csv(outDir+"/csv/raw.csv")
    data["percentAway"] = data.apply(lambda row: (row["tpm"]/dictBase[row["ID"]])*100 if row["ID"] in dictBase else np.nan,axis=1)
    print("< Done calculating percent away from the original tpm estimation")
    print("> Begin calculating quartiles, mean and whiskers")
    dataSF["q25"] = pd.DataFrame(data.groupby(["sf"])["percentAway"].quantile(0.25)).reset_index()["percentAway"]
    dataSF["median"] = pd.DataFrame(data.groupby(["sf"])["percentAway"].quantile(0.50)).reset_index()["percentAway"]
    dataSF["q75"] = pd.DataFrame(data.groupby(["sf"])["percentAway"].quantile(0.75)).reset_index()["percentAway"]
    # get average tpm
    dataSF["mean"] = pd.DataFrame(data.groupby(["sf"])["percentAway"].mean()).reset_index()["percentAway"]

    dataSF[['whiskLow','whiskHigh','extremes']] = pd.DataFrame([x for x in dataSF.apply(lambda row: calcWhisk(row,data),axis=1)])
    print("< Done calculating quartiles, mean and whiskers")
    # now it's time to count the number of outliers for each sf
    print("> Begin counting, weighing and normalizing counts for outliers")
    dataSF["numExtremes"] = dataSF.apply(lambda row: len(row["extremes"]),axis=1)
    # now we shall weigh the number of extremes by the mean distance from the median point at each sf
    # the reason mean is chosen (and not the median) is because we want to assign weight based on how
    # badly the distribution is kewed and thus mean is a better option since it is biased heavily by the outliers
    dataSF["weightedNumExtremes"] = dataSF.apply(lambda row: 0 if row["sf"] == 1.0 else abs(np.array(row["extremes"])-row["median"]).mean()*row["numExtremes"],axis=1)
    # weighted extremes normalized by number of transcripts in the bin
    dataSF["weightedNormalizedNumExtremes"] = dataSF.apply(lambda row: row["weightedNumExtremes"]/row["NumTranscripts"],axis=1)
    print("< Done counting, weighing and normalizing counts for outliers")
    # add standard deviation and coefficient of variation
    print("> Begin calculating std and cv")
    dataSF["std"] = pd.DataFrame(data.groupby(["sf"])["percentAway"].std()).reset_index()["percentAway"]
    dataSF["cv"] = dataSF.apply(lambda row: (row["std"]/row['mean'])*100,axis=1)
    print("< Done calculating std and cv")
    # STRATEGY 3: Time to include kendal tau ranking correlation results for each sf-base pair
    # create an ID for distinguishing between samples
    print("> Begin Ranking and Tau coeeficient calculation")
    data["RankSampleID"] = data["ID"]+":"+data["sample"].astype(str)

    dataSF["tauFull"] = dataSF.apply(lambda row: KendalTau(data[data["sf"] == row["sf"]][["RankSampleID","tpm"]],data[data["sf"] == 1.0][["RankSampleID","tpm"]],1.0,True),axis=1)
    dataSF["tauTop10"] = dataSF.apply(lambda row: KendalTau(data[data["sf"] == row["sf"]][["RankSampleID","tpm"]],data[data["sf"] == 1.0][["RankSampleID","tpm"]],0.1,True),axis=1)
    dataSF["tauTop20"] = dataSF.apply(lambda row: KendalTau(data[data["sf"] == row["sf"]][["RankSampleID","tpm"]],data[data["sf"] == 1.0][["RankSampleID","tpm"]],0.2,True),axis=1)
    dataSF["tauTop50"] = dataSF.apply(lambda row: KendalTau(data[data["sf"] == row["sf"]][["RankSampleID","tpm"]],data[data["sf"] == 1.0][["RankSampleID","tpm"]],0.5,True),axis=1)
    dataSF["tauBottom10"] = dataSF.apply(lambda row: KendalTau(data[data["sf"] == row["sf"]][["RankSampleID","tpm"]],data[data["sf"] == 1.0][["RankSampleID","tpm"]],0.1,False),axis=1)
    dataSF["tauBottom20"] = dataSF.apply(lambda row: KendalTau(data[data["sf"] == row["sf"]][["RankSampleID","tpm"]],data[data["sf"] == 1.0][["RankSampleID","tpm"]],0.2,False),axis=1)
    dataSF["tauBottom50"] = dataSF.apply(lambda row: KendalTau(data[data["sf"] == row["sf"]][["RankSampleID","tpm"]],data[data["sf"] == 1.0][["RankSampleID","tpm"]],0.5,False),axis=1)
    print("< Done Ranking and Tau coeeficient calculation")
    dataSF[["sf",
            "cov",
            "tpm",
            "falsePositives",
            "falseNegatives",
            "falseNegativesFull",
            "NumTranscripts",
            "q25",
            "median",
            "q75",
            "mean",
            "whiskLow",
            "whiskHigh",
            "weightedNumExtremes",
            "weightedNormalizedNumExtremes",
            "std",
            "cv",
            "tauFull",
            "tauTop10",
            "tauTop20",
            "tauTop50",
            "tauBottom10",
            "tauBottom20",
            "tauBottom50"]].to_csv(outDir+"/csv/groupedBySF.csv")
    
    del dataSF
    # del data
    print("< Done grouping transcripts by downsampling factor")

def readStatsID(path,outDir):
    print("> Begin grouping transcripts by unique ID")
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

    dataPrime = pd.read_csv(path,sep="\t",skiprows=1,names=names,dtype=dtypeC)
    dataPrime["ID"] = dataPrime["tissue"]+":"+dataPrime["chr"]+":"+dataPrime['refID']+":"+dataPrime["tID"]+":"+dataPrime["start"].astype(str)+"-"+dataPrime["end"].astype(str)
    # Aggregate a column for uniqueID
    # ID column will be used for quantification of results and building the results DF
    unique=dataPrime["sf"].unique().tolist()
    frames = []
    for sf in unique[:-1]:
        print(sf)
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
        dfTPM[["ID","tpmSTD"]] = pd.DataFrame(df.groupby(by=["ID"])["percentAway"].std()).reset_index()[["ID","percentAway"]]
        dfG = pd.merge(dfG,dfTPM[["ID","tpmSTD"]],on=["ID"],how='left')
        dfTPM[["ID","tpmQ25"]] = pd.DataFrame(df.groupby(by=["ID"])["percentAway"].quantile(0.25)).reset_index()[["ID","percentAway"]]
        dfG = pd.merge(dfG,dfTPM[["ID","tpmQ25"]],on=["ID"],how='left')
        dfTPM[["ID","tpmQ50"]] = pd.DataFrame(df.groupby(by=["ID"])["percentAway"].quantile(0.50)).reset_index()[["ID","percentAway"]]
        dfG = pd.merge(dfG,dfTPM[["ID","tpmQ50"]],on=["ID"],how='left')
        dfTPM[["ID","tpmQ75"]] = pd.DataFrame(df.groupby(by=["ID"])["percentAway"].quantile(0.75)).reset_index()[["ID","percentAway"]]
        dfG = pd.merge(dfG,dfTPM[["ID","tpmQ75"]],on=["ID"],how='left')
        dfG["tpmCV"] = dfG["tpmSTD"]/dfG["tpmMEAN"]*100
        dfG["tpmIQR"] = dfG["tpmQ75"]-dfG["tpmQ25"]
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
            'sf']])
    data=pd.concat(frames)
    data=data.reset_index(drop=True)
    data.to_csv(outDir+"/csv/groupedByID.csv")
    print("< Done grouping transcripts by unique ID")

def main(args):
    if not os.path.exists(os.path.abspath(args.out)):
        os.makedirs(os.path.abspath(args.out))
        os.makedirs(os.path.abspath(args.out)+'/csv')

    readStatsSF(os.path.abspath(args.input),os.path.abspath(args.out))
    readStatsID(os.path.abspath(args.input),os.path.abspath(args.out))