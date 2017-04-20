# This script can be invoked from the main EVA application
# This script is designed to optimize identification of the most suitable alignments for the analysis of variation
# The script requires alignment logs as produced by bowtie2
# The script will extract information about the total number of pairs/reads used for each alignment as well as the percent of pairs/reads aligned
# The script will then use this data in order to identify the optimal bin of specified size based on density within the bin and the mean number of aligned pairs/reads

import os
import argparse
import sys
import subprocess
import numpy as np
import pandas as pd

def isfloat(string):
    try:
        float(string)
        return True
    except ValueError:
        return False

def locateBestCrams(args):

    def alInt(df,num,inputDir,output,numBins):
        df=df[df['nPairs'].apply(lambda x: x.isdigit())]
        df=df[df['percentAlign'].apply(lambda x: isfloat(x))]
        df["nPairs"]=df['nPairs'].astype(int)
        df["percentAlign"]=df['percentAlign'].astype(float)
        df["pairsAlign"]=((df["percentAlign"]/100.0)*df["nPairs"].astype(float)).astype(int)

        bins = np.histogram(df["pairsAlign"], bins=numBins)[1] # returns the bins automatically calculated
        binIdx = np.digitize(df["pairsAlign"],bins)
        for i in np.unique(binIdx)[::-1]:
            idx = np.where(binIdx == i)[0].tolist()
            subDF = df.iloc[idx]
            if(len(subDF) >= num):
                submit = df[np.logical_or.reduce([df["pairsAlign"] for x in subDF["pairsAlign"].tolist()])]
                break
        outDF=subDF.sort_values(by="percentAlign",ascending=False).reset_index().drop("index",axis=1).iloc[0:num]
        outDF["path"]=outDF.apply(lambda row: os.path.normpath(inputDir+row["tissue"]+"/")+"/"+row["indiv"]+".cram",axis=1)
        outDF.to_csv(output+".top")
        # here is a boxplot of the number of aligned pairs within the 12 cram files chosen for the analysis
        ax=outDF["pairsAlign"].plot(kind="box",title="Distribution of the number of aligned pairs from select alignments",grid=None,xticks=None)
        plt.savefig(output+".top.png")
        print(", ".join(outDF["path"].tolist()))
        del subDF

    def alRange(df,num,inputDir,output,minBound,maxBound):
        df=df[df['nPairs'].apply(lambda x: x.isdigit())]
        df=df[df['percentAlign'].apply(lambda x: isfloat(x))]
        df["nPairs"]=df['nPairs'].astype(int)
        df["percentAlign"]=df['percentAlign'].astype(float)
        df["pairsAlign"]=((df["percentAlign"]/100.0)*df["nPairs"].astype(float)).astype(int)

        outDF=df[(df["pairsAlign"]>minBound)&(df["pairsAlign"]<maxBound)].reset_index().drop("index",axis=1)
        outDF["path"]=outDF.apply(lambda row: os.path.normpath(inputDir+row["tissue"]+"/")+"/"+row["indiv"]+".cram",axis=1)
        outDF.to_csv(output+".top")
        # here is a boxplot of the number of aligned pairs within the 12 cram files chosen for the analysis
        ax=outDF["pairsAlign"].plot(kind="box",title="Distribution of the number of aligned pairs from select alignments",grid=None,xticks=None)
        plt.savefig(output+".top.png")
        print(", ".join(outDF["path"].tolist()))

    def alAuto(df,num,inputDir,output):
        df=df[df['nPairs'].apply(lambda x: x.isdigit())]
        df=df[df['percentAlign'].apply(lambda x: isfloat(x))]
        df["nPairs"]=df['nPairs'].astype(int)
        df["percentAlign"]=df['percentAlign'].astype(float)
        df["pairsAlign"]=((df["percentAlign"]/100.0)*df["nPairs"].astype(float)).astype(int)

        bins = np.histogram(df["pairsAlign"], bins="auto")[1] # returns the bins automatically calculated
        binIdx = np.digitize(df["pairsAlign"],bins)
        for i in np.unique(binIdx)[::-1]:
            idx = np.where(binIdx == i)[0].tolist()
            subDF = df.iloc[idx]
            if(len(subDF) >= num):
                submit = df[np.logical_or.reduce([df["pairsAlign"] for x in subDF["pairsAlign"].tolist()])]
                break
        outDF=subDF.sort_values(by="percentAlign",ascending=False).reset_index().drop("index",axis=1).iloc[0:num]
        outDF["path"]=outDF.apply(lambda row: os.path.normpath(inputDir+row["tissue"]+"/")+"/"+row["indiv"]+".cram",axis=1)
        outDF.to_csv(output+".top")
        # here is a boxplot of the number of aligned pairs within the 12 cram files chosen for the analysis
        ax=outDF["pairsAlign"].plot(kind="box",title="Distribution of the number of aligned pairs from select alignments",grid=None,xticks=None)
        plt.savefig(output+".top.png")
        print(", ".join(outDF["path"].tolist()))
        del subDF

    header = ['tissue',
              'indiv',
              'nPairs',
              'percentAlign']
    dtypeC={'tissue':'object',
            'indiv':'object',
            'nPairs':'int',
            'percentAlign':'float'}
    alignmentLogDF = pd.read_csv(args.out)
    alignmentLogDF.columns = header

    inputDir = args.input+"/"

    try:
        alInt(alignmentLogDF,args.number,inputDir,args.out,int(args.range))
    except:
        if args.range == "auto":
            alAuto(alignmentLogDF,args.number,inputDir,args.out)
        else:
            try:
                bounds = args.range.split(":")
                alRange(alignmentLogDF,args.number,inputDir,args.out,int(bounds[0]),int(bounds[1]))
            except:
                print("Seems the range parameter is specified incorectly")
    del alignmentLogDF

def main(args):
    scriptCMD="./readLog.sh "+os.path.abspath(args.input)+"/ "+args.out
    os.system(scriptCMD) # run the alignment log parsing and extraction script
    locateBestCrams(args)