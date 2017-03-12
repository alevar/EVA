#!/usr/bin/env python

import os
import subprocess
import numpy as np
import random

def xfrange(start, stop, step):
    i = 0
    while start + i * step < stop:
        yield start + i * step
        i += 1

def parseDir(path,sample,params):
#parameters are given as a list
# [statPath,scalingFactor,repetition,baseCov,baseTPM]
    assemPath=path+"/transcripts.gtf"
    assemData=np.genfromtxt(assemPath,comments="#",delimiter="\t",dtype={'names':['seqname','source','feature','start','end','score','strand','frame','geneID','transcriptID','refID','refGeneID','refGeneName','cov','fpkm','tpm'],'formats':['object','object','object','object','object','object','object','object','object','object','object','object','object','object','object','object']})

    #assemData = np.hstack((assemData["f0"], np.char.split(str(assemData["f8"]),';'), assemData["f8"]))
    #assemData=np.char.split(str(assemData["f8"]), ';')
    # and this is where we shall find all the stats and then output them into the correct file
    #print(len(assemData))
    print(assemData["cov"])
    #print(len(np.unique(assemData["f10"])))

    statsOut = open(params[0],'a')
    for referenceID in np.unique(assemData["refGeneName"]):
        #outLine = referenceID+"\t"+params[3]+"\t"+sample+"\t"+params[2]+"\t"+params[1]+"\t"+"COVERAGE"+"\t"+"TPM"+"\t"+params[4]+"\n"
        #statsOut.write(outLine)
        for point in assemData[assemData["refGeneName"] == referenceID]:
            outLine = referenceID+"\t"+str(params[3])+"\t"+sample+"\t"+str(params[2])+"\t"+str(params[1])+"\t"+str(point['cov'])+"\t"+str(point['tpm'])+"\t"+str(params[4])+"\n"
            statsOut.write(outLine)
    statsOut.close()

def readStat(path):
    statData = np.genfromtxt(path,skip_header=1, delimiter='\t',dtype='S5,i8,S7,i8,f8,i8,i8,i8')

    uniqueREGIONS = np.unique(statData["f0"])
    uniqueFACTORS = np.unique(statData["f4"])
    uniqueSAMPLES = np.unique(statData["f3"])
    uniqueTISSUES = np.unique(statData["f2"])
    uniqueCOVERAG = np.unique(statData["f5"])

    # could group values of coverage by closest

    # Just for this experiment we shall do the grouping as follows:
    # 1. Group by gene
    # 2. group by coverage
    # 3. Collect quartile information about the TPM for each of these point

    # Also add ability to report errors on the graph:
    # Report them as dots on the graph - the closer to the middle line - the better

    #plt.clf()
    #plt.hold(1)


    # perhaps what we could also do is to plot the graph of not the actual values and their quartiles and means
    # but rather use percentages and deviations from the original value
    # That way we will begin at a zero point


    # this line selects all that match
    # print(statData[(statData["f0"] == "GENE0") & (statData["f4"] == 0.1)])


    # begins a loop for each unique region identifier in the 0th column
    xs = []
    ys = []
    for region in uniqueREGIONS:
        xs = []
        ys = []
        baseTPM = 0
        percentile25 = np.array([])
        percentile50 = np.array([])
        percentile75 = np.array([])
        percentileMin = np.array([])
        percentileMax = np.array([])

        # Now just need to convert the data into percentages

        regionDF = statData[statData["f0"] == region]
        for coveragePoint in uniqueCOVERAG:
            region_covDF = regionDF[regionDF["f5"] == coveragePoint]
            newAr = np.array(region_covDF["f6"].tolist())
            # baseTPM = region_covDF["f6"].tolist()
            if not len(newAr) == 0:
                percentileMin = np.append(percentileMin,newAr.min())
                percentile25 = np.append(percentile25,np.percentile(newAr, 25))
                percentile50 = np.append(percentile50,np.percentile(newAr, 50))
                percentile75 = np.append(percentile75,np.percentile(newAr, 75))
                percentileMax = np.append(percentileMax,newAr.max())
                baseTPM = np.unique(region_covDF["f7"])
                #TPM_quartiles = (newAr.min(),np.percentile(newAr, 25),np.percentile(newAr, 50),np.percentile(newAr, 75),newAr.max())
                xs.append(coveragePoint)

    percentileMin = (percentileMin*100)/baseTPM
    percentile25 = (percentile25*100)/baseTPM
    percentile50 = (percentile50*100)/baseTPM
    percentile75 = (percentile75*100)/baseTPM
    percentileMax = (percentileMax*100)/baseTPM

    #plt.plot(np.array(xs), percentile50,'k',color='#CC4F1B')
    #plt.fill_between(np.array(xs), percentile25, percentile75,
    #    alpha=0.5, edgecolor='#CC4F1B', facecolor='#FF9848')
    #plt.fill_between(np.array(xs), percentileMin, percentileMax,
    #    alpha=0.5, edgecolor='#CC4F1B', facecolor='#FF9848')

    #plt.title("Change in TPM deviation from the mean as a function of coverage")
    #plt.xlabel("Coverage")
    #plt.ylabel("TPM percentage")

    #plt.show()

def main():
    path = "/scratch0/avaraby/replica/protocolCov/test3"
    statFilePath = path+"/stats.log"
    statFile=open(statFilePath,'w+')
    statFile.write("@GENE BASE_COV TISSUE SAMPLE FACTOR COV TPM BASE_TPM\n")
    statFile.close()
    samples = ["ERR188044_chrX"]
    reps = 2
    # [statPath,scalingFactor,repetition,baseCov,baseTPM]
    for sample in samples:
        for sf in xfrange(0.1,1.1,0.3):
            for rep in range(reps):
                params=[statFilePath,sf,rep,1000,1000]
                endDir = path+"/"+sample+"/"+sample+"_F:"+str(sf)+"_R:"+str(rep)
                parseDir(endDir,sample,params)

if __name__=="__main__":
    main()
