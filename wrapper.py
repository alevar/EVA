#!/usr/bin/env python

# ./wrapper.py -i ./data/cram/ -r 0.05:1.05:0.05 -a /scratch0/genomes/hg38/annotation/hg38_p8.biotype_flt.cls.gff3 -e /scratch0/genomes/hg38/hg38.fa -t 2 -o /scratch0/avaraby/replica/protocolCov/test -w 2 -f 10

import os
import argparse
import sys
import subprocess
import random
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import multiprocessing
import signal


def signal_handler(signal, frame,forks):
    print('You pressed Ctrl+C!')
    for fork in forks:
        fork.terminate()
    sys.exit(0)
# The following function is meant to parse the data and produce an output for plotting
def parseDir(path,sample,params):
    assemPath=path+"/transcripts.gtf"
    assemData=np.genfromtxt(assemPath,comments="#",delimiter="\t",dtype={'names':['seqname','source','feature','start','end','score','strand','frame','geneID','transcriptID','refID','refGeneID','refGeneName','cov','fpkm','tpm'],'formats':['object','object','object','object','object','object','object','object','object','object','object','object','object','object','object','object']})

    statsOut = open(params[0],'a')
    for referenceID in np.unique(assemData["refGeneName"]):
        for point in assemData[assemData["refGeneName"] == referenceID]:
            outLine = referenceID+"\t"+sample+"\t"+str(params[2])+"\t"+str(params[1])+"\t"+str(point['cov'])+"\t"+str(point['tpm'])+"\n"
            statsOut.write(outLine)
    statsOut.close()

# the followoing function will be used in order to find optimal genes for analysis
def findGenes(path):
    covData = np.genfromtxt("assembly.coverage",dtype='i8')
    covData.sort()
    minCov = np.percentile(covData[covData >= 100],90) # find the 90 percentile of what is higher than a 100 coverage
    if not len(minCov) == 0:
        pass

# Also need to write a method which will take the final .stat file and produce a sorted one
def sortStat():
    pass

# currently analyzes the top 5 by coverage referenced transcripts
def readStat(path,regionNum,percent):
    def plotIndiv(outPath,statData,binCov):
        plt.clf()
        plt.hold(1)

        uniqueFACTORS = np.unique(statData["sf"])
        uniqueREGIONS = np.unique(statData["region"])
        uniqueTISSUES = np.unique(statData["tissue"])

        # second method of selecting the transcripts to use
        baseStatData = statData[statData["sf"] == 1.0]

        percentile25 = np.array([])
        percentile50 = np.array([])
        percentile75 = np.array([])
        percentileMin = np.array([])
        percentileMax = np.array([])

        covTPM = dict.fromkeys(uniqueFACTORS.tolist(),()) # dictionary for storing lists of normalized TPM from tissues corresponding to keys being scale factors
        baseCovRange = [statData[(statData["sf"] == 1.0) & (statData["sample"] == 0)]["cov"].min(),statData[(statData["sf"] == 1.0) & (statData["sample"] == 0)]["cov"].max()]

        # this one first groups by tissue
        for tissue in uniqueTISSUES:
            for region in uniqueREGIONS:
                currentFrame = statData[(statData["region"] == region) & (statData["tissue"] == tissue)]
                baseTPM = currentFrame[currentFrame["sf"] == 1.0]["tpm"][0]
                for coveragePoint in uniqueFACTORS:
                    region_covDF = currentFrame[currentFrame["sf"] == coveragePoint]
                    newAr = (np.array(region_covDF["tpm"].tolist()).astype(float)*100.0)/float(baseTPM)
                    if not len(newAr) == 0:
                        # if len(covTPM[coveragePoint]) > 2:
                            # if max(covTPM[coveragePoint]) < float(np.average(newAr)):
                                # print("The new largest value is: "+str(coveragePoint)+" "+region+" ",float(np.average(newAr)))

                        # this calculation of percent-leveled variation is very crude
                        covTPM[coveragePoint] = covTPM[coveragePoint]+(float(np.average(newAr)),)

        for factor in uniqueFACTORS:
            # percentileMin = np.append(percentileMin,np.array(list(covTPM[factor])).min())
            percentileMin = np.append(percentileMin,np.percentile(np.array(list(covTPM[factor])), 5))
            percentile25 = np.append(percentile25,np.percentile(np.array(list(covTPM[factor])), 25))
            percentile50 = np.append(percentile50,np.percentile(np.array(list(covTPM[factor])), 50))
            percentile75 = np.append(percentile75,np.percentile(np.array(list(covTPM[factor])), 75))
            percentileMax = np.append(percentileMax,np.percentile(np.array(list(covTPM[factor])), 95))
            # percentileMax = np.append(percentileMin,np.array(list(covTPM[factor])).max())

        plt.plot(uniqueFACTORS, percentile50,'k',color='#CC4F1B')
        plt.fill_between(np.array(uniqueFACTORS), percentile25, percentile75,
            alpha=0.5, edgecolor='#CC4F1B', facecolor='#FF9848')
        plt.fill_between(np.array(uniqueFACTORS), percentileMin, percentileMax,
            alpha=0.5, edgecolor='#CC4F1B', facecolor='#FF9848')

        plt.suptitle("Change in TPM deviation from the mean as a function of coverage")
        plt.title("#REG: "+str(len(uniqueREGIONS))+" BASE_COV_RANGE: "+str(baseCovRange[0])+"-"+str(baseCovRange[1]))
        plt.xlabel("Coverage")
        plt.ylabel("TPM percentage")

        plt.savefig(outPath+"/"+str(binCov)+".png")

        # Here we will attempt to calculate Coefficient of variation
        # And use it in the plots

        cvL = np.array([])
        for factor in uniqueFACTORS:
            std = np.std(covTPM[factor])
            mean = np.mean(covTPM[factor])
            cv = (std/mean)*100
            cvL = np.append(cvL,cv)

        plt.clf()

        plt.plot(uniqueFACTORS, cvL,'k',color='#CC4F1B')
        plt.scatter(uniqueFACTORS, cvL)
        plt.title("TPM variation as a function of coverage")
        plt.xlabel("Coverage")
        plt.ylabel("CV")

        plt.savefig(outPath+"/"+str(binCov)+"_CV.png")

    if not os.path.exists(path+"/pngs"):
        os.makedirs(path+"/pngs")

    statPath = path+"/stats.log"
    statData = np.genfromtxt(statPath,skip_header=1, delimiter='\t',dtype={'names':['region','tissue','sample','sf','cov','tpm'],'formats':['object','object','int','float','float','float']})
    uniqueFACTORS = np.unique(statData["sf"])
    uniqueREGIONS = np.unique(statData["region"])
    uniqueTISSUES = np.unique(statData["tissue"])

    setRegions = set(uniqueREGIONS)

    u, indices = np.unique(statData["region"], return_inverse=True)
    freq=u[np.argmax(np.bincount(indices))] # selecting highest frequency transcript

    # begins a loop for each unique region identifier in the 0th column
    xs = []
    ys = []

    # this has to be done for all regions
    uniqueREGIONS = np.unique(statData[(statData["sf"] == 1.0)]["region"])
    #verify that each selected tissue has all coverages represented by scaling factor

    # Removes every region that has incomplete coverage or tissue data and reloads the data
    setRegions = set(uniqueREGIONS)
    for tissue in uniqueTISSUES:
        for factor in uniqueFACTORS:
            setRegions.intersection_update(set(list(np.unique(statData[(statData["tissue"] == tissue) & (statData["sf"] == factor)]["region"]))))
    # slow approach
    # indeces = np.where((statData["region"][:,None]==np.array(list(setRegions))[:,None][:,None]).all(-1))[1]
    # statData = statData[indeces]

    # here we shall attempt to remove all the duplicate regions
    expectedCount = len(uniqueFACTORS)*len(np.unique(statData["sample"]))*len(uniqueTISSUES) # counts the expected number of occurences in the stats file
    unique, counts = np.unique(statData["region"], return_counts=True)
    bins = dict(zip(unique, counts))
    keysExpected = dict((k, v) for k, v in bins.items() if v == expectedCount).keys()
    setRegions.intersection_update(set(keysExpected))
    # use the output in the awk command below
    # Certainly faster
    cmpCMD = ""
    newPath = path+"/statsNew.log"
    for region in setRegions:
        cmpCMD = cmpCMD+"""$1==\""""+region+"""\" || """
    cmpCMD = cmpCMD[:-4]
    awkCMD = """awk -F "\\t" '"""+cmpCMD+" {print}' "+statPath+" > "+newPath
    os.system(awkCMD)

    # Next we hsall break the information down into bins
    # and plot several graphs based on base coverage bins

    # OK. Here it goes. The grand plan:
    #. Below we shall use the above function to create one graph

    statData = np.genfromtxt(newPath,skip_header=1, delimiter='\t',dtype={'names':['region','tissue','sample','sf','cov','tpm'],'formats':['object','object','int','float','float','float']})
    unique, counts = np.unique(statData["region"], return_counts=True)
    bins = dict(zip(unique, counts))
    uniqueFACTORS = np.unique(statData["sf"])
    uniqueREGIONS = np.unique(statData["region"])
    uniqueTISSUES = np.unique(statData["tissue"])

    # second method of selecting the transcripts to use
    baseStatData = statData[statData["sf"] == 1.0]
    baseStatData[::-1].sort(order=["cov"])
    lenBase = int(len(baseStatData)*0.05)
    baseStatData = baseStatData[baseStatData["sample"] == 0][0:lenBase:1]
    bestTissuesByPercent = np.unique(baseStatData["region"])

    # now we need to identify those regions that have sufficient base coverage
    ind = np.argpartition(statData[statData["sf"] == 1.0]["cov"], -6)[-6:] # selects top 5 (multiplied by number of replicates) candidates for analysis based on the base coverage
    bestTissuesByNum = np.unique(statData[statData["sf"] == 1.0][ind]["region"]) # finds out the tissue names of those

    bestRegions = []

    if percent:
        bestRegions = bestTissuesByPercent
    else:
        bestRegions = bestTissuesByNum
    statData = statData[np.logical_or.reduce([statData["region"] == x for x in bestRegions.tolist()])]

    plotIndiv(path+"/pngs",statData,"BASE"+str(len(np.unique(statData["region"]))))

    # Below we shall use the above plotting function to create binned graphs
    # reports the total number of regions; mean base coverage
    # the excerpt bins the coverage data and groups it together into separate clusters
    baseDF = statData[(statData["sf"] == 1.0) & (statData["sample"] == 0)]
    intData = baseDF["cov"].astype(int)
    bins = np.histogram(intData, bins="auto")[1] # returns the bins automatically calculated
    binIdx = np.digitize(intData,bins)
    
    for i in np.unique(binIdx):
        idx = np.where(binIdx == i)[0].tolist()
        subDF = baseDF[idx]
        if(len(subDF) > 2):
            submit = statData[np.logical_or.reduce([statData["region"] == x for x in subDF["region"].tolist()])]
            plotIndiv(path+"/pngs",submit,bins[i-1])

            # def piecewise_linear(x, x0, y0, k1, k2):
            #     return np.piecewise(x, [x < x0], [lambda x:k1*x + y0-k1*x0, lambda x:k2*x + y0-k2*x0])

            # p , e = optimize.curve_fit(piecewise_linear, x, y)
            # xd = np.linspace(0, 15, 100)
            # pl.plot(x, y, "o")
            # pl.plot(xd, piecewise_linear(xd, *p))

    # Below we shall use the above plotting function to create individual graphs of all the regions
    # Reports all the referenceIDs; mean base coverage



    # lastly we will use segmented regression to plot and analyze the data
    # Reports the referenceID of the gene; base coverage; base TPM

# Also need to include a logger so that i can tell what is happening while the application is running
def logger():
    pass

#The following function will be used in order to randomly generate a sample .stat file similar to those outputed by the application
# the funtion will create a sample.stat file in the current directory
def randStat():
    sampleStat = open('sample.stat','w')
    sampleStat.write("@GENE BASE_COV TISSUE SAMPLE FACTOR COV TPM BASE_TPM\n")
    for gene in range(2):
        base_cov = random.randint(1,100)*100
        base_TPM = random.randint(10,100)*20+random.randint(1,99)
        HD = "GENE"+str(gene)+"\t"+str(base_cov)
        for factor in xfrange(0.1,1.0,0.1):
            coverageM = int(base_cov*factor)
            coverage = random.randint(coverageM-int(0.01*coverageM),coverageM+int(0.01*coverageM))
            for tissue in range(50):
                for sample in range(50):
                    TPM = random.randint(int(base_TPM*factor),base_TPM+(base_TPM-int(base_TPM*factor)))
                    LN = HD+"\tTISSUE"+str(tissue)+"\t"+str(sample)+"\t"+str(factor)+"\t"+str(coverage)+"\t"+str(TPM)+"\t"+str(base_TPM)+"\n"
                    sampleStat.write(LN)
    sampleStat.close()

# The following function will be used to verify the necessary software is installed and configured
def verifySoftware():
    pass

# The function below will be used to parse the alignment and identify regions of high coverage which correspond to exons
# The information will then be used to downsample and calculate the statistics
def isolateHighCoverage(path):
    pass

# The function below will be used to calculate the statistics
def calcStats():
    pass

# This function shall calculate the correct number of:
# 1. forks to maintain simultaneously
# 2. threads to allocate to each fork
def CalcNumThreads(requestedTreads):
    cpuinfo = open('/proc/cpuinfo','r')
    forkNum = 0
    for line in cpuinfo:
        try:
            if line.split()[0] == "siblings":
                forkNum = line.split()[2]/2/requestedTreads
                if forkNum == 0:
                    forkNum = 1
        except:
            pass
    cpuinfo.close()
    cpu_count = multiprocessing.cpu_count()
    forkNum = int(cpu_count/(requestedTreads*2))
    return forkNum # needs to be something else

def xfrange(start, stop, step):
    i = 0
    while start + i * step < stop:
        yield start + i * step
        i += 1

# This function will be used in order to verify whether the input information is submitted in bulk or as a single piece and process the information accordingly
def checkDirFile(inStr):
    toAnalyze = []
    rootDir = os.popen("pwd").read()[:-1]

    if os.path.isdir(inStr):
        for d in os.listdir(inStr):
            if not os.path.isdir(d):
                # is directory. compile a list of all files with full paths
                toAnalyze.append(os.path.abspath(inStr+"/"+d))
            else:
                # Needs to raise an exception
                return "ERROR"
        return toAnalyze
    else:
        return [inStr]

def main(argv):
    signal.signal(signal.SIGINT, signal_handler)

    curPath = os.path.dirname(os.path.realpath(__file__))
    outDir = curPath

    threads = 1

    parser = argparse.ArgumentParser(description='''Help Page''')
    parser.add_argument('-i','--input',required=True,type=str,help="path to the CRAM allignment which is to be used in coverage efficacy analysis by downsampling")
    parser.add_argument('-s','--stats',type=str,help="Path to the location of the informational log about the simulation")
    parser.add_argument('-r','--range',required=True,type=str,help="Specify downsampling range and step where 1.0 is full and 0.0 is none. Format: [start]:[stop]:[step]. Stop point is not included")
    parser.add_argument('-a','--annotation',required=True,type=str,help="Provide the path to the reference annotation")
    parser.add_argument('-e','--reference',required=True,type=str,help="Provide path to the referernce sequence in fasta format which will be used in order to transalate CRAM into SAM and feed it into the stringtie for assembly")
    parser.add_argument('-t','--threads',type=int,help="Indicate the maximum number of threads to be used by the pipeline")
    parser.add_argument('-f','--forks',type=int,help="Indicate the maximum numberof forks to maintain at the same time. Each fork will receive either precalculated or predefined by -t number of threads")
    parser.add_argument('-o','--out',type=str,help="Directory where all the output will be saved")
    parser.add_argument('-w','--random',required=True,type=int,help="The number of random sampling repetitions to perform for each genomic region at each scaling factor")
    parser.add_argument('-n','--num',type=int,help="The number of genomic regions to be included in the analysis")

    args=parser.parse_args()
    inputs = args.input
    assert len(args.range.split(":"))==3
    covRange = [float(args.range.split(":")[0]),float(args.range.split(":")[1]),float(args.range.split(":")[2])]

    numReps = args.random
    sequenceRef = args.reference
    annotationRef = args.annotation
    assert os.path.isfile(annotationRef)
    assert os.path.isfile(sequenceRef)

    if args.stats is not  None:
        outLog = args.stats

    if args.threads is not None:
        threads = args.threads
    else:
        threads = 1

    regionNum = 6
    percent = True
    if args.num is not None:
        regionNum = args.num
        percent = False

    forkNum = CalcNumThreads(threads)

    if args.stats is not None:
        open(outLog,'a').close()
    else:
        open(outDir+"/stats.log",'a').close()
    if args.out is not None:
        outDir = args.out
        os.system("mkdir -p "+outDir)

    def child(path):
        baseDirName = path.split("/")[-1].split(".")[:-1][0]
        finDir = outDir+"/"+baseDirName+"/"+baseDirName
        # If this method with the script works - consider writing the config file from here rather than passing parameters
        for scalefactor in xfrange(covRange[0],covRange[1],covRange[2]):
            randSeed = random.sample(range(1,10000),numReps) # Need to think where to place this line
            for rep in range(numReps):
                inFile = path
                scriptCMD = "./rnaseq_al_pipe.sh "+inFile+" "+finDir+" "+str(randSeed[rep]+scalefactor)+" "+str(rep)+" "+sequenceRef+" "+annotationRef
                os.system(scriptCMD)
        os.remove(finDir+".sam")
        os._exit(0)

    def parent(inputs):

        theorizedForkNum = CalcNumThreads(threads) # Calculates the number of forks to make
        if args.forks is not None and args.forks <= theorizedForkNum*2: # safeguards against excessive forking
            forkNum = args.forks
        else:
            forkNum = theorizedForkNum

        inputAls = checkDirFile(inputs)
        childPIDS = [] # list to be populated with the fork PIDS
        # while ( we parse through the list of tissue alignments to analyze )
        for inputFile in inputAls:
            path = inputFile
            if len(childPIDS) >= forkNum:
                    childPIDS[0].join()
                    childPIDS.remove(childPIDS[0])
            else:
                p = multiprocessing.Process(target=child, args=(os.path.abspath(path),))
                childPIDS.append(p)
                p.start()

        while(len(childPIDS) > 0):
            jobs.join()
            childPIDS.remove(job)

    parent(inputs)

    # here we begin parsing the assembled output .gtf
    statFilePath = outDir+"/stats.log"
    statFile=open(statFilePath,'w+')
    statFile.write("@GENE TISSUE SAMPLE FACTOR COV TPM\n")
    statFile.close()

    inputAls = checkDirFile(inputs)
    # [statPath,scalingFactor,repetition,baseCov,baseTPM]
    for inputFile in inputAls:
        path = os.path.abspath(inputFile)
        sample = path.split("/")[-1].split(".")[:-1][0]
        for sf in xfrange(covRange[0],covRange[1],covRange[2]):
            for rep in range(numReps):
                params=[statFilePath,sf,rep]
                endDir = outDir+"/"+sample+"/"+sample+"_F:"+str(sf)+"_R:"+str(rep)
                parseDir(endDir,sample,params)

    # here we shall parse the stat file and produce the graph
    readStat(outDir,regionNum,percent)

if __name__=="__main__":
    main(sys.argv[1:])
