#!/usr/bin/env python

# ./wrapper.py -i ./data/ERR188044_chrX.sam -r 0.1:1.1:0.3 -a /scratch0/genomes/hg38/annotation/hg38_p8.biotype_flt.cls.gff3 -e /scratch0/genomes/hg38/hg38.fa -t 8 -o /scratch0/avaraby/replica/protocolCov/test3 -w 2 -f 2

# another useful command for selecting protein coding genes from the annotation
# awk -F "\t" '$3 == "gene" {print $9}' hg38_p8.biotype_flt.cls.gff3 | awk -F ";" '$7 == "gene_biotype=protein_coding"'

# Another useful command which allows outputting only the reads that were mapped
# samtools view -F 4 -C -t /scratch0/genomes/hg38/hg38.fa SRR821573.cram > SRR821573.MAPPED.cram

# Another command for generating an output of only high-quality (>50) information as a pileup with a reference
# samtools view -h SRR821573.cram | samtools mpileup -C50 -f  /scratch0/genomes/hg38/hg38.fa -s - > pileup

# perhaps the previous 2 commands could be combined with a pipe for a slight increase in performance - likely to be negligible

# In order to make a tangible statistic we could:
# 1. Write stat log individually for each tissue sample (since they will have different starting coverage)
# 2. For each scaling factor (coverage point) calculate the box-plot quartile groups
# The most important point is to somehow normalize the coverage and the TPM information (although TPM should not have to be normalized since it is a normalized value on its own)

# The ideal graph for the analysis will be:
# 1. mean line of the averages of tpm values for the current coverage
# 2. darker-shade colored 2-nd&3-rd quartile groups of the TPM range
# 3. lighter-shade colored 1-st&4-th quartile groups of the TPM range

# Also do not forget about repeting the pipeline for each tissue-gene-etc run multiple times since the downsampling will be random

# Useful command fr extracting transcripts with a referenceID from assembly gtf files
# awk -F "\t|;" '$3 == "transcript" && match($11,/reference_id*/) {print} ERR18XXXXX.gtf'
# Better command. This command will yield everything with a coverage higher than 1000
# awk -F"\t|;" '$3 == "transcript" && match($11,/reference_id*/) {if(int(substr($14,7,length($14-1)) > 1000) print} ERR18XXXXX.gtf'

# Idea. The top n-th percentile genes can be specified as a parameter

# Another parameter to specify could be to tell the application whether we are interested in information
# regarding genes/transcripts. Basically, what feature are we looking for in the gtf file?

import os
import argparse
import sys
import subprocess
import random
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

# The following function is meant to parse the data and produce an output for plotting
def parseDir(path,sample,params):
    assemPath=path+"/transcripts.gtf"
    assemData=np.genfromtxt(assemPath,comments="#",delimiter="\t",dtype={'names':['seqname','source','feature','start','end','score','strand','frame','geneID','transcriptID','refID','refGeneID','refGeneName','cov','fpkm','tpm'],'formats':['object','object','object','object','object','object','object','object','object','object','object','object','object','object','object','object']})

    statsOut = open(params[0],'a')
    for referenceID in np.unique(assemData["refGeneName"]):
        for point in assemData[assemData["refGeneName"] == referenceID]:
            outLine = referenceID+"\t"+str(params[3])+"\t"+sample+"\t"+str(params[2])+"\t"+str(params[1])+"\t"+str(point['cov'])+"\t"+str(point['tpm'])+"\t"+str(params[4])+"\n"
            statsOut.write(outLine)
    statsOut.close()

# the followoing function will be used in order to find optimal genes for analysis
def findGenes(path):
    covData = np.genfromtxt("assembly.coverage",dtype='i8')
    covData.sort()
    minCov = np.percentile(covData[covData >= 100],90) # find the 90 percentile of what is higher than a 100 coverage
    if not len(minCov) == 0:

        # use this number in awk command in order to identify suitable transcripts/regions and their positions
        # awk -F"\t|;" '$3 == "transcript" && match($11,/reference_id*/) {if(int(substr($14,7,length($14-1)) > 1000) print} ERR18XXXXX.gtf'
        pass

# Also need to write a method which will take the final .stat file and produce a sorted one
def sortStat():
    pass

def readStat(path):
    statPath=path+"/stats.log"
    statData = np.genfromtxt(path,skip_header=1, delimiter='\t',dtype={'names':['region','baseCOV','tissue','sample','sf','cov','tpm','baseTPM'],'formats':['object','i8','object','i8','f8','f8','f8','i8']})
    uniqueFACTORS = np.unique(statData["sf"])
    uniqueTISSUES = np.unique(statData["tissue"])

    u, indices = np.unique(statData["tissue"], return_inverse=True)
    freq=u[np.argmax(np.bincount(indices))] # selecting highest frequency transcript

    plt.clf()
    plt.hold(1)

    # begins a loop for each unique region identifier in the 0th column
    xs = []
    ys = []

    # this has to be done for all regions
    uniqueTISSUES = np.unique(statData[statData["sf"] == 1.0]["tissue"])
    #verify that each selected tissue has all coverages represented by scaling factor

    # Removes everything that has incomplete coverage data
    # highly inefficient - might need to think of a better way of doing this
    tissueBASE_TPM = {}
    for tissue in uniqueTISSUES:
        for factor in uniqueFACTORS:
            if len(statData[(statData["sf"] == factor) & (statData["tissue"] == tissue)]) == 0:
                statData = statData[statData["tissue"]!=tissue]
    
    # now we need to identify those regions that have sufficient base coverage
    ind = np.argpartition(statData[statData["sf"] == 1.0]["cov"], -6)[-6:] # selects top 5 (multiplied by number of replicates) candidates for analysis based on the base tpm
    bestTissues = np.unique(statData[statData["sf"] == 1.0][ind]["tissue"]) # finds out the tissue names of those
    print(bestTissues)
    # print(statData[(statData["tissue"] == bestTissues[0]) & (statData["sf"] == 1.0)])
    
    percentile25 = np.array([])
    percentile50 = np.array([])
    percentile75 = np.array([])
    percentileMin = np.array([])
    percentileMax = np.array([])

    covTPM = dict.fromkeys(uniqueFACTORS.tolist(),()) # dictionary for storing lists of normalized TPM from tissues corresponding to keys being scale factors
    print(covTPM)
    # this one first groups by tissue
    for tissue in bestTissues:
        currentFrame = statData[(statData["tissue"] == tissue)]
        baseTPM = currentFrame[currentFrame["sf"] == 1.0]["tpm"]
        for coveragePoint in uniqueFACTORS:
            region_covDF = currentFrame[currentFrame["sf"] == coveragePoint]
            newAr = (np.array(region_covDF["tpm"].tolist())*100)/float(baseTPM[0])
            if not len(newAr) == 0:
                covTPM[coveragePoint] = covTPM[coveragePoint]+(float(np.average(newAr)),)

    for factor in uniqueFACTORS:
        percentileMin = np.append(percentileMin,np.array(list(covTPM[factor])).min())
        percentile25 = np.append(percentile25,np.percentile(np.array(list(covTPM[factor])), 25))
        percentile50 = np.append(percentile50,np.percentile(np.array(list(covTPM[factor])), 50))
        percentile75 = np.append(percentile75,np.percentile(np.array(list(covTPM[factor])), 75))
        percentileMax = np.append(percentileMax,np.array(list(covTPM[factor])).max())

    print(uniqueFACTORS)
    print(percentile25)
    print(percentile50)
    print(percentile75)
    print(percentileMax)
    print(percentileMin)
    
    plt.plot(uniqueFACTORS, percentile50,'k',color='#CC4F1B')
    plt.fill_between(np.array(uniqueFACTORS), percentile25, percentile75,
        alpha=0.5, edgecolor='#CC4F1B', facecolor='#FF9848')
    plt.fill_between(np.array(uniqueFACTORS), percentileMin, percentileMax,
        alpha=0.5, edgecolor='#CC4F1B', facecolor='#FF9848')

    plt.title("Change in TPM deviation from the mean as a function of coverage")
    plt.xlabel("Coverage")
    plt.ylabel("TPM percentage")

    plt.savefig(path+"covTMP.png")

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

# The followingi tab delimited format will be used to save tpm information
# @GENE BASE_COV TISSUE SAMPLE FACTOR COV TPM BASE_TPM
# GENE0   7600    TISSUE0 0   0.1 753 205 788
# GENE0   7600    TISSUE0 1   0.1 753 663 788
# GENE0   7600    TISSUE0 2   0.1 753 1360    788
# GENE0   7600    TISSUE0 3   0.1 753 1227    788

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

    curPath = os.path.dirname(os.path.realpath(__file__))
    outDir = curPath

    threads = 1

    # Arguments to consider adding in future:
    # 1. Directory with all input alignments
    # 2. Directory with all input reference annotations
    # 3. Directory with all input reference sequences for samtools to be able to view the CRAM files correctly

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

    # Not yet implemented
    parser.add_argument('-n','-num',type=int,help="The number of genomic regions to be included in the analysis")

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

    forkNum = CalcNumThreads(threads)

    if args.stats is not None:
        open(outLog,'a').close()
        print("Your stats log is created in"+outLog)
    else:
        open(outDir+"/stats.log",'a').close()
        print("The default log was created in"+outDir+"/stats.log")

    if args.out is not None:
        outDir = args.out
        os.system("mkdir -p "+outDir)

    # For MARCC need to add parallelization within this wrapper:
    # this level of paralellization will allow running the wrapper simultaneously on multiple tissue samples (read as input files).
    # perhaps do it simply by system forking the process - is very likely to be the most efficient way of doing it

    def child(path):
    # Child should receive its own pseudo-terminal for proper threading of the application
        # Add parsing of the filename
        baseDirName = path.split("/")[-1].split(".")[:-1][0]

        # if(not os.path.exists(outDir+"/downsamp/"+baseDirName)):
        #     os.makedirs(outDir+"/downsamp/"+baseDirName)
        # for i in xfrange(covRange[0],covRange[1],covRange[2]):
        #     if(os.path.exists(outDir+"/downsamp/"+baseDirName+"/"+baseDirName+str(i)+".cram")):
        #         os.system("rm -r "+outDir+"/downsamp/"+baseDirName+"/"+baseDirName+str(i)+".cram")

#!

        # if(not os.path.exists(outDir+"/downsamp/"+baseDirName)):
        #     os.makedirs(outDir+"/downsamp/"+baseDirName)
        # for i in xfrange(covRange[0],covRange[1],covRange[2]):
        #     if(os.path.exists(outDir+"/downsamp/"+baseDirName+"/"+baseDirName+str(i))):
        #         os.system("rm -r "+outDir+"/downsamp/"+baseDirName+"/"+baseDirName+str(i))
        #     os.makedirs(outDir+"/downsamp/"+baseDirName+"/"+baseDirName+str(i))

        # if(not os.path.exists(outDir+"/assembly/"+baseDirName)):
        #     os.makedirs(outDir+"/assembly/"+baseDirName)
        # for i in xfrange(covRange[0],covRange[1],covRange[2]):
        #     if(os.path.exists(outDir+"/assembly/"+baseDirName+"/"+baseDirName+str(i))):
        #         os.system("rm -r "+outDir+"/assembly/"+baseDirName+"/"+baseDirName+str(i))
        #     os.makedirs(outDir+"/assembly/"+baseDirName+"/"+baseDirName+str(i))

        # if(not os.path.exists(outDir+"/statsAl/"+baseDirName)):
        #     os.makedirs(outDir+"/statsAl/"+baseDirName)
        # for i in xfrange(covRange[0],covRange[1],covRange[2]):
        #     if(os.path.exists(outDir+"/statsAl/"+baseDirName+"/"+baseDirName+str(i))):
        #         os.system("rm -r "+outDir+"/statsAl/"+baseDirName+"/"+baseDirName+str(i))
        #     os.makedirs(outDir+"/statsAl/"+baseDirName+"/"+baseDirName+str(i))

#!

        # Perhaps also add a step to break the alignment into functional pieces based on the reference, so that separate directories/files may be created for different genes. Analysis will be much faster and perhaps easier

        # also do not forget to create a separate directory for a base assembly
        # perhaps just save it in the root of each directory where all the subdirs for downsampled assemblies are

        #==============================================
        #first run the base assembly as described above
        #==============================================

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        # print("==================================")
        # print(sequenceRef)
        # print("==================================")

        #assembleFullCov = "samtools view -h -T "+os.path.abspath(sequenceRef)+" "+path+" | stringtie -p "+str(threads)+" -m 150 -G "+os.path.abspath(annotationRef)+" -o "+outDir+"/assembly/"+baseDirName+"/"+baseDirName+".gtf -A "+outDir+"/assembly/"+baseDirName+"/Genes.tab -"
        #os.system(assembleFullCov)

        # Now lets write the first downsampling method using samtools
        #for i in xfrange(covRange[0],covRange[1],covRange[2]):
        #    downSampleCmd = "samtools view -h -s "+str(i)+" -C "+path+" > "+outDir+"/downsamp/"+baseDirName+"/"+baseDirName+str(i)+".cram"
        #    os.system(downSampleCmd)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# The snippet below is for current testing only
# the final version will be slightly more elaborate
        gtf = "/scratch0/RNAseq_protocol/chrX_data/genes/chrX.gtf"

        # If this method with the script works - consider writing the config file from here rather than passing parameters

        for scalefactor in xfrange(covRange[0],covRange[1],covRange[2]):
            for rep in range(numReps):
                inFile = path
                finDir = outDir+"/"+baseDirName+"/"+baseDirName
                scriptCMD = "./rnaseq_al_pipe.sh "+inFile+" "+finDir+" "+str(scalefactor)+" "+str(rep)

                os.system(scriptCMD)

        # for scalefactor in xfrange(covRange[0],covRange[1],covRange[2]):
        #     for rep in range(numReps):

        #         finSampDir = outDir+"/downsamp/"+baseDirName+"/"+baseDirName+str(scalefactor)+"/"
        #         finAssemDir = outDir+"/assembly/"+baseDirName+"/"+baseDirName+str(scalefactor)+"/"
        #         finStatDir = outDir+"/statsAl/"+baseDirName+"/"+baseDirName+str(scalefactor)+"/"

        #         outSampName = finSampDir+baseDirName+"_F:"+str(scalefactor)+"_R:"+str(rep)
        #         outAssemName = finAssemDir+baseDirName+"_F:"+str(scalefactor)+"_R:"+str(rep)

        #         downSampleCmd = "samtools view -h -s "+str(scalefactor)+" "+path+" > "+outSampName+".sam"
        #         sortCommand = "samtools view -S -b "+outSampName+".sam | samtools sort -@ "+str(threads)+" -o "+outSampName+".bam -"
        #         assemCommand = "stringtie -p "+str(threads)+" -G "+gtf+" -o "+outAssemName+".gtf -l "+outAssemName+" "+outSampName+".bam"

        #         mergeCommand0 = "ls -l "+outAssemName+".gtf > "+finStatDir+"mergelist_F:"+str(scalefactor)+"_R:"+str(rep)+".txt"
        #         mergeCommand1 = "stringtie --merge -p "+str(threads)+" -G "+gtf+" -o "+finStatDir+"stringtie_merged_F:"+str(scalefactor)+"_R:"+str(rep)+".gtf "+finStatDir+"mergelist_F:"+str(scalefactor)+"_R:"+str(rep)+".txt"
        #         mergeCommand2 = "stringtie -e -B -p "+str(threads)+" -G "+finStatDir+"stringtie_merged_F:"+str(scalefactor)+"_R:"+str(rep)+".gtf -o "+outAssemName+".FIN.gtf "+outSampName+".bam"

        #         print("BEGIN DOWNSAMPLING REP: "+str(rep))
        #         os.system(downSampleCmd)
        #         print("BEGIN SORTING")
        #         print(sortCommand)
        #         os.system(sortCommand)
        #         print("BEGIN ASSEM")
        #         os.system(assemCommand)
        #         os.system(mergeCommand0)
        #         os.system(mergeCommand1)
        #         os.system(mergeCommand2)

        # for each child Fork create a temporary log in its location so that information can be saved sequentially
        # At the end of all the child processes unify all the temp logs into a single comprehensive under the parent thread

        # Exit after the process is over

        # for scalefactor in xfrange(covRange[0],covRange[1],covRange[2]):
        #     for rep in range(numReps):
        #         finSampDir = outDir+"/downsamp/"+baseDirName+"/"+baseDirName+str(scalefactor)+"/"
        #         finAssemDir = outDir+"/assembly/"+baseDirName+"/"+baseDirName+str(scalefactor)+"/"
        #         finStatDir = outDir+"/statsAl/"+baseDirName+"/"+baseDirName+str(scalefactor)+"/"

        os._exit(0)

    def parent(inputs):
        # the block below is responsible for creating the base directories
#?
#        if(not os.path.exists(outDir+"/downsamp")):
#            os.makedirs(outDir+"/downsamp/")
#        if(not os.path.exists(outDir+"/assembly")):
#            os.makedirs(outDir+"/assembly/")
#        if(not os.path.exists(outDir+"/statsAl")):
#            os.makedirs(outDir+"/statsAl/")

        theorizedForkNum = CalcNumThreads(threads) # Calculates the number of forks to make
        if args.forks is not None and args.forks <= theorizedForkNum*2: # safeguards against excessive forking
            forkNum = args.forks
        else:
            forkNum = theorizedForkNum

        inputAls = checkDirFile(inputs)
        childPIDS = [] # list to be populated with the fork PIDS
        # while ( we parse through the list of tissue alignments to analyze )
        for inputFile in inputAls:
            newpid = os.fork()
            path = inputFile
            if newpid == 0:
                child(os.path.abspath(path))
            else:
                childPIDS.append((os.getpid(),newpid))
                print("PARENT: NEW CHILD CREATED: ",childPIDS)
                if len(childPIDS) >= forkNum:
                    expiredPID = os.wait()
                    print("EXPIRED: ",expiredPID)
                    childPIDS.remove((os.getpid(),expiredPID[0]))

            while(len(childPIDS) > 0):
                expiredPID = os.wait()
                print("FINAL EXPIRATION: ",expiredPID)
                childPIDS.remove((os.getpid(),expiredPID[0]))

    parent(inputs)

    # here we begin parsing the assembled output .gtf

    statFilePath = outDir+"/stats.log"
    statFile=open(statFilePath,'w+')
    statFile.write("@GENE BASE_COV TISSUE SAMPLE FACTOR COV TPM BASE_TPM\n")
    statFile.close()

    # [statPath,scalingFactor,repetition,baseCov,baseTPM]
    for inputFile in inputAls:
        path = os.path.abspath(inputFile)
        sample = path.split("/")[-1].split(".")[:-1][0]
        for sf in xfrange(covRange[0],covRange[1],covRange[2]):
            for rep in range(numReps):
                params=[statFilePath,sf,rep,1000,1000]
                endDir = path+"/"+sample+"/"+sample+"_F:"+str(sf)+"_R:"+str(rep)
                parseDir(endDir,sample,params)

    # here we shall parse the stat file and produce the graph
    readStat(outDir)

if __name__=="__main__":
    main(sys.argv[1:])
