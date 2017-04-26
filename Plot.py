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

spotsOriginal = 92884446.8 # eventuall will be suplied by the program

def plotBoxSF(data,outDir,res):
    # First we plotted the median, q25 and q75
    plt.close('all')
    plt.clf()
    plt.figure(figsize=(int(res[0]), int(res[1])))
    plt.title('Change in median, 2nd and 3rd quartiles of transcript expression levels(TPM)')
    plt.ylabel('Deviation of expression estimate from control (% TPM)')
    plt.xlabel("Portion of aligned spots used in assembly. Full alignment contained "+str(spotsOriginal)+" spots")
    plt.plot(data["sf"], data["median"],'k',color='#CC4F1B')
    plt.fill_between(data["sf"], data["q25"], data["q75"],
        alpha=0.5, edgecolor='#CC4F1B', facecolor='#FF9848')
    plt.fill_between(data["sf"], data["whiskLow"], data["whiskHigh"],
        alpha=0.5, edgecolor='#CC4F1B', facecolor='#FF9848')

    caption = "Figure. Plot shows the change in median, 2nd and 3rd quartiles for the deviation of transcript expression levels (%TPM) as a function of the number of aligned reads used in the assembly."
    plt.figtext(.05, .05, caption,wrap=True,fontsize=18)
    plt.savefig(outDir+"/png/boxSF.png")
    plt.close('all')
    plt.clf()

def plotTauSF(data,outDir,res):
    # Next we plotted all kendal tau coefficients
    title = "Portion of aligned spots versus Kendall's Tau ranking of transcripts by expression levels"
    ax=data[["sf","tauFull","tauTop10","tauTop20","tauTop50","tauBottom50","tauBottom20","tauBottom10"]].plot(x="sf",subplots=False, figsize=(int(res[0]), int(res[1])), title=title);
    ax.set_xlabel("Portion of aligned spots used in assembly. Full alignment contained "+str(spotsOriginal)+" spots")
    ax.set_ylabel("Kendall's Tau ranking correlation coefficient")
    numTransc = len(data["ID"].unique())/len(data["tissue"].unique())
    caption = "Figure. Plot shows the change in Kendall's tau ranking correlation coefficient between the ranking of transcripts by expression levels(TPM) using all reads in the alignment and the ranking of same transcripts using a portion of aligned reads. Assembly and expression estimation was performed using stringtie. Original number of aligned reads is "+str(spotsOriginal)+" Number of transcripts is "+str(numTransc)+" etc"
    plt.figtext(.05, .05, caption,wrap=True,fontsize=18)
    plt.savefig(outDir+"/png/tauSF.png")

def plotScattermatrixSF(data,outDir,res):
    # Then we created a scatter matrix of each quantitative measure with density plots on the diagonal
    # The data was visually inspected to determine any signs of linearity
    ax=sbn.pairplot(data[['sf','falsePositives','falseNegatives','falseNegativesFull','std','tauFull','median']])
    caption = "Figure. Pairplot of ... The diagonal shows distibutions of ..."
    plt.figtext(.05, .05, caption,wrap=True,fontsize=18)
    plt.savefig(outDir+"/png/scatterMatrixSF.png")

def plotPCA(data,outDir,res):
    Y = data["sf"]
    X = data[["falsePositives","falseNegatives","falseNegativesFull","median","weightedNormalizedNumExtremes","std","tauFull"]]
    X_std = StandardScaler().fit_transform(X)

    sklearn_pca = sklearnPCA(n_components=2)
    Y_sklearn = sklearn_pca.fit_transform(X_std)*-1

    with plt.style.context('seaborn-whitegrid'):
        plt.figure(figsize=(6, 4))
        plt.scatter(Y.tolist(),Y_sklearn[:,0].tolist())
        plt.xlabel('Sampling Factor')
        plt.ylabel('Principal Component')
        plt.legend(loc='lower center')
        plt.tight_layout()

    xs=Y_sklearn[:,0]
    ys=Y_sklearn[:,1]*-1
    n=sklearn_pca.components_.shape[1]
    labels = ["falsePositives","falseNegatives","falseNegativesFull","median","weightedNormalizedNumExtremes","std","tauFull"]
    plt.figure(figsize=(30,30))
    plt.scatter(xs,ys)
    for i in range(n):
        plt.arrow(0, 0, sklearn_pca.components_[0,i]*-1, sklearn_pca.components_[1,i],color='r',alpha=0.5)
        plt.text(sklearn_pca.components_[0,i]*-1.15, sklearn_pca.components_[1,i] * 1.15, labels[i], color='g', ha='center', va='center',size=26)
    plt.savefig(outDir+"/png/pcaSF.png")

def plotAll(data,outDir,res,iqrCoefficient,gif):
    print(len(data))
    try:
        iqrC = float(iqrCoefficient)
        # lets try identifying upper outliers in covBase
        q25,q50,q75 = data['covBase'].quantile([0.25,0.5,0.75])
        iqr = q75-q25
        thw = q75+iqrC*iqr
        tlw = q25-iqrC*iqr
        ahw = data[data["covBase"]<thw]["covBase"].max()
        alw = data[data["covBase"]>tlw]["covBase"].min()
        data = data[(data['covBase']<ahw)&(data['covBase']>alw)]
        data.reset_index(inplace=True)
        data = data.drop("index",axis=1)
    except:
        if iqrCoefficient == "full":
            pass
        else:
            try:
                bounds = iqrCoefficient.split(":")
                data = data[(data['covBase']<float(bounds[1]))&(data['covBase']>float(bounds[0]))]
                data.reset_index(inplace=True)
                data = data.drop("index",axis=1)
            except:
                print("Seems the coverage parameter is specified incorectly: ", sys.exc_info())

    print(len(data))
    eDF = pd.DataFrame(data["ID"])
    for cl in list(data)[1:-1]:
        eDF[str(cl)] = pd.DataFrame(preprocessing.scale(boxcox(data[cl]+1)[0]))
    eDF["sf"]=data['sf']
    eDF.to_csv(outDir+"/csv/groupedByIDTransformed.csv")

    fig = plt.figure()
    ims = []
    line, = plt.plot([], [], lw=0.25)
    plt.xlim(eDF["covBase"].min(), eDF["covBase"].max())
    plt.ylim(eDF["tpmQ50"].min(), eDF["tpmQ50"].max())
    def update_line(i,sfs):
        line.set_xdata(ims[i][0])
        line.set_ydata(ims[i][1])
        ax = line.get_axes()
        spotsRetained = spotsOriginal*sfs[i]
        ax.set_xlabel("Change in median of transcript expression levels(TPM) at "+str(sfs[i])+" or mean of "+str(int(spotsRetained))+" spots across 12 alignments")
        ax.set_ylabel('Deviation of expression estimate from control (% TPM)')
        return line,

    fig2 = plt.figure()
    ims2 = []
    line1, = plt.plot([], [], lw=0.25)
    plt.xlim(data["covBase"].min(), data["covBase"].max())
    plt.ylim(data["tpmQ50"].min(), data["tpmQ50"].max())
    def update_line2(i,sfs):
        line1.set_xdata(ims2[i][0])
        line1.set_ydata(ims2[i][1])
        ax = line1.get_axes()
        spotsRetained = spotsOriginal*sfs[i]
        ax.set_xlabel("Change in median of transcript expression levels(TPM) at "+str(sfs[i])+" or mean of "+str(int(spotsRetained))+" spots across 12 alignments")
        ax.set_ylabel('Deviation of expression estimate from control (% TPM)')
        return line1,

    unique=data["sf"].unique().tolist()
    for sf in unique:
        print(sf)
        dataTMP = data[data["sf"]==sf]
        print(len(dataTMP))

        plt.close("all")
        ax=sbn.pairplot(dataTMP[['covBase','falseNegative','tpmMEAN','tpmQ50','tpmIQR','tpmSTD']])
        ax.savefig(outDir+"/png/"+str(sf)+"scatterMatrixByID.png")

        plt.close("all")
        ax=sbn.pairplot(eDF[eDF["sf"]==sf][['covBase','falseNegative','tpmMEAN','tpmQ50','tpmIQR','tpmSTD']])
        ax.savefig(outDir+"/png/"+str(sf)+"scatterMatrixByIDTransformed.png")

        plt.close('all')
        plt.clf()
        plt.figure(figsize=(int(res[0]), int(res[1])))
        plt.title('Normalized TPM Deviation')
        plt.plot(dataTMP["covBase"], dataTMP["tpmQ50"],'k',color='#CC4F1B',lw=0.25)
        plt.fill_between(dataTMP["covBase"], dataTMP["tpmQ25"], dataTMP["tpmQ75"],
            alpha=0.5, edgecolor='#CC4F1B', facecolor='#FF9848')
        spotsRetained = spotsOriginal*sf
        caption = "Figure. Change in median of transcript expression levels(TPM) at "+str(sf)+" or mean of "+str(int(spotsRetained))+" spots across 12 alignments"
        plt.figtext(.02, .02, caption,wrap=True)
        plt.savefig(outDir+"/png/"+str(sf)+"boxID.png")
        if gif == True:
            ims2.append((dataTMP["covBase"],dataTMP["tpmQ50"]))

        plt.close('all')
        plt.clf()
        plt.figure(figsize=(int(res[0]), int(res[1])))
        plt.title('Change in median, 2nd and 3rd quartiles of transcript expression levels(TPM)')
        plt.ylabel('Deviation of expression estimate from control (% TPM)')
        plt.xlabel("Portion of aligned spots")
        plt.plot(eDF[eDF["sf"]==sf]["covBase"], eDF[eDF["sf"]==sf]["tpmQ50"],'k',color='#CC4F1B',lw=0.25)
        spotsRetained = spotsOriginal*sf
        caption = "Figure. Change in median of transcript expression levels(TPM) at "+str(sf)+" or mean of "+str(int(spotsRetained))+" spots across 12 alignments. A BoxCoxTransformation was applied to coverage and tpm values for the purpose of normalization"
        plt.figtext(.02, .02, caption,wrap=True)
        plt.savefig(outDir+"/png/"+str(sf)+"boxIDTransformed.png")
        if gif == True:
            ims.append((eDF[eDF["sf"]==sf]["covBase"],eDF[eDF["sf"]==sf]["tpmQ50"]))
            
        # plt.close("all")
        # ax=dataTMP2[['covBase','tpmBase','falseNegative','tpmMEAN','tpmSTD','tpmQ25','tpmQ50','tpmQ75','tpmCV','tpmIQR']].plot(x=dataTMP2["covBase"],subplots=True, layout=(4, 4), figsize=(int(res[0]), int(res[1])), sharex=False)
        # plt.savefig(outDir+"/png/"+str(sf)+"groupedID.png")

        # plt.close("all")
        # ax=eDF[eDF["sf"]==sf][['covBase','falseNegative','tpmMEAN','tpmQ50','tpmSTD','tpmIQR']].plot(x=eDF[eDF["sf"]==sf]["covBase"],subplots=True, layout=(4, 4), figsize=(int(res[0]), int(res[1])), sharex=False)
        # plt.savefig(outDir+"/png/"+str(sf)+"groupedIDTransformed.png")

        del ax
        del dataTMP

    if gif == True:
        ims.reverse()
        im_ani = animation.FuncAnimation(fig, update_line, len(unique), fargs=(list(reversed(unique)),), interval=600, blit=True)
        im_ani.save(outDir+'/png/boxIDTransformed.gif',writer="imagemagick",dpi=200)

        ims2.reverse()
        im_ani2 = animation.FuncAnimation(fig2, update_line2, len(unique),fargs=(list(reversed(unique)),), interval=600, blit=True)
        im_ani2.save(outDir+'/png/boxID.gif',writer="imagemagick",dpi=200)

def main(args):
    if not os.path.exists(os.path.abspath(args.out)):
        os.makedirs(os.path.abspath(args.out))
    if not os.path.exists(os.path.abspath(args.out)+"/csv"):
        os.makedirs(os.path.abspath(args.out)+"/csv")
    if not os.path.exists(os.path.abspath(args.out)+"/png"):
        os.makedirs(os.path.abspath(args.out)+'/png')

    if not args.sf == None:
        headersSF = ["sf",
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
                    "tauBottom50"]

        dataSF = pd.read_csv(os.path.abspath(args.sf)).drop("Unnamed: 0",axis=1)
        dataSF.columns = headersSF
        plotBoxSF(dataSF,os.path.abspath(args.out),args.resolution.split(":"))
        plotTauSF(dataSF,os.path.abspath(args.out),args.resolution.split(":"))
        plotScattermatrixSF(dataSF,os.path.abspath(args.out),args.resolution.split(":"))
        plotPCA(dataSF,os.path.abspath(args.out),args.resolution.split(":"))

    if not args.id == None:
        headersID = ['ID',
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
                    'sf']
        dataID = pd.read_csv(os.path.abspath(args.id)).drop("Unnamed: 0",axis=1)
        dataID.columns = headersID
        plotAll(dataID,os.path.abspath(args.out),args.resolution.split(":"),args.coverage,args.gif)