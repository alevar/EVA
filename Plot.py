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
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA as sklearnPCA
mpl.rcParams['agg.path.chunksize'] = 10000

spotsOriginal = 92884447 # Mean number of spots in the selected alignments. eventually will be suplied by the program

def plotBoxSF(data,outDir,res):
    # First we plotted the median and quartiles
    plt.close('all')
    plt.clf()
    fig = plt.figure(figsize=(int(res[0]),int(res[1])))
    ax1 = fig.add_axes((0.1,0.25,0.8,0.7))
    ax1.set_title('Change in variation of transcript expression (%TPM) estimation as a function of the portion of aligned spots')
    ax1.set_ylabel('Deviation of expression estimate from control (% TPM)')
    ax1.set_xlabel("Portion of aligned spots")
    ax1.plot(data["sf"], data["median"],'k',color='#CC4F1B')
    ax1.set_xlim(data["sf"].min(),data["sf"].max())
    ax1.set_xticks(data["sf"])
    ax1.fill_between(data["sf"], data["q25"], data["q75"],
        alpha=0.5, edgecolor='#CC4F1B', facecolor='#FF9848')
    ax1.fill_between(data["sf"], data["whiskLow"], data["whiskHigh"],
        alpha=0.5, edgecolor='#CC4F1B', facecolor='#FF9848')

    ax2 = fig.add_axes((0.1,0.2,0.8,0.0))
    ax2.set_xlim(data["sf"].min(),data["sf"].max())
    ax2.yaxis.set_visible(False)
    ax2.set_xticks(data["sf"])
    ax2.set_xlabel("# aligned spots")
    ticks2F = [str('%.2f' %((elem*spotsOriginal)/1000000))+"M" for elem in data["sf"].tolist()]
    ax2.set_xticklabels(ticks2F)

    caption = "Figure. Plot shows the change in variation of estimated transcript expression levels (%TPM) as a function of the number of aligned reads used in the assembly. The plot shows the median and four quartiles of the distribution of estimated expression levels (%TPM) for each downsampling."
    # plt.figtext(.05, .05, caption,wrap=True,fontsize=18)
    plt.savefig(outDir+"/png/boxSF.png")
    plt.close('all')
    plt.clf()

def plotTauSF(data,outDir,res):
    # Next we plotted all kendal tau coefficients
    plt.close('all')
    plt.clf()
    fig = plt.figure(figsize=(int(res[0]),int(res[1])))
    ax1 = fig.add_axes((0.1,0.25,0.8,0.7))
    title = "Portion of aligned spots versus Kendall's Tau ranking of transcripts by expression levels"
    ax1.set_title(title)
    ax1.set_ylabel("Kendall's Tau ranking correlation coefficient")
    ax1.set_xlabel("Portion of aligned spots")
    for tauSeries in data[["tauFull","tauTop10","tauTop20","tauTop50","tauBottom50","tauBottom20","tauBottom10"]]:
        ax1.plot(data["sf"], data[tauSeries],label=tauSeries)
    ax1.legend()
    ax1.set_xlim(data["sf"].min(),data["sf"].max())
    ax1.set_xticks(data["sf"])

    ax2 = fig.add_axes((0.1,0.2,0.8,0.0))
    ax2.set_xlim(data["sf"].min(),data["sf"].max())
    ax2.yaxis.set_visible(False)
    ax2.set_xticks(data["sf"])
    ax2.set_xlabel("# aligned spots")
    ticks2F = [str('%.2f' %((elem*spotsOriginal)/1000000))+"M" for elem in data["sf"].tolist()]
    ax2.set_xticklabels(ticks2F)
    
    numTransc = int(data[data["sf"]==1.0]["NumTranscripts"]/120)
    caption = "Figure. Plot shows the change in Kendall's tau ranking correlation coefficient between the ranking of transcripts by expression levels(TPM) using all reads in the alignment and the ranking of same transcripts using a portion of aligned reads. Assembly and expression estimation was performed using stringtie. Original number of aligned reads is "+str(spotsOriginal)+". Number of transcripts is "+str(numTransc)+"."
    # plt.figtext(.05, .05, caption,wrap=True,fontsize=18)
    plt.savefig(outDir+"/png/tauSF.png")

def plotScattermatrixSFFull(data,outDir,res):
    # Then we created a scatter matrix of each quantitative measure with density plots on the diagonal
    # The data was visually inspected to determine any signs of linearity
    ax=sbn.pairplot(data[['sf','falsePositives','falseNegatives','falseNegativesFull','std','tauFull','median']])
    caption = "Figure. Pairplot of ... The diagonal shows distibutions of ..."
    # plt.figtext(.05, .05, caption,wrap=True,fontsize=18)
    plt.savefig(outDir+"/png/scatterMatrixSF.png")

def plotScattermatrixSFRange(data,outDir,res):
    # Then we created a scatter matrix of each quantitative measure with density plots on the diagonal
    # The data was visually inspected to determine any signs of linearity
    ax=sbn.pairplot(data[['sf','falseNegatives','falseNegativesFull','std','tauFull','median']])
    caption = "Figure. Pairplot of ... The diagonal shows distibutions of ..."
    # plt.figtext(.05, .05, caption,wrap=True,fontsize=18)
    plt.savefig(outDir+"/png/scatterMatrixSF.png")

def plotPCAFull(data,outDir,res):
    plt.close("all")
    plt.clf()
    Y = data["sf"]
    X = data[["falsePositives","falseNegatives","falseNegativesFull","median","weightedNormalizedNumExtremes","std","tauFull"]]
    X_std = StandardScaler().fit_transform(X)

    sklearn_pca = sklearnPCA(n_components=2)
    Y_sklearn = sklearn_pca.fit_transform(X_std)*-1

    xs=Y_sklearn[:,0]
    ys=Y_sklearn[:,1]*-1
    labels = ["falsePositives","falseNegatives","falseNegativesFull","median","weightedNormalizedNumExtremes","std","tauFull"]
    fig = plt.figure(figsize=(int(res[0]),int(res[1])))
    ax1 = fig.add_axes((0.1,0.1,0.8,0.85))
    ax1.set_title("PCA1 versus PCA2")
    ax1.set_ylabel("PCA2")
    ax1.set_xlabel("PCA1")
    ax1.scatter(xs,ys)

    n=sklearn_pca.components_.shape[1]
    for i in range(n):
        ax1.arrow(0, 0, sklearn_pca.components_[0,i]*-1, sklearn_pca.components_[1,i],color='r',alpha=1.0)
        ax1.annotate(labels[i],(sklearn_pca.components_[0,i]*-1.15, sklearn_pca.components_[1,i] * 1.15))
    
    ticks2F = [str('%.2f' %((elem*spotsOriginal)/1000000))+"M" for elem in data["sf"].tolist()]
    for i, txt in enumerate(data["sf"].tolist()):
        ax1.annotate(str(txt)+" ("+ticks2F[i]+")", (xs[i],ys[i]))


    plt.savefig(outDir+"/png/pcaSF.png")

def plotPCARange(data,outDir,res):
    plt.close("all")
    plt.clf()
    Y = data["sf"]
    X = data[["falseNegatives","falseNegativesFull","median","weightedNormalizedNumExtremes","std","tauFull"]]
    X_std = StandardScaler().fit_transform(X)

    sklearn_pca = sklearnPCA(n_components=2)
    Y_sklearn = sklearn_pca.fit_transform(X_std)*-1

    xs=Y_sklearn[:,0]
    ys=Y_sklearn[:,1]*-1
    labels = ["falseNegatives","falseNegativesFull","median","weightedNormalizedNumExtremes","std","tauFull"]
    fig = plt.figure(figsize=(int(res[0]),int(res[1])))
    ax1 = fig.add_axes((0.1,0.1,0.8,0.85))
    ax1.set_title("PCA1 versus PCA2")
    ax1.set_ylabel("PCA2")
    ax1.set_xlabel("PCA1")
    ax1.scatter(xs,ys)

    n=sklearn_pca.components_.shape[1]
    for i in range(n):
        ax1.arrow(0, 0, sklearn_pca.components_[0,i]*-1, sklearn_pca.components_[1,i],color='r',alpha=1.0)
        ax1.annotate(labels[i],(sklearn_pca.components_[0,i]*-1.15, sklearn_pca.components_[1,i] * 1.15))
    
    ticks2F = [str('%.2f' %((elem*spotsOriginal)/1000000))+"M" for elem in data["sf"].tolist()]
    for i, txt in enumerate(data["sf"].tolist()):
        ax1.annotate(str(txt)+" ("+ticks2F[i]+")", (xs[i],ys[i]))


    plt.savefig(outDir+"/png/pcaSF.png")

def plotPrecision(data,outDir,res):
    plt.close('all')
    plt.clf()
    fig = plt.figure(figsize=(int(res[0]),int(res[1])))
    ax1 = fig.add_axes((0.1,0.25,0.8,0.7))
    title = "Precision of Assembly"
    ax1.set_title(title)
    ax1.set_ylabel("Precision")
    ax1.set_xlabel("Portion of aligned spots")
    ax1.scatter(data["sf"], data["precision"])
    ax1.set_xlim(data["sf"].min(),data["sf"].max())
    ax1.set_xticks(data["sf"])

    ax2 = fig.add_axes((0.1,0.2,0.8,0.0))
    ax2.set_xlim(data["sf"].min(),data["sf"].max())
    ax2.yaxis.set_visible(False)
    ax2.set_xticks(data["sf"])
    ax2.set_xlabel("# aligned spots")
    ticks2F = [str('%.2f' %((elem*spotsOriginal)/1000000))+"M" for elem in data["sf"].tolist()]
    ax2.set_xticklabels(ticks2F)

    plt.savefig(outDir+"/png/precision.png")

def plotPrecision_VS_Recall(data,outDir,res):
    plt.close("all")
    plt.clf()
    plt.figure(figsize=(int(res[0]),int(res[1])))
    plt.title("Recall vs Precision")
    fig, ax = plt.subplots()
    ax.scatter(data["recall"], data["precision"])
    ax.set_title("Recall versus Precision")
    ax.set_xlabel("Recall")
    ax.set_ylabel('Precision')
    minXTick = data["recall"].min()-0.0001*(data["recall"].min())
    maxXTick = data["recall"].max()+0.0001*(data["recall"].min())
    ax.set_xlim(minXTick,maxXTick)
    minYTick = data["precision"].min()-0.01*(data["precision"].min())
    maxYTick = data["precision"].max()+0.01*(data["precision"].min())
    ax.set_ylim(minYTick,maxYTick)

    ticks2F = [str('%.2f' %((elem*spotsOriginal)/1000000))+"M" for elem in data["sf"].tolist()]
    for i, txt in enumerate(data["sf"].tolist()):
        ax.annotate(str(txt)+" ("+ticks2F[i]+")", (data["recall"].tolist()[i],data["precision"].tolist()[i]))

    plt.savefig(outDir+"/png/recallPrecision.png")

def plotAll(data,outDir,res,iqrCoefficient,gif):
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

    eDF = pd.DataFrame(data["ID"])
    for cl in list(data)[1:-1]:
        if cl == "covBase":
            eDF[str(cl)] = pd.DataFrame(preprocessing.scale(boxcox(data[cl]+1)[0]))
        else:
            eDF[str(cl)] = data[str(cl)]
    eDF["sf"]=data['sf']
    eDF.to_csv(outDir+"/csv/groupedByIDTransformed.csv")

    distXCovBase = (max(data["covBase"])-min(data["covBase"]))
    distXCovBaseT = (max(eDF["covBase"])-min(eDF["covBase"]))
    distYExp = (max(data["tpmQ75"])-min(data["tpmQ25"]))
    distYExpT = (max(data["tpmQ50"])-min(data["tpmQ50"]))
    tickXCovBase = np.arange(min(data["covBase"]), max(data["covBase"])+distXCovBase*0.01, distXCovBase/20).tolist()
    tickXCovBaseT = np.arange(min(eDF["covBase"]), max(eDF["covBase"])+distXCovBaseT*0.01, distXCovBaseT/20).tolist()
    tickYExp = np.arange(min(data["tpmQ25"]), max(data["tpmQ75"])+distYExp*0.01, distYExp/20).tolist()
    tickYExpT = np.arange(min(data["tpmQ50"]), max(data["tpmQ50"])+distYExpT*0.01, distYExpT/20).tolist()

    plt.close("all")
    plt.clf()

    fig = plt.figure(figsize=(int(res[0]),int(res[1])))
    ims = []
    line, = plt.plot([], [], lw=0.25)
    plt.xlim(eDF["covBase"].min(), eDF["covBase"].max())
    plt.ylim(data["tpmQ50"].min(), data["tpmQ50"].max())
    def update_line(i,sfs):
        line.set_xdata(ims[i][0])
        line.set_ydata(ims[i][1])
        ax = line.properties()['axes']
        spotsRetained = spotsOriginal*sfs[i]
        ax.set_xlabel("Change in median of transcript expression levels(TPM) at "+str(sfs[i]*100)+"%% of original numner of spots or mean of "+str(int(spotsRetained))+" spots across 12 alignments")
        ax.set_ylabel('Deviation of expression estimate from control (%TPM)')
        ax.set_xticks(tickXCovBaseT)
        ax.set_yticks(tickYExpT)
        return line,

    fig2 = plt.figure(figsize=(int(res[0]),int(res[1])))
    ims2 = []
    line2, = plt.plot([], [], lw=0.25)
    plt.xlim(data["covBase"].min(), data["covBase"].max())
    plt.ylim(data["tpmQ50"].min(), data["tpmQ50"].max())
    def update_line2(i,sfs):
        line2.set_xdata(ims2[i][0])
        line2.set_ydata(ims2[i][1])
        ax = line2.properties()['axes']
        spotsRetained = spotsOriginal*sfs[i]
        ax.set_xlabel("Change in median of transcript expression levels(TPM) at "+str(sfs[i]*100)+"%% of original numner of spots or mean of "+str(int(spotsRetained))+" spots across 12 alignments")
        ax.set_ylabel('Deviation of expression estimate from control (%TPM)')
        ax.set_xticks(tickXCovBase)
        ax.set_yticks(tickYExpT)
        return line2,

    fig3 = plt.figure(figsize=(int(res[0]),int(res[1])))
    ims3 = []
    line3, = plt.plot([], [], lw=0.25)
    plt.xlim(data["covBase"].min(), data["covBase"].max())
    plt.ylim(data["falseNegative"].min(), data["falseNegative"].max())
    def update_line3(i,sfs):
        line3.set_xdata(ims3[i][0])
        line3.set_ydata(ims3[i][1])
        ax = line3.properties()['axes']
        spotsRetained = spotsOriginal*sfs[i]
        caption = "Figure. Number of false negatives reported across all random downsamplings to "+str(sfs[i]*100)+"%% spots retained for each transcript."
        ax.set_xlabel(caption)
        ax.set_ylabel('Number of false negatives')
        ax.set_xticks(tickXCovBase)
        return line3,

    unique=data["sf"].unique().tolist()
    for sf in unique:
        print("Plotting: ",sf)
        dataTMP = data[data["sf"]==sf]

        plt.close("all")
        ax=sbn.pairplot(dataTMP[['covBase','falseNegative','tpmMEAN','tpmQ50','tpmIQR','tpmSTD']],size=2,aspect=1)
        ax.savefig(outDir+"/png/"+str(sf)+"scatterMatrixByID.png")

        plt.close("all")
        ax=sbn.pairplot(eDF[eDF["sf"]==sf][['covBase','falseNegative','tpmMEAN','tpmQ50','tpmIQR','tpmSTD']],size=2,aspect=1)
        ax.savefig(outDir+"/png/"+str(sf)+"scatterMatrixByIDTransformed.png")

        plt.close('all')
        plt.clf()
        plt.figure(figsize=(int(res[0]), int(res[1])))
        plt.title('Normalized TPM Deviation')
        plt.xticks(np.arange(min(dataTMP["covBase"]), max(dataTMP["covBase"])+1, (max(dataTMP["covBase"])-min(dataTMP["covBase"]))/20).tolist())
        plt.yticks(np.arange(min(dataTMP["tpmQ25"]), max(dataTMP["tpmQ75"])+1, (max(dataTMP["tpmQ75"])-min(dataTMP["tpmQ25"]))/20).tolist())
        plt.ylabel('Deviation of expression estimate from control (% TPM)')
        plt.xlabel("Transcript Coverage")
        plt.plot(dataTMP["covBase"], dataTMP["tpmQ50"],'k',color='#CC4F1B',lw=0.25)
        plt.fill_between(dataTMP["covBase"], dataTMP["tpmQ25"], dataTMP["tpmQ75"],
            alpha=0.5, edgecolor='#CC4F1B', facecolor='#FF9848')
        spotsRetained = spotsOriginal*sf
        caption = "Figure. Change in median of transcript expression levels(TPM) at "+str(sf*100)+"%% of original numner of spots or mean of "+str(int(spotsRetained))+" spots across 12 alignments"
        # plt.figtext(.02, .02, caption,wrap=True)
        plt.savefig(outDir+"/png/"+str(sf)+"boxID.png")
        if gif == True:
            ims2.append((dataTMP["covBase"],dataTMP["tpmQ50"]))

        plt.close('all')
        plt.clf()
        plt.figure(figsize=(int(res[0]), int(res[1])))
        plt.title('Change in median, 2nd and 3rd quartiles of transcript expression levels(TPM)')
        plt.xticks(np.arange(min(eDF["covBase"]), max(eDF["covBase"])+1, (max(eDF["covBase"])-min(eDF["covBase"]))/20).tolist())
        plt.yticks(np.arange(min(eDF["tpmQ50"]), max(eDF["tpmQ50"])+1, (max(eDF["tpmQ50"])-min(eDF["tpmQ50"]))/20).tolist())
        plt.ylabel('Deviation of expression estimate from control (% TPM)')
        plt.xlabel("Transcript Coverage")
        plt.plot(eDF[eDF["sf"]==sf]["covBase"], eDF[eDF["sf"]==sf]["tpmQ50"],'k',color='#CC4F1B',lw=0.25)
        spotsRetained = spotsOriginal*sf
        caption = "Figure. Change in median of transcript expression levels(TPM) at "+str(sf*100)+"%% of original numner of spots or mean of "+str(int(spotsRetained))+" spots across 12 alignments. A BoxCoxTransformation was applied to coverage values to normalize the distribution."
        # plt.figtext(.02, .02, caption,wrap=True)
        plt.savefig(outDir+"/png/"+str(sf)+"boxIDTransformed.png")
        if gif == True:
            ims.append((eDF[eDF["sf"]==sf]["covBase"],eDF[eDF["sf"]==sf]["tpmQ50"]))
            
        plt.close("all")
        plt.clf()
        plt.figure(figsize=(int(res[0]),int(res[1])))
        plt.title("Number of false negatives versus coverage")
        plt.ylabel('Number of false negatives')
        plt.xlabel("Coverage")
        plt.xticks(np.arange(min(dataTMP["covBase"]), max(dataTMP["covBase"])+1, (max(dataTMP["covBase"])-min(dataTMP["covBase"]))/20).tolist())
        plt.plot(dataTMP["covBase"], dataTMP["falseNegative"],'k',color='#CC4F1B',lw=0.25)
        spotsRetained = spotsOriginal*sf
        caption = "Figure. Number of false negatives reported across all random downsamplings to "+str(sf*100)+"%% spots retained for each transcript."
        # plt.figtext(.02, .02, caption,wrap=True)
        plt.savefig(outDir+"/png/"+str(sf)+"falseNegative.png")
        if gif == True:
            ims3.append((dataTMP["covBase"],dataTMP["falseNegative"]))

        # plt.close("all")
        # plt.clf()
        # plt.figure(figsize=(int(res[0]),int(res[1])))
        # plt.title("Number of false negatives reported for each "++" number of samples per downsampling")
        # ax=eDF[eDF["sf"]==sf]['falseNegative'].plot(x=eDF[eDF["sf"]==sf]["covBase"],subplots=True, layout=(4, 4), figsize=(int(res[0]), int(res[1])), sharex=False)
        # plt.savefig(outDir+"/png/"+str(sf)+"groupedIDTransformed.png")

        plt.close("all")
        plt.clf()

        del ax
        del dataTMP

    if gif == True:
        ims.reverse()
        im_ani = animation.FuncAnimation(fig, update_line, len(unique), fargs=(list(reversed(unique)),), interval=600, blit=True)
        im_ani.save(outDir+'/png/boxIDTransformed.gif',writer="imagemagick",dpi=50)

        ims2.reverse()
        im_ani2 = animation.FuncAnimation(fig2, update_line2, len(unique),fargs=(list(reversed(unique)),), interval=600, blit=True)
        im_ani2.save(outDir+'/png/boxID.gif',writer="imagemagick",dpi=50)

        ims3.reverse()
        im_ani3 = animation.FuncAnimation(fig3, update_line3, len(unique),fargs=(list(reversed(unique)),), interval=600, blit=True)
        im_ani3.save(outDir+'/png/falseNegatives.gif',writer="imagemagick",dpi=50)

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
                    "tauBottom50",
                    "recall",
                    "precision"]

        dataSF = pd.read_csv(os.path.abspath(args.sf)).drop("Unnamed: 0",axis=1)
        dataSF.columns = headersSF
        plotBoxSF(dataSF,os.path.abspath(args.out),args.resolution.split(":"))
        plotTauSF(dataSF,os.path.abspath(args.out),args.resolution.split(":"))
        if (len(dataSF["recall"].unique()) == 1) and np.isnan(dataSF["recall"].unique().tolist()[0]):
            plotScattermatrixSFRange(dataSF,os.path.abspath(args.out),args.resolution.split(":"))
            plotPCARange(dataSF,os.path.abspath(args.out),args.resolution.split(":"))
            plotPrecision(dataSF,os.path.abspath(args.out),args.resolution.split(":"))
        else:
            plotScattermatrixSFFull(dataSF,os.path.abspath(args.out),args.resolution.split(":"))
            plotPCAFull(dataSF,os.path.abspath(args.out),args.resolution.split(":"))
            plotPrecision_VS_Recall(dataSF,os.path.abspath(args.out),args.resolution.split(":"))

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