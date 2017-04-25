# EVA - Expression Variation Analysis

(https://github.com/alevar/EVA/blob/master/sampleOut/png/boxID.gif)

EVA is a tool designed to facilitate better understanding of the effects sequence deth and coverage have on the estimation of transcription levels (TPM) in a transcriptome.

The end goal of running EVA is to gain sufficient understanding of what levels of variations in transcript expression estimates to expect based on the number of reads sequenced and particular research interests (such as but not limited to: interest in transcripts with low or high expressions levels). We hope this tool will aid researchers in estimating the costs for their RNA-seq experiments.

As an example EVA was run on the full GTEX dataset in order to infer the estimation precision and cost tradeoffs of lower sequencing coverage for human transcriptomes.

Workflow:
1. "eva select -i <path to the gtex data alignments made with hisat2> -n 20 -r auto"
	- this command parses through the alignment logs and builds information about each tissue and individual sample in the dataset and the total number of reads and percent of reads aligned.
	- the constructed dataset is saved locally as the output of the command
	- the auto range parameter binned the data and outputed the bin with the greatest mean number of aligned reads and saved this subset of the original dataset locally as well as a box plot for the number of aligned reads
	- after close inspection of both dataframes we came to the following conclusion:
		- samples suggested by the auto range parameter had very low standard deviation
		- the mean of 68M reads was too low as compared to many other samples in the original dataset which were distributed a little further from their means when binned
	- We ran the same command with custom ranges ("eva select -i <path to the gtex data alignments made with hisat2> -n 20 -r minbound:maxbound") to establish the best bin for the task.
		- in our case we've got 12 samples in the range 90M-96M
2. "eva assemble -i <paths to the samples as outputed by eva select> -r 0.1:1.0:0.1 -o ./test -a <path to the annotation> -e <path to the reference> -w 10 -f 10 -t 2"
	- Explanation (more info in eva assemble --help):
		-i - perform for each of the given alignments
		-r - uses samtools hidden -S option to randomly downsample each alignment starting at 0.1 (10% retention of original reads) to 1.0 (all original reads retained - no downsampling) in 0.1 (10%) increments
		-o - output directory
		-w - randomness of each downsampling is based on a randomly generated seed. this option specifies how many times each downsampling with a different random seed should be done
		-f - split eva assemble into 10 forks each carrying out the pipeline for one assigned input sample
		-t - the number of cpu threads each fork is allowed to use. This gets passed to stringtie and samtools.
	- For each assembly eva assemble will also extract information about transcripts and will produce a single dataset for all assemblies (for each input: for each downsampling: for each repetition). The datase contains the following information:
		- TISSUE CHROMOSOME REFERENCE_ID TRANSCRIPT_NAME START END SAMPLE FACTOR COVERAGE TPM
3. "eva analyze -i ./test/stats.log"
	- This step produces several types of grouped datasets from the input and constructs a multitude of plots for close inspection of the results
	- First the TPM values at each downsampling are normalized against their original value at no sampling. The resulting value is expressed as a percent away from the original. The rest of the statistical computations are based off of this value not the actual TPM.
	- Secondly the dataset is grouped by the sampling factor (% retained reads)
		- such yields the number of rows equal to the number of sampling factors. in our case of 0.1:1.0:0.1 we had 10 output rows.
		- The following statistical measures are computed for each row/sampling factor:
			- falsePositives - number of transcripts identified at downsampling but not at the original number of reads
			- falseNegatives - number of reads which were present in the assembly of the original alignment but are no longer detected after downsampling.
			- falseNegativesFull - number of reads which were present in the assembly of the original alignment but are no longer detected after downsampling in any of the repetitions of that particular downsampling (completely lost).
			- Quartiles, median,mean and whiskers of the percentAway value
				- low and high whiskers were defined as the lowest percentAway value above (quartile25-1.5*InterQuartileRange) and the greatest percentAway value under (quartile75+1.5*InterQuartileRange) respectively
			- Extreme outliers were defined as all the values that lie beyond the whiskers. Extremes were temporarily saved in a dictionary but are not reported in the dataset. However the number of extremes was reported as numExtremes
			- weightedNumExtremes - each number of extremes was weighted by the mean of the distance of all extremes for that sf from their original values
			- Standard Deviation
			- coefficient of variation
			- kendal's tau coefficient - for each downsampling the transcripts were ranked by their expression levels (TPM) and kendal's tau was used to compare the ranking against the ranking of the assembly with no downsampling.
		- The following plots are produced by eva analyze for the data grouped by sampling factor:
			- tauSF - kendal's tau coeeficient plots with a shared y-axis (tau coeeficient). Plots include:
				- all transcripts
				- top 10% of transcripts by percentAway value
				- top 20% of transcripts by percentAway value
				- top 50% of transcripts by percentAway value
				- bottom 50% of transcripts by percentAway value
				- bottom 20% of transcripts by percentAway value
				- bottom 10% of transcripts by percentAway value
			- boxSF - median and quartiles
			- PCA - a biplot of principal component 1 vs principal component 2. The biplot also shows the direction and weights of the vectors for each of the following input features:
				- falsePositives
				- falseNegatives
				- falseNegativesFull
				- median
				- weightedNormalizedNumExtremes
				- std
				- tauFull
	- Thirdly the dataset is grouped by the uniqueID of each transcript. Such grouping aggregates information from all random downsamplings for a particular sampling factor for each transcript.
		- The following statistical measures are computed for each row/transcript:
			- falseNegative
            - tpmMEAN
            - tpmSTD
            - tpmQ25
            - tpmQ50
            - tpmQ75
            - tpmCV
            - tpmIQR
        - the following plots are produced for each sampling factor:
        	- boxID - median and quartiles. X-axis: coverage; y-axis: percentAway
        	- groupedID - each statistical measure plotted separately against the coverage (as reported for unsampled assembly)
        	- scatterMatrixByID - scatter matrix of the measures against each other. Diagonal shows histograms
        	- Additionally the each measure was transformed via BoxCox transformation and the same plots were produced for each downsampling from the transformed data.

4. By analysing the results we came to the following conclusions:
	- To be updated in a few days. The final run is still working.

Please make sure the following packages/libraries are installed: 
1. Python:
	- appdirs
	- cycler
	- functools32
	- packaging
	- pyparsing
	- python-dateutil
	- pytz
	- six
	- pandas
	- matplotlib
	- seaborn
	- numpy
	- scipy
	- scikit-learn
	- subprocess
2. R:
	- ggplot2
	- ggbiplot
3. Other:
	- stringtie
	- hisat2
	- samtools