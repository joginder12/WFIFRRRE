# Integrating Fuzzy Rough Set Based Entropies for Identifying Drug Resistant miRNAs in Cancer (WFIFRRRE)
                               Authors: Joginder Singh and Shubhra Sankar Ray
                             E-Mail: joginder265@gmail.com and shubhra@isical.ac.in

## Download drug-resistant miRNA expression data from:
a)    For esophageal cancer data: Esophageal_Cancer.csv from https://drive.google.com/file/d/15bkTE8p5gpJkQmvlbcmhExbBi7ohHaPW/view.  

b)        For lung cancer data: Lung_cancer.csv from https://drive.google.com/file/d/1dIWvaRnXesxZU7STZ_zOvMZJmZKt4mBj/view.
2.        Open python and install the packages numpy, math, csv, pandas, sklearn, matplotlib, time, scipy.
        ( Use command 'pip install package_name' e.g., 'pip install pandas'. In higher versions of python use pip3 in place of pip. )
*         In windows environment if Spyder is used for python then one has to install pip package first using command "python get-pip.py"
3.        Download the code for WFIFRRRE from:WFIFRRRE.py
4.        Keep the code and the datasets in the same folder, otherwise change the folder path along with the name of the dataset in the code (Line number 18).
5.        Run 'WFIFRRRE.py' to produce 'WFIFRRRE_Result.csv' and 'WFIFRRRE_performance.csv'. 'WFIFRRRE_Result.csv' contains the miRNA names. 'WFIFRRRE_performance.csv' contains the classification performance achieved by the method considering the selected 1% of miRNAs.
6.        The performance of WFIFRRRE is provided using one percent of miRNAs. For checking the performance using higher number of miRNAs, provide the percentage of miRNAs (e.g.: 0.01 for 1%, 0.05 for 5%) to be selected in Line Number 374, in 'WFIFRRRE.py' code.
