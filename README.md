# Integrating Fuzzy Rough Set Based Entropies for Identifying Drug Resistant miRNAs in Cancer (WFIFRRRE)

## Abstract
MicroRNAs (miRNAs) are key biomarkers in cancer diagnosis and treatment. Identification of drug-resistant
miRNAs may help in effective treatment of cancer. Two new z score based fuzzy rough relevance and
redundancy entropies are developed and then a weighted framework is introduced to integrate the entropies
for ranking and selecting miRNAs in classifying control and drug resistant patients. Here, two key components
of soft computing, fuzzy set and rough set are utilized. The methodology is called a weighted framework
for integrating fuzzy rough set-based relevance and redundancy entropies (WFIFRRRE). The z score is used
to compute the fuzzy membership of expression values required for both entropies. Fuzziness deals with the
overlapping nature of miRNA expression profiles and rough set helps in determining the exact class size.
The weights in WFIFRRRE, assigned to relevance and redundancy entropies, are determined in a supervised
manner to maximize the 𝐹 score used for validating the classification performance in discriminating the
control and drug-resistant patients. The weights are varied from 0 to 1 in steps of 0.01 which enables an
integration between relevance and redundancy entropies. A subset of miRNAs is selected from the ranked list
and the performance is evaluated using three benchmark classifiers on eight drug-resistant cancer datasets.
Experimental results show that WFIFRRRE provides better prediction accuracy than the popular methods used
for comparison. The classification accuracy in terms of 𝐹 score, achieved by WFIFRRRE, ranges from 0.74 to
1.0, 0.75 to 1.0, and 0.73 to 1.0 using random forest, Naive Bayes, and linear SVM classifiers, respectively.
The resultant set of miRNAs obtained using WFIFRRRE is also verified with the help of existing biological
studies. The source code of WFIFRRRE is available at https://www.isical.ac.in/~shubhra/WFIFRRRE.html.


Download drug-resistant miRNA expression data from:
a)        For esophageal cancer data: Esophageal_Cancer.csv from https://drive.google.com/file/d/15bkTE8p5gpJkQmvlbcmhExbBi7ohHaPW/view.
b)        For lung cancer data: Lung_cancer.csv from https://drive.google.com/file/d/1dIWvaRnXesxZU7STZ_zOvMZJmZKt4mBj/view.
2.        Open python and install the packages numpy, math, csv, pandas, sklearn, matplotlib, time, scipy.
        ( Use command 'pip install package_name' e.g., 'pip install pandas'. In higher versions of python use pip3 in place of pip. )
*         In windows environment if Spyder is used for python then one has to install pip package first using command "python get-pip.py"
3.        Download the code for WFIFRRRE from:WFIFRRRE.py
4.        Keep the code and the datasets in the same folder, otherwise change the folder path along with the name of the dataset in the code (Line number 18).
5.        Run 'WFIFRRRE.py' to produce 'WFIFRRRE_Result.csv' and 'WFIFRRRE_performance.csv'. 'WFIFRRRE_Result.csv' contains the miRNA names. 'WFIFRRRE_performance.csv' contains the classification performance achieved by the method considering the selected 1% of miRNAs.
6.        The performance of WFIFRRRE is provided using one percent of miRNAs. For checking the performance using higher number of miRNAs, provide the percentage of miRNAs (e.g.: 0.01 for 1%, 0.05 for 5%) to be selected in Line Number 374, in 'WFIFRRRE.py' code.
