
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
from scipy.special import logsumexp
from sklearn.model_selection import train_test_split, LeaveOneOut, cross_val_score, cross_val_predict
from sklearn import svm
from sklearn.metrics import classification_report
from sklearn.pipeline import Pipeline
from sklearn.naive_bayes import GaussianNB
from sklearn.metrics import f1_score
import time
import os


# Read the CSV file
df2 = pd.read_csv("Esophageal_cancer.csv", delimiter = ',')
df2.columns = df2.columns.str.split('.').str[0]

# Separate columns with the same name into different DataFrames
print(df2)

df3 = df2['miRNA']
#print(df)

df = df2.filter(like='Control')

#print(df1)

df1 = df2.filter(like='Resistant')
deta_s = np.asarray(df).tolist()
deta_r = np.asarray(df1).tolist()
miRNA_list = df3.values.tolist()

data_s = np.array(deta_s)
data_r = np.array(deta_r)
data_f = np.concatenate((data_s, data_r), axis = 1)
a,b  = data_s.shape
c, d = data_r.shape
e,f = data_f.shape 
#print(data_f)
print (a,b,c,d,e,f)
############################################################################################################################

mean_of_each = np.mean(data_s, axis = 1)

mean_of_each_r = np.mean(data_r, axis = 1)

#print(len(mean_of_each), len(mean_of_each_r))

#print(mean_of_each, mean_of_each_r)
sensitive_array = np.ndarray((e, f))
resistive_array = np.ndarray((e, f))

print(df3)
print (a,b,c,d,e,f)
############################################################################################################################
mean_g = np.mean(data_f, axis = 1)

print(mean_g.shape)
h = e-1
g = 2*f
membership = np.ndarray((e, e, f))
membership_compliment = np.ndarray((e, e, f))
for i in range(e):
	for j in range(e):
		for k in range(f):
			if (i == j):
				membership[i][j][k] = 0
				membership_compliment[i][j][k] = 1
			else:
				if (mean_g[i] == data_f[i][k] and mean_g[j] == data_f[i][k]):
					membership[i][j][k] = 1/(1+1)
					membership_compliment[i][j][k] = 1-membership[i][j][k]
				elif(mean_g[j] == data_f[i][k]):
					membership[i][j][k] = 1/(1+1)
					membership_compliment[i][j][k] = 1-membership[i][j][k]
				else:
					membership[i][j][k] = 1/(1 + ((mean_g[i]-data_f[i][k])/(mean_g[j]-data_f[i][k]))**2)
					membership_compliment[i][j][k] = 1-membership[i][j][k]

#print(membership)
#print(membership_compliment)

#alpha =membership[ ~np.all(membership == 0, axis = 2)]

#bita = alpha.reshape(4,3,6)
#gammaa =membership_compliment[ ~np.all(membership_compliment == 1, axis = 2)]
#jammaa = gammaa.reshape(4,3,6)
#print(bita)#, jammaa)


membership_arrange = np.ndarray((e, e, f))
for i in range(e):
	for j in range(e):
		for k in range(f):
			membership_arrange[i][j][k] = membership_compliment[j][i][k]
			
#print(membership_arrange)
#gammaa =membership_arrange[ ~np.all(membership_arrange == 1, axis = 2)]
#jammaa = gammaa.reshape(4,3,6)

#print(jammaa)

real_membership = np.concatenate((membership, membership_arrange), axis = 2)

#print(real_membership)

Save_ones = (e,e,f)

array_ones = np.ones(Save_ones)
array_zeros = np.zeros(Save_ones)

max_min_rough_array = np.concatenate((array_ones, array_zeros), axis = 2)

#print(max_min_rough_array)

#upper_rough = np.max((max_min_rough_array, real_membership), axis = 2)     ###########this is not working

#print(upper_rough)   ######################this is also not working

upper_rough_array = np.ndarray((e, e, g))
lower_rough_array = np.ndarray((e, e, g))
for i in range(e):
	for j in range(e):
		for k in range(f):
			upper_rough_array[i][j][k] = max(max_min_rough_array[i][j][k], real_membership[i][j][k]) 
			lower_rough_array[i][j][k] = min(max_min_rough_array[i][j][k], real_membership[i][j][k]) 
#			

#print(upper_rough_array)
#print(lower_rough_array)

relative_fre_lower = np.sum(lower_rough_array, axis = 2)

#print(relative_fre_lower)
			
mean_rela_fre = np.divide(relative_fre_lower, f)

avg_rela_fre = np.ndarray((e,e))

for i in range(e):
	for j in range(e):
		avg_rela_fre[i][j] = 0.5 * (mean_rela_fre[i][j] + mean_rela_fre[j][i])	

#print(avg_rela_fre)
print(np.sum(avg_rela_fre, axis =1))
print(np.sum(avg_rela_fre, axis =0))
rel_fre_sum = np.sum(avg_rela_fre, axis =1)
###############################################this segment represent code to check different combination##################################################################
entropy = np.ndarray((e,e))
for i in range(e):
	for j in range(e):
		if (avg_rela_fre[i][j] == 0):
			entropy[i][j] = 1
		else:
			entropy[i][j] = -(avg_rela_fre[i][j] * math.log(avg_rela_fre[i][j]) + (1 - avg_rela_fre[i][j])* math.log(1 - avg_rela_fre[i][j]))


entropy2 = np.sum(entropy, axis =1)
entropy3 = []
for i in range(e):
	value = entropy2[i] - 1
	value1 = value/h
	entropy3.append(value1)

print(entropy3)
			
check_rel_fre = np.divide(rel_fre_sum, h)	
entropy1 = []
for i in range(e):
	value = -(check_rel_fre[i] * math.log(check_rel_fre[i]) + (1 - check_rel_fre[i])* math.log(1 - check_rel_fre[i]))
	entropy1.append(value)
	
print(entropy1)


#np.savetxt('entropy1.csv', entropy1, delimiter = '\t')
np.savetxt('Zredundancy.csv', entropy3, delimiter = '\t')
#################################################END#####################################################

mean_of_each = np.mean(data_s, axis = 1)

mean_of_each_r = np.mean(data_r, axis = 1)
#print(mean_of_each, mean_of_each_r)
sensitive_array = np.ndarray((e, f))
resistive_array = np.ndarray((e, f))


for i in range(e):
	for j in range(f):
		if(mean_of_each_r[i] == data_f[i][j] and mean_of_each[i] == data_f[i][j]):
			sensitive_array[i][j] = 1/(1 + 1)
			resistive_array[i][j] = 1/(1 + 1)		
		elif(mean_of_each[i] == data_f[i][j]):				
			sensitive_array[i][j] = 1/(1 + ((mean_of_each[i] - data_f[i][j])/(mean_of_each_r[i] - data_f[i][j]))**2)
			resistive_array[i][j] = 1/(1 + 1)
		elif(mean_of_each_r[i] == data_f[i][j]):				
			sensitive_array[i][j] = 1/(1 + 1)
			resistive_array[i][j] = 1/(1 + ((mean_of_each_r[i] - data_f[i][j])/(mean_of_each[i] - data_f[i][j]))**2)	
		else:				
			sensitive_array[i][j] = 1/(1 + ((mean_of_each[i] - data_f[i][j])/(mean_of_each_r[i] - data_f[i][j]))**2)
			resistive_array[i][j] = 1/(1 + ((mean_of_each_r[i] - data_f[i][j])/(mean_of_each[i] - data_f[i][j]))**2)

			
			

def assignment_fun1():
	for i in range(e):
		for j in range(f):
			sensitive_array[i][j] = float("{0:.8f}".format(sensitive_array[i][j]))
			resistive_array[i][j] = float("{0:.8f}".format(resistive_array[i][j]))
	return sensitive_array, resistive_array

new_arra1 = assignment_fun1()

#print(sensitive_array + resistive_array)	
###############################just to check##########################################################3

data_s_zero1 = np.zeros((a, b))
data_s_one1 = np.ones((a, b))
data_r_zero1 = np.zeros((c, d))
data_r_one1 = np.ones((c, d))

#print(data_s_zero.shape, data_s_one.shape, data_r_zero.shape, data_r_one.shape)

##########################################################fuction to caclculate upper and lower approximation of set########################################3
sensitive_n1 = np.append(data_s_one1, data_r_zero1, axis = 1)
resistive_c1 = np.append(data_s_zero1, data_r_one1, axis = 1)
#print(sensitive_n.sum())
#print(resistive_c.sum())
upper_sensi1 = np.ndarray((e, f))
lower_sensi1 = np.ndarray((e, f))
upper_resis1 = np.ndarray((e, f))
lower_resis1 = np.ndarray((e, f))
######################################this whole function is for mj lower amd upper###################
def upper_lower_rough1():
	for i in range(e):
		for j in range(f):
			upper_sensi1[i][j] = max(sensitive_array[i][j], sensitive_n1[i][j])
			lower_sensi1[i][j] = min(sensitive_array[i][j], sensitive_n1[i][j])
			upper_resis1[i][j] = max(resistive_array[i][j], resistive_c1[i][j])
			lower_resis1[i][j] = min(resistive_array[i][j], resistive_c1[i][j])
	return upper_sensi1, lower_sensi1, upper_resis1, lower_resis1

upper_low1 = upper_lower_rough1()

#print upper_sensi,lower_sensi, upper_resis, lower_resis 
fre_upper_s1 = upper_sensi1.sum(axis = 1)

fre_lower_s1 = lower_sensi1.sum(axis = 1)

fre_upper_r1 = upper_resis1.sum(axis = 1)

fre_lower_r1 = lower_resis1.sum(axis = 1)

#print(fre_lower_r, fre_lower_s)
#@$$@%%@@@@@@@@@@@@@@###############################

#print fre_upper_s, fre_lower_s, fre_upper_r, fre_lower_r 
rel_fre_us1 = np.ndarray(e)
rel_fre_ls1 = np.ndarray(e)
rel_fre_ur1 = np.ndarray(e)
rel_fre_lr1 = np.ndarray(e)

def relative_frequency1():
	for i in range(e):
		rel_fre_us1[i] = fre_upper_s1[i]/b
		rel_fre_ls1[i] = fre_lower_s1[i]/b
		rel_fre_ur1[i] = fre_upper_r1[i]/d		
		rel_fre_lr1[i] = fre_lower_r1[i]/d
	return rel_fre_us1, rel_fre_ls1, rel_fre_ur1, rel_fre_lr1

frequency_rough1 = relative_frequency1()
#print(rel_fre_ls, rel_fre_lr)
average_frequency_lower1 = np.ndarray(e)
average_frequency_overlapping1 =np.ndarray(e)

for i in range(e):
	average_frequency_lower1[i] = (rel_fre_ls1[i] + rel_fre_lr1[i])/2
	average_frequency_overlapping1[i] = ((1-rel_fre_ls1[i]) + (1- rel_fre_lr1[i]))/2
	#average_frequency_lower.append(value), average_frequency_overlapping.append(value1)

#print (average_frequency_lower, average_frequency_overlapping)
#####################################to check the properties##############################################################


################################################entropy########################################################################3
entropy_array2 = []

def entropy_func1():
	for i in range(e):
		value = -(average_frequency_lower1[i]*math.log(average_frequency_lower1[i]) + average_frequency_overlapping1[i]*math.log(average_frequency_overlapping1[i]))
		entropy_array2.append(value)
	return entropy_array2

entropy_f1 = entropy_func1()

#print(entropy_array)
for i in range(e):
	global newss1
	
	if (average_frequency_overlapping1[i] + average_frequency_lower1[i] > 1):
		newss1 = average_frequency_overlapping1[i] + average_frequency_lower1[i] -1
		print(newss1)

print(entropy_array2)
print(max(entropy_array2), min(entropy_array2))

np.savetxt('Zrelevance.csv', entropy_array2, delimiter = "\t")

#######################################################################################

#print(data_s)
df4 = pd.read_csv("Zrelevance.csv", delimiter = '\t', header = None)
df5 = pd.read_csv("Zredundancy.csv", delimiter = '\t', header = None)
deta_s = np.asarray(df)
deta_r = np.asarray(df1)
#deta_f = np.asarray(df2).tolist()
#print (deta_s[1998][4])
miRNA_list = df3.values.tolist()

data_s = np.array(deta_s)
data_r = np.array(deta_r)
data_f = np.concatenate((data_s, data_r), axis = 1)
a,b  = data_s.shape
c, d = data_r.shape
e,f = data_f.shape 
#print(data_f)
print (a,b,c,d,e,f)
print(data_s)
############################################################################################################################

mean_of_each = np.mean(data_s, axis = 1)

mean_of_each_r = np.mean(data_r, axis = 1)

#print(len(mean_of_each), len(mean_of_each_r))
entropy1 =  np.array(df4)
entropy2 =  np.array(df5)

#print(entropy1, entropy2)

weighted_array = np.ndarray(a)
sensitivity_array = []
specificity_array = []
fscore = []
Accuracy = []
MCC_array = []

#input("Enter the value of weight", weight)
weight = 0.01
while(weight <1):
    for i in range(a):
        weighted_array[i] = weight*entropy1[i] + (1-weight)*(1-entropy2[i])
    sorted_entropy = []
    for i in range(a):
        value = weighted_array[i]
        sorted_entropy.append(value)
    sorted_entropy.sort(reverse = False)
    array_value = []
    def sorted_algo():
        for i in range(a):
            for j in range(a):
                if (sorted_entropy[i] == weighted_array[j]):
                    value = j
                    array_value.append(value)
                else:
                    continue
        return array_value
    algo_sorted = sorted_algo()
    #selected_percentage_miRNA = float(input("Enter the number of miRNAs in percentage: "))
    selection_len = math.floor(0.01 * a)
    classifier_array = np.ndarray((selection_len, b))
    miRNA_print = []
    for i in range(selection_len):
        value = array_value[i]
        value1 = miRNA_list[value]
        miRNA_print.append(value1)
    for i in range(selection_len):
        for j in range(b):
            classifier_array[i][j] = data_s[array_value[i]][j]
    classifier_array1 = np.ndarray((selection_len, d))
    for i in range(selection_len):
        for j in range(d):
            classifier_array1[i][j] = data_r[array_value[i]][j]
    data_f12 = np.concatenate((classifier_array, classifier_array1), axis = 1)
    data_f11 = data_f12.T
    alphaa, bita = data_f11.shape
    join11 = np.zeros(b)
    join12 = np.ones(d)
    join13 = np.concatenate((join11, join12), axis = 0)
    y = join13.T
    z = data_f11
    array_fine = np.asarray(z)
    z_array = np.ndarray((bita,alphaa))
    acc_array = np.ndarray((bita,alphaa))
    X = np.ndarray(alphaa)
    for i in range(bita):
        new_alphaa = np.array(array_fine[:, i])
        X = np.array(new_alphaa).reshape(-1, 1)
        pipe = Pipeline([
			('model', svm.SVC(kernel='linear', C=5))
		])
        start = time.time()
        result_clf = cross_val_score(estimator = pipe, X = X, y = y, scoring = 'accuracy', cv = LeaveOneOut())
        result_clf1 = cross_val_predict(estimator = pipe, X = X, y = y, cv = LeaveOneOut())
        acc_array[i] = result_clf
        z_array[i] = result_clf1
    np.array(y)
    repetitions = bita
    accuracy_array = np.tile(y, (repetitions, 1))
    TP = []
    TN = []
    FP = []
    FN = []
    def confusion_matrixes():
        for i in range(bita):
            for j in range(alphaa):
                value = z_array[i][j]
                if(value == 0 and accuracy_array[i][j] == 0):
                    TN.append(value)
                elif(value ==0 and accuracy_array[i][j] == 1):
                    FN.append(value)
                elif(value == 1 and accuracy_array[i][j] == 1):
                    TP.append(value)
                elif(value == 1 and accuracy_array[i][j] == 0):
                    FP.append(value)
                else:
                    break
        return TP,TN,FP,FN
    confusion_matrixess = confusion_matrixes()
    f11 = len(TP)/(len(TP) + (1/2)*(len(FP) + len(FN)))
    s = len(TP)
    t = len(FP)
    u = len(TN)
    v = len(FN)
    sensitivity = s/(s + v)
    specificity = u/(u + t)
    MCC = (s*u - t*v)/(math.sqrt((s+t)*(s+v)*(u+t)*(u+v)))
    sensitivity_array.append(sensitivity)
    specificity_array.append(specificity)
    fscore.append(f11)
    Accuracy.append(np.mean(acc_array)*100)
    MCC_array.append(MCC)
    print(sensitivity, specificity, f11, np.mean(acc_array)*100, MCC)
    weight += 0.01
    
filename = os.path.splitext(os.path.basename('Book1.csv'))[0]
f1score = max(fscore)
max_index = fscore.index(f1score)
sensitivity = sensitivity_array[max_index]
specificity = specificity_array[max_index]
MCC = MCC_array[max_index]
Accuracy = Accuracy[max_index]

data = {
    'Measures': ['sensitivity', 'specificity', 'f1score', 'Accuracy', 'MCC'],
    filename: [sensitivity, specificity, f1score, Accuracy, MCC]
}

# Create a DataFrame
result_df = pd.DataFrame(data)

print(result_df)
result_df.to_csv('WFIFRRRE_Performance.csv', index=False)

weight = (max_index+2)/100
for i in range(a):
    weighted_array[i] = weight*entropy1[i] + (1-weight)*(1-entropy2[i])
sorted_entropy = []
for i in range(a):
    value = weighted_array[i]
    sorted_entropy.append(value)
sorted_entropy.sort(reverse = False)
array_value = []
def sorted_algo():
    for i in range(a):
        for j in range(a):
            if (sorted_entropy[i] == weighted_array[j]):
                value = j
                array_value.append(value)
            else:
                continue
    return array_value
algo_sorted = sorted_algo()
selection_len = math.floor(0.01 * a)
classifier_array = np.ndarray((selection_len, b))
miRNA_print = []
for i in range(selection_len):
    value = array_value[i]
    value1 = miRNA_list[value]
    miRNA_print.append(value1)

print(miRNA_print)