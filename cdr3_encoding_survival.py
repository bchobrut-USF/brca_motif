import pandas as pd
from localcider.sequenceParameters import SequenceParameters
import numpy as np
import itertools
from lifelines.statistics import logrank_test
from lifelines import KaplanMeierFitter
import os
from multiprocessing import Pool, cpu_count
import time
start_time = time.time()
from tqdm import tqdm

li = []
for i in itertools.product(["N", "P", "H", "A", "X"], repeat=6):
    li.append(''.join(map(str, i)))
 
for i in itertools.product(["N", "P", "H", "A", "X"], repeat=5):
    li.append(''.join(map(str, i)))
    
for i in itertools.product(["N", "P", "H", "A", "X"], repeat=4):
    li.append(''.join(map(str, i)))



def getdata(cancer, receptor, sample, maindir):
    cdr3 = pd.read_hdf(maindir+"VDJ Recoveries/all_vdj_encoded.hdf", "with_aromatic")
    cdr3 = cdr3[cdr3["Cancer"] == cancer]
    cdr3 = cdr3.loc[cdr3["Receptor"].str.contains(receptor)]
    cdr3 = cdr3.loc[cdr3["Sample"].str.contains(sample)]
    clinical = pd.read_csv(maindir+"VDJ Recoveries/%s_Results/clinical.csv"%cancer)
    clinical["OS_STATUS"] = np.where(clinical["OS_STATUS"]=="LIVING", 0, 1)
    clinical["DFS_STATUS"] = np.where(clinical["DFS_STATUS"]=="DiseaseFree", 0, 1)
    clinical["Filename"] = clinical["PATIENT_ID"]
    cdr3 = cdr3.merge(clinical[["Filename", "OS_MONTHS", "OS_STATUS", "DFS_MONTHS", "DFS_STATUS"]], on='Filename', how='inner')
    cdr3["encoded"] = cdr3["encoded"].fillna("Z")
    
    cdr3 = cdr3.reset_index(drop=True)
    cdr3[["OS_MONTHS", "OS_STATUS", "DFS_MONTHS", "DFS_STATUS"]] = cdr3[["OS_MONTHS", "OS_STATUS", "DFS_MONTHS", "DFS_STATUS"]].apply(pd.to_numeric, errors='coerce')
    return cdr3
    
def medians(time_A, time_B, event_observed_A, event_observed_B):
    kmfa = KaplanMeierFitter()
    kmfb = KaplanMeierFitter()
    kmfa.fit(time_A, event_observed=event_observed_A)
    kmfb.fit(time_B, event_observed=event_observed_B)
    return kmfa.median_survival_time_, kmfb.median_survival_time_
    
def statscalc(i):
    global cdr3
    try:
        cdr3_os = cdr3.dropna(subset=["OS_STATUS", "OS_MONTHS"])
        cdr3_os = cdr3_os.reset_index(drop = True)
        cdr3_dfs = cdr3.dropna(subset=["DFS_STATUS", "DFS_MONTHS"])
        cdr3_dfs = cdr3_dfs.reset_index(drop = True)

        df1_os = cdr3_os[cdr3_os["encoded"].str.contains(i)]
        df1_os = df1_os.drop_duplicates("Filename")
        df1_os = df1_os.reset_index(drop = True)

        n_1_os = len(df1_os)


        if n_1_os < 50:
            return (np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan)
        else:

            df2_os = cdr3_os[~cdr3_os["Filename"].isin(df1_os["Filename"])]
            df2_os = df2_os.drop_duplicates("Filename")
            df2_os = df2_os.reset_index(drop = True)


            df1_dfs = cdr3_dfs[cdr3_dfs["encoded"].str.contains(i)]
            df1_dfs = df1_dfs.drop_duplicates("Filename")
            df1_dfs = df1_dfs.reset_index(drop = True)
            df2_dfs = cdr3_dfs[~cdr3_dfs["Filename"].isin(df1_dfs["Filename"])]
            df2_dfs = df2_dfs.drop_duplicates("Filename")
            df2_dfs = df2_dfs.reset_index(drop = True)

            logrank_os = logrank_test(df1_os["OS_MONTHS"], df2_os["OS_MONTHS"], event_observed_A = df1_os["OS_STATUS"], event_observed_B = df2_os["OS_STATUS"])

            logrank_dfs = logrank_test(df1_dfs["DFS_MONTHS"], df2_dfs["DFS_MONTHS"], event_observed_A = df1_dfs["DFS_STATUS"], event_observed_B = df2_dfs["DFS_STATUS"])

            median_1_os, median_2_os = medians(df1_os["OS_MONTHS"], df2_os["OS_MONTHS"], event_observed_A = df1_os["OS_STATUS"], event_observed_B = df2_os["OS_STATUS"])

            median_1_dfs, median_2_dfs = medians(df1_dfs["DFS_MONTHS"], df2_dfs["DFS_MONTHS"], event_observed_A = df1_dfs["DFS_STATUS"], event_observed_B = df2_dfs["DFS_STATUS"])


            n_2_os = len(df2_os)
            n_1_dfs = len(df1_dfs)
            n_2_dfs = len(df2_dfs)
            return (logrank_os.p_value, logrank_dfs.p_value, median_1_os, median_2_os, median_1_dfs, median_2_dfs, n_1_os, n_2_os, n_1_dfs, n_2_dfs, i)
    except:
        return (np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan)
        
    

def kmcurve(cdr3):
    results = pd.DataFrame()
    substring_list = []
    logrank_os_p_value_list = []
    logrank_dfs_p_value_list = []
    median_1_os_list = []
    median_2_os_list = []
    median_1_dfs_list = []
    median_2_dfs_list = []
    n_1_os_list = []
    n_2_os_list = []
    n_1_dfs_list = []
    n_2_dfs_list = []

    p = Pool(cpu_count())
    source = p.map(statscalc, li)
    for result in source:
        # for x in li: 
        #     result = statscalc(x) 
        logrank_os_p_value_list.append(result[0])
        logrank_dfs_p_value_list.append(result[1])
        median_1_os_list.append(result[2])
        median_2_os_list.append(result[3])
        median_1_dfs_list.append(result[4])
        median_2_dfs_list.append(result[5])
        n_1_os_list.append(result[6])
        n_2_os_list.append(result[7])
        n_1_dfs_list.append(result[8])
        n_2_dfs_list.append(result[9])
        substring_list.append(result[10])
    p.close()
    p.join()
    
    results["substring"] = substring_list
    results["km_p_os"] = logrank_os_p_value_list
    results["km_p_dfs"] = logrank_dfs_p_value_list
    results["median_substring_os"] = median_1_os_list
    results["median_no_substring_os"] = median_2_os_list
    results["median_substring_dfs"] = median_1_dfs_list
    results["median_no_substring_dfs"] = median_2_dfs_list
    results["n_substring_os"] = n_1_os_list
    results["n_no_substring_os"] = n_2_os_list
    results["n_substring_dfs"] = n_1_dfs_list
    results["n_no_substring_dfs"] = n_2_dfs_list
    
    print(results)
    results = results.dropna()
    results = results.reset_index(drop=True)
    print(results)
    return results


allresults = pd.DataFrame()
maindir = "/mnt/usfboxsync/Boris Blanck Lab Work/"
cancers = [ name for name in os.listdir(maindir + "VDJ Recoveries/") if os.path.isdir(os.path.join(maindir, "VDJ Recoveries/", name))]
for cancer in tqdm(cancers):
    cancer = cancer.replace("_Results", "")
    print (cancer)
    for receptor in ["TRG"]:#["TRA|TRB", "TRA", "TRB", "IGH", "IGH|IGK|IGL", "IGH|IGK", "IGH|IGL", "TRG", "TRD", "TRG|TRD"]:
        print (receptor)
        for sample in ["01|06|10", "01|06", "10"]:
            print (sample)
            try:
                cdr3 = getdata(cancer, receptor, sample, maindir)
                print ("Patients: " + str(len(set(cdr3["Filename"]))))

                results = kmcurve(cdr3)
                results.insert(loc = 0, column = "Sample", value = sample)
                results.insert(loc = 0, column = "Receptor", value = receptor)
                results.insert(loc = 0, column = "Cancer", value = cancer)
                allresults = pd.concat([allresults, results], ignore_index=True)
                allresults = allresults.sort_values("km_p_os", ascending = True)
                allresults = allresults.reset_index(drop=True)
                allresults.to_csv("/mnt/usfboxsync/Andrea/Blanck/BreastCancer/all_cancer_motif_survival.csv", index=False)
            except:
                print ("ERROR")
                pass

        
            
