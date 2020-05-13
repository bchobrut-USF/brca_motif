import pandas as pd
from lifelines.statistics import logrank_test
from os import listdir

def get_data(path):
    df_master = pd.read_csv(path)
    
    p_values_list = []
    n_a_list = []
    n_b_list = []
    
    for i in range(5):
        df = df_master.sample(n = int(0.5 * len(df_master)))

        list_df_cols = list(df)
        filename_col = list_df_cols[0]
        months_col = list_df_cols[1]
        group_a_name = list_df_cols[2]
        group_b_name = list_df_cols[3]
        # only for 2 variable km curves
        df[months_col] = pd.to_numeric(df[months_col], errors = "coerce")
        df[group_a_name] = pd.to_numeric(df[group_a_name], errors = "coerce")
        df[group_b_name] = pd.to_numeric(df[group_b_name], errors = "coerce")
        
        group_a_cols = [filename_col, months_col, group_a_name]
        group_b_cols = [filename_col, months_col, group_b_name]

        group_a_df = df[df[group_a_name].notnull()]
        group_b_df = df[df[group_b_name].notnull()]

        group_a_df = group_a_df[group_a_cols]
        group_b_df = group_b_df[group_b_cols]

        group_a_df_rand = group_a_df
        group_b_df_rand = group_b_df
        km = logrank_test(group_a_df_rand[months_col], group_b_df_rand[months_col], event_observed_A = group_a_df_rand[group_a_name], event_observed_B = group_b_df_rand[group_b_name])
        p_values_list.append(km.p_value)
        
        n_one = len(group_a_df_rand)
        n_two = len(group_b_df_rand)
        
        n_a_list.append(n_one)
        n_b_list.append(n_two)
        
    return p_values_list, n_a_list, n_b_list


all_results = pd.DataFrame()

for file in listdir("/mnt/usfboxsync/Andrea/Blanck/BreastCancer/survival_randomization"):
    if ".csv" in file and "results" not in file:      
        results = get_data("/mnt/usfboxsync/Andrea/Blanck/BreastCancer/survival_randomization/" + file)   
        
        p_values_list = results[0]
        n_a_list = results[1]
        n_b_list = results[2]
        
        mini_results = pd.DataFrame()
        
        mini_results["p_value"] = p_values_list
        mini_results["file"] = file
        mini_results["n_a"] = n_a_list
        mini_results["n_b"] = n_b_list
        
        all_results = pd.concat([all_results, mini_results], ignore_index=True)

all_results = all_results[["file", "p_value", "n_a", "n_b"]]

all_results.to_csv("/mnt/usfboxsync/Andrea/Blanck/BreastCancer/survival_randomization/results_whole_rand_v5.csv", index=False)