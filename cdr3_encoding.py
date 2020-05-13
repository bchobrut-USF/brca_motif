import pandas as pd
from localcider.sequenceParameters import SequenceParameters
import numpy as np
from tqdm import tqdm
import itertools
from multiprocessing import Pool, cpu_count

li = []
for i in itertools.product(["N", "P", "A", "H", "X"], repeat=6):
    li.append(''.join(map(str, i)))

for i in itertools.product(["N", "P", "A", "H", "X"], repeat=5):
    li.append(''.join(map(str, i)))
    
for i in itertools.product(["N", "P", "A", "H", "X"], repeat=4):
    li.append(''.join(map(str, i)))

print len(li)
##motif vs all remaining
##blood and tumor separate and together
##OS only
##put aromatic back and combined
##with and without clipped ends
##6-mers, 5-mers
##n occurences


masterdf = pd.read_hdf("/home/boris/emcloud/usfboxsync/Boris Blanck Lab Work/VDJ Recoveries/all_vdj_recoveries.hdf", "physicochem")


grps = [["E", "D"], ["R", "K"], ["F", "Y", "W"], ["A", "V", "I","L", "M"]]

def encode(cdr3):
    #print cdr3
    cdr3 = cdr3.replace("X", "V")
    cdr3 = cdr3[1:-3]
    #print cdr3
    cidercdr3 = SequenceParameters(str(cdr3))
    a = cidercdr3.get_linear_sequence_composition(blobLen=1, grps=grps)
    
    res = np.where(a[1][0]==1, "N",a[1][0]) 
    res = np.where(a[1][1]==1, "P", res) 
    res = np.where(a[1][2]==1, "H", res)
    res = np.where(a[1][3]==1, "H", res)
    res = np.where(res=="0.0", "X", res)
    res = "".join(res)
    
    return res



encoded = []

for i in tqdm(range(len(masterdf))):
    encoded.append(encode(masterdf["CDR3"].iloc[i]))

masterdf["encoded"] = encoded
masterdf.to_hdf("/home/boris/emcloud/usfboxsync/Boris Blanck Lab Work/VDJ Recoveries/all_vdj_encoded.hdf", "no_aromatic")

# masterdf = pd.read_hdf("/home/boris/emcloud/usfboxsync/Boris Blanck Lab Work/VDJ Recoveries/all_vdj_encoded.hdf", "with_aromatic")
# results = pd.DataFrame()
# 
# def subsearch(i):
#     global masterdf
#     df = masterdf[masterdf["encoded"].str.contains(i)]
#     df["substring"] = i
#     df = df.groupby(["substring", "Cancer", "Receptor"]).nunique()
#     df = df[["encoded","CDR3", "Filename"]] 
#     df = df.reset_index()
#     return df
#    
# 
# p = Pool(cpu_count())
# source = p.map(subsearch, li)
# for result in source:
#     results = pd.concat([results, result], ignore_index=True)
# results.to_hdf("/home/boris/emcloud/usfboxsync/Boris Blanck Lab Work/VDJ Recoveries/all_vdj_encoded.hdf", "with_aromatic_counts")
# 
# results = results.rename(columns={"encoded": "NumEncoded", "CDR3": "NumCDR3", "Filename":"NumFilename"})
"""



masterdf.to_csv("/Users/boris/Desktop/pan_cdr3_endocded_raw.csv", index=True)

df = masterdf.groupby(["encoded", "Receptor", "Cancer"]).nunique()

df.to_csv("/Users/boris/Desktop/pan_cdr3_counts_endocded.csv", index=True)

df = df.sort_values("Cancer", ascending=False)

df = df[["CDR3", "Cancer", "Filename"]]

df = df.reset_index()

df = df[["encoded", "CDR3", "Receptor", "Cancer", "Filename"]] 

df.to_csv("/Users/boris/Desktop/pan_cdr3_counts_endocded.csv", index=True)

df = df.merge(masterdf[["CDR3", "ncpr", "uversky_hydropathy", "aromaticity", "delta"]], how='left', on = 'CDR3')

vdjdb = pd.read_csv("/Users/boris/Box Sync/human_CDR3.csv")

df = df.merge(vdjdb, on="CDR3", how="left")
#once = df[df["Cancer"] == 1]
# 
df = df[["encoded", "CDR3", "Receptor", "Cancer", "Filename", "Epitope species", "Epitope gene", "ncpr", "uversky_hydropathy", "aromaticity", "delta"]] 



df1 = df[df["Cancer"] >1]

df2 = df[df["Cancer"] ==1]

df2 = df2.merge(masterdf[["CDR3", "Cancer"]], on="CDR3", how='left')

df2["Cancer"] = df2["Cancer_y"]

df2 = df2[["encoded", "CDR3", "Receptor", "Cancer", "Filename", "Epitope species", "Epitope gene", "ncpr", "uversky_hydropathy", "aromaticity", "delta"]] 


df = pd.concat([df1,df2], ignore_index=True)

df = df.rename(columns={"Cancer": "NumCancers", "Filename": "NumFilenames", "Epitope species":"VDJDB_Epitope_Species", "Epitope gene":"VDJDB_Epitope_Gene"})

df = df.drop_duplicates("CDR3")

df.to_csv("/Users/boris/Desktop/pan_cdr3_counts_endocded.csv", index=False)

"""
