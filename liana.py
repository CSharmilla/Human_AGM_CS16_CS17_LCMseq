# -*- coding: utf-8 -*-
import pandas as pd
import os
import liana as ln
import decoupler as dc

folder = "C:\\Users\\schandr3\\Documents\\R_scripts\\AGM_CS16_CS17_LCMseq\\python_input"

expr_dict = {}

for file in os.listdir(folder):
    #print(file)
    #file = "V_D.csv"
    if file.endswith(".csv"):
        name = file.replace(".csv", "")
        
        df = pd.read_csv(os.path.join(folder, file), index_col=0)
        
        # clean
        df.index = df.index.str.upper()
        df = df[~df.index.duplicated()]
        df = df[df.index.notna()]
        
        expr_dict[name] = df

print(expr_dict.keys())


liana_lr = ln.resource.select_resource()
#liana_lr = ln.resource.explode_complexes(liana_lr)

# Create two new DataFrames, each containing one of the pairs of columns to be concatenated
#df1 = liana_lr[['interaction', 'ligand']]
#df2 = liana_lr[['interaction', 'receptor']]

# Rename the columns in each new DataFrame
#df1.columns = ['interaction', 'genes']
#df2.columns = ['interaction', 'genes']

liana_lr.columns = ['source', 'target']

# Concatenate the two new DataFrames
#liana_lr = pd.concat([df1, df2], axis=0)
liana_lr['weight'] = 1

# Find duplicated rows
duplicates = liana_lr.duplicated()

# Remove duplicated rows
liana_lr = liana_lr[~duplicates]

liana_lr

# run for all datasets
lr_results = {}

for name, expr in expr_dict.items():
    
    print(f"Running LR for {name}...")
    
    lr_score, lr_pval = dc.mt.ulm(
        data=expr,
        net=liana_lr,
        tmin=0,
        verbose=False
    )
    
    lr_results[name] = lr_score

overlap = set(expr.index) & set(liana_lr['target'])
print("Overlap:", len(overlap))
