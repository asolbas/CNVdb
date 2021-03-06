import pandas as pd
import numpy as np

df = pd.read_csv('MAF_CNV_Database_PC.tsv', sep="\t", low_memory=False)


#remove X and Y chromosomes
df = df.loc[~((df['0']=="X") | (df['0']=="Y"))]

df=df.rename({'0':'chr', '1':'start', '2':'end'}, axis=1) #rename columns

#duplications
df_DUP=df[df['AC_DUP']>=1].filter(regex= "chr|start|end|^AN$|DUP")
df_DUP.insert(3,"SV_type","DUP")

#deletions
df_DEL=df[df['AC_DEL']>=1].filter(regex= "chr|start|end|^AN$|DEL")
df_DEL.insert(3,"SV_type","DEL")

DUP_file = 'DB_DUP_PC.tsv'
DEL_file = 'DB_DEL_PC.tsv'


df_DUP.to_csv(DUP_file, sep='\t', index=False)
df_DEL.to_csv(DEL_file, sep='\t', index=False)
