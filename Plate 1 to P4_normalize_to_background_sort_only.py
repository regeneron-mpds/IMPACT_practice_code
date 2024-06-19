#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib


# In[2]:


pd.set_option('display.max_rows', 200)


# In[3]:


sample_list= []


# In[4]:


sample_list_2=[]


# In[5]:


path= "(insert file path)"
sample='(insert sample name'


# In[6]:


counts_df=pd.read_csv(rf'file.txt', sep='\t', header=0, index_col=0, skiprows=1)


# In[7]:


counts_df=counts_df.drop(['Chr','Start', 'End', 'Strand', 'Length'], axis=1)


# In[8]:


counts_df['Total_Counts']=counts_df[f'(insert file path and name)'] + counts_df[f'(insert file path and name)']


# In[9]:


counts_df.drop([f'(insert file name)' , f'(insert file name)'], axis=1)


# In[10]:


List_of_samples_mapped_to_hACVR1=[]

for sample in sample_list:
    counts_df=pd.read_csv(rf'(insert text file)', sep='\t', header=0, index_col=0, skiprows=1)
    counts_df_2=counts_df.drop([add columns to drop], axis=1)
    counts_df_2[f'Total_Counts{sample}']=counts_df_2[f'(insert file name)'] + counts_df_2[f'(insert file name)']
    counts_df_3=counts_df_2.drop([f'(insert file name)' , f'(insert file name)'], axis=1)
    Total_Count=int(counts_df_3[f'Total_Counts{sample}'])
    List_of_samples_mapped_to_hACVR1.append(f'{sample}' + ',' + str(Total_Count))
    



# In[11]:


Mapped_S=pd.Series(List_of_samples_mapped_to_hACVR1)


# In[12]:


S_mapped_sep = Mapped_S.str.split(',', 1, expand=False)


# In[13]:


Sample_Name=S_mapped_sep.str[0]
Count_Mapped=S_mapped_sep.str[1]


# In[14]:


Sample_Counts_df=pd.DataFrame({'Sample Name':Sample_Name , 'Counts':Count_Mapped})


# In[15]:


Sample_Counts_df.set_index('Sample Name', drop=True)


# In[16]:


Sample_Counts_df['Counts'].loc[0]


# In[17]:


Sample_Counts_df['Counts']=Sample_Counts_df['Counts'].astype(int)


# In[18]:


Sample_Counts_df['Subtract from Raw Matches']=(20 * Sample_Counts_df['Counts'])/Sample_Counts_df['Counts'].loc[0]
Sample_Counts_df


# In[19]:


List_of_samples_mapped_to_hACVR1_2=[]

for sample in sample_list_2:
    counts_df_4=pd.read_csv(rf'(insert file name)', sep='\t', header=0, index_col=0, skiprows=1)
    counts_df_5=counts_df_4.drop([insert columns to drop], axis=1)
    counts_df_5[f'Total_Counts{sample}']=counts_df_5[f'(insert file name)'] + counts_df_5[f'(insert file name)']
    counts_df_6=counts_df_5.drop([f'(insert file name)' , f'(insert file name)'], axis=1)
    Total_Count=int(counts_df_6[f'Total_Counts{sample}'])
    List_of_samples_mapped_to_hACVR1_2.append(f'{sample}' + ',' + str(Total_Count))
    

# In[20]:


Mapped_S_2=pd.Series(List_of_samples_mapped_to_hACVR1_2)


# In[21]:


S_mapped_sep_2 = Mapped_S_2.str.split(',', 1, expand=False)


# In[22]:


Sample_Name_2=S_mapped_sep_2.str[0]
Count_Mapped_2=S_mapped_sep_2.str[1]


# In[23]:


Sample_Counts_df_2=pd.DataFrame({'Sample Name':Sample_Name_2 , 'Counts':Count_Mapped_2})


# In[24]:


Sample_Counts_df_2.set_index('Sample Name', drop=True)


# In[25]:


Sample_Counts_df_2['Counts']=Sample_Counts_df_2['Counts'].astype(int)


# In[26]:


Sample_Counts_df_2['Subtract from Raw Matches']=(20 * Sample_Counts_df_2['Counts'])/Sample_Counts_df_2['Counts'].loc[0]


# In[27]:


reference = SeqIO.read("(insert fast file name), "fasta")


# In[28]:


codontab = {
    'TCA': 'S',    # Serine
    'TCC': 'S',    # Serine
    'TCG': 'S',    # Serine
    'TCT': 'S',    # Serine
    'TTC': 'F',    # Phenylalanine
    'TTT': 'F',    # Phenylalanine
    'TTA': 'L',    # Leucine
    'TTG': 'L',    # Leucine
    'TAC': 'Y',    # Tyrosine
    'TAT': 'Y',    # Tyrosine
    'TAA': '*',    # Stop
    'TAG': '*',    # Stop
    'TGC': 'C',    # Cysteine
    'TGT': 'C',    # Cysteine
    'TGA': '*',    # Stop
    'TGG': 'W',    # Tryptophan
    'CTA': 'L',    # Leucine
    'CTC': 'L',    # Leucine
    'CTG': 'L',    # Leucine
    'CTT': 'L',    # Leucine
    'CCA': 'P',    # Proline
    'CCC': 'P',    # Proline
    'CCG': 'P',    # Proline
    'CCT': 'P',    # Proline
    'CAC': 'H',    # Histidine
    'CAT': 'H',    # Histidine
    'CAA': 'Q',    # Glutamine
    'CAG': 'Q',    # Glutamine
    'CGA': 'R',    # Arginine
    'CGC': 'R',    # Arginine
    'CGG': 'R',    # Arginine
    'CGT': 'R',    # Arginine
    'ATA': 'I',    # Isoleucine
    'ATC': 'I',    # Isoleucine
    'ATT': 'I',    # Isoleucine
    'ATG': 'M',    # Methionine
    'ACA': 'T',    # Threonine
    'ACC': 'T',    # Threonine
    'ACG': 'T',    # Threonine
    'ACT': 'T',    # Threonine
    'AAC': 'N',    # Asparagine
    'AAT': 'N',    # Asparagine
    'AAA': 'K',    # Lysine
    'AAG': 'K',    # Lysine
    'AGC': 'S',    # Serine
    'AGT': 'S',    # Serine
    'AGA': 'R',    # Arginine
    'AGG': 'R',    # Arginine
    'GTA': 'V',    # Valine
    'GTC': 'V',    # Valine
    'GTG': 'V',    # Valine
    'GTT': 'V',    # Valine
    'GCA': 'A',    # Alanine
    'GCC': 'A',    # Alanine
    'GCG': 'A',    # Alanine
    'GCT': 'A',    # Alanine
    'GAC': 'D',    # Aspartic acid
    'GAT': 'D',    # Aspartic acid
    'GAA': 'E',    # Aspartic acid
    'GAG': 'E',    # Aspartic acid
    'GGA': 'G',    # Glycine
    'GGC': 'G',    # Glycine
    'GGG': 'G',    # Glycine
    'GGT': 'G'     # Glycine
}


# In[29]:


WT_seqs_tuples_join = []
seq = str(reference.seq)
for position in range(438, len(seq)-3,3):
    bp = seq[position:position+3]
    WT_position = str(position)
    WT_codon= bp
    WT_AA=codontab[WT_codon]
    WT_seqs_tuples_add=(WT_position + WT_codon)
    WT_seqs_tuples_join.append(WT_seqs_tuples_add)


# In[30]:


rows_to_remove = WT_seqs_tuples_join


# In[31]:


S_rows_to_remove=pd.Series(rows_to_remove)
S_rows_to_remove_sep = S_rows_to_remove.str.split('(\d+)', 1, expand=False)
Positions_from_index_remove=S_rows_to_remove_sep.str[1]
Codons_from_index_remove=S_rows_to_remove_sep.str[2]


# In[32]:


list_of_WT_AA= []
position_WT=[]
AA_for_WT=[]
Codon_for_WT=[]
for i in range(len(S_rows_to_remove_sep)):
    codon = Codons_from_index_remove[i]
    list_value_i=(Positions_from_index_remove[i], codontab[codon])
    list_value_i_2=Positions_from_index_remove[i]
    list_value_i_3=codontab[codon]
    list_value_i_4=codon
    list_of_WT_AA.append(list_value_i)
    position_WT.append(list_value_i_2)
    AA_for_WT.append(list_value_i_3)
    Codon_for_WT.append(list_value_i_4)
list_of_WT_AA[0:5]


# ## Sample Title 1

# In[33]:


sample = '(insert sample name)'


# In[34]:


perfect_match_df = pd.read_csv("(insert file name)", sep=',', header=0, index_col=0, skiprows=0)
perfect_match_df.rename(columns = {f'0':'Matches_R1'}, inplace=True)
perfect_match_df['Matches_R2'] = pd.read_csv("(insert file name)", sep=',', header=0, index_col=0, skiprows=0)
perfect_match_df['Total_Matches'] = perfect_match_df['Matches_R1'] + perfect_match_df['Matches_R2']
perfect_match_df


# In[35]:


new_index_pre=perfect_match_df.index.difference(rows_to_remove)
perfect_match_without_WT_df=perfect_match_df.loc[new_index_pre]


# In[36]:


S_matches = pd.Series(perfect_match_without_WT_df.index, dtype="string")
S_matches_sep = S_matches.str.split('(\d+)', 1, expand=False)
Positions_from_index=S_matches_sep.str[1]
Codons_from_index=S_matches_sep.str[2]


# In[37]:


list_of_tuples_for_multiindex = []
for i in range(len(S_matches_sep)):
    codon = Codons_from_index[i]
    multi_index_value_i=(Positions_from_index[i], codon, codontab[codon])
    list_of_tuples_for_multiindex.append(multi_index_value_i)
new_index = pd.MultiIndex.from_tuples(list_of_tuples_for_multiindex)
perfect_match_without_WT_df.index = new_index
perfect_match_without_WT_df.drop(['Matches_R1','Matches_R2'], axis=1, inplace=True)


# In[38]:


positions_to_drop=perfect_match_without_WT_df.loc[(insert range to remove)]
L_positions_to_drop=list(positions_to_drop.index)


# In[39]:


new_index_pre_2=perfect_match_without_WT_df.index.difference(L_positions_to_drop)
Plate1_without_WT_df=perfect_match_without_WT_df.loc[new_index_pre_2]


# In[40]:


S_Total_Matches_pre_Plate1=pd.Series(Plate1_without_WT_df['Total_Matches'])


# In[41]:


Count_mutant=int(Plate1_without_WT_df['Total_Matches'].sum())
Count_mutant


# In[42]:


Plate1_without_WT_df['Total_Matches'].sum()


# In[43]:


Plate1_without_WT_df['Normalized Total_Matches_Pre-Sort']=Plate1_without_WT_df['Total_Matches']/Count_mutant
perfect_match_without_WT_normalized_df_1=Plate1_without_WT_df.drop(['Total_Matches'], axis=1)


# ## Sample Title 2

# In[44]:


sample_2 = '(insert sample name)'


# In[45]:


perfect_match_df_2 = pd.read_csv("(insert file name), sep=',', header=0, index_col=0, skiprows=0)
perfect_match_df_2.rename(columns = {f'0':'Matches_R1'}, inplace=True)
perfect_match_df_2['Matches_R2'] = pd.read_csv("(insert file name)", sep=',', header=0, index_col=0, skiprows=0)
perfect_match_df_2['Total_Matches'] = perfect_match_df_2['Matches_R1'] + perfect_match_df_2['Matches_R2']



# In[46]:


new_index_sort_No_AA=perfect_match_df_2.index.difference(rows_to_remove)
perfect_match_without_WT_df_2=perfect_match_df_2.loc[new_index_sort_No_AA]


# In[47]:


S_matches_2 = pd.Series(perfect_match_without_WT_df_2.index, dtype="string")
S_matches_sep_2 = S_matches_2.str.split('(\d+)', 1, expand=False)
Positions_from_index_2=S_matches_sep_2.str[1]
Codons_from_index_2=S_matches_sep_2.str[2]


# In[48]:


list_of_tuples_for_multiindex_2 = []
for i in range(len(S_matches_sep_2)):
    codon = Codons_from_index_2[i]
    multi_index_value_i_2=(Positions_from_index_2[i], codon, codontab[codon])
    list_of_tuples_for_multiindex_2.append(multi_index_value_i_2)
new_index_2 = pd.MultiIndex.from_tuples(list_of_tuples_for_multiindex_2)
perfect_match_without_WT_df_2.index = new_index_2
perfect_match_without_WT_df_2.drop(['Matches_R1','Matches_R2'], axis=1, inplace=True)


# In[49]:


new_index_sort_3=perfect_match_without_WT_df_2.index.difference(L_positions_to_drop)
Plate1_without_WT_df_2=perfect_match_without_WT_df_2.loc[new_index_sort_3]


# In[50]:


Plate1_without_WT_df_2['Normalize Raw Counts']=Plate1_without_WT_df_2['Total_Matches'] - Sample_Counts_df['Subtract from Raw Matches'].loc[2]


# In[51]:


Plate1_without_WT_df_2[Plate1_without_WT_df_2 < 0] = 0   


# In[52]:


S_Total_Matches_sort_Plate1=pd.Series(Plate1_without_WT_df_2['Normalize Raw Counts'])


# In[53]:


S_Total_Matches_sort_Plate1_2=pd.Series(Plate1_without_WT_df_2['Total_Matches'])


# In[54]:


S_Total_Matches_sort_Plate1_2


# In[55]:


Raw_Totals_Plate1_df=pd.DataFrame({'Pre-Sort':S_Total_Matches_pre_Plate1, 'Sort AA': S_Total_Matches_sort_Plate1})

# In[56]:


Plate1_without_WT_df_2.value_counts()


# In[57]:


Plate1_without_WT_df_2.sum()


# In[58]:


int(Plate1_without_WT_df_2['Normalize Raw Counts'].sum())


# In[59]:


Count_mutant_2=int(Plate1_without_WT_df_2['Normalize Raw Counts'].sum())


# In[60]:


Plate1_without_WT_df_2.sum()


# In[61]:


Plate1_without_WT_df_2['Normalized Total_Matches_Sort']=Plate1_without_WT_df_2['Normalize Raw Counts']/Count_mutant_2
perfect_match_without_WT_normalized_df_2=Plate1_without_WT_df_2.drop(['Total_Matches', 'Normalize Raw Counts'], axis=1)
perfect_match_without_WT_normalized_df_2


# ## Sample 3 Title

# In[62]:


sample_3='(insert sample name)'


# In[63]:


perfect_match_df_3 = pd.read_csv("(insert file name)", sep=',', header=0, index_col=0, skiprows=0)
perfect_match_df_3.rename(columns = {f'0':'Matches_R1'}, inplace=True)
perfect_match_df_3['Matches_R2'] = pd.read_csv("(insert sample name)", sep=',', header=0, index_col=0, skiprows=0)
perfect_match_df_3['Total_Matches'] = perfect_match_df_3['Matches_R1'] + perfect_match_df_3['Matches_R2']


# In[64]:


new_index_pre_3=perfect_match_df_3.index.difference(rows_to_remove)
perfect_match_without_WT_df_3=perfect_match_df_3.loc[new_index_pre_3]


# In[65]:


S_matches_3 = pd.Series(perfect_match_without_WT_df_3.index, dtype="string")
S_matches_sep_3 = S_matches_3.str.split('(\d+)', 1, expand=False)
Positions_from_index_3=S_matches_sep_3.str[1]
Codons_from_index_3=S_matches_sep_3.str[2]


# In[66]:


list_of_tuples_for_multiindex_3 = []
for i in range(len(S_matches_sep_3)):
    codon = Codons_from_index_3[i]
    multi_index_value_i_3=(Positions_from_index_3[i], codon, codontab[codon])
    list_of_tuples_for_multiindex_3.append(multi_index_value_i_3)
new_index_3 = pd.MultiIndex.from_tuples(list_of_tuples_for_multiindex_3)
perfect_match_without_WT_df_3.index = new_index_3
perfect_match_without_WT_df_3.drop(['Matches_R1','Matches_R2'], axis=1, inplace=True)


# In[67]:


positions_to_drop_2=perfect_match_without_WT_df_3.loc['(insert range of positions)']
positions_to_drop_3=perfect_match_without_WT_df_3.loc['(insert range of positions)']


# In[68]:


L_positions_to_drop_2=list(positions_to_drop_2.index)


# In[69]:


L_positions_to_drop_3=list(positions_to_drop_3.index)


# In[70]:


new_index_pre_3=perfect_match_without_WT_df_3.index.difference(L_positions_to_drop_2)
Plate2_without_WT_df=perfect_match_without_WT_df_3.loc[new_index_pre_3]


# In[71]:


new_index_pre_4=Plate2_without_WT_df.index.difference(L_positions_to_drop_3)
Plate2_without_WT_df=Plate2_without_WT_df.loc[new_index_pre_4]


# In[72]:


S_Total_Matches_pre_Plate2=pd.Series(Plate2_without_WT_df['Total_Matches'])


# In[73]:


Count_mutant_3=int(Plate2_without_WT_df['Total_Matches'].sum())


# In[74]:


Plate2_without_WT_df['Normalized Total_Matches_Pre-Sort']=Plate2_without_WT_df['Total_Matches']/Count_mutant_3
perfect_match_without_WT_normalized_df_3=Plate2_without_WT_df.drop(['Total_Matches'], axis=1)


# In[75]:


perfect_match_without_WT_normalized_df_3.value_counts()


# ## Sample Title 4

# In[76]:


sample_4='(Insert Sample Name)'


# In[77]:


perfect_match_df_4 = pd.read_csv("(insert file name)", sep=',', header=0, index_col=0, skiprows=0)
perfect_match_df_4.rename(columns = {f'0':'Matches_R1'}, inplace=True)
perfect_match_df_4['Matches_R2'] = pd.read_csv("(insert file name)", sep=',', header=0, index_col=0, skiprows=0)
perfect_match_df_4['Total_Matches'] = perfect_match_df_4['Matches_R1'] + perfect_match_df_4['Matches_R2']



# In[78]:


new_index_sort_4=perfect_match_df_4.index.difference(rows_to_remove)
perfect_match_without_WT_df_4=perfect_match_df_4.loc[new_index_sort_4]


# In[79]:


S_matches_4 = pd.Series(perfect_match_without_WT_df_4.index, dtype="string")
S_matches_sep_4 = S_matches_4.str.split('(\d+)', 1, expand=False)
Positions_from_index_4=S_matches_sep_4.str[1]
Codons_from_index_4=S_matches_sep_4.str[2]


# In[80]:


list_of_tuples_for_multiindex_4 = []
for i in range(len(S_matches_sep_4)):
    codon = Codons_from_index_4[i]
    multi_index_value_i_4=(Positions_from_index_4[i], codon, codontab[codon])
    list_of_tuples_for_multiindex_4.append(multi_index_value_i_4)
new_index_4 = pd.MultiIndex.from_tuples(list_of_tuples_for_multiindex_4)
perfect_match_without_WT_df_4.index = new_index_4
perfect_match_without_WT_df_4.drop(['Matches_R1','Matches_R2'], axis=1, inplace=True)


# In[81]:


new_index_sort_4=perfect_match_without_WT_df_4.index.difference(L_positions_to_drop_2)
Plate2_without_WT_df_2=perfect_match_without_WT_df_4.loc[new_index_sort_4]


# In[82]:


new_index_sort_5=Plate2_without_WT_df_2.index.difference(L_positions_to_drop_3)
Plate2_without_WT_df_2=Plate2_without_WT_df_2.loc[new_index_sort_5]


# In[83]:


Sample_Counts_df_2.loc[2]


# In[84]:


Plate2_without_WT_df_2['Normalize Raw Counts']=Plate2_without_WT_df_2['Total_Matches'] - Sample_Counts_df_2['Subtract from Raw Matches'].loc[2]


# In[85]:


Plate2_without_WT_df_2[Plate2_without_WT_df_2 < 0] = 0   


# In[86]:


S_Total_Matches_sort_Plate2=pd.Series(Plate2_without_WT_df_2['Normalize Raw Counts'])


# In[87]:


S_Total_Matches_sort_Plate2_2=pd.Series(Plate2_without_WT_df_2['Total_Matches'])


# In[88]:


Raw_Totals_Plate2_df=pd.DataFrame({'Pre-Sort':S_Total_Matches_pre_Plate2, 'Sort AA': S_Total_Matches_sort_Plate2})


# In[89]:


Raw_Totals_Plate2_df


# In[90]:


Plate2_without_WT_df_2.sum()


# In[91]:


Count_mutant_4=int(Plate2_without_WT_df_2['Normalize Raw Counts'].sum())


# In[92]:


Plate2_without_WT_df_2['Normalized Total_Matches_Sort']=Plate2_without_WT_df_2['Normalize Raw Counts']/Count_mutant_4
perfect_match_without_WT_normalized_df_4=Plate2_without_WT_df_2.drop(['Total_Matches', 'Normalize Raw Counts'], axis=1)


# ## Sample 5 Title

# In[93]:


sample_5='(insert sample)'


# In[94]:


perfect_match_df_5 = pd.read_csv("(insert file name)", sep=',', header=0, index_col=0, skiprows=0)
perfect_match_df_5.rename(columns = {f'0':'Matches_R1'}, inplace=True)
perfect_match_df_5['Matches_R2'] = pd.read_csv("(insert file name)", sep=',', header=0, index_col=0, skiprows=0)
perfect_match_df_5['Total_Matches'] = perfect_match_df_5['Matches_R1'] + perfect_match_df_5['Matches_R2']

# In[95]:

new_index_pre_5=perfect_match_df_5.index.difference(rows_to_remove)
perfect_match_without_WT_df_5=perfect_match_df_5.loc[new_index_pre_5]

# In[96]:


S_matches_5 = pd.Series(perfect_match_without_WT_df_5.index, dtype="string")
S_matches_sep_5 = S_matches_5.str.split('(\d+)', 1, expand=False)
Positions_from_index_5=S_matches_sep_5.str[1]
Codons_from_index_5=S_matches_sep_5.str[2]


# In[97]:


list_of_tuples_for_multiindex_5 = []
for i in range(len(S_matches_sep_5)):
    codon = Codons_from_index_5[i]
    multi_index_value_i_5=(Positions_from_index_5[i], codon, codontab[codon])
    list_of_tuples_for_multiindex_5.append(multi_index_value_i_5)
new_index_5 = pd.MultiIndex.from_tuples(list_of_tuples_for_multiindex_5)
perfect_match_without_WT_df_5.index = new_index_5
perfect_match_without_WT_df_5.drop(['Matches_R1','Matches_R2'], axis=1, inplace=True)


# In[98]:


positions_to_drop_4=perfect_match_without_WT_df_5.loc['(insert range)']
positions_to_drop_5=perfect_match_without_WT_df_5.loc['(insert range)']
positions_to_drop_6=perfect_match_without_WT_df_5.loc['(insert range)']


# In[99]:


L_positions_to_drop_4=list(positions_to_drop_4.index)
L_positions_to_drop_5=list(positions_to_drop_5.index)
L_positions_to_drop_6=list(positions_to_drop_6.index)


# In[100]:


new_index_pre_6=perfect_match_without_WT_df_5.index.difference(L_positions_to_drop_4)
Plate3_without_WT_df=perfect_match_without_WT_df_5.loc[new_index_pre_6]

# In[101]:


new_index_pre_7=Plate3_without_WT_df.index.difference(L_positions_to_drop_5)
Plate3_without_WT_df=Plate3_without_WT_df.loc[new_index_pre_7]

# In[102]:


new_index_pre_8=Plate3_without_WT_df.index.difference(L_positions_to_drop_6)
Plate3_without_WT_df=Plate3_without_WT_df.loc[new_index_pre_8]

# In[103]:

S_Total_Matches_pre_Plate3=pd.Series(Plate3_without_WT_df['Total_Matches'])

# In[104]:

Count_mutant_5=int(Plate3_without_WT_df['Total_Matches'].sum())

# In[105]:

Plate3_without_WT_df['Normalized Total_Matches_Pre-Sort']=Plate3_without_WT_df['Total_Matches']/Count_mutant_5
perfect_match_without_WT_normalized_df_5=Plate3_without_WT_df.drop(['Total_Matches'], axis=1)

# ## Sample 6 Title

# In[106]:


sample_6='(insert sample)'


# In[107]:


perfect_match_df_6 = pd.read_csv("(insert file name)", sep=',', header=0, index_col=0, skiprows=0)
perfect_match_df_6.rename(columns = {f'0':'Matches_R1'}, inplace=True)
perfect_match_df_6['Matches_R2'] = pd.read_csv("/d(insert file name)", sep=',', header=0, index_col=0, skiprows=0)
perfect_match_df_6['Total_Matches'] = perfect_match_df_6['Matches_R1'] + perfect_match_df_6['Matches_R2']

# In[108]:

new_index_sort_6=perfect_match_df_6.index.difference(rows_to_remove)
perfect_match_without_WT_df_6=perfect_match_df_6.loc[new_index_sort_6]

# In[109]:

S_matches_6 = pd.Series(perfect_match_without_WT_df_6.index, dtype="string")
S_matches_sep_6 = S_matches_6.str.split('(\d+)', 1, expand=False)
Positions_from_index_6=S_matches_sep_6.str[1]
Codons_from_index_6=S_matches_sep_6.str[2]


# In[110]:


list_of_tuples_for_multiindex_6 = []
for i in range(len(S_matches_sep_6)):
    codon = Codons_from_index_6[i]
    multi_index_value_i_6=(Positions_from_index_6[i], codon, codontab[codon])
    list_of_tuples_for_multiindex_6.append(multi_index_value_i_6)
new_index_6 = pd.MultiIndex.from_tuples(list_of_tuples_for_multiindex_6)
perfect_match_without_WT_df_6.index = new_index_6
perfect_match_without_WT_df_6.drop(['Matches_R1','Matches_R2'], axis=1, inplace=True)


# In[111]:


new_index_sort_7=perfect_match_without_WT_df_6.index.difference(L_positions_to_drop_4)
Plate3_without_WT_df_2=perfect_match_without_WT_df_6.loc[new_index_sort_7]


# In[112]:


new_index_sort_8=Plate3_without_WT_df_2.index.difference(L_positions_to_drop_5)
Plate3_without_WT_df_2=Plate3_without_WT_df_2.loc[new_index_sort_8]


# In[113]:


new_index_sort_9=Plate3_without_WT_df_2.index.difference(L_positions_to_drop_6)
Plate3_without_WT_df_2=Plate3_without_WT_df_2.loc[new_index_sort_9]


# In[114]:


Sample_Counts_df_2.loc[4]


# In[115]:


Plate3_without_WT_df_2['Normalize Raw Counts']=Plate3_without_WT_df_2['Total_Matches'] - Sample_Counts_df_2['Subtract from Raw Matches'].loc[4]

# In[116]:


Plate3_without_WT_df_2[Plate3_without_WT_df_2 < 0] = 0    


# In[117]:


S_Total_Matches_sort_Plate3=pd.Series(Plate3_without_WT_df_2['Normalize Raw Counts'])


# In[118]:


S_Total_Matches_sort_Plate3_2=pd.Series(Plate3_without_WT_df_2['Total_Matches'])


# In[119]:


Raw_Totals_Plate3_df=pd.DataFrame({'Pre-Sort':S_Total_Matches_pre_Plate3, 'Sort AA': S_Total_Matches_sort_Plate3})


# In[120]:


Plate3_without_WT_df_2.sum()


# In[121]:


Count_mutant_6=int(Plate3_without_WT_df_2['Normalize Raw Counts'].sum())


# In[122]:


Plate3_without_WT_df_2['Normalized Total_Matches_Sort']=Plate3_without_WT_df_2['Normalize Raw Counts']/Count_mutant_6
perfect_match_without_WT_normalized_df_6=Plate3_without_WT_df_2.drop(['Total_Matches', 'Normalize Raw Counts'], axis=1)


# ## Sample 7 Title

# In[123]:


sample_7='(insert sample name)'


# In[124]:


perfect_match_df_7 = pd.read_csv("(insert file name)", sep=',', header=0, index_col=0, skiprows=0)
perfect_match_df_7.rename(columns = {f'0':'Matches_R1'}, inplace=True)
perfect_match_df_7['Matches_R2'] = pd.read_csv("(insert sample name)", sep=',', header=0, index_col=0, skiprows=0)
perfect_match_df_7['Total_Matches'] = perfect_match_df_7['Matches_R1'] + perfect_match_df_7['Matches_R2']


# In[125]:


new_index_pre_9=perfect_match_df_7.index.difference(rows_to_remove)
perfect_match_without_WT_df_7=perfect_match_df_7.loc[new_index_pre_9]


# In[126]:


S_matches_7 = pd.Series(perfect_match_without_WT_df_7.index, dtype="string")
S_matches_sep_7 = S_matches_7.str.split('(\d+)', 1, expand=False)
Positions_from_index_7=S_matches_sep_7.str[1]
Codons_from_index_7=S_matches_sep_7.str[2]


# In[127]:


Codons_from_index_7


# In[128]:


list_of_tuples_for_multiindex_7 = []
for i in range(len(S_matches_sep_7)):
    codon = Codons_from_index_7[i]
    multi_index_value_i_7=(Positions_from_index_7[i], codon, codontab[codon])
    list_of_tuples_for_multiindex_7.append(multi_index_value_i_7)
new_index_7 = pd.MultiIndex.from_tuples(list_of_tuples_for_multiindex_7)
perfect_match_without_WT_df_7.index = new_index_7
perfect_match_without_WT_df_7.drop(['Matches_R1','Matches_R2'], axis=1, inplace=True)


# In[129]:


positions_to_drop_7=perfect_match_without_WT_df_7.loc['(insert range)']
positions_to_drop_8=perfect_match_without_WT_df_7.loc['(insert range)']
positions_to_drop_9=perfect_match_without_WT_df_7.loc['(insert range)']


# In[130]:


L_positions_to_drop_7=list(positions_to_drop_7.index)
L_positions_to_drop_8=list(positions_to_drop_8.index)
L_positions_to_drop_9=list(positions_to_drop_9.index)


# In[131]:


new_index_pre_10=perfect_match_without_WT_df_7.index.difference(L_positions_to_drop_7)
Plate4_without_WT_df=perfect_match_without_WT_df_7.loc[new_index_pre_10]


# In[132]:


new_index_pre_11=Plate4_without_WT_df.index.difference(L_positions_to_drop_8)
Plate4_without_WT_df=Plate4_without_WT_df.loc[new_index_pre_11]


# In[133]:


new_index_pre_12=Plate4_without_WT_df.index.difference(L_positions_to_drop_9)
Plate4_without_WT_df=Plate4_without_WT_df.loc[new_index_pre_12]


# In[134]:


S_Total_Matches_pre_Plate4=pd.Series(Plate4_without_WT_df['Total_Matches'])


# In[135]:


Count_mutant_7=int(Plate4_without_WT_df['Total_Matches'].sum())


# In[136]:


Plate4_without_WT_df['Normalized Total_Matches_Pre-Sort']=Plate4_without_WT_df['Total_Matches']/Count_mutant_7
perfect_match_without_WT_normalized_df_7=Plate4_without_WT_df.drop(['Total_Matches'], axis=1)


# ## Sample 8 Title

# In[137]:


sample_8='(insert sample name)'


# In[138]:


perfect_match_df_8 = pd.read_csv("(insert file name)", sep=',', header=0, index_col=0, skiprows=0)
perfect_match_df_8.rename(columns = {f'0':'Matches_R1'}, inplace=True)
perfect_match_df_8['Matches_R2'] = pd.read_csv("(insert file name)", sep=',', header=0, index_col=0, skiprows=0)
perfect_match_df_8['Total_Matches'] = perfect_match_df_8['Matches_R1'] + perfect_match_df_8['Matches_R2']


# In[139]:


new_index_sort_10=perfect_match_df_8.index.difference(rows_to_remove)
perfect_match_without_WT_df_8=perfect_match_df_8.loc[new_index_sort_10]


# In[140]:


S_matches_8 = pd.Series(perfect_match_without_WT_df_8.index, dtype="string")
S_matches_sep_8 = S_matches_8.str.split('(\d+)', 1, expand=False)
Positions_from_index_8=S_matches_sep_8.str[1]
Codons_from_index_8=S_matches_sep_8.str[2]


# In[141]:


list_of_tuples_for_multiindex_8 = []
for i in range(len(S_matches_sep_8)):
    codon = Codons_from_index_8[i]
    multi_index_value_i_8=(Positions_from_index_8[i], codon, codontab[codon])
    list_of_tuples_for_multiindex_8.append(multi_index_value_i_8)
new_index_8 = pd.MultiIndex.from_tuples(list_of_tuples_for_multiindex_8)
perfect_match_without_WT_df_8.index = new_index_8
perfect_match_without_WT_df_8.drop(['Matches_R1','Matches_R2'], axis=1, inplace=True)


# In[142]:


new_index_sort_11=perfect_match_without_WT_df_8.index.difference(L_positions_to_drop_7)
Plate4_without_WT_df_2=perfect_match_without_WT_df_8.loc[new_index_sort_11]


# In[143]:


new_index_sort_12=Plate4_without_WT_df_2.index.difference(L_positions_to_drop_8)
Plate4_without_WT_df_2=Plate4_without_WT_df_2.loc[new_index_sort_12]


# In[144]:


new_index_sort_13=Plate4_without_WT_df_2.index.difference(L_positions_to_drop_9)
Plate4_without_WT_df_2=Plate4_without_WT_df_2.loc[new_index_sort_13]


# In[145]:


Sample_Counts_df_2.loc[6]


# In[146]:


Plate4_without_WT_df_2['Normalize Raw Counts']=Plate4_without_WT_df_2['Total_Matches'] - Sample_Counts_df_2['Subtract from Raw Matches'].loc[6]


# In[147]:


Plate4_without_WT_df_2[Plate4_without_WT_df_2 < 0] = 0   


# In[148]:


S_Total_Matches_sort_Plate4=pd.Series(Plate4_without_WT_df_2['Normalize Raw Counts'])


# In[149]:


S_Total_Matches_sort_Plate4_2=pd.Series(Plate4_without_WT_df_2['Total_Matches'])


# In[150]:


Raw_Totals_Plate4_df=pd.DataFrame({'Pre-Sort':S_Total_Matches_pre_Plate4,'Sort AA': S_Total_Matches_sort_Plate4})


# In[151]:


Plate4_without_WT_df_2.sum()


# In[152]:


Count_mutant_8=int(Plate4_without_WT_df_2['Normalize Raw Counts'].sum())


# In[153]:


Plate4_without_WT_df_2['Normalized Total_Matches_Sort']=Plate4_without_WT_df_2['Normalize Raw Counts']/Count_mutant_8
perfect_match_without_WT_normalized_df_8=Plate4_without_WT_df_2.drop(['Total_Matches', 'Normalize Raw Counts'], axis=1)

# ## Analysis

# In[154]:


Raw_Totals_df=pd.concat([Raw_Totals_Plate1_df, Raw_Totals_Plate2_df, Raw_Totals_Plate3_df, Raw_Totals_Plate4_df], ignore_index=False, axis=0)


# In[155]:


Raw_Totals_df['Sort AA'].value_counts()


# In[156]:


Raw_Totals_by_AA=Raw_Totals_df.groupby(level=[0,2]).sum()
Raw_Totals_by_AA['Pre-Sort'].value_counts()


# In[157]:


Raw_Totals_by_AA['Sort AA'].value_counts()


# In[158]:


int(Raw_Totals_by_AA['Sort AA'].sum())


# In[159]:


Raw_Totals_by_AA['Pre-Sort'].sum()


# In[160]:


plt.figure(figsize=(8,8))
plot=sns.histplot(data=Raw_Totals_df,binwidth=1, kde=True, element='step', stat='density')
plot.set_xlabel('Raw Totals', fontsize='x-large')
plot.set_ylabel('Density', fontsize='x-large')
plt.ylim(0,0.01)
plt.xlim(-1,100)


# In[161]:


Raw_Totals_df.replace(0, 0.5, inplace=True)


# In[162]:


log_Raw_Totals_df=Raw_Totals_df.applymap(np.log2)


# In[163]:


log_Raw_Totals_df['Position']=log_Raw_Totals_df.index.get_level_values(0)
log_Raw_Totals_df


# In[ ]:





# In[164]:


log_Raw_Totals_df['Position']=log_Raw_Totals_df['Position'].astype(int)
log_Raw_Totals_df.sort_values(by=['Position'], ascending=True, inplace=True)


# In[165]:


log_Raw_Totals_df['Position']=log_Raw_Totals_df['Position'].astype(str)


# In[166]:


Plate1_df=perfect_match_without_WT_normalized_df_1.join(perfect_match_without_WT_normalized_df_2)


# In[167]:


Plate2_df=perfect_match_without_WT_normalized_df_3.join(perfect_match_without_WT_normalized_df_4)


# In[168]:


Plate3_df=perfect_match_without_WT_normalized_df_5.join(perfect_match_without_WT_normalized_df_6)


# In[169]:


Plate4_df=perfect_match_without_WT_normalized_df_7.join(perfect_match_without_WT_normalized_df_8)


# In[170]:


Analysis_df=pd.concat([Plate1_df, Plate2_df, Plate3_df, Plate4_df], ignore_index=False, axis=0)


# In[171]:


Analysis_df['Normalized Total_Matches_Pre-Sort'].value_counts()


# In[172]:


Analysis_df['Normalized Total_Matches_Sort'].value_counts()


# In[173]:


Analysis_df_for_log=Analysis_df.replace(0, 0.000001, inplace=False)
Analysis_df_for_log.replace(np.inf, 0.3, inplace=True) 
log2_Analysis_df=Analysis_df_for_log.applymap(np.log2)


# In[174]:


log2_Analysis_df['Position']=log2_Analysis_df.index.get_level_values(0)


# In[175]:


log2_Analysis_df['Position']=log2_Analysis_df['Position'].astype(int)
log2_Analysis_df.sort_values(by=['Position'], ascending=True, inplace=True)


# In[176]:


log2_Analysis_df['Position']=log2_Analysis_df['Position'].astype(str)


# In[177]:


Analysis_df['Fold Change Sort AA Matches'] = Analysis_df['Normalized Total_Matches_Sort'] /Analysis_df['Normalized Total_Matches_Pre-Sort']
Analysis_df_fold_change=Analysis_df.drop(['Normalized Total_Matches_Sort','Normalized Total_Matches_Pre-Sort'], axis=1, inplace=False)


# In[178]:


Analysis_df_fold_change['Fold Change Sort AA Matches'].value_counts()


# In[179]:


Analysis_df_fold_change.isna().sum()


# In[180]:


Analysis_df_fold_change.replace(0, 0.001, inplace=True)
Analysis_df_fold_change.replace(np.inf, 1000, inplace=True)
log2_Analysis_df_fold_change=Analysis_df_fold_change.applymap(np.log2)


# In[181]:


log2_Analysis_df_fold_change['Position']=log2_Analysis_df_fold_change.index.get_level_values(0)
log2_Analysis_df_fold_change


# In[182]:


log2_Analysis_df_fold_change['Position']=log2_Analysis_df_fold_change['Position'].astype(int)
log2_Analysis_df_fold_change.sort_values(by=['Position'], ascending=True, inplace=True)


# In[183]:


log2_Analysis_df_fold_change['Position']=log2_Analysis_df_fold_change['Position'].astype(str)


# In[184]:


Analysis_df_index_list=list(log2_Analysis_df_fold_change.index)


# In[185]:


List_Analysis=[]
S_Analysis_df_index=pd.Series(Analysis_df_index_list)
position_from_index=S_Analysis_df_index.str[0]
AA_from_index=S_Analysis_df_index.str[2]

for p in range(len(S_Analysis_df_index)):
    value_p=(position_from_index[p], AA_from_index[p])
    List_Analysis.append(value_p)


# In[186]:


WT_AA_loop= []
for x in list_of_WT_AA:
    WT_AA_loop.append(x)


# In[187]:


log2_Analysis_df


# In[188]:


legend2=[]
for i in List_Analysis:
    if i  == ('585', 'P'):
        legend2.append('Known mutant')
    elif i == ('603', 'I'):
        legend2.append('Known mutant')
    elif i == ('615', 'H'):
        legend2.append('Known mutant')
    elif i == ('618', 'E'):
        legend2.append('Known mutant')
    elif i == ('618', 'D'):
        legend2.append('Known mutant')
    elif i  == ('972', 'A'):
        legend2.append('Known mutant')
    elif i == ('981', 'E'):
        legend2.append('Known mutant')
    elif i == ('981', 'W'):
        legend2.append('Known mutant')
    elif i == ('981', 'R'):
        legend2.append('Known mutant')
    elif i == ('771', 'S'):
        legend2.append('Known mutant') 
    elif i == ('771', 'G'):
        legend2.append('Known mutant')
    elif i == ('1122', 'P'):
        legend2.append('Known mutant') 
    elif i == ('1065', 'D'):
        legend2.append('Known mutant')
    elif i[1] == ('*'):
        legend2.append('Stop_Codon')
    elif i in WT_AA_loop:
        legend2.append('WT_AA')
    else:
        legend2.append('unknown')


# In[189]:


legend2.count('WT_AA')


# In[190]:


size_legend_2=[]
for legend in legend2:
    if legend == ('Known mutant'):
        size_legend_2.append(40)
    if legend == ('Stop_Codon'):
        size_legend_2.append(10)
    if legend == ('WT_AA'):
        size_legend_2.append(10)
    if legend == ('unknown'):
        size_legend_2.append(5)


# In[191]:


log_Raw_Totals_df.loc['insert range']


# In[192]:


Remove_from_Analysis=log_Raw_Totals_df.loc[(log_Raw_Totals_df['Pre-Sort'] <=2) & (log_Raw_Totals_df['Pre-Sort'] >-1)]


# In[193]:


len(Remove_from_Analysis.loc[(Remove_from_Analysis['Sort AA'] > -1)])


# In[194]:


palette={'unknown':'#999999', 'WT_AA':'#377eb8', 'Stop_Codon':'#e41a1c' , 'Known mutant':'#984ea3'}
plot=sns.JointGrid(data=log_Raw_Totals_df, x='Pre-Sort', y='Sort AA', height=10, ratio=2)
plot.plot_joint(sns.scatterplot, hue=legend2, palette=palette, size=size_legend_2, sizes=(30,100))
plot.plot_marginals(sns.histplot, hue=legend2, palette=palette, multiple='layer',hue_order=['Known mutant','Stop_Codon', 'WT_AA', 'unknown'] ,binwidth=1, kde=True, element='step', stat='density')
plot.ax_joint.axvline(x=-1, color='red', linestyle='dashed')
plot.ax_joint.axhline(y=-1, color='red', linestyle='dashed')
#plot.ax_joint.axvline(x=2, color='black', linestyle='dashed')
plot.ax_joint.axhline(y=6.6, color='black', linestyle='dashed')
#plot.ax_joint.set_xlim(-22,0)
#plot.ax_joint.set_ylim(-22,0)
plot.ax_joint.plot([-1.5,13], [-1.5,13], 'b-', linewidth = 2)
plot.ax_joint.set_xlabel('Log2 Raw Total_Matches_Pre-Sort', fontsize='x-large')
plot.ax_joint.set_ylabel('Log2 Raw Total_Matches_Sort_AA', fontsize='x-large')
plot.ax_joint.annotate('insert name', (coordinates), xytext=(text coordinates))

handles=[]
for label, color in palette.items():
    handle=matplotlib.patches.Patch(color=color, label=label)
    handles.append(handle)
plot.ax_joint.legend(handles=handles)
sns.move_legend(plot.ax_joint, loc=(1.2,1))


# In[195]:



# In[196]:


len(log2_Analysis_df.loc[(log2_Analysis_df['Normalized Total_Matches_Pre-Sort'] <=-16) & (log2_Analysis_df['Normalized Total_Matches_Pre-Sort'] >-19.9)& (log2_Analysis_df['Normalized Total_Matches_Sort'] >-19.9)])


# In[197]:


palette={'unknown':'#999999', 'WT_AA':'#377eb8', 'Stop_Codon':'#e41a1c' , 'Known mutant':'#984ea3'}
plot=sns.JointGrid(data=log2_Analysis_df, x='Normalized Total_Matches_Pre-Sort', y='Normalized Total_Matches_Sort', height=10, ratio=2)
plot.plot_joint(sns.scatterplot, hue=legend2, palette=palette, size=size_legend_2, sizes=(30,100))
plot.plot_marginals(sns.histplot, hue=legend2, palette=palette, multiple='layer',hue_order=['Known mutant','Stop_Codon', 'WT_AA', 'unknown'] ,binwidth=1, kde=True, element='step', stat='density')
plot.ax_joint.axvline(x=-19.9, color='red', linestyle='dashed')
plot.ax_joint.axhline(y=-19.9, color='red', linestyle='dashed')
#plot.ax_joint.ayvline(x=-16, color='black', linestyle='dashed')
plot.ax_joint.set_xlim(-22,0)
plot.ax_joint.set_ylim(-22,0)
plot.ax_joint.plot([-22,0], [-22,-0], 'b-', linewidth = 2)
plot.ax_joint.set_xlabel('Log2 Normalized Total_Matches_Pre-Sort', fontsize='x-large')
plot.ax_joint.set_ylabel('Log2 Normalized Total_Matches_Sort_AA', fontsize='x-large')
handles=[]
for label, color in palette.items():
    handle=matplotlib.patches.Patch(color=color, label=label)
    handles.append(handle)
plot.ax_joint.legend(handles=handles)
sns.move_legend(plot.ax_joint, loc=(1.2,1))
plot.ax_joint.annotate


# In[198]:


# Category=[]
# for i in log2_Analysis_df_fold_change.index:
#     if i  == ('585', 'P'):
#         Category.append('Known mutant')
#     elif i == ('603', 'I'):
#         Category.append('Known mutant')
#     elif i == ('615', 'H'):
#         Category.append('Known mutant')
#     elif i == ('618', 'E'):
#         Category.append('Known mutant')
#     elif i == ('618', 'D'):
#         Category.append('Known mutant')
#     elif i  == ('972', 'A'):
#         Category.append('Known mutant')
#     elif i == ('981', 'E'):
#         Category.append('Known mutant')
#     elif i == ('981', 'W'):
#         Category.append('Known mutant')
#     elif i == ('981', 'R'):
#         Category.append('Known mutant')
#     elif i == ('771', 'S'):
#         Category.append('Known mutant') 
#     elif i == ('771', 'G'):
#         Category.append('Known mutant')
#     elif i == ('1122', 'P'):
#         Category.append('Known mutant') 
#     elif i == ('1065', 'D'):
#         Category.append('Known mutant')
#     elif i[1] == ('*'):
#         Category.append('Stop_Codon')
#     elif i in WT_AA_loop:
#         Category.append('WT_AA')
#     else:
#         Category.append('unknown')


# In[199]:

# In[200]:

# In[201]:


Category_S=pd.Series(legend2)


# In[202]:


Category_S.index = log2_Analysis_df_fold_change.index


# In[203]:


log2_Analysis_df_fold_change.index


# In[204]:


# In[205]:

# In[206]:


log2_Analysis_df_fold_change_category_df=pd.concat([log2_Analysis_df_fold_change, Category_S], axis=1)


# In[207]:


log2_Analysis_df_fold_change_category_df.rename(columns = {0:'Category'}, inplace=True)


# In[208]:


# In[209]:


unknown_df=log2_Analysis_df_fold_change_category_df[log2_Analysis_df_fold_change_category_df['Category'].str.contains('unknown')]


# In[210]:


WT_AA_df=log2_Analysis_df_fold_change_category_df[log2_Analysis_df_fold_change_category_df['Category'].str.contains('WT_AA')]


# In[211]:

# In[212]:


len(WT_AA_df.loc[WT_AA_df['Fold Change Sort AA Matches'] >= 2])


# In[213]:


len(WT_AA_df.loc[WT_AA_df['Fold Change Sort AA Matches'] >= 4.5])


# In[214]:


Known_mutants_df=log2_Analysis_df_fold_change_category_df[log2_Analysis_df_fold_change_category_df['Category'].str.contains('Known mutant')]


# In[215]:


Known_mutants_df.describe()


# In[216]:


stop_codons_df=log2_Analysis_df_fold_change_category_df[log2_Analysis_df_fold_change_category_df['Category'].str.contains('Stop_Codon')]
len(stop_codons_df)


# In[217]:


palette6={'WT_AA':'#377eb8' , 'Known mutant':'#984ea3', 'unknown': '#999999'}
A=Known_mutants_df
B=stop_codons_df
C=WT_AA_df
D=unknown_df
plt.figure(figsize=(12,8))
plot=sns.histplot(data=Known_mutants_df, x='Fold Change Sort AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density', legend=True, color='#984ea3')
#plot=sns.histplot(data=B, x='Fold Change Sort AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density',legend=True, color='#e41a1c')
plot=sns.histplot(data=C, x='Fold Change Sort AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density',legend=True, color='#377eb8')
plot=sns.histplot(data=D, x='Fold Change Sort AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density',legend=True, color='#999999')
plot.set_ylabel('Density', fontsize='x-large')
plot.set_xlabel('Log2 Fold Change Sort/Pre-sort', fontsize='x-large')
plt.ylim(0,0.4)
plt.xlim(-11,11)
plot.axvline(x=4, color='black', linestyle='dashed')
plot.axvline(x=2, color='black')
handles=[]
for label, color in palette6.items():
    handle=matplotlib.patches.Patch(color=color, label=label)
    handles.append(handle)

plot.legend(handles=handles,fontsize='xx-large')
sns.move_legend(plot, "best")


# In[218]:

# In[219]:


len(unknown_df.loc[unknown_df['Fold Change Sort AA Matches'] >= 4.5])


# In[220]:


len(unknown_df.loc[(unknown_df['Fold Change Sort AA Matches'] >= 4) & (unknown_df['Fold Change Sort AA Matches'] <9.96)])


# In[221]:


len(unknown_df.loc[(unknown_df['Fold Change Sort AA Matches'] >= 2) & (unknown_df['Fold Change Sort AA Matches'] <9.96)])


# In[222]:


len(unknown_df.loc[(unknown_df['Fold Change Sort AA Matches'] >= 9.96)])


# In[223]:


len(unknown_df.loc[unknown_df['Fold Change Sort AA Matches'] >= 4])


# In[224]:


unknown_4_and_above=unknown_df.loc[unknown_df['Fold Change Sort AA Matches'] >= 4]


# In[225]:


unknown_4_to_9=unknown_df.loc[(unknown_df['Fold Change Sort AA Matches'] >= 4) & (unknown_df['Fold Change Sort AA Matches'] <9.96)]


# In[226]:


unknown_4_to_9.sort_values(by='Fold Change Sort AA Matches', ascending=False)


# In[227]:


palette2={'WT_AA':'#377eb8' , 'Known mutant':'#984ea3'}
A=Known_mutants_df
B=stop_codons_df
C=WT_AA_df
D=unknown_df
plt.figure(figsize=(12,8))
plot=sns.histplot(data=Known_mutants_df, x='Fold Change Sort AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density', legend=True, color='#984ea3')
#plot=sns.histplot(data=B, x='Fold Change Sort AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density',legend=True, color='#e41a1c')
plot=sns.histplot(data=C, x='Fold Change Sort AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density',legend=True, color='#377eb8')
#plot=sns.histplot(data=D, x='Fold Change Sort AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density',legend=True, color='#999999')
plot.set_ylabel('Density', fontsize='x-large')
plot.set_xlabel('Log2 Fold Change Sort/Pre-sort', fontsize='x-large')
plt.ylim(0,0.4)
plt.xlim(-11,11)
plot.axvline(x=4, color='black', linestyle='dashed')
#plot.axvline(x=2, color='black')
handles=[]
for label, color in palette2.items():
    handle=matplotlib.patches.Patch(color=color, label=label)
    handles.append(handle)

plot.legend(handles=handles,fontsize='xx-large')
sns.move_legend(plot, "best")


# In[228]:


palette3={'Known mutant':'#984ea3'}
A=Known_mutants_df
B=stop_codons_df
C=WT_AA_df
D=unknown_df
plt.figure(figsize=(12,8))
plot=sns.histplot(data=Known_mutants_df, x='Fold Change Sort AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density', legend=True, color='#984ea3')
#plot=sns.histplot(data=B, x='Fold Change Sort AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density',legend=True, color='#e41a1c')
#plot=sns.histplot(data=C, x='Fold Change Sort AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density',legend=True, color='#377eb8')
#plot=sns.histplot(data=D, x='Fold Change Sort AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density',legend=True, color='#999999')
plot.set_ylabel('Density', fontsize='x-large')
plot.set_xlabel('Log2 Fold Change Sort/Pre-sort', fontsize='x-large')
plt.ylim(0,0.4)
plt.xlim(-11,10)

handles=[]
for label, color in palette3.items():
    handle=matplotlib.patches.Patch(color=color, label=label)
    handles.append(handle)

plot.legend(handles=handles,fontsize='xx-large')
sns.move_legend(plot, "best")


# In[229]:


palette4={'WT_AA':'#377eb8'}
A=Known_mutants_df
B=stop_codons_df
C=WT_AA_df
D=unknown_df
plt.figure(figsize=(12,8))
#plot=sns.histplot(data=Known_mutants_df, x='Fold Change Sort AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density', legend=True, color='#984ea3')
#plot=sns.histplot(data=B, x='Fold Change Sort AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density',legend=True, color='#e41a1c')
plot=sns.histplot(data=C, x='Fold Change Sort AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density',legend=True, color='#377eb8')
#plot=sns.histplot(data=D, x='Fold Change Sort AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density',legend=True, color='#999999')
plot.set_ylabel('Density', fontsize='x-large')
plot.set_xlabel('Log2 Fold Change Sort/Pre-sort', fontsize='x-large')
plt.ylim(0,0.4)
plt.xlim(-11,10)

handles=[]
for label, color in palette4.items():
    handle=matplotlib.patches.Patch(color=color, label=label)
    handles.append(handle)

plot.legend(handles=handles,fontsize='xx-large')
sns.move_legend(plot, "best")


# In[230]:

# In[231]:


palette5={'unknown': '#999999'}
A=Known_mutants_df
B=stop_codons_df
C=WT_AA_df
D=unknown_df
plt.figure(figsize=(12,8))
#plot=sns.histplot(data=Known_mutants_df, x='Fold Change Sort AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density', legend=True, color='#984ea3')
#plot=sns.histplot(data=B, x='Fold Change Sort AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density',legend=True, color='#e41a1c')
#plot=sns.histplot(data=C, x='Fold Change Sort AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density',legend=True, color='#377eb8')
plot=sns.histplot(data=D, x='Fold Change Sort AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density',legend=True, color='#999999')
plot.set_ylabel('Density', fontsize='x-large')
plot.set_xlabel('Log2 Fold Change Sort/Pre-sort', fontsize='x-large')
plt.ylim(0,0.4)
plt.xlim(-11,11)

handles=[]
for label, color in palette5.items():
    handle=matplotlib.patches.Patch(color=color, label=label)
    handles.append(handle)

plot.legend(handles=handles,fontsize='xx-large')
sns.move_legend(plot, "best")


# In[232]:


A=Known_mutants_df
B=stop_codons_df
C=WT_AA_df
D=unknown_df
plt.figure(figsize=(12,8))
#plot=sns.histplot(data=Known_mutants_df, x='Fold Change Sort AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density', legend=True, color='#984ea3')
plot=sns.histplot(data=B, x='Fold Change Sort AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density',legend=True, color='#e41a1c')
#plot=sns.histplot(data=C, x='Fold Change Sort AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density',legend=True, color='#377eb8')
#plot=sns.histplot(data=D, x='Fold Change Sort AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density',legend=True, color='#999999')
plot.set_ylabel('Density', fontsize='x-large')
plot.set_xlabel('Log2 Fold Change Sort/Pre-sort Stop Codons', fontsize='x-large')
plt.ylim(0,0.4)
plt.xlim(-11,11)



# In[233]:


A=Known_mutants_df
B=stop_codons_df
C=WT_AA_df
D=unknown_df
plt.figure(figsize=(8,8))
plot=sns.histplot(data=Known_mutants_df, x='Fold Change Sort AA Matches' , binwidth=1,kde=True, stat='density', legend=True, color='#984ea3')
plot=sns.histplot(data=B, x='Fold Change Sort AA Matches' , binwidth=1,kde=True, stat='density',legend=True, color='#e41a1c')
plot=sns.histplot(data=C, x='Fold Change Sort AA Matches' , binwidth=1,  kde=True, stat='density',legend=True, color='#377eb8')
plot=sns.histplot(data=D, x='Fold Change Sort AA Matches' , binwidth=1,  kde=True, stat='density',legend=True, color='#999999')
plot.set_ylabel('Density', fontsize='x-large')
plot.set_xlabel('Log2 Fold Change', fontsize='x-large')
plt.ylim(0,0.4)
plt.xlim(-11,10)

handles=[]
for label, color in palette.items():
    handle=matplotlib.patches.Patch(color=color, label=label)
    handles.append(handle)

plot.legend(handles=handles,fontsize='xx-large')
sns.move_legend(plot, "best")


# In[234]:


len(log2_Analysis_df_fold_change.loc[(log2_Analysis_df_fold_change['Fold Change Sort AA Matches'] >= 2)& (log2_Analysis_df_fold_change['Fold Change Sort AA Matches'] < 9.96)])


# In[235]:


len(log2_Analysis_df_fold_change.loc[(log2_Analysis_df_fold_change['Fold Change Sort AA Matches'] >=9.96)])


# In[236]:


log2_4_and_up=log2_Analysis_df_fold_change.loc[(log2_Analysis_df_fold_change['Fold Change Sort AA Matches'] >= 4)& (log2_Analysis_df_fold_change['Fold Change Sort AA Matches'] < 9.96)]


# In[237]:


len(log2_4_and_up)


# In[238]:


log2_4_and_up.sort_values(by='Fold Change Sort AA Matches', ascending=False)


# In[239]:


log_Raw_Totals_df.loc[(log_Raw_Totals_df['Pre-Sort'] <=2) & (log_Raw_Totals_df['Pre-Sort'] >-1)]


# In[240]:


log2_Analysis_df_fold_change.loc['618',:,'E']


# In[298]:


plt.figure(figsize=(100,40))
#plt.ylim(0, 30)
plt.xlim(-0.5,356.5)
#plt.xlabel('Position', fontsize=1)
#plt.ylabel('Fold Change Sort AA Matches', fontsize=1)
plt.xticks(rotation=90, fontsize=50)
plt.yticks(fontsize=50)
#sns.set(font_scale=10)
plot=sns.scatterplot(data=log2_Analysis_df_fold_change, x='Position', y='Fold Change Sort AA Matches', hue=legend2, palette=palette, size=size_legend_2, sizes=(400,800))
plot.set_xlabel('Position', fontsize=60)
plot.set_ylabel('Log2 Fold Change Sort AA Matches/Pre-sort', fontsize=60)
plt.axhline(y=0, color='black', linestyle='dashed')
plot.annotate
for label in plot.xaxis.get_ticklabels()[::2]:
    label.set_visible(False)

handles=[]
for label, color in palette.items():
    handle=matplotlib.patches.Patch(color=color, label=label)
    handles.append(handle)

plot.legend(handles=handles,fontsize=200)
sns.move_legend(plot, loc=(1,1))


# In[302]:


plt.figure(figsize=(100,20))
#plt.ylim(0, 30)
plt.xlim(-0.5,178.5)
#plt.xlabel('Position', fontsize=1)
#plt.ylabel('Fold Change Sort AA Matches', fontsize=1)
plt.xticks(rotation=90, fontsize=40)

plt.yticks(fontsize=50)
#sns.set(font_scale=10)
plot=sns.scatterplot(data=log2_Analysis_df_fold_change, x='Position', y='Fold Change Sort AA Matches', hue=legend2, palette=palette, size=size_legend_2, sizes=(400,800))
plot.set_xlabel('Position', fontsize=60)
plot.set_ylabel('Log2 Fold Change Sort AA Matches/Pre-sort', fontsize=50)
plt.axhline(y=0, color='black', linestyle='dashed')
plot.annotate
plot.axhline(y=4, color='black', linestyle='dashed', linewidth=10)
#plot.axhline(y=2, color='black', linestyle='dashed', linewidth=5)
for label in plot.xaxis.get_ticklabels()[::2]:
    label.set_visible(False)
handles=[]
for label, color in palette.items():
    handle=matplotlib.patches.Patch(color=color, label=label)
    handles.append(handle)

plot.legend(handles=handles,fontsize=200)
sns.move_legend(plot, loc=(1,1))


# In[304]:


plt.figure(figsize=(100,20))
plt.ylim(3.5, 10.5)
plt.xlim(-0.5,178.5)
#plt.xlabel('Position', fontsize=1)
#plt.ylabel('Fold Change Sort AA Matches', fontsize=1)
plt.xticks(rotation=90, fontsize=40)

plt.yticks(fontsize=50)
#sns.set(font_scale=10)
plot=sns.scatterplot(data=log2_Analysis_df_fold_change, x='Position', y='Fold Change Sort AA Matches', hue=legend2, palette=palette, size=size_legend_2, sizes=(400,800))
plot.set_xlabel('Position', fontsize=60)
plot.set_ylabel('Log2 Fold Change Sort AA Matches/Pre-sort', fontsize=50)
plt.axhline(y=0, color='black', linestyle='dashed')
plot.annotate
plot.axhline(y=4, color='black', linestyle='dashed', linewidth=10)
#plot.axhline(y=2, color='black', linestyle='dashed', linewidth=5)
for label in plot.xaxis.get_ticklabels()[::2]:
    label.set_visible(False)
handles=[]
for label, color in palette.items():
    handle=matplotlib.patches.Patch(color=color, label=label)
    handles.append(handle)

plot.legend(handles=handles,fontsize=200)
sns.move_legend(plot, loc=(1,1))


# In[303]:


plt.figure(figsize=(100,20))
#plt.ylim(0, 30)
plt.xlim(178.5,356.5)
#plt.xlabel('Position', fontsize=1)
#plt.ylabel('Fold Change Sort AA Matches', fontsize=1)
plt.xticks(rotation=90, fontsize=40)

plt.yticks(fontsize=50)
#sns.set(font_scale=10)
plot=sns.scatterplot(data=log2_Analysis_df_fold_change, x='Position', y='Fold Change Sort AA Matches', hue=legend2, palette=palette, size=size_legend_2, sizes=(400,800))
plot.set_xlabel('Position', fontsize=60)
plot.set_ylabel('Log2 Fold Change Sort AA Matches/Pre-sort', fontsize=50)
plt.axhline(y=0, color='black', linestyle='dashed')
plot.annotate

plot.axhline(y=4, color='black', linestyle='dashed', linewidth=10)
#plot.axhline(y=2, color='black', linestyle='dashed', linewidth=5)
for label in plot.xaxis.get_ticklabels()[::2]:
    label.set_visible(False)
handles=[]
for label, color in palette.items():
    handle=matplotlib.patches.Patch(color=color, label=label)
    handles.append(handle)

plot.legend(handles=handles,fontsize=200)
sns.move_legend(plot, loc=(1,1))


# In[305]:


plt.figure(figsize=(100,20))
plt.ylim(3.5, 10.5)
plt.xlim(178.5,356.5)
#plt.xlabel('Position', fontsize=1)
#plt.ylabel('Fold Change Sort AA Matches', fontsize=1)
plt.xticks(rotation=90, fontsize=40)

plt.yticks(fontsize=50)
#sns.set(font_scale=10)
plot=sns.scatterplot(data=log2_Analysis_df_fold_change, x='Position', y='Fold Change Sort AA Matches', hue=legend2, palette=palette, size=size_legend_2, sizes=(400,800))
plot.set_xlabel('Position', fontsize=60)
plot.set_ylabel('Log2 Fold Change Sort AA Matches/Pre-sort', fontsize=50)
plt.axhline(y=0, color='black', linestyle='dashed')
plot.annotate

plot.axhline(y=4, color='black', linestyle='dashed', linewidth=10)
#plot.axhline(y=2, color='black', linestyle='dashed', linewidth=5)
for label in plot.xaxis.get_ticklabels()[::2]:
    label.set_visible(False)
handles=[]
for label, color in palette.items():
    handle=matplotlib.patches.Patch(color=color, label=label)
    handles.append(handle)

plot.legend(handles=handles,fontsize=200)
sns.move_legend(plot, loc=(1,1))


# In[243]:


plt.figure(figsize=(100,40))
plt.ylim(1.5, 10.5)
plt.xlim(-0.5,356.5)
#plt.xlabel('Position', fontsize=1)
#plt.ylabel('Fold Change Sort AA Matches', fontsize=1)
plt.xticks(rotation=90, fontsize=50)
plt.yticks(fontsize=50)
#sns.set(font_scale=10)
plot=sns.scatterplot(data=log2_Analysis_df_fold_change, x='Position', y='Fold Change Sort AA Matches', hue=legend2, palette=palette, size=size_legend_2, sizes=(400,800))
plot.set_xlabel('Position', fontsize=60)
plot.set_ylabel('Log2 Fold Change Sort AA Matches/Pre-sort', fontsize=60)
plt.axhline(y=0, color='black', linestyle='dashed')
plot.annotate
plot.axhline(y=4, color='black', linestyle='dashed', linewidth=10)
plot.axhline(y=2, color='black', linestyle='dashed', linewidth=5)
handles=[]
for label, color in palette.items():
    handle=matplotlib.patches.Patch(color=color, label=label)
    handles.append(handle)

plot.legend(handles=handles,fontsize=200)
sns.move_legend(plot, loc=(1,1))


# In[244]:


Sort_Raw_Totals_Plate1_df=pd.DataFrame({'Raw Sort Count': S_Total_Matches_sort_Plate1_2})

# In[245]:


Sort_Raw_Totals_Plate2_df=pd.DataFrame({'Raw Sort Count': S_Total_Matches_sort_Plate2_2})


# In[246]:


Sort_Raw_Totals_Plate3_df=pd.DataFrame({'Raw Sort Count': S_Total_Matches_sort_Plate3_2})


# In[247]:


Sort_Raw_Totals_Plate4_df=pd.DataFrame({'Raw Sort Count': S_Total_Matches_sort_Plate4_2})


# In[248]:


Sort_Raw_Totals_df=pd.concat([Sort_Raw_Totals_Plate1_df,Sort_Raw_Totals_Plate2_df, Sort_Raw_Totals_Plate3_df, Sort_Raw_Totals_Plate4_df])


# In[249]:


Summary_Plate_Sort_df=Raw_Totals_df.join(Sort_Raw_Totals_df, how='right')
Summary_Plate_Sort_df.replace(0.5, 0, inplace=True)
Summary_Plate_Sort_df=Summary_Plate_Sort_df.join(Analysis_df, how='right')
Summary_Plate_Sort_df.rename(columns={'Normalized Total_Matches_Pre-Sort': 'Raw Normalized Total_Matches_Pre-Sort', 'Normalized Total_Matches_Sort': 'Raw Normalized Total_Matches_Sort', 'Fold Change Sort AA Matches':'Raw Fold Change Sort AA Matches'}, inplace=True)


# In[250]:


Summary_Plate_Sort_df=Summary_Plate_Sort_df.join(log2_Analysis_df, how='right')
Summary_Plate_Sort_df.rename(columns={'Normalized Total_Matches_Pre-Sort':'Log2 Normalized Total_Matches_Pre-Sort', 'Normalized Total_Matches_Sort':'Log2 Normalized Total_Matches_Sort'}, inplace=True)


# In[251]:


Summary_Plate_Sort_df.drop(['Position'], axis=1, inplace=True)


# In[252]:

# In[253]:


Summary_Plate_Sort_df=Summary_Plate_Sort_df.join(log2_Analysis_df_fold_change, how='right')
Summary_Plate_Sort_df.rename(columns={'Fold Change Sort AA Matches': 'Log2 Fold Change Sort AA Matches'}, inplace=True)

# In[254]:

Summary_Plate_Sort_df.drop(['Position'], axis=1, inplace=True)
Summary_Plate_Sort_df=pd.concat([Summary_Plate_Sort_df, Category_S], axis=1)
Summary_Plate_Sort_df.rename(columns={f'0':'Category'}, inplace=True)

# In[255]:


Summary_Plate_Sort_df.rename(columns={0:'Category'}, inplace=True)

# In[256]:


Summary_Plate_Sort_df.loc[Summary_Plate_Sort_df['Raw Fold Change Sort AA Matches'] > 1]


# In[257]:


# In[258]:


Analysis_df_by_AA=Analysis_df.groupby(level=[0,2]).sum()


# In[259]:

# In[260]:


Analysis_df_by_AA.drop(['Fold Change Sort AA Matches'],axis=1, inplace=True)
Analysis_df_by_AA['Fold Change Sort AA Matches'] = Analysis_df_by_AA['Normalized Total_Matches_Sort'] /Analysis_df_by_AA['Normalized Total_Matches_Pre-Sort']
Analysis_df_by_AA_fold_change=Analysis_df_by_AA.drop(['Normalized Total_Matches_Sort', 'Normalized Total_Matches_Pre-Sort'], axis=1, inplace=False)
Analysis_df_by_AA_fold_change['Position']=Analysis_df_by_AA_fold_change.index.get_level_values(0)
Analysis_df_by_AA_fold_change['Amino Acid']=Analysis_df_by_AA_fold_change.index.get_level_values(1)

# In[261]:


Analysis_df_by_AA_fold_change['Fold Change Sort AA Matches'].value_counts()


# In[262]:


Analysis_df_by_AA_fold_change.isna().sum()


# In[263]:


Analysis_df_by_AA_fold_change['Position']=Analysis_df_by_AA_fold_change['Position'].astype(int)


# In[264]:


Analysis_df_by_AA_fold_change.sort_values(by=['Position'], ascending=True, inplace=True)


# In[265]:


Analysis_df_by_AA_fold_change['Position']=Analysis_df_by_AA_fold_change['Position'].astype(str)


# In[266]:


Analysis_df_fold_change_by_AA_for_log=Analysis_df_by_AA_fold_change.drop(["Position"], axis=1, inplace=False)
Analysis_df_fold_change_by_AA_for_log=Analysis_df_fold_change_by_AA_for_log.drop(["Amino Acid"], axis=1)


# In[267]:


log2_Analysis_df_fold_change_by_AA=Analysis_df_fold_change_by_AA_for_log.replace(0,0.001)
log2_Analysis_df_fold_change_by_AA.replace(np.inf, 1000, inplace=True)
log2_Analysis_df_fold_change_by_AA=log2_Analysis_df_fold_change_by_AA.applymap(np.log2)


# In[268]:


log2_Analysis_df_fold_change_by_AA['Position']=log2_Analysis_df_fold_change_by_AA.index.get_level_values(0)
log2_Analysis_df_fold_change_by_AA['Amino Acid']=log2_Analysis_df_fold_change_by_AA.index.get_level_values(1)


# In[269]:


Analysis_df_index_list_2=list(log2_Analysis_df_fold_change_by_AA.index)


# In[270]:


List_Analysis_2=[]
S_Analysis_df_index_2=pd.Series(Analysis_df_index_list_2)
position_from_index_2=S_Analysis_df_index_2.str[0]
AA_from_index_2=S_Analysis_df_index_2.str[1]

for p2 in range(len(S_Analysis_df_index_2)):
    value_p_2=(position_from_index_2[p2], AA_from_index_2[p2])
    List_Analysis_2.append(value_p_2)

# In[271]:


legend3=[]
for i in List_Analysis_2:
    if i  == ('585', 'P'):
        legend3.append('Known mutant')
    elif i == ('603', 'I'):
        legend3.append('Known mutant')
    elif i == ('615', 'H'):
        legend3.append('Known mutant')
    elif i == ('618', 'E'):
        legend3.append('Known mutant')
    elif i == ('618', 'D'):
        legend3.append('Known mutant')    
    elif i  == ('972', 'A'):
        legend3.append('Known mutant')
    elif i == ('981', 'E'):
        legend3.append('Known mutant')
    elif i == ('981', 'W'):
        legend3.append('Known mutant')
    elif i == ('981', 'R'):
        legend3.append('Known mutant')
    elif i == ('771', 'S'):
        legend3.append('Known mutant') 
    elif i == ('771', 'G'):
        legend3.append('Known mutant')
    elif i == ('1122', 'P'):
        legend3.append('Known mutant') 
    elif i == ('1065', 'D'):
        legend3.append('Known mutant')
    elif i[1] == ('*'):
        legend3.append('Stop_Codon')
    elif i in WT_AA_loop:
        legend3.append('WT_AA')
    else:
        legend3.append('unknown')


# In[272]:


size_legend_3=[]
for legend in legend3:
    if legend == ('Known mutant'):
        size_legend_3.append(40)
    if legend == ('Stop_Codon'):
        size_legend_3.append(10)
    if legend == ('WT_AA'):
        size_legend_3.append(10)
    if legend == ('unknown'):
        size_legend_3.append(5)


# In[273]:


Category_S_2=pd.Series(legend3)


# In[274]:


Category_S_2.index = log2_Analysis_df_fold_change_by_AA.index


# In[275]:


log2_Analysis_df_fold_change_by_AA_category_df=pd.concat([log2_Analysis_df_fold_change_by_AA, Category_S_2], axis=1)

# In[276]:


log2_Analysis_df_fold_change_by_AA_category_df.rename(columns = {0:'Category'}, inplace=True)


# In[277]:


unknown_df_2=log2_Analysis_df_fold_change_by_AA_category_df[log2_Analysis_df_fold_change_by_AA_category_df['Category'].str.contains('unknown')]


# In[278]:


known_mutants_df_2=log2_Analysis_df_fold_change_by_AA_category_df[log2_Analysis_df_fold_change_by_AA_category_df['Category'].str.contains('Known mutant')]


# In[279]:


WT_AA_df_2=log2_Analysis_df_fold_change_by_AA_category_df[log2_Analysis_df_fold_change_by_AA_category_df['Category'].str.contains('WT_AA')]


# In[280]:


stop_codons_df_2=log2_Analysis_df_fold_change_by_AA_category_df[log2_Analysis_df_fold_change_by_AA_category_df['Category'].str.contains('Stop_Codon')]

# In[281]:


A=known_mutants_df_2
B=stop_codons_df_2
C=WT_AA_df_2
D=unknown_df_2
plt.figure(figsize=(12,8))
plot=sns.histplot(data=A, x='Fold Change Sort AA Matches' , binwidth=0.5,kde=True, stat='density', legend=True, color='#984ea3')
plot=sns.histplot(data=B, x='Fold Change Sort AA Matches' , binwidth=0.5,kde=True, stat='density',legend=True, color='#e41a1c')
plot=sns.histplot(data=C, x='Fold Change Sort AA Matches' , binwidth=0.5,  kde=True, stat='density',legend=True, color='#377eb8')
plot=sns.histplot(data=D, x='Fold Change Sort AA Matches' , binwidth=0.5,  kde=True, stat='density',legend=True, color='#999999')
plot.set_ylabel('Density', fontsize='x-large')
plot.set_xlabel('Log2 Fold Change Sort/Pre-sort', fontsize='x-large')
plt.ylim(0,0.7)
plt.xlim(-11,10)

handles=[]
for label, color in palette.items():
    handle=matplotlib.patches.Patch(color=color, label=label)
    handles.append(handle)

plot.legend(handles=handles,fontsize='xx-large')
sns.move_legend(plot, "best")


# In[282]:


palette6={'WT_AA':'#377eb8' , 'Known mutant':'#984ea3', 'unknown': '#999999'}
A=known_mutants_df_2
B=stop_codons_df_2
C=WT_AA_df_2
D=unknown_df_2
plt.figure(figsize=(12,8))
plot=sns.histplot(data=A, x='Fold Change Sort AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density', legend=True, color='#984ea3')
#plot=sns.histplot(data=B, x='Fold Change Sort AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density',legend=True, color='#e41a1c')
plot=sns.histplot(data=C, x='Fold Change Sort AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density',legend=True, color='#377eb8')
plot=sns.histplot(data=D, x='Fold Change Sort AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density',legend=True, color='#999999')
plot.set_ylabel('Density', fontsize='x-large')
plot.set_xlabel('Log2 Fold Change Sort/Pre-sort', fontsize='x-large')
plt.ylim(0,0.7)
plt.xlim(-11,11)
#plot.axvline(x=4.5, color='black', linestyle='dashed')
#plot.axvline(x=2, color='black')
handles=[]
for label, color in palette6.items():
    handle=matplotlib.patches.Patch(color=color, label=label)
    handles.append(handle)

plot.legend(handles=handles,fontsize='xx-large')
sns.move_legend(plot, "best")


# In[283]:


palette2={'WT_AA':'#377eb8' , 'Known mutant':'#984ea3'}
A=known_mutants_df_2
B=stop_codons_df_2
C=WT_AA_df_2
D=unknown_df_2
plt.figure(figsize=(12,8))
plot=sns.histplot(data=A, x='Fold Change Sort AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density', legend=True, color='#984ea3')
#plot=sns.histplot(data=B, x='Fold Change Sort AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density',legend=True, color='#e41a1c')
plot=sns.histplot(data=C, x='Fold Change Sort AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density',legend=True, color='#377eb8')
#plot=sns.histplot(data=D, x='Fold Change Sort AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density',legend=True, color='#999999')
plot.set_ylabel('Density', fontsize='x-large')
plot.set_xlabel('Log2 Fold Change Sort/Pre-sort', fontsize='x-large')
plt.ylim(0,0.7)
plt.xlim(-11,10)
plot.axvline(x=4, color='black', linestyle='dashed')
plot.axvline(x=2, color='black')
handles=[]
for label, color in palette2.items():
    handle=matplotlib.patches.Patch(color=color, label=label)
    handles.append(handle)

plot.legend(handles=handles,fontsize='xx-large')
sns.move_legend(plot, "best")


# In[284]:


A=known_mutants_df_2
B=stop_codons_df_2
C=WT_AA_df_2
D=unknown_df_2
plt.figure(figsize=(12,8))
plot=sns.histplot(data=A, x='Fold Change Sort AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density', legend=True, color='#984ea3')
#plot=sns.histplot(data=B, x='Fold Change Sort AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density',legend=True, color='#e41a1c')
#plot=sns.histplot(data=C, x='Fold Change Sort AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density',legend=True, color='#377eb8')
#plot=sns.histplot(data=D, x='Fold Change Sort AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density',legend=True, color='#999999')
plot.set_ylabel('Density', fontsize='x-large')
plot.set_xlabel('Log2 Fold Change Sort/Pre-sort', fontsize='x-large')
plt.ylim(0,0.7)
plt.xlim(-11,10)
#plot.axvline(x=4, color='black', linestyle='dashed')
#plot.axvline(x=2, color='black')
handles=[]
for label, color in palette3.items():
    handle=matplotlib.patches.Patch(color=color, label=label)
    handles.append(handle)

plot.legend(handles=handles,fontsize='xx-large')
sns.move_legend(plot, "best")


# In[285]:


A=known_mutants_df_2
B=stop_codons_df_2
C=WT_AA_df_2
D=unknown_df_2
plt.figure(figsize=(12,8))
#plot=sns.histplot(data=A, x='Fold Change Sort AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density', legend=True, color='#984ea3')
#plot=sns.histplot(data=B, x='Fold Change Sort AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density',legend=True, color='#e41a1c')
plot=sns.histplot(data=C, x='Fold Change Sort AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density',legend=True, color='#377eb8')
#plot=sns.histplot(data=D, x='Fold Change Sort AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density',legend=True, color='#999999')
plot.set_ylabel('Density', fontsize='x-large')
plot.set_xlabel('Log2 Fold Change Sort/Pre-sort', fontsize='x-large')
plt.ylim(0,0.7)
plt.xlim(-11,10)
#plot.axvline(x=4, color='black', linestyle='dashed')
#plot.axvline(x=2, color='black')
handles=[]
for label, color in palette4.items():
    handle=matplotlib.patches.Patch(color=color, label=label)
    handles.append(handle)

plot.legend(handles=handles,fontsize='xx-large')
sns.move_legend(plot, "best")


# In[286]:

# In[287]:


A=known_mutants_df_2
B=stop_codons_df_2
C=WT_AA_df_2
D=unknown_df_2
plt.figure(figsize=(12,8))
#plot=sns.histplot(data=A, x='Fold Change Sort AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density', legend=True, color='#984ea3')
#plot=sns.histplot(data=B, x='Fold Change Sort AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density',legend=True, color='#e41a1c')
#plot=sns.histplot(data=C, x='Fold Change Sort AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density',legend=True, color='#377eb8')
plot=sns.histplot(data=D, x='Fold Change Sort AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density',legend=True, color='#999999')
plot.set_ylabel('Density', fontsize='x-large')
plot.set_xlabel('Log2 Fold Change Sort/Pre-sort', fontsize='x-large')
plt.ylim(0,0.7)
plt.xlim(-11,10)
#plot.axvline(x=4, color='black', linestyle='dashed')
#plot.axvline(x=2, color='black')
handles=[]
for label, color in palette5.items():
    handle=matplotlib.patches.Patch(color=color, label=label)
    handles.append(handle)

plot.legend(handles=handles,fontsize='xx-large')
sns.move_legend(plot, "best")


# In[288]:


A=known_mutants_df_2
B=stop_codons_df_2
C=WT_AA_df_2
D=unknown_df_2
plt.figure(figsize=(12,8))
#plot=sns.histplot(data=A, x='Fold Change Sort AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density', legend=True, color='#984ea3')
plot=sns.histplot(data=B, x='Fold Change Sort AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density',legend=True, color='#e41a1c')
#plot=sns.histplot(data=C, x='Fold Change Sort AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density',legend=True, color='#377eb8')
#plot=sns.histplot(data=D, x='Fold Change Sort AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density',legend=True, color='#999999')
plot.set_ylabel('Density', fontsize='x-large')
plot.set_xlabel('Log2 Fold Change Sort/Pre-sort', fontsize='x-large')
plt.ylim(0,0.7)
plt.xlim(-11,10)
#plot.axvline(x=4, color='black', linestyle='dashed')
#plot.axvline(x=2, color='black')
#handles=[]
#for label, color in palette5.items():
    #handle=matplotlib.patches.Patch(color=color, label=label)
    #handles.append(handle)

#plot.legend(handles=handles,fontsize='xx-large')
#sns.move_legend(plot, "best")


# In[289]:
# In[290]:


plt.figure(figsize=(100,40))
#plt.ylim(0, 30)
plt.xlim(-0.5,356.5)
#plt.xlabel('Position', fontsize=1)
#plt.ylabel('Fold Change Sort AA Matches', fontsize=1)
plt.xticks(rotation=90, fontsize=50)
plt.yticks(fontsize=50)
#sns.set(font_scale=10)
plot=sns.scatterplot(data=log2_Analysis_df_fold_change_by_AA, x='Position', y='Fold Change Sort AA Matches', hue=legend3, palette=palette, size=size_legend_3, sizes=(400,800))
plot.set_xlabel('Position', fontsize=60)
plot.set_ylabel('Log2 Fold Change Sort AA Matches/Pre-sort', fontsize=60)
plt.axhline(y=0, color='black', linestyle='dashed')
plot.annotate(add info)
plot.axhline(y=2, color='black', linewidth=5)
plot.axhline(y=4, color='black', linestyle='dashed', linewidth=10)
handles=[]
for label, color in palette.items():
    handle=matplotlib.patches.Patch(color=color, label=label)
    handles.append(handle)
    
#for label in plot.xaxis.get_ticklabels()[::2]:
#    label.set_visible(False)

plot.legend(handles=handles,fontsize='xx-large')
sns.move_legend(plot, "best")


# In[309]:


plt.figure(figsize=(100,20))
#plt.ylim(0, 30)
plt.xlim(-0.5,178.5)
#plt.xlabel('Position', fontsize=1)
#plt.ylabel('Fold Change Sort AA Matches', fontsize=1)
plt.xticks(rotation=90, fontsize=50)
plt.yticks(fontsize=50)
#sns.set(font_scale=10)
plot=sns.scatterplot(data=log2_Analysis_df_fold_change_by_AA, x='Position', y='Fold Change Sort AA Matches', hue=legend3, palette=palette, size=size_legend_3, sizes=(400,800))
plot.set_xlabel('Position', fontsize=60)
plot.set_ylabel('Log2 Fold Change Sort AA Matches/Pre-sort', fontsize=50)
plt.axhline(y=0, color='black', linestyle='dashed')
plot.annotate(add info)
#plot.axhline(y=2, color='black', linewidth=5)
plot.axhline(y=4, color='black', linestyle='dashed', linewidth=10)
handles=[]
for label, color in palette.items():
    handle=matplotlib.patches.Patch(color=color, label=label)
    handles.append(handle)
    
for label in plot.xaxis.get_ticklabels()[::2]:
    label.set_visible(False)

plot.legend(handles=handles,fontsize='xx-large')
sns.move_legend(plot, "best")


# In[310]:


plt.figure(figsize=(100,20))
#plt.ylim(0, 30)
plt.xlim(178.5,356.5)
#plt.xlabel('Position', fontsize=1)
#plt.ylabel('Fold Change Sort AA Matches', fontsize=1)
plt.xticks(rotation=90, fontsize=50)
plt.yticks(fontsize=50)
#sns.set(font_scale=10)
plot=sns.scatterplot(data=log2_Analysis_df_fold_change_by_AA, x='Position', y='Fold Change Sort AA Matches', hue=legend3, palette=palette, size=size_legend_3, sizes=(400,800))
plot.set_xlabel('Position', fontsize=60)
plot.set_ylabel('Log2 Fold Change Sort AA Matches/Pre-sort', fontsize=50)
plt.axhline(y=0, color='black', linestyle='dashed')
plot.annotate(add info)
#plot.axhline(y=2, color='black', linewidth=5)
plot.axhline(y=4, color='black', linestyle='dashed', linewidth=10)
handles=[]
for label, color in palette.items():
    handle=matplotlib.patches.Patch(color=color, label=label)
    handles.append(handle)
    
for label in plot.xaxis.get_ticklabels()[::2]:
    label.set_visible(False)

plot.legend(handles=handles,fontsize='xx-large')
sns.move_legend(plot, "best")

