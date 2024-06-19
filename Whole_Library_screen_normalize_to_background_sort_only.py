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


pd.set_option('display.max_rows', 100)


# In[3]:


sample_list= [(insert samples)]


# In[4]:


sample_list_2=[(insert samples)]


# In[5]:


path= "(insert path name)"
sample='(insert sample)'


# In[6]:


List_of_samples_mapped_to_hACVR1=[]

for sample in sample_list:
    counts_df=pd.read_csv(rf'(insert file name)', sep='\t', header=0, index_col=0, skiprows=1)
    counts_df_2=counts_df.drop([(insert columns to drop)], axis=1)
    counts_df_2[f'Total_Counts{sample}']=counts_df_2[f'(insert file name)']
    counts_df_3=counts_df_2.drop([f'(insert file name)'], axis=1)
    Total_Count=int(counts_df_3[f'Total_Counts{sample}'])
    List_of_samples_mapped_to_hACVR1.append(f'{sample}' + ',' + str(Total_Count))
    
# In[7]:


Mapped_S=pd.Series(List_of_samples_mapped_to_hACVR1)


# In[8]:


S_mapped_sep = Mapped_S.str.split(',', 1, expand=False)


# In[9]:


Sample_Name=S_mapped_sep.str[0]
Count_Mapped=S_mapped_sep.str[1]


# In[10]:


Sample_Counts_df=pd.DataFrame({'Sample Name':Sample_Name , 'Counts':Count_Mapped})


# In[11]:


Sample_Counts_df.set_index('Sample Name', drop=True)


# In[12]:


Sample_Counts_df['Counts']=Sample_Counts_df['Counts'].astype(int)


# ##### Sample_Counts_df['Subtract from Raw Matches']=(20 * Sample_Counts_df['Counts'])/Sample_Counts_df['Counts'].loc[0]
# Sample_Counts_df

# In[14]:


List_of_samples_mapped_to_hACVR1_2=[]

for sample in sample_list_2:
    counts_df_4=pd.read_csv(rf'(insert file name)', sep='\t', header=0, index_col=0, skiprows=1)
    counts_df_5=counts_df_4.drop([(insert column names)], axis=1)
    counts_df_5[f'Total_Counts{sample}']=counts_df_5[f'(insert file name)'] + counts_df_5[f'(insert sample name)']
    counts_df_6=counts_df_5.drop([f'(insert column names)'], axis=1)
    Total_Count=int(counts_df_6[f'Total_Counts{sample}'])
    List_of_samples_mapped_to_hACVR1_2.append(f'{sample}' + ',' + str(Total_Count))
    
# In[15]:


Mapped_S_2=pd.Series(List_of_samples_mapped_to_hACVR1_2)


# In[16]:


S_mapped_sep_2 = Mapped_S_2.str.split(',', 1, expand=False)


# In[17]:


Sample_Name_2=S_mapped_sep_2.str[0]
Count_Mapped_2=S_mapped_sep_2.str[1]


# In[18]:


Sample_Counts_df_2=pd.DataFrame({'Sample Name':Sample_Name_2 , 'Counts':Count_Mapped_2})


# In[19]:


Sample_Counts_df_2.set_index('Sample Name', drop=True)


# In[20]:


Sample_Counts_df_2['Counts']=Sample_Counts_df_2['Counts'].astype(int)


# In[21]:


Sample_Counts_df_2['Subtract from Raw Matches']=(20 * Sample_Counts_df_2['Counts'])/Sample_Counts_df_2['Counts'].loc[0]


# In[22]:


reference = SeqIO.read("(insert fast file name)", "fasta")


# In[23]:


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


# In[24]:


WT_seqs_tuples_join = []
seq = str(reference.seq)
for position in range(438, len(seq)-3,3):
    bp = seq[position:position+3]
    WT_position = str(position)
    WT_codon= bp
    WT_AA=codontab[WT_codon]
    WT_seqs_tuples_add=(WT_position + WT_codon)
    WT_seqs_tuples_join.append(WT_seqs_tuples_add)


# In[25]:


rows_to_remove = WT_seqs_tuples_join


# In[26]:


S_rows_to_remove=pd.Series(rows_to_remove)
S_rows_to_remove_sep = S_rows_to_remove.str.split('(\d+)', 1, expand=False)
Positions_from_index_remove=S_rows_to_remove_sep.str[1]
Codons_from_index_remove=S_rows_to_remove_sep.str[2]


# In[27]:


list_of_WT_AA= []
position_WT=[]
AA_for_WT=[]
Codon_for_WT=[]
for i in range(len(S_rows_to_remove_sep)):
    codon = Codons_from_index_remove[i]
    list_value_i=(Positions_from_index_remove[i]+'P1', codontab[codon])
    list_value_i_2=(Positions_from_index_remove[i]+'P1')
    list_value_i_3=codontab[codon]
    list_value_i_4=codon
    list_value_i_5=(Positions_from_index_remove[i]+'P2', codontab[codon])
    list_value_i_6=(Positions_from_index_remove[i]+'P2')
    list_of_WT_AA.append(list_value_i)
    list_of_WT_AA.append(list_value_i_5)
    position_WT.append(list_value_i_2)
    position_WT.append(list_value_i_6)
    AA_for_WT.append(list_value_i_3)
    Codon_for_WT.append(list_value_i_4)
list_of_WT_AA[0:5]


# In[28]:


S_WT_AA_2=pd.Series(list_of_WT_AA)
position_WT_2=S_WT_AA_2.str[0]
amino_acid_2=S_WT_AA_2.str[1]
S_WT_AA_2


# In[29]:


WT_dict={}
for i in range(len(S_WT_AA_2)):
    WT_dict_AA = amino_acid_2[i]
    WT_dict_pos = position_WT_2[i]
    WT_dict[WT_dict_pos] = WT_dict_AA  
WT_dict

           


# ## sample 1 name

# In[30]:


sample = '(insert sample name)'


# In[31]:


perfect_match_df = pd.read_csv("(insert file name)", sep=',', header=0, index_col=0, skiprows=0)
perfect_match_df.rename(columns = {f'0':'Matches_R1'}, inplace=True)
perfect_match_df['Matches_R2'] = pd.read_csv("(insert file name)", sep=',', header=0, index_col=0, skiprows=0)
perfect_match_df['Total_Matches'] = perfect_match_df['Matches_R1'] + perfect_match_df['Matches_R2']

# In[32]:


new_index_pre=perfect_match_df.index.difference(rows_to_remove)
perfect_match_without_WT_df=perfect_match_df.loc[new_index_pre]

# In[33]:


S_matches = pd.Series(perfect_match_without_WT_df.index, dtype="string")
S_matches_sep = S_matches.str.split('(\d+)', 1, expand=False)
Positions_from_index=S_matches_sep.str[1]
Codons_from_index=S_matches_sep.str[2]


# In[34]:


list_of_tuples_for_multiindex = []
for i in range(len(S_matches_sep)):
    codon = Codons_from_index[i]
    multi_index_value_i=(Positions_from_index[i]+'P1', codon, codontab[codon])
    list_of_tuples_for_multiindex.append(multi_index_value_i)
new_index = pd.MultiIndex.from_tuples(list_of_tuples_for_multiindex)
perfect_match_without_WT_df.index = new_index
perfect_match_without_WT_df.drop(['Matches_R1','Matches_R2'], axis=1, inplace=True)


# In[35]:


positions_to_drop=perfect_match_without_WT_df.loc['804P1':'822P1']


# In[36]:


L_positions_to_drop=list(positions_to_drop.index)


# In[37]:


new_index_pre_2=perfect_match_without_WT_df.index.difference(L_positions_to_drop)
perfect_match_without_WT_df=perfect_match_without_WT_df.loc[new_index_pre_2]

# In[38]:


S_Total_Matches_pre_P1=pd.Series(perfect_match_without_WT_df['Total_Matches'])

# In[39]:


perfect_match_without_WT_df['Total_Matches'].value_counts()


# In[40]:


Count_mutant=int(perfect_match_without_WT_df.sum())

# In[41]:


perfect_match_without_WT_df['Normalized Total_Matches_Pre-Sort']=perfect_match_without_WT_df['Total_Matches']/Count_mutant
perfect_match_without_WT_normalized_df_1=perfect_match_without_WT_df.drop(['Total_Matches'], axis=1)


# ## Sample 2 name

# In[42]:


sample_2 = '(insert sample name)'


# In[43]:


perfect_match_df_2 = pd.read_csv("(insert file name)", sep=',', header=0, index_col=0, skiprows=0)
perfect_match_df_2.rename(columns = {f'0':'Matches_R1'}, inplace=True)
perfect_match_df_2['Matches_R2'] = pd.read_csv("(insert sample name)", sep=',', header=0, index_col=0, skiprows=0)
perfect_match_df_2['Total_Matches'] = perfect_match_df_2['Matches_R1'] + perfect_match_df_2['Matches_R2']

# In[44]:


new_index_sort_No_AA=perfect_match_df_2.index.difference(rows_to_remove)
perfect_match_without_WT_df_2=perfect_match_df_2.loc[new_index_sort_No_AA]

# In[45]:


S_matches_2 = pd.Series(perfect_match_without_WT_df_2.index, dtype="string")
S_matches_sep_2 = S_matches_2.str.split('(\d+)', 1, expand=False)
Positions_from_index_2=S_matches_sep_2.str[1]
Codons_from_index_2=S_matches_sep_2.str[2]


# In[46]:


list_of_tuples_for_multiindex_2 = []
for i in range(len(S_matches_sep_2)):
    codon = Codons_from_index_2[i]
    multi_index_value_i_2=(Positions_from_index_2[i]+'P1', codon, codontab[codon])
    list_of_tuples_for_multiindex_2.append(multi_index_value_i_2)
new_index_2 = pd.MultiIndex.from_tuples(list_of_tuples_for_multiindex_2)
perfect_match_without_WT_df_2.index = new_index_2
perfect_match_without_WT_df_2.drop(['Matches_R1','Matches_R2'], axis=1, inplace=True)


# In[47]:


new_index_sort_2=perfect_match_without_WT_df_2.index.difference(L_positions_to_drop)
perfect_match_without_WT_df_2=perfect_match_without_WT_df_2.loc[new_index_sort_2]

# In[48]:


S_Total_Matches_sort_No_AA=pd.Series(perfect_match_without_WT_df_2['Total_Matches'])


# In[49]:


Sample_Counts_df.loc[1]


# In[50]:


perfect_match_without_WT_df_2['Normalize Raw Counts']=perfect_match_without_WT_df_2['Total_Matches'] - Sample_Counts_df['Subtract from Raw Matches'].loc[1]

# In[51]:


perfect_match_without_WT_df_2[perfect_match_without_WT_df_2 < 0] = 0   

# In[52]:


S_Total_Matches_sort_No_AA_1=pd.Series(perfect_match_without_WT_df_2['Normalize Raw Counts'])

# In[53]:


perfect_match_without_WT_df_2['Normalize Raw Counts'].value_counts()


# In[54]:


perfect_match_without_WT_df_2.loc[perfect_match_without_WT_df_2['Normalize Raw Counts'] == 0]


# In[55]:


Count_mutant_2=int(perfect_match_without_WT_df_2['Normalize Raw Counts'].sum())


# In[56]:


perfect_match_without_WT_df_2['Normalized Total_Matches_Sort_No_AA']=perfect_match_without_WT_df_2['Normalize Raw Counts']/Count_mutant_2
perfect_match_without_WT_normalized_df_2=perfect_match_without_WT_df_2.drop(['Total_Matches', 'Normalize Raw Counts'], axis=1)


# ## Sample 3 Name

# In[57]:


sample_3 = '(insert sample name)'


# In[58]:


perfect_match_df_3 = pd.read_csv("(insert file name)", sep=',', header=0, index_col=0, skiprows=0)
perfect_match_df_3.rename(columns = {f'0':'Matches_R1'}, inplace=True)
perfect_match_df_3['Matches_R2'] = pd.read_csv("(insert file name)", sep=',', header=0, index_col=0, skiprows=0)
perfect_match_df_3['Total_Matches'] = perfect_match_df_3['Matches_R1'] + perfect_match_df_3['Matches_R2']

# In[59]:


new_index_sort_AA=perfect_match_df_3.index.difference(rows_to_remove)
perfect_match_without_WT_df_3=perfect_match_df_3.loc[new_index_sort_AA]


# In[60]:


S_matches_3 = pd.Series(perfect_match_without_WT_df_3.index, dtype="string")
S_matches_sep_3 = S_matches_3.str.split('(\d+)', 1, expand=False)
Positions_from_index_3=S_matches_sep_3.str[1]
Codons_from_index_3=S_matches_sep_3.str[2]


# In[61]:


list_of_tuples_for_multiindex_3 = []
for i in range(len(S_matches_sep_3)):
    codon = Codons_from_index_3[i]
    multi_index_value_i_3=(Positions_from_index_3[i]+'P1', codon, codontab[codon])
    list_of_tuples_for_multiindex_3.append(multi_index_value_i_3)
new_index_3 = pd.MultiIndex.from_tuples(list_of_tuples_for_multiindex_3)
perfect_match_without_WT_df_3.index = new_index_3
perfect_match_without_WT_df_3.drop(['Matches_R1','Matches_R2'], axis=1, inplace=True)


# In[62]:


new_index_sort_3=perfect_match_without_WT_df_3.index.difference(L_positions_to_drop)
perfect_match_without_WT_df_3=perfect_match_without_WT_df_3.loc[new_index_sort_3]

# In[63]:


Sample_Counts_df.loc[2]


# In[64]:


perfect_match_without_WT_df_3['Normalize Raw Counts']=perfect_match_without_WT_df_3['Total_Matches'] - Sample_Counts_df['Subtract from Raw Matches'].loc[2]


# In[65]:


perfect_match_without_WT_df_3[perfect_match_without_WT_df_3 < 0] = 0   


# In[66]:


S_Total_Matches_sort_AA_P1=pd.Series(perfect_match_without_WT_df_3['Total_Matches'])


# In[67]:


S_Total_Matches_sort_AA_1=pd.Series(perfect_match_without_WT_df_3['Normalize Raw Counts'])


# In[68]:


len(perfect_match_without_WT_df_3.loc[perfect_match_without_WT_df_3['Normalize Raw Counts'] == 0])


# In[69]:


Count_mutant_3=int(perfect_match_without_WT_df_3['Normalize Raw Counts'].sum())


# In[70]:


perfect_match_without_WT_df_3['Normalize Raw Counts'].value_counts()


# In[71]:


perfect_match_without_WT_df_3['Normalized Total_Matches_Sort_AA']=perfect_match_without_WT_df_3['Normalize Raw Counts']/Count_mutant_3
perfect_match_without_WT_normalized_df_3=perfect_match_without_WT_df_3.drop(['Total_Matches', 'Normalize Raw Counts'], axis=1)


# ## Sample 4 Name

# In[72]:


sample4='(insert sample name)'


# In[73]:


perfect_match_df_4 = pd.read_csv("(insert file name)", sep=',', header=0, index_col=0, skiprows=0)
perfect_match_df_4.rename(columns = {f'0':'Matches_R1'}, inplace=True)
perfect_match_df_4['Matches_R2'] = pd.read_csv("(insert file name)", sep=',', header=0, index_col=0, skiprows=0)
perfect_match_df_4['Total_Matches'] = perfect_match_df_4['Matches_R1'] + perfect_match_df_4['Matches_R2']


# In[74]:


new_index_pre_3=perfect_match_df_4.index.difference(rows_to_remove)
perfect_match_without_WT_df_4=perfect_match_df_4.loc[new_index_pre_3]


# In[75]:


S_matches_4 = pd.Series(perfect_match_without_WT_df_4.index, dtype="string")
S_matches_sep_4 = S_matches_4.str.split('(\d+)', 1, expand=False)
Positions_from_index_4=S_matches_sep_4.str[1]
Codons_from_index_4=S_matches_sep_4.str[2]


# In[76]:


list_of_tuples_for_multiindex_4 = []
for i in range(len(S_matches_sep_4)):
    codon = Codons_from_index_4[i]
    multi_index_value_i_4=(Positions_from_index_4[i]+'P2', codon, codontab[codon])
    list_of_tuples_for_multiindex_4.append(multi_index_value_i_4)
new_index_4 = pd.MultiIndex.from_tuples(list_of_tuples_for_multiindex_4)
perfect_match_without_WT_df_4.index = new_index_4
perfect_match_without_WT_df_4.drop(['Matches_R1','Matches_R2'], axis=1, inplace=True)


# In[77]:


positions_to_drop_2=perfect_match_without_WT_df_4.loc[(insert range)]
positions_to_drop_3=perfect_match_without_WT_df_4.loc[(insert range)]



# In[78]:


L_positions_to_drop_2=list(positions_to_drop_2.index)


# In[79]:


L_positions_to_drop_3=list(positions_to_drop_3.index)


# In[80]:


new_index_pre_4=perfect_match_without_WT_df_4.index.difference(L_positions_to_drop_2)
perfect_match_without_WT_df_4=perfect_match_without_WT_df_4.loc[new_index_pre_4]


# In[81]:


new_index_pre_5=perfect_match_without_WT_df_4.index.difference(L_positions_to_drop_3)
perfect_match_without_WT_df_4=perfect_match_without_WT_df_4.loc[new_index_pre_5]

# In[82]:


S_Total_Matches_pre_P2=pd.Series(perfect_match_without_WT_df_4['Total_Matches'])


# In[83]:


perfect_match_without_WT_df_4.value_counts()


# In[84]:


Count_mutant_4=int(perfect_match_without_WT_df_4.sum())


# In[85]:


perfect_match_without_WT_df_4['Normalized Total_Matches_Pre-Sort']=perfect_match_without_WT_df_4['Total_Matches']/Count_mutant_4
perfect_match_without_WT_normalized_df_4=perfect_match_without_WT_df_4.drop(['Total_Matches'], axis=1)


# ## Sample 5 Name

# In[86]:


sample_5='(insert name)'


# In[87]:


perfect_match_df_5 = pd.read_csv("(insert file name)", sep=',', header=0, index_col=0, skiprows=0)
perfect_match_df_5.rename(columns = {f'0':'Matches_R1'}, inplace=True)
perfect_match_df_5['Matches_R2'] = pd.read_csv("(insert file name)", sep=',', header=0, index_col=0, skiprows=0)
perfect_match_df_5['Total_Matches'] = perfect_match_df_5['Matches_R1'] + perfect_match_df_5['Matches_R2']

# In[88]:


new_index_sort_No_AA_2=perfect_match_df_5.index.difference(rows_to_remove)
perfect_match_without_WT_df_5=perfect_match_df_5.loc[new_index_sort_No_AA_2]

# In[89]:


S_matches_5 = pd.Series(perfect_match_without_WT_df_5.index, dtype="string")
S_matches_sep_5 = S_matches_5.str.split('(\d+)', 1, expand=False)
Positions_from_index_5=S_matches_sep_5.str[1]
Codons_from_index_5=S_matches_sep_5.str[2]


# In[90]:


list_of_tuples_for_multiindex_5 = []
for i in range(len(S_matches_sep_5)):
    codon = Codons_from_index_5[i]
    multi_index_value_i_5=(Positions_from_index_5[i]+'P2', codon, codontab[codon])
    list_of_tuples_for_multiindex_5.append(multi_index_value_i_5)
new_index_5 = pd.MultiIndex.from_tuples(list_of_tuples_for_multiindex_5)
perfect_match_without_WT_df_5.index = new_index_5
perfect_match_without_WT_df_5.drop(['Matches_R1','Matches_R2'], axis=1, inplace=True)


# In[91]:


new_index_sort_5=perfect_match_without_WT_df_5.index.difference(L_positions_to_drop_2)
perfect_match_without_WT_df_5=perfect_match_without_WT_df_5.loc[new_index_sort_5]

# In[92]:


new_index_sort_6=perfect_match_without_WT_df_5.index.difference(L_positions_to_drop_3)
perfect_match_without_WT_df_5=perfect_match_without_WT_df_5.loc[new_index_sort_6]

# In[93]:


Sample_Counts_df_2.loc[1]


# In[94]:


perfect_match_without_WT_df_5['Normalize Raw Counts']=perfect_match_without_WT_df_5['Total_Matches'] - Sample_Counts_df_2['Subtract from Raw Matches'].loc[1]


# In[95]:


perfect_match_without_WT_df_5[perfect_match_without_WT_df_5 < 0] = 0   


# In[96]:


S_Total_Matches_sort_No_AA_3=pd.Series(perfect_match_without_WT_df_5['Total_Matches'])

# In[97]:


S_Total_Matches_sort_No_AA_2=pd.Series(perfect_match_without_WT_df_5['Normalize Raw Counts'])


# In[98]:


perfect_match_without_WT_df_5['Normalize Raw Counts'].value_counts()


# In[99]:


Count_mutant_5=int(perfect_match_without_WT_df_5['Normalize Raw Counts'].sum())


# In[100]:


perfect_match_without_WT_df_5['Normalized Total_Matches_Sort_No_AA']=perfect_match_without_WT_df_5['Normalize Raw Counts']/Count_mutant_5
perfect_match_without_WT_normalized_df_5=perfect_match_without_WT_df_5.drop(['Total_Matches', 'Normalize Raw Counts'], axis=1)


# ## Sample 6 name

# In[101]:


sample_6= '(insert name)'


# In[102]:


perfect_match_df_6 = pd.read_csv("(insert file name)", sep=',', header=0, index_col=0, skiprows=0)
perfect_match_df_6.rename(columns = {f'0':'Matches_R1'}, inplace=True)
perfect_match_df_6['Matches_R2'] = pd.read_csv("(insert file name)", sep=',', header=0, index_col=0, skiprows=0)
perfect_match_df_6['Total_Matches'] = perfect_match_df_6['Matches_R1'] + perfect_match_df_6['Matches_R2']

# In[103]:


new_index_sort_AA_2=perfect_match_df_6.index.difference(rows_to_remove)
perfect_match_without_WT_df_6=perfect_match_df_6.loc[new_index_sort_AA_2]


# In[104]:


S_matches_6 = pd.Series(perfect_match_without_WT_df_6.index, dtype="string")
S_matches_sep_6 = S_matches_6.str.split('(\d+)', 1, expand=False)
Positions_from_index_6=S_matches_sep_6.str[1]
Codons_from_index_6=S_matches_sep_6.str[2]


# In[105]:


list_of_tuples_for_multiindex_6 = []
for i in range(len(S_matches_sep_6)):
    codon = Codons_from_index_6[i]
    multi_index_value_i_6=(Positions_from_index_6[i]+'P2', codon, codontab[codon])
    list_of_tuples_for_multiindex_6.append(multi_index_value_i_6)
new_index_6 = pd.MultiIndex.from_tuples(list_of_tuples_for_multiindex_6)
perfect_match_without_WT_df_6.index = new_index_6
perfect_match_without_WT_df_6.drop(['Matches_R1','Matches_R2'], axis=1, inplace=True)


# In[106]:


new_index_sort_7=perfect_match_without_WT_df_6.index.difference(L_positions_to_drop_2)
perfect_match_without_WT_df_6=perfect_match_without_WT_df_6.loc[new_index_sort_7]

# In[107]:


new_index_sort_8=perfect_match_without_WT_df_6.index.difference(L_positions_to_drop_3)
perfect_match_without_WT_df_6=perfect_match_without_WT_df_6.loc[new_index_sort_8]

# In[108]:


Sample_Counts_df_2.loc[2]


# In[109]:


perfect_match_without_WT_df_6['Normalize Raw Counts']=perfect_match_without_WT_df_6['Total_Matches'] - Sample_Counts_df_2['Subtract from Raw Matches'].loc[2]

# In[110]:


perfect_match_without_WT_df_6[perfect_match_without_WT_df_6 < 0] = 0   


# In[111]:


S_Total_Matches_sort_AA_2=pd.Series(perfect_match_without_WT_df_6['Normalize Raw Counts'])


# In[112]:


S_Total_Matches_sort_AA_3=pd.Series(perfect_match_without_WT_df_6['Total_Matches'])


# In[113]:


perfect_match_without_WT_df_6['Normalize Raw Counts'].value_counts()


# In[114]:


Count_mutant_6=int(perfect_match_without_WT_df_6['Normalize Raw Counts'].sum())


# In[115]:


perfect_match_without_WT_df_6['Normalized Total_Matches_Sort_AA']=perfect_match_without_WT_df_6['Normalize Raw Counts']/Count_mutant_6
perfect_match_without_WT_normalized_df_6=perfect_match_without_WT_df_6.drop(['Total_Matches', 'Normalize Raw Counts'], axis=1)


# ## Analysis

# In[116]:


# In[117]:


Raw_Totals_P1_df=pd.DataFrame({'Pre-Sort':S_Total_Matches_pre_P1, 'Sort No AA':S_Total_Matches_sort_No_AA_1,'Sort AA': S_Total_Matches_sort_AA_1})


# In[118]:


plt.figure(figsize=(8,8))
plot=sns.histplot(data=Raw_Totals_P1_df,binwidth=1, kde=True, element='step', stat='density')
plot.set_xlabel('Raw Totals P1', fontsize='x-large')
plot.set_ylabel('Density', fontsize='x-large')
plt.ylim(0,0.01)
plt.xlim(-1,100)


# In[119]:


Raw_Totals_P2_df=pd.DataFrame({'Pre-Sort':S_Total_Matches_pre_P2, 'Sort No AA':S_Total_Matches_sort_No_AA_2,'Sort AA': S_Total_Matches_sort_AA_2})


# In[120]:


plt.figure(figsize=(8,8))
plot=sns.histplot(data=Raw_Totals_P2_df,binwidth=1, kde=True, element='step', stat='density')
plot.set_xlabel('Raw Totals P2', fontsize='x-large')
plot.set_ylabel('Density', fontsize='x-large')
plt.ylim(0,0.01)
plt.xlim(-1,100)


# In[121]:


Raw_Totals_df=pd.concat([Raw_Totals_P1_df, Raw_Totals_P2_df], ignore_index=False, axis=0)


# In[122]:


plt.figure(figsize=(8,8))
plot=sns.histplot(data=Raw_Totals_df,binwidth=1, kde=True, element='step', stat='density')
plot.set_xlabel('Raw Totals including overlap', fontsize='x-large')
plot.set_ylabel('Density', fontsize='x-large')
plt.ylim(0,0.01)
plt.xlim(-1,100)


# In[123]:


Raw_Totals_df.replace(0, 0.5, inplace=True)


# In[124]:


log_Raw_Totals_df=Raw_Totals_df.applymap(np.log2)


# In[125]:


log_Raw_Totals_df['Position']=log_Raw_Totals_df.index.get_level_values(0)


# In[126]:


S_Position = pd.Series(log_Raw_Totals_df['Position'], dtype="string")
S_Position_sep = S_Position.str.split('(\d+)', 1, expand=False)
Positions_from_index=S_Position_sep.str[1]
PCR_Product_from_index=S_Position_sep.str[2]


# In[127]:


log_Raw_Totals_df['Position Number']=Positions_from_index


# In[128]:


log_Raw_Totals_df['Position Number']=log_Raw_Totals_df['Position Number'].astype(int)
log_Raw_Totals_df.sort_values(by=['Position Number'], ascending=True, inplace=True)


# In[129]:

# In[130]:


#log_Raw_Totals_df['Position Number']=log_Raw_Totals_df['Position Number'].astype(str)


# In[131]:


overlap_df=Raw_Totals_P1_df.loc[(insert range)]


# In[132]:


overlap_df_2=Raw_Totals_P2_df.loc[(insert range)]


# In[133]:


plt.figure(figsize=(8,8))
plot=sns.histplot(data=overlap_df,binwidth=1, kde=True, element='step', stat='density')
plot.set_xlabel('Raw Totals overlap P1', fontsize='x-large')
plot.set_ylabel('Density', fontsize='x-large')
plt.ylim(0,0.01)
plt.xlim(-1,100)


# In[134]:


plt.figure(figsize=(8,8))
plot=sns.histplot(data=overlap_df_2,binwidth=1, kde=True, element='step', stat='density')
plot.set_xlabel('Raw Totals overlap P2', fontsize='x-large')
plot.set_ylabel('Density', fontsize='x-large')
plt.ylim(0,0.01)
plt.xlim(-1,100)


# In[135]:


Raw_Totals_P1_df.replace(0, 0.5, inplace=True)


# In[136]:


Raw_Totals_P2_df.replace(0, 0.5, inplace=True)


# In[137]:


Analysis_df_P1=pd.concat([perfect_match_without_WT_normalized_df_1,perfect_match_without_WT_normalized_df_2], axis=1)


# In[138]:


Analysis_df_P1=pd.concat([Analysis_df_P1, perfect_match_without_WT_normalized_df_3], axis=1)


# In[139]:


Analysis_df_P1['Position']=Analysis_df_P1.index.get_level_values(0)


# In[140]:


S_Position_3 = pd.Series(Analysis_df_P1['Position'], dtype="string")
S_Position_sep_3 = S_Position_3.str.split('(\d+)', 1, expand=False)
Positions_from_index_3=S_Position_sep_3.str[1]
PCR_Product_from_index_3=S_Position_sep_3.str[2]


# In[141]:


Analysis_df_P1['Position Number']=Positions_from_index_3


# In[142]:

# In[143]:


Analysis_df_P2=pd.concat([perfect_match_without_WT_normalized_df_4, perfect_match_without_WT_normalized_df_5, perfect_match_without_WT_normalized_df_6], axis=1)


# In[144]:

# In[145]:


Analysis_df_P2['Position']=Analysis_df_P2.index.get_level_values(0)


# In[146]:


S_Position_2 = pd.Series(Analysis_df_P2['Position'], dtype="string")
S_Position_sep_2 = S_Position_2.str.split('(\d+)', 1, expand=False)
Positions_from_index_2=S_Position_sep_2.str[1]
PCR_Product_from_index_2=S_Position_sep_2.str[2]


# In[147]:


Analysis_df_P2['Position Number']=Positions_from_index_2


# In[148]:


# Analysis_df_P2['Position Number']=Analysis_df_P2['Position Number'].astype(int)
# Analysis_df_P2.sort_values(by=['Position Number'], ascending=True, inplace=True)

# In[149]:


Analysis_df=pd.concat([Analysis_df_P1, Analysis_df_P2])

# In[150]:


Analysis_df['Position Number']=Analysis_df['Position Number'].astype(int)
Analysis_df.sort_values(by=['Position Number'], ascending=True, inplace=True)


# In[151]:


Analysis_df['Normalized Total_Matches_Pre-Sort'].value_counts()


# In[152]:


Analysis_df['Normalized Total_Matches_Sort_No_AA'].value_counts()


# In[153]:


Analysis_df['Normalized Total_Matches_Sort_AA'].value_counts()


# In[154]:


Analysis_df_for_log=Analysis_df.replace(0, 0.000001, inplace=False)
Analysis_df_for_log.replace(np.inf, 0.3, inplace=True) 
Analysis_df_for_log.drop(['Position', 'Position Number'], axis=1, inplace=True)
log2_Analysis_df=Analysis_df_for_log.applymap(np.log2)


# In[155]:


Analysis_df['Fold Change Sort AA Matches'] = Analysis_df['Normalized Total_Matches_Sort_AA'] /Analysis_df['Normalized Total_Matches_Pre-Sort']
Analysis_df['Fold Change Sort No AA Matches']= Analysis_df['Normalized Total_Matches_Sort_No_AA'] /Analysis_df['Normalized Total_Matches_Pre-Sort']
Analysis_df_fold_change=Analysis_df.drop(['Normalized Total_Matches_Sort_AA','Normalized Total_Matches_Sort_No_AA','Normalized Total_Matches_Pre-Sort', 'Position', 'Position Number'], axis=1, inplace=False)

# In[156]:


Analysis_df_fold_change['Fold Change Sort AA Matches'].value_counts()


# In[157]:


Analysis_df_fold_change['Fold Change Sort No AA Matches'].value_counts()


# In[158]:


Analysis_df_fold_change.isna().sum()


# In[159]:


Analysis_df_fold_change.replace(0, 0.001, inplace=True)
Analysis_df_fold_change.replace(np.inf, 700, inplace=True)
log2_Analysis_df_fold_change=Analysis_df_fold_change.applymap(np.log2)


# In[160]:


log2_Analysis_df_fold_change['Position']=log2_Analysis_df_fold_change.index.get_level_values(0)


# In[161]:


S_Position_4 = pd.Series(log2_Analysis_df_fold_change['Position'], dtype="string")
S_Position_sep_4 = S_Position_4.str.split('(\d+)', 1, expand=False)
Positions_from_index_4=S_Position_sep_4.str[1]
PCR_Product_from_index_4=S_Position_sep_4.str[2]

# In[162]:


log2_Analysis_df_fold_change['Position Number']=Positions_from_index_4
log2_Analysis_df_fold_change['PCR Product']=PCR_Product_from_index_4
log2_Analysis_df_fold_change['Codon']=log2_Analysis_df_fold_change.index.get_level_values(1)
log2_Analysis_df_fold_change['Amino']=log2_Analysis_df_fold_change.index.get_level_values(2)


# In[163]:


PCR_Product_1_df_log2_fold_change=log2_Analysis_df_fold_change[log2_Analysis_df_fold_change['PCR Product'].str.contains('P1')]
len(PCR_Product_1_df_log2_fold_change)


# In[164]:


PCR_Product_2_df_log2_fold_change=log2_Analysis_df_fold_change[log2_Analysis_df_fold_change['PCR Product'].str.contains('P2')]
len(PCR_Product_2_df_log2_fold_change)


# In[165]:

# In[166]:


PCR_Overlap_P2_log2_fold_change=PCR_Product_2_df_log2_fold_change.head(3591)


# In[167]:


PCR_Overlap_P1_log2_fold_change=PCR_Product_1_df_log2_fold_change.loc[(insert range)]


# In[168]:


len(PCR_Overlap_P1_log2_fold_change)


# In[169]:


#P1_Sort_AA_overlap_S=pd.Series(PCR_Overlap_P1_log2_fold_change['Fold Change Sort AA Matches'])


# In[170]:


#P2_Sort_AA_overlap_S=pd.Series(PCR_Overlap_P2_log2_fold_change['Fold Change Sort AA Matches'])


# In[171]:


##do this for the other Product 2

PCR_Overlap_P1_log2_fold_change.set_index(['Position Number','Codon','Amino'], inplace=True)
PCR_Overlap_P1_log2_fold_change.rename({'Fold Change Sort AA Matches':'P1 Fold Change Sort AA Matches', 'Fold Change Sort No AA Matches':'P1 Fold Change Sort No AA Matches'}, axis=1, inplace=True)


# In[172]:


PCR_Overlap_P1_log2_fold_change.drop(['Position', 'PCR Product'], axis=1, inplace=True)


# In[173]:


PCR_Overlap_P2_log2_fold_change.set_index(['Position Number','Codon','Amino'], inplace=True)
PCR_Overlap_P2_log2_fold_change.rename({'Fold Change Sort AA Matches':'P2 Fold Change Sort AA Matches', 'Fold Change Sort No AA Matches':'P2 Fold Change Sort No AA Matches'}, axis=1, inplace=True)


# In[174]:


PCR_Overlap_P2_log2_fold_change.drop(['Position', 'PCR Product'], axis=1, inplace=True)


# In[175]:


Overlap_log2_fold_change=pd.concat([PCR_Overlap_P1_log2_fold_change,PCR_Overlap_P2_log2_fold_change], axis=1)


# In[176]:

# In[177]:


Overlap_df_index_list=list(Overlap_log2_fold_change.index)
List_Analysis_5=[]
S_Analysis_df_index_5=pd.Series(Overlap_df_index_list)
position_from_index_5=S_Analysis_df_index_5.str[0]
AA_from_index_5=S_Analysis_df_index_5.str[2]

for p in range(len(S_Analysis_df_index_5)):
    value_p_5=(position_from_index_5[p], AA_from_index_5[p])
    List_Analysis_5.append(value_p_5)

# In[178]:

# In[179]:


list_of_WT_AA_2= []
for i in range(len(S_rows_to_remove_sep)):
    codon = Codons_from_index_remove[i]
    list_value_i_7=(Positions_from_index_remove[i], codontab[codon])
    list_of_WT_AA_2.append(list_value_i_7)


# In[180]:


WT_AA_loop_2= []
for x in list_of_WT_AA_2:
    WT_AA_loop_2.append(x)


# In[181]:


legend6=[]
for i in List_Analysis_5:
    if i == ('771', 'S'):
        legend6.append('Known mutant') 
    elif i == ('771', 'G'):
        legend6.append('Known mutant')
    elif i[1] == ('*'):
        legend6.append('Stop_Codon')
    elif i in WT_AA_loop_2:
        legend6.append('WT_AA')
    else:
        legend6.append('unknown')


# In[182]:

# In[183]:


size_legend_6=[]
for legend in legend6:
    if legend == ('Known mutant'):
        size_legend_6.append(40)
    if legend == ('Stop_Codon'):
        size_legend_6.append(10)
    if legend == ('WT_AA'):
        size_legend_6.append(10)
    if legend == ('unknown'):
        size_legend_6.append(5)


# In[184]:

# In[185]:


palette={'unknown':'#999999', 'WT_AA':'#377eb8', 'Stop_Codon':'#e41a1c' , 'Known mutant':'#984ea3'}
plot=sns.JointGrid(data=Overlap_log2_fold_change ,x='P1 Fold Change Sort AA Matches', y='P2 Fold Change Sort AA Matches', height=10, ratio=2)
plot.plot_joint(sns.scatterplot, hue=legend6, palette=palette, size=size_legend_6, sizes=(30,100))
plot.plot_marginals(sns.histplot, hue= legend6, palette=palette, multiple='layer',hue_order=['Known mutant','Stop_Codon', 'WT_AA', 'unknown'] , binwidth=1, kde=True, element='step', stat='density')
plot.ax_joint.axvline(x=-9.96, color='red', linestyle='dashed')
plot.ax_joint.axhline(y=-9.96, color='red', linestyle='dashed')
#plot.ax_joint.axhline(y=6.6, color='black', linestyle='dashed')
#plot.ax_joint.axhline(y=3.3, color='black', linestyle='dashed')
#plot.ax_joint.set_xlim(-22,0)
#plot.ax_joint.set_ylim(-2,14.5)
plot.ax_joint.plot([-10,10], [-10,10], 'b-', linewidth = 2)
plot.ax_joint.set_xlabel('Log2 P1 Fold Change Sort AA', fontsize='x-large')
plot.ax_joint.set_ylabel('Log2 P2 Fold Change Sort AA', fontsize='x-large')
#x=Known_mutant_log2_Raw_Totals_category_df['Pre-Sort']
#y=Known_mutant_log2_Raw_Totals_category_df['Sort AA']
plot.ax_joint.annotate(insert annotation)

handles=[]
for label, color in palette.items():
    handle=matplotlib.patches.Patch(color=color, label=label )
    handles.append(handle)
plot.ax_joint.legend(handles=handles)
sns.move_legend(plot.ax_joint, loc=(1.2,1))


# In[186]:

# In[187]:


palette={'unknown':'#999999', 'WT_AA':'#377eb8', 'Stop_Codon':'#e41a1c' , 'Known mutant':'#984ea3'}
plot=sns.JointGrid(data=Overlap_log2_fold_change ,x='P1 Fold Change Sort No AA Matches', y='P2 Fold Change Sort No AA Matches', height=10, ratio=2)
plot.plot_joint(sns.scatterplot, hue=legend6, palette=palette, size=size_legend_6, sizes=(30,100))
plot.plot_marginals(sns.histplot, hue= legend6, palette=palette, multiple='layer',hue_order=['Known mutant','Stop_Codon', 'WT_AA', 'unknown'] , binwidth=1, kde=True, element='step', stat='density')
plot.ax_joint.axvline(x=-9.96, color='red', linestyle='dashed')
plot.ax_joint.axhline(y=-9.96, color='red', linestyle='dashed')
#plot.ax_joint.axhline(y=6.6, color='black', linestyle='dashed')
#plot.ax_joint.axhline(y=3.3, color='black', linestyle='dashed')
#plot.ax_joint.set_xlim(-22,0)
#plot.ax_joint.set_ylim(-2,14.5)
plot.ax_joint.plot([-10,10], [-10,10], 'b-', linewidth = 2)
plot.ax_joint.set_xlabel('Log2 P1 Fold Change Sort No AA', fontsize='x-large')
plot.ax_joint.set_ylabel('Log2 P2 Fold Change Sort No AA', fontsize='x-large')
#x=Known_mutant_log2_Raw_Totals_category_df['Pre-Sort']
#y=Known_mutant_log2_Raw_Totals_category_df['Sort AA']
plot.ax_joint.annotate(insert annotation)

handles=[]
for label, color in palette.items():
    handle=matplotlib.patches.Patch(color=color, label=label )
    handles.append(handle)
plot.ax_joint.legend(handles=handles)
sns.move_legend(plot.ax_joint, loc=(1.2,1))


# In[188]:


log2_Analysis_df_fold_change['Position Number']=log2_Analysis_df_fold_change['Position Number'].astype(int)
log2_Analysis_df_fold_change.sort_values(by=['Position Number'], ascending=True, inplace=True)


# In[189]:


Analysis_df_index_list=list(log2_Analysis_df_fold_change.index)

# In[190]:


List_Analysis=[]
S_Analysis_df_index=pd.Series(Analysis_df_index_list)
position_from_index=S_Analysis_df_index.str[0]
AA_from_index=S_Analysis_df_index.str[2]

for p in range(len(S_Analysis_df_index)):
    value_p=(position_from_index[p], AA_from_index[p])
    List_Analysis.append(value_p)


# In[191]:


WT_AA_loop= []
for x in list_of_WT_AA:
    WT_AA_loop.append(x)


# In[192]:

# In[193]:


legend2=[]
for i in List_Analysis:
    if i  == ('585P1', 'P'):
        legend2.append('Known mutant')
    elif i == ('603P1', 'I'):
        legend2.append('Known mutant')
    elif i == ('615P1', 'H'):
        legend2.append('Known mutant')
    elif i == ('618P1', 'E'):
        legend2.append('Known mutant')
    elif i == ('618P1', 'D'):
        legend2.append('Known mutant')
    elif i  == ('972P2', 'A'):
        legend2.append('Known mutant')
    elif i == ('981P2', 'E'):
        legend2.append('Known mutant')
    elif i == ('981P2', 'W'):
        legend2.append('Known mutant')
    elif i == ('981P2', 'R'):
        legend2.append('Known mutant')
    elif i == ('771P1', 'S'):
        legend2.append('Known mutant') 
    elif i == ('771P2', 'S'):
        legend2.append('Known mutant')
    elif i == ('771P1', 'G'):
        legend2.append('Known mutant')
    elif i == ('771P2', 'G'):
        legend2.append('Known mutant')
    elif i == ('1122P2', 'P'):
        legend2.append('Known mutant') 
    elif i == ('1065P2', 'D'):
        legend2.append('Known mutant')
    elif i[1] == ('*'):
        legend2.append('Stop_Codon')
    elif i in WT_AA_loop:
        legend2.append('WT_AA')
    else:
        legend2.append('unknown')


# In[194]:


legend2.count('WT_AA')


# In[195]:


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


# In[196]:

# In[197]:


Category_S=pd.Series(legend2)


# In[198]:


Category_S.index=log2_Analysis_df_fold_change.index


# In[199]:


log_Raw_Totals_df_index_list=list(log_Raw_Totals_df.index)


# In[200]:


List_Analysis_2=[]
S_Raw_Totals_index_2=pd.Series(log_Raw_Totals_df_index_list)
position_from_index_2=S_Raw_Totals_index_2.str[0]
AA_from_index_2=S_Raw_Totals_index_2.str[2]

for p in range(len(S_Raw_Totals_index_2)):
    value_p=(position_from_index_2[p], AA_from_index_2[p])
    List_Analysis_2.append(value_p)

# In[201]:


legend3=[]
for i in List_Analysis_2:
    if i  == ('585P1', 'P'):
        legend3.append('Known mutant')
    elif i == ('603P1', 'I'):
        legend3.append('Known mutant')
    elif i == ('615P1', 'H'):
        legend3.append('Known mutant')
    elif i == ('618P1', 'E'):
        legend3.append('Known mutant')
    elif i == ('618P1', 'D'):
        legend3.append('Known mutant')
    elif i  == ('972P2', 'A'):
        legend3.append('Known mutant')
    elif i == ('981P2', 'E'):
        legend3.append('Known mutant')
    elif i == ('981P2', 'W'):
        legend3.append('Known mutant')
    elif i == ('981P2', 'R'):
        legend3.append('Known mutant')
    elif i == ('771P1', 'S'):
        legend3.append('Known mutant') 
    elif i == ('771P2', 'S'):
        legend3.append('Known mutant')
    elif i == ('771P1', 'G'):
        legend3.append('Known mutant')
    elif i == ('771P2', 'G'):
        legend3.append('Known mutant')
    elif i == ('1122P2', 'P'):
        legend3.append('Known mutant') 
    elif i == ('1065P2', 'D'):
        legend3.append('Known mutant')
    elif i[1] == ('*'):
        legend3.append('Stop_Codon')
    elif i in WT_AA_loop:
        legend3.append('WT_AA')
    else:
        legend3.append('unknown')


# In[202]:


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


# In[203]:


Category_S_2=pd.Series(legend3)


# In[204]:


Category_S_2.index=log_Raw_Totals_df.index


# In[205]:


log2_Raw_Totals_category_df=pd.concat([log_Raw_Totals_df, Category_S_2], axis=1)
log2_Raw_Totals_category_df.rename(columns = {0:'Category'}, inplace=True)


# In[206]:


Amino_Acid_Position=list((log2_Raw_Totals_category_df['Position Number']/3)+1)


# In[207]:


log2_Raw_Totals_category_df.index.get_level_values(2)


# In[208]:


S_WT_AA=pd.Series(list_of_WT_AA)
WT_position=S_WT_AA.str[0]
WT_AA=S_WT_AA.str[1]


# In[209]:


log2_Raw_Totals_category_df['Mutant Name']= (log2_Raw_Totals_category_df['Position Number']/3)+1
log2_Raw_Totals_category_df['Mutant Name']=log2_Raw_Totals_category_df['Mutant Name'].astype(int)
log2_Raw_Totals_category_df['Mutant Name']=log2_Raw_Totals_category_df['Mutant Name'].astype(str)
log2_Raw_Totals_category_df
log2_Raw_Totals_category_df['Codon']=log2_Raw_Totals_category_df.index.get_level_values(1)
log2_Raw_Totals_category_df['Amino']=log2_Raw_Totals_category_df.index.get_level_values(2)
log2_Raw_Totals_category_df['Full Mutant Name']=log2_Raw_Totals_category_df['Mutant Name'] + log2_Raw_Totals_category_df['Amino'] + '_' +log2_Raw_Totals_category_df['Codon']
log2_Raw_Totals_category_df['Position']=log2_Raw_Totals_category_df['Position'].astype(str)

# In[210]:


S_Positions = pd.Series(log2_Raw_Totals_category_df.index.get_level_values(0))
S_Positions.reindex(index=None)
#S_Positions_sep = S_Positions.str.split('(\d+)', 1, expand=False)
#Positions_1=Known_mutant_log2_Raw_Totals_category_df['Position']
#Positions_only=Positions_1.str[0]

# In[211]:


WT_List=[]
for i in range(len(S_Positions)):
    Positions=S_Positions[i]
    Amino_Acid_WT= WT_dict[Positions]
    WT_List.append(Amino_Acid_WT)


# In[212]:


S_WT_List=pd.Series(WT_List)


# In[213]:


S_WT_List.index=log2_Raw_Totals_category_df.index


# In[214]:


log2_Raw_Totals_category_df=pd.concat([log2_Raw_Totals_category_df, S_WT_List], axis=1)
log2_Raw_Totals_category_df.rename(columns = {0:'WT_Amino'}, inplace=True)


# In[215]:


log2_Raw_Totals_category_df['Full Mutant Name 2']= log2_Raw_Totals_category_df['WT_Amino'] + log2_Raw_Totals_category_df['Full Mutant Name']


# In[216]:


Known_mutant_log2_Raw_Totals_category_df=log2_Raw_Totals_category_df.loc[log2_Raw_Totals_category_df['Category'] == 'Known mutant']


# In[217]:


Mutant_List=list(Known_mutant_log2_Raw_Totals_category_df['Full Mutant Name 2'])


# In[218]:

# In[219]:

# In[221]:


palette={'unknown':'#999999', 'WT_AA':'#377eb8', 'Stop_Codon':'#e41a1c' , 'Known mutant':'#984ea3'}
plot=sns.JointGrid(data=log_Raw_Totals_df, x='Pre-Sort', y='Sort AA', height=10, ratio=2)
plot.plot_joint(sns.scatterplot, hue=legend3, palette=palette, size=size_legend_3, sizes=(30,100))
plot.plot_marginals(sns.histplot, hue=legend3, palette=palette, multiple='layer',hue_order=['Known mutant','Stop_Codon', 'WT_AA', 'unknown'] ,binwidth=1, kde=True, element='step', stat='density')
plot.ax_joint.axvline(x=-1, color='red', linestyle='dashed')
plot.ax_joint.axhline(y=-1, color='red', linestyle='dashed')
plot.ax_joint.axhline(y=6.6, color='black', linestyle='dashed')
#plot.ax_joint.axhline(y=3.3, color='black', linestyle='dashed')
#plot.ax_joint.set_xlim(-22,0)
plot.ax_joint.set_ylim(-2,14.5)
plot.ax_joint.plot([-1.5,13], [-1.5,13], 'b-', linewidth = 2)
plot.ax_joint.set_xlabel('Log2 Raw Total_Matches_Pre-Sort', fontsize='x-large')
plot.ax_joint.set_ylabel('Log2 Raw Total_Matches_Sort_AA', fontsize='x-large')
x=Known_mutant_log2_Raw_Totals_category_df['Pre-Sort']
y=Known_mutant_log2_Raw_Totals_category_df['Sort AA']
for i, txt in enumerate(Mutant_List):
    plot.ax_joint.annotate(txt, ((x[i]+0.3), y[i]))
handles=[]
for label, color in palette.items():
    handle=matplotlib.patches.Patch(color=color, label=label, )
    handles.append(handle)
plot.ax_joint.legend(handles=handles)
sns.move_legend(plot.ax_joint, loc=(1.2,1))


# In[222]:


palette={'unknown':'#999999', 'WT_AA':'#377eb8', 'Stop_Codon':'#e41a1c' , 'Known mutant':'#984ea3'}
plot=sns.JointGrid(data=log_Raw_Totals_df, x='Pre-Sort', y='Sort No AA', height=10, ratio=2)
plot.plot_joint(sns.scatterplot, hue=legend3, palette=palette, size=size_legend_3, sizes=(30,100))
plot.plot_marginals(sns.histplot, hue=legend3, palette=palette, multiple='layer',hue_order=['Known mutant','Stop_Codon', 'WT_AA', 'unknown'] ,binwidth=1, kde=True, element='step', stat='density')
plot.ax_joint.axvline(x=-1, color='red', linestyle='dashed')
plot.ax_joint.axhline(y=-1, color='red', linestyle='dashed')
plot.ax_joint.axhline(y=6.6, color='black', linestyle='dashed')
#plot.ax_joint.axhline(y=3.3, color='black', linestyle='dashed')
#plot.ax_joint.set_xlim(-22,0)
#plot.ax_joint.set_ylim(-22,0)
plot.ax_joint.plot([-1.5,13], [-1.5,13], 'b-', linewidth = 2)
plot.ax_joint.set_xlabel('Log2 Raw Total_Matches_Pre-Sort', fontsize='x-large')
plot.ax_joint.set_ylabel('Log2 Raw Total_Matches_Sort_No_AA', fontsize='x-large')
x=Known_mutant_log2_Raw_Totals_category_df['Pre-Sort']
y=Known_mutant_log2_Raw_Totals_category_df['Sort No AA']
for i, txt in enumerate(Mutant_List):
    plot.ax_joint.annotate(txt, ((x[i]+0.3), y[i]))
handles=[]
for label, color in palette.items():
    handle=matplotlib.patches.Patch(color=color, label=label, )
    handles.append(handle)
plot.ax_joint.legend(handles=handles)
sns.move_legend(plot.ax_joint, loc=(1.2,1))


# In[223]:

# In[224]:


log2_Analysis_df_index_list=list(log2_Analysis_df.index)


# In[225]:


List_Analysis_3=[]
S_Raw_Totals_index_3=pd.Series(log2_Analysis_df_index_list)
position_from_index_3=S_Raw_Totals_index_3.str[0]
AA_from_index_3=S_Raw_Totals_index_3.str[2]

for p in range(len(S_Raw_Totals_index_3)):
    value_p=(position_from_index_3[p], AA_from_index_3[p])
    List_Analysis_3.append(value_p)

# In[226]:


legend4=[]
for i in List_Analysis_3:
    if i  == ('585P1', 'P'):
        legend4.append('Known mutant')
    elif i == ('603P1', 'I'):
        legend4.append('Known mutant')
    elif i == ('615P1', 'H'):
        legend4.append('Known mutant')
    elif i == ('618P1', 'E'):
        legend4.append('Known mutant')
    elif i == ('618P1', 'D'):
        legend4.append('Known mutant')
    elif i  == ('972P2', 'A'):
        legend4.append('Known mutant')
    elif i == ('981P2', 'E'):
        legend4.append('Known mutant')
    elif i == ('981P2', 'W'):
        legend4.append('Known mutant')
    elif i == ('981P2', 'R'):
        legend4.append('Known mutant')
    elif i == ('771P1', 'S'):
        legend4.append('Known mutant') 
    elif i == ('771P2', 'S'):
        legend4.append('Known mutant')
    elif i == ('771P1', 'G'):
        legend4.append('Known mutant')
    elif i == ('771P2', 'G'):
        legend4.append('Known mutant')
    elif i == ('1122P2', 'P'):
        legend4.append('Known mutant') 
    elif i == ('1065P2', 'D'):
        legend4.append('Known mutant')
    elif i[1] == ('*'):
        legend4.append('Stop_Codon')
    elif i in WT_AA_loop:
        legend4.append('WT_AA')
    else:
        legend4.append('unknown')


# In[227]:


size_legend_4=[]
for legend in legend4:
    if legend == ('Known mutant'):
        size_legend_4.append(40)
    if legend == ('Stop_Codon'):
        size_legend_4.append(10)
    if legend == ('WT_AA'):
        size_legend_4.append(10)
    if legend == ('unknown'):
        size_legend_4.append(5)


# In[228]:


Category_S_3=pd.Series(legend4)


# In[229]:


Category_S_3.index=log2_Analysis_df.index


# In[230]:


log2_Analysis_category_df=pd.concat([log2_Analysis_df, Category_S_3], axis=1)
log2_Analysis_category_df.rename(columns = {0:'Category'}, inplace=True)


# In[231]:


log2_Analysis_category_df['Position']=log2_Analysis_category_df.index.get_level_values(0)


# In[232]:


S_Position_5 = pd.Series(log2_Analysis_category_df['Position'], dtype="string")
S_Position_sep_5 = S_Position_5.str.split('(\d+)', 1, expand=False)
Positions_from_index_5=S_Position_sep_5.str[1]
PCR_Product_from_index_5=S_Position_sep_5.str[2]


# In[233]:


log2_Analysis_category_df['Position Number']=Positions_from_index_5
log2_Analysis_category_df['Position Number']=log2_Analysis_category_df['Position Number'].astype(int)


# In[234]:


log2_Analysis_category_df['Mutant Name']= (log2_Analysis_category_df['Position Number']/3)+1
log2_Analysis_category_df['Mutant Name']=log2_Analysis_category_df['Mutant Name'].astype(int)
log2_Analysis_category_df['Mutant Name']=log2_Analysis_category_df['Mutant Name'].astype(str)

log2_Analysis_category_df['Codon']=log2_Analysis_category_df.index.get_level_values(1)
log2_Analysis_category_df['Amino']=log2_Analysis_category_df.index.get_level_values(2)
log2_Analysis_category_df['Full Mutant Name']=log2_Analysis_category_df['Mutant Name'] + log2_Analysis_category_df['Amino'] + '_' + log2_Analysis_category_df['Codon']


# In[235]:


S_Positions_2 = pd.Series(log2_Analysis_category_df.index.get_level_values(0))
S_Positions_2.reindex(index=None)
#S_Positions_sep = S_Positions.str.split('(\d+)', 1, expand=False)
#Positions_1=Known_mutant_log2_Raw_Totals_category_df['Position']
#Positions_1
#Positions_only=Positions_1.str[0]


# In[236]:


WT_List_2=[]
for i in range(len(S_Positions_2)):
    Positions_2=S_Positions_2[i]
    Amino_Acid_WT_2= WT_dict[Positions_2]
    WT_List_2.append(Amino_Acid_WT_2)

# In[237]:


S_WT_List_2=pd.Series(WT_List_2)


# In[238]:


S_WT_List_2.index=log2_Analysis_category_df.index


# In[239]:


log2_Analysis_category_df=pd.concat([log2_Analysis_category_df, S_WT_List_2], axis=1)
log2_Analysis_category_df.rename(columns = {0:'WT_Amino'}, inplace=True)


# In[240]:


log2_Analysis_category_df['Full Mutant Name 2']= log2_Analysis_category_df['WT_Amino'] + log2_Analysis_category_df['Full Mutant Name'


# In[241]:


Known_mutant_log2_Analysis_category_df=log2_Analysis_category_df.loc[log2_Analysis_category_df['Category'] == 'Known mutant']


# In[242]:


Mutant_List_2=list(Known_mutant_log2_Analysis_category_df['Full Mutant Name 2'])


# In[243]:


palette={'unknown':'#999999', 'WT_AA':'#377eb8', 'Stop_Codon':'#e41a1c' , 'Known mutant':'#984ea3'}
plot=sns.JointGrid(data=log2_Analysis_df, x='Normalized Total_Matches_Pre-Sort', y='Normalized Total_Matches_Sort_AA', height=10, ratio=2)
plot.plot_joint(sns.scatterplot, hue=legend4, palette=palette, size=size_legend_4, sizes=(30,100))
plot.plot_marginals(sns.histplot, hue=legend4, palette=palette, multiple='layer',hue_order=['Known mutant','Stop_Codon', 'WT_AA', 'unknown'] ,binwidth=1, kde=True, element='step', stat='density')
plot.ax_joint.axvline(x=-19.9, color='red', linestyle='dashed')
plot.ax_joint.axhline(y=-19.9, color='red', linestyle='dashed')
#plot.ax_joint.axvline(x=-16, color='black', linestyle='dashed')
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
x=Known_mutant_log2_Analysis_category_df['Normalized Total_Matches_Pre-Sort']
y=Known_mutant_log2_Analysis_category_df['Normalized Total_Matches_Sort_AA']
for i, txt in enumerate(Mutant_List_2):
    plot.ax_joint.annotate(txt, ((x[i]+0.3), y[i]))


# In[244]:


palette={'unknown':'#999999', 'WT_AA':'#377eb8', 'Stop_Codon':'#e41a1c' , 'Known mutant':'#984ea3'}
plot=sns.JointGrid(data=log2_Analysis_df, x='Normalized Total_Matches_Pre-Sort', y='Normalized Total_Matches_Sort_No_AA', height=10, ratio=2)
plot.plot_joint(sns.scatterplot, hue=legend4, palette=palette, size=size_legend_4, sizes=(30,100))
plot.plot_marginals(sns.histplot, hue=legend4, palette=palette, multiple='layer',hue_order=['Known mutant','Stop_Codon', 'WT_AA', 'unknown'] ,binwidth=1, kde=True, element='step', stat='density')
plot.ax_joint.axvline(x=-19.9, color='red', linestyle='dashed')
plot.ax_joint.axhline(y=-19.9, color='red', linestyle='dashed')
#plot.ax_joint.axvline(x=-16, color='black', linestyle='dashed')
plot.ax_joint.set_xlim(-22,0)
plot.ax_joint.set_ylim(-22,0)
plot.ax_joint.plot([-22,0], [-22,-0], 'b-', linewidth = 2)
plot.ax_joint.set_xlabel('Log2 Normalized Total_Matches_Pre-Sort', fontsize='x-large')
plot.ax_joint.set_ylabel('Log2 Normalized Total_Matches_Sort_No_AA', fontsize='x-large')
handles=[]
for label, color in palette.items():
    handle=matplotlib.patches.Patch(color=color, label=label)
    handles.append(handle)
plot.ax_joint.legend(handles=handles)
sns.move_legend(plot.ax_joint, loc=(1.2,1))
x=Known_mutant_log2_Analysis_category_df['Normalized Total_Matches_Pre-Sort']
y=Known_mutant_log2_Analysis_category_df['Normalized Total_Matches_Sort_No_AA']
for i, txt in enumerate(Mutant_List_2):
    plot.ax_joint.annotate(txt, ((x[i]+0.3), y[i]))


# In[245]:

# In[246]:


log2_Analysis_fold_change_category_df=pd.concat([log2_Analysis_df_fold_change, Category_S], axis=1)
log2_Analysis_fold_change_category_df.rename(columns = {0:'Category'}, inplace=True)


# In[247]:


log2_Analysis_fold_change_category_df['Mutant Name']= (log2_Analysis_fold_change_category_df['Position Number']/3)+1
log2_Analysis_fold_change_category_df['Mutant Name']=log2_Analysis_fold_change_category_df['Mutant Name'].astype(int)
log2_Analysis_fold_change_category_df['Mutant Name']=log2_Analysis_fold_change_category_df['Mutant Name'].astype(str)

log2_Analysis_fold_change_category_df['Codon']=log2_Analysis_fold_change_category_df.index.get_level_values(1)
log2_Analysis_fold_change_category_df['Amino']=log2_Analysis_fold_change_category_df.index.get_level_values(2)
log2_Analysis_fold_change_category_df['Full Mutant Name']=log2_Analysis_fold_change_category_df['Mutant Name'] + log2_Analysis_fold_change_category_df['Amino'] + '_' + log2_Analysis_fold_change_category_df['Codon']


# In[248]:


S_Positions_3 = pd.Series(log2_Analysis_fold_change_category_df.index.get_level_values(0))
S_Positions_3.reindex(index=None)


# In[249]:


WT_List_3=[]
for i in range(len(S_Positions_3)):
    Positions_3=S_Positions_3[i]
    Amino_Acid_WT_3= WT_dict[Positions_3]
    WT_List_3.append(Amino_Acid_WT_3)

# In[250]:


S_WT_List_3=pd.Series(WT_List_3)


# In[251]:


S_WT_List_3.index=log2_Analysis_fold_change_category_df.index


# In[252]:


log2_Analysis_fold_change_category_df=pd.concat([log2_Analysis_fold_change_category_df, S_WT_List_3], axis=1)
log2_Analysis_fold_change_category_df.rename(columns = {0:'WT_Amino'}, inplace=True)

# In[253]:


log2_Analysis_fold_change_category_df['Full Mutant Name 2']= log2_Analysis_fold_change_category_df['WT_Amino'] + log2_Analysis_fold_change_category_df['Full Mutant Name']

# In[254]:


Known_mutant_log2_Analysis_fold_change_category_df=log2_Analysis_fold_change_category_df.loc[log2_Analysis_fold_change_category_df['Category'] == 'Known mutant']


# In[255]:


Mutant_List_3=list(Known_mutant_log2_Analysis_fold_change_category_df['Full Mutant Name 2'])


# In[256]:

# In[257]:


unknown_df=log2_Analysis_fold_change_category_df[log2_Analysis_fold_change_category_df['Category'].str.contains('unknown')]


# In[258]:


len(unknown_df.loc[(unknown_df['Fold Change Sort AA Matches'] >= 2)& (unknown_df['Fold Change Sort AA Matches'] < 9.45)])


# In[259]:


len(unknown_df.loc[(unknown_df['Fold Change Sort No AA Matches'] >= 2)& (unknown_df['Fold Change Sort No AA Matches'] < 9.45)])


# In[260]:


len(unknown_df.loc[(unknown_df['Fold Change Sort No AA Matches'] >= 4)& (unknown_df['Fold Change Sort No AA Matches'] < 9.45)])


# In[261]:


unknown_df.loc[(unknown_df['Fold Change Sort AA Matches'] >= 4)& (unknown_df['Fold Change Sort AA Matches'] < 9.45)]


# In[262]:


len(unknown_df.loc[(unknown_df['Fold Change Sort AA Matches'] >= 9.45)])


# In[263]:


len(unknown_df.loc[(unknown_df['Fold Change Sort No AA Matches'] >= 9.45)])


# In[264]:


WT_AA_df=log2_Analysis_fold_change_category_df[log2_Analysis_fold_change_category_df['Category'].str.contains('WT_AA')]


# In[265]:


Known_mutants_df=log2_Analysis_fold_change_category_df[log2_Analysis_fold_change_category_df['Category'].str.contains('Known mutant')]

# In[266]:


stop_codons_df=log2_Analysis_fold_change_category_df[log2_Analysis_fold_change_category_df['Category'].str.contains('Stop_Codon')]


# In[267]:


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
plot.set_xlabel('Log2 Fold Change Sort AA/Pre-sort', fontsize='x-large')
plt.ylim(0,0.4)
plt.xlim(-11,10)
plot.axvline(x=4.5, color='black', linestyle='dashed')
plot.axvline(x=2, color='black')
handles=[]
for label, color in palette6.items():
    handle=matplotlib.patches.Patch(color=color, label=label)
    handles.append(handle)

plot.legend(handles=handles,fontsize='xx-large')
sns.move_legend(plot, "best")


# In[268]:


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
plot.set_xlabel('Log2 Fold Change Sort AA/Pre-sort', fontsize='x-large')
plt.ylim(0,0.4)
plt.xlim(-11,10)
plot.axvline(x=4, color='black', linestyle='dashed')
plot.axvline(x=2, color='black')
handles=[]
for label, color in palette2.items():
    handle=matplotlib.patches.Patch(color=color, label=label)
    handles.append(handle)

plot.legend(handles=handles,fontsize='xx-large')
sns.move_legend(plot, "best")


# In[269]:


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


# In[270]:


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


# In[271]:


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
plot.set_xlabel('Log2 Fold Change Sort AA/Pre-sort', fontsize='x-large')
plt.ylim(0,0.4)
plt.xlim(-11,10)

handles=[]
for label, color in palette5.items():
    handle=matplotlib.patches.Patch(color=color, label=label)
    handles.append(handle)

plot.legend(handles=handles,fontsize='xx-large')
sns.move_legend(plot, "best")


# In[272]:


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
plot.set_xlabel('Log2 Fold Change Sort AA/Pre-sort Stop Codons', fontsize='x-large')
plt.ylim(0,0.4)
plt.xlim(-11,10)


# In[273]:


palette6={'WT_AA':'#377eb8' , 'Known mutant':'#984ea3', 'unknown': '#999999'}
A=Known_mutants_df
B=stop_codons_df
C=WT_AA_df
D=unknown_df
plt.figure(figsize=(12,8))
plot=sns.histplot(data=Known_mutants_df, x='Fold Change Sort No AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density', legend=True, color='#984ea3')
#plot=sns.histplot(data=B, x='Fold Change Sort AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density',legend=True, color='#e41a1c')
plot=sns.histplot(data=C, x='Fold Change Sort No AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density',legend=True, color='#377eb8')
plot=sns.histplot(data=D, x='Fold Change Sort No AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density',legend=True, color='#999999')
plot.set_ylabel('Density', fontsize='x-large')
plot.set_xlabel('Log2 Fold Change Sort No AA/Pre-sort', fontsize='x-large')
plt.ylim(0,0.4)
plt.xlim(-11,10)
plot.axvline(x=4, color='black', linestyle='dashed')
plot.axvline(x=2, color='black')
handles=[]
for label, color in palette6.items():
    handle=matplotlib.patches.Patch(color=color, label=label)
    handles.append(handle)

plot.legend(handles=handles,fontsize='xx-large')
sns.move_legend(plot, "best")


# In[274]:


palette2={'WT_AA':'#377eb8' , 'Known mutant':'#984ea3'}
A=Known_mutants_df
B=stop_codons_df
C=WT_AA_df
D=unknown_df
plt.figure(figsize=(12,8))
plot=sns.histplot(data=Known_mutants_df, x='Fold Change Sort No AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density', legend=True, color='#984ea3')
#plot=sns.histplot(data=B, x='Fold Change Sort AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density',legend=True, color='#e41a1c')
plot=sns.histplot(data=C, x='Fold Change Sort No AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density',legend=True, color='#377eb8')
#plot=sns.histplot(data=D, x='Fold Change Sort AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density',legend=True, color='#999999')
plot.set_ylabel('Density', fontsize='x-large')
plot.set_xlabel('Log2 Fold Change Sort No AA/Pre-sort', fontsize='x-large')
plt.ylim(0,0.4)
plt.xlim(-11,10)
plot.axvline(x=4, color='black', linestyle='dashed')
plot.axvline(x=2, color='black')
handles=[]
for label, color in palette2.items():
    handle=matplotlib.patches.Patch(color=color, label=label)
    handles.append(handle)

plot.legend(handles=handles,fontsize='xx-large')
sns.move_legend(plot, "best")


# In[275]:


palette3={'Known mutant':'#984ea3'}
A=Known_mutants_df
B=stop_codons_df
C=WT_AA_df
D=unknown_df
plt.figure(figsize=(12,8))
plot=sns.histplot(data=Known_mutants_df, x='Fold Change Sort No AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density', legend=True, color='#984ea3')
#plot=sns.histplot(data=B, x='Fold Change Sort AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density',legend=True, color='#e41a1c')
#plot=sns.histplot(data=C, x='Fold Change Sort AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density',legend=True, color='#377eb8')
#plot=sns.histplot(data=D, x='Fold Change Sort AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density',legend=True, color='#999999')
plot.set_ylabel('Density', fontsize='x-large')
plot.set_xlabel('Log2 Fold Change Sort No AA/Pre-sort', fontsize='x-large')
plt.ylim(0,0.4)
plt.xlim(-11,10)

handles=[]
for label, color in palette3.items():
    handle=matplotlib.patches.Patch(color=color, label=label)
    handles.append(handle)

plot.legend(handles=handles,fontsize='xx-large')
sns.move_legend(plot, "best")


# In[276]:


palette4={'WT_AA':'#377eb8'}
A=Known_mutants_df
B=stop_codons_df
C=WT_AA_df
D=unknown_df
plt.figure(figsize=(12,8))
#plot=sns.histplot(data=Known_mutants_df, x='Fold Change Sort AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density', legend=True, color='#984ea3')
#plot=sns.histplot(data=B, x='Fold Change Sort AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density',legend=True, color='#e41a1c')
plot=sns.histplot(data=C, x='Fold Change Sort No AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density',legend=True, color='#377eb8')
#plot=sns.histplot(data=D, x='Fold Change Sort AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density',legend=True, color='#999999')
plot.set_ylabel('Density', fontsize='x-large')
plot.set_xlabel('Log2 Fold Change Sort No AA/Pre-sort', fontsize='x-large')
plt.ylim(0,0.4)
plt.xlim(-11,10)

handles=[]
for label, color in palette4.items():
    handle=matplotlib.patches.Patch(color=color, label=label)
    handles.append(handle)

plot.legend(handles=handles,fontsize='xx-large')
sns.move_legend(plot, "best")


# In[277]:


palette5={'unknown': '#999999'}
A=Known_mutants_df
B=stop_codons_df
C=WT_AA_df
D=unknown_df
plt.figure(figsize=(12,8))
#plot=sns.histplot(data=Known_mutants_df, x='Fold Change Sort AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density', legend=True, color='#984ea3')
#plot=sns.histplot(data=B, x='Fold Change Sort AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density',legend=True, color='#e41a1c')
#plot=sns.histplot(data=C, x='Fold Change Sort AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density',legend=True, color='#377eb8')
plot=sns.histplot(data=D, x='Fold Change Sort No AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density',legend=True, color='#999999')
plot.set_ylabel('Density', fontsize='x-large')
plot.set_xlabel('Log2 Fold Change Sort No AA/Pre-sort', fontsize='x-large')
plt.ylim(0,0.4)
plt.xlim(-11,10)

handles=[]
for label, color in palette5.items():
    handle=matplotlib.patches.Patch(color=color, label=label)
    handles.append(handle)

plot.legend(handles=handles,fontsize='xx-large')
sns.move_legend(plot, "best")


# In[278]:


A=Known_mutants_df
B=stop_codons_df
C=WT_AA_df
D=unknown_df
plt.figure(figsize=(12,8))
#plot=sns.histplot(data=Known_mutants_df, x='Fold Change Sort AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density', legend=True, color='#984ea3')
plot=sns.histplot(data=B, x='Fold Change Sort No AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density',legend=True, color='#e41a1c')
#plot=sns.histplot(data=C, x='Fold Change Sort AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density',legend=True, color='#377eb8')
#plot=sns.histplot(data=D, x='Fold Change Sort AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density',legend=True, color='#999999')
plot.set_ylabel('Density', fontsize='x-large')
plot.set_xlabel('Log2 Fold Change Sort No AA/Pre-sort Stop Codons', fontsize='x-large')
plt.ylim(0,0.4)
plt.xlim(-11,10)


# In[279]:


palette7={'Known mutant_AA':'#984ea3', 'Known mutant_No_AA':'black'}
A=Known_mutants_df
B=stop_codons_df
C=WT_AA_df
D=unknown_df
plt.figure(figsize=(12,8))
plot=sns.histplot(data=Known_mutants_df, x='Fold Change Sort No AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density', legend=True, color='black')
plot=sns.histplot(data=Known_mutants_df, x='Fold Change Sort AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density', legend=True, color='#984ea3')
#plot=sns.histplot(data=B, x='Fold Change Sort AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density',legend=True, color='#e41a1c')
#plot=sns.histplot(data=C, x='Fold Change Sort AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density',legend=True, color='#377eb8')
#plot=sns.histplot(data=D, x='Fold Change Sort AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density',legend=True, color='#999999')
plot.set_ylabel('Density', fontsize='x-large')
plot.set_xlabel('Log2 Fold Change Sort No AA/Pre-sort', fontsize='x-large')
plt.ylim(0,0.4)
plt.xlim(-11,10)

handles=[]
for label, color in palette7.items():
    handle=matplotlib.patches.Patch(color=color, label=label)
    handles.append(handle)

plot.legend(handles=handles,fontsize='xx-large')
sns.move_legend(plot, "best")


# In[280]:


palette8={'unknown_AA': '#999999','unknown No AA': 'black'}
A=Known_mutants_df
B=stop_codons_df
C=WT_AA_df
D=unknown_df
plt.figure(figsize=(12,8))
#plot=sns.histplot(data=Known_mutants_df, x='Fold Change Sort AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density', legend=True, color='#984ea3')
#plot=sns.histplot(data=B, x='Fold Change Sort AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density',legend=True, color='#e41a1c')
#plot=sns.histplot(data=C, x='Fold Change Sort AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density',legend=True, color='#377eb8')
plot=sns.histplot(data=D, x='Fold Change Sort AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density',legend=True, color='#999999')
plot=sns.histplot(data=D, x='Fold Change Sort No AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density',legend=True, color='black')
plot.set_ylabel('Density', fontsize='x-large')
plot.set_xlabel('Log2 Fold Change Sort No AA/Pre-sort', fontsize='x-large')
plt.ylim(0,0.4)
plt.xlim(-11,10)

handles=[]
for label, color in palette8.items():
    handle=matplotlib.patches.Patch(color=color, label=label)
    handles.append(handle)

plot.legend(handles=handles,fontsize='xx-large')
sns.move_legend(plot, "best")


# In[281]:

# In[282]:


# In[283]:


log2_Analysis_df_fold_change['Position Number']=log2_Analysis_df_fold_change['Position Number'].astype(int)


# In[284]:


log2_Analysis_df_fold_change['Mutant Name']= (log2_Analysis_df_fold_change['Position Number']/3)+1
log2_Analysis_df_fold_change['Mutant Name']=log2_Analysis_df_fold_change['Mutant Name'].astype(int)
log2_Analysis_df_fold_change['Mutant Name']=log2_Analysis_df_fold_change['Mutant Name'].astype(str)

log2_Analysis_df_fold_change['Codon']=log2_Analysis_df_fold_change.index.get_level_values(1)
log2_Analysis_df_fold_change['Amino']=log2_Analysis_df_fold_change.index.get_level_values(2)
log2_Analysis_df_fold_change['Full Mutant Name']=log2_Analysis_df_fold_change['Mutant Name'] + log2_Analysis_df_fold_change['Amino'] + '_' + log2_Analysis_df_fold_change['Codon']


# In[285]:


plt.figure(figsize=(170,50))
#plt.ylim(0, 30)
plt.xlim(-0.5,413.5)
#plt.xlabel('Position', fontsize=1)
#plt.ylabel('Fold Change Sort AA Matches', fontsize=1)
plt.xticks(rotation=90, fontsize=60)
plt.yticks(fontsize=70)
#sns.set(font_scale=10)
plot=sns.scatterplot(data=log2_Analysis_df_fold_change, x='Position', y='Fold Change Sort AA Matches', hue=legend2, palette=palette, size=size_legend_2, sizes=(600,1000))
plot.set_xlabel('Position', fontsize=70)
plot.set_ylabel('Log2 Fold Change Sort AA Matches/Pre-sort', fontsize=70)
plt.axhline(y=0, color='black', linestyle='dashed')

#plot.axhline(y=4, color='black', linestyle='dashed', linewidth=10)
#plot.axhline(y=2, color='black', linestyle='dashed', linewidth=5)
handles=[]
for label, color in palette.items():
    handle=matplotlib.patches.Patch(color=color, label=label)
    handles.append(handle)

plot.legend(handles=handles,fontsize=200)
sns.move_legend(plot, "best")

x=Known_mutant_log2_Analysis_fold_change_category_df['Position']
y=Known_mutant_log2_Analysis_fold_change_category_df['Fold Change Sort AA Matches']
for i, txt in enumerate(Mutant_List_3):
    plot.annotate(txt, (x[i], (y[i])), fontsize=60)
    
plot.axhline(y=4, color='black', linestyle='dashed', linewidth=10)
plot.axhline(y=2, color='black', linestyle='dashed', linewidth=5)


# In[286]:


plt.figure(figsize=(100,20))
#plt.ylim(0, 30)
plt.xlim(206.5,413.5)
#plt.xlabel('Position', fontsize=1)
#plt.ylabel('Fold Change Sort AA Matches', fontsize=1)
plt.xticks(rotation=90, fontsize=60)
plt.yticks(fontsize=60)
#sns.set(font_scale=10)
plot=sns.scatterplot(data=log2_Analysis_df_fold_change, x='Position', y='Fold Change Sort AA Matches', hue=legend2, palette=palette, size=size_legend_2, sizes=(600,1000))
plot.set_xlabel('Position', fontsize=60)
plot.set_ylabel('Log2 Fold Change Sort AA Matches/Pre-sort', fontsize=50)
plt.axhline(y=0, color='black', linestyle='dashed')

#plot.axhline(y=4, color='black', linestyle='dashed', linewidth=10)
#plot.axhline(y=2, color='black', linestyle='dashed', linewidth=5)
handles=[]
for label, color in palette.items():
    handle=matplotlib.patches.Patch(color=color, label=label)
    handles.append(handle)

plot.legend(handles=handles,fontsize=200)
sns.move_legend(plot, "best")

x=Known_mutant_log2_Analysis_fold_change_category_df['Position']
y=Known_mutant_log2_Analysis_fold_change_category_df['Fold Change Sort AA Matches']
for i, txt in enumerate(Mutant_List_3):
    plot.annotate(txt, (x[i], (y[i])), fontsize=60)
for label in plot.xaxis.get_ticklabels()[::2]:
    label.set_visible(False)    
plot.axhline(y=4, color='black', linestyle='dashed', linewidth=10)
#plot.axhline(y=2, color='black', linestyle='dashed', linewidth=5)


# In[287]:


plt.figure(figsize=(100,20))
#plt.ylim(1.5, 10.5)
plt.xlim(-0.5,206.5)
#plt.xlabel('Position', fontsize=1)
#plt.ylabel('Fold Change Sort AA Matches', fontsize=1)
plt.xticks(rotation=90, fontsize=50)
plt.yticks(fontsize=50)
#sns.set(font_scale=10)
plot=sns.scatterplot(data=log2_Analysis_df_fold_change, x='Position', y='Fold Change Sort AA Matches', hue=legend2, palette=palette, size=size_legend_2, sizes=(600,1000))
plot.set_xlabel('Position', fontsize=70)
plot.set_ylabel('Log2 Fold Change Sort AA Matches/Pre-sort', fontsize=50)
plt.axhline(y=0, color='black', linestyle='dashed')

handles=[]
for label, color in palette.items():
    handle=matplotlib.patches.Patch(color=color, label=label)
    handles.append(handle)

plot.legend(handles=handles,fontsize=200)
sns.move_legend(plot, "best")

for label in plot.xaxis.get_ticklabels()[::2]:
    label.set_visible(False)
x=Known_mutant_log2_Analysis_fold_change_category_df['Position']
y=Known_mutant_log2_Analysis_fold_change_category_df['Fold Change Sort AA Matches']
for i, txt in enumerate(Mutant_List_3):
    plot.annotate(txt, (x[i], (y[i])), fontsize=60)
    
plot.axhline(y=4, color='black', linestyle='dashed', linewidth=10)
#plot.axhline(y=2, color='black', linestyle='dashed', linewidth=5)


# In[288]:


plt.figure(figsize=(170,50))
#plt.ylim(1.5, 10.5)
plt.xlim(64.5,178.5)
#plt.xlabel('Position', fontsize=1)
#plt.ylabel('Fold Change Sort AA Matches', fontsize=1)
plt.xticks(rotation=90, fontsize=60)
plt.yticks(fontsize=70)
#sns.set(font_scale=10)
plot=sns.scatterplot(data=log2_Analysis_df_fold_change, x='Position', y='Fold Change Sort AA Matches', hue=legend2, palette=palette, size=size_legend_2, sizes=(600,1000))
plot.set_xlabel('Position', fontsize=70)
plot.set_ylabel('Log2 Fold Change Sort AA Matches/Pre-sort', fontsize=70)
plt.axhline(y=0, color='black', linestyle='dashed')

handles=[]
for label, color in palette.items():
    handle=matplotlib.patches.Patch(color=color, label=label)
    handles.append(handle)

plot.legend(handles=handles,fontsize=200)
sns.move_legend(plot, "best")

x=Known_mutant_log2_Analysis_fold_change_category_df['Position']
y=Known_mutant_log2_Analysis_fold_change_category_df['Fold Change Sort AA Matches']
for i, txt in enumerate(Mutant_List_3):
    plot.annotate(txt, (x[i], (y[i])), fontsize=60)
    
plot.axhline(y=4, color='black', linestyle='dashed', linewidth=10)
plot.axhline(y=2, color='black', linestyle='dashed', linewidth=5)


# In[289]:


plt.figure(figsize=(12,6))
#plt.ylim(1.5, 10.5)
plt.xlim(154.5,162.5)
#plt.xlabel('Position', fontsize=1)
#plt.ylabel('Fold Change Sort AA Matches', fontsize=1)
plt.xticks(rotation=90, fontsize=10)
plt.yticks(fontsize=10)
#sns.set(font_scale=10)
plot=sns.scatterplot(data=log2_Analysis_df_fold_change, x='Position', y='Fold Change Sort AA Matches', hue=legend2, palette=palette, size=size_legend_2, sizes=(50,100))
plot.set_xlabel('Position', fontsize=10)
plot.set_ylabel('Log2 Fold Change Sort AA Matches/Pre-sort', fontsize=10)
plt.axhline(y=0, color='black', linestyle='dashed')

handles=[]
for label, color in palette.items():
    handle=matplotlib.patches.Patch(color=color, label=label)
    handles.append(handle)

plot.legend(handles=handles,fontsize=200)
sns.move_legend(plot, "best")

x=Known_mutant_log2_Analysis_fold_change_category_df['Position']
y=Known_mutant_log2_Analysis_fold_change_category_df['Fold Change Sort AA Matches']
for i, txt in enumerate(Mutant_List_3):
    plot.annotate(txt, (x[i], (y[i]+0.3)), fontsize=10)
    
plot.axhline(y=4, color='black', linestyle='dashed', linewidth=5)
plot.axhline(y=2, color='black', linestyle='dashed', linewidth=2)


# In[290]:


plt.figure(figsize=(170,50))
#plt.ylim(0, 30)
plt.xlim(-0.5,413.5)
#plt.xlabel('Position', fontsize=1)
#plt.ylabel('Fold Change Sort AA Matches', fontsize=1)
plt.xticks(rotation=90, fontsize=60)
plt.yticks(fontsize=70)
#sns.set(font_scale=10)
plot=sns.scatterplot(data=log2_Analysis_df_fold_change, x='Position', y='Fold Change Sort No AA Matches', hue=legend2, palette=palette, size=size_legend_2, sizes=(600,1000))
plot.set_xlabel('Position', fontsize=70)
plot.set_ylabel('Log2 Fold Change Sort No AA Matches/Pre-sort', fontsize=70)
plt.axhline(y=0, color='black', linestyle='dashed')

#plot.axhline(y=4, color='black', linestyle='dashed', linewidth=10)
#plot.axhline(y=2, color='black', linestyle='dashed', linewidth=5)
handles=[]
for label, color in palette.items():
    handle=matplotlib.patches.Patch(color=color, label=label)
    handles.append(handle)

plot.legend(handles=handles,fontsize=200)
sns.move_legend(plot, "best")

x=Known_mutant_log2_Analysis_fold_change_category_df['Position']
y=Known_mutant_log2_Analysis_fold_change_category_df['Fold Change Sort No AA Matches']
for i, txt in enumerate(Mutant_List_3):
    plot.annotate(txt, (x[i], (y[i])), fontsize=60)
    
plot.axhline(y=4, color='black', linestyle='dashed', linewidth=10)
plot.axhline(y=2, color='black', linestyle='dashed', linewidth=5)


# In[363]:


plt.figure(figsize=(100,20))
#plt.ylim(0, 30)
plt.xlim(-0.5,206.5)
#plt.xlabel('Position', fontsize=1)
#plt.ylabel('Fold Change Sort AA Matches', fontsize=1)
plt.xticks(rotation=90, fontsize=40)
plt.yticks(fontsize=40)
#sns.set(font_scale=10)
plot=sns.scatterplot(data=log2_Analysis_df_fold_change, x='Position', y='Fold Change Sort No AA Matches', hue=legend2, palette=palette, size=size_legend_2, sizes=(600,1000))
plot.set_xlabel('Position', fontsize=60)
plot.set_ylabel('Log2 Fold Change Sort No AA Matches/Pre-sort', fontsize=40)
plt.axhline(y=0, color='black', linestyle='dashed')

#plot.axhline(y=4, color='black', linestyle='dashed', linewidth=10)
#plot.axhline(y=2, color='black', linestyle='dashed', linewidth=5)
handles=[]
for label, color in palette.items():
    handle=matplotlib.patches.Patch(color=color, label=label)
    handles.append(handle)

plot.legend(handles=handles,fontsize=200)
sns.move_legend(plot, "best")
for label in plot.xaxis.get_ticklabels()[::2]:
    label.set_visible(False)
x=Known_mutant_log2_Analysis_fold_change_category_df['Position']
y=Known_mutant_log2_Analysis_fold_change_category_df['Fold Change Sort No AA Matches']
for i, txt in enumerate(Mutant_List_3):
    plot.annotate(txt, (x[i], (y[i])), fontsize=60)
    
#plot.axhline(y=4, color='black', linestyle='dashed', linewidth=10)
#plot.axhline(y=2, color='black', linestyle='dashed', linewidth=5)


# In[364]:


plt.figure(figsize=(100,20))
#plt.ylim(0, 30)
plt.xlim(206.5,413.5)
#plt.xlabel('Position', fontsize=1)
#plt.ylabel('Fold Change Sort AA Matches', fontsize=1)
plt.xticks(rotation=90, fontsize=40)
plt.yticks(fontsize=40)
#sns.set(font_scale=10)
plot=sns.scatterplot(data=log2_Analysis_df_fold_change, x='Position', y='Fold Change Sort No AA Matches', hue=legend2, palette=palette, size=size_legend_2, sizes=(600,1000))
plot.set_xlabel('Position', fontsize=60)
plot.set_ylabel('Log2 Fold Change Sort No AA Matches/Pre-sort', fontsize=40)
plt.axhline(y=0, color='black', linestyle='dashed')

#plot.axhline(y=4, color='black', linestyle='dashed', linewidth=10)
#plot.axhline(y=2, color='black', linestyle='dashed', linewidth=5)
handles=[]
for label, color in palette.items():
    handle=matplotlib.patches.Patch(color=color, label=label)
    handles.append(handle)

plot.legend(handles=handles,fontsize=200)
sns.move_legend(plot, "best")
for label in plot.xaxis.get_ticklabels()[::2]:
    label.set_visible(False)
x=Known_mutant_log2_Analysis_fold_change_category_df['Position']
y=Known_mutant_log2_Analysis_fold_change_category_df['Fold Change Sort No AA Matches']
for i, txt in enumerate(Mutant_List_3):
    plot.annotate(txt, (x[i], (y[i])), fontsize=60)
    
#plot.axhline(y=4, color='black', linestyle='dashed', linewidth=10)
#plot.axhline(y=2, color='black', linestyle='dashed', linewidth=5)


# In[291]:


plt.figure(figsize=(170,50))
plt.ylim(1.5, 10.5)
plt.xlim(-0.5,413.5)
#plt.xlabel('Position', fontsize=1)
#plt.ylabel('Fold Change Sort AA Matches', fontsize=1)
plt.xticks(rotation=90, fontsize=60)
plt.yticks(fontsize=70)
#sns.set(font_scale=10)
plot=sns.scatterplot(data=log2_Analysis_df_fold_change, x='Position', y='Fold Change Sort No AA Matches', hue=legend2, palette=palette, size=size_legend_2, sizes=(600,1000))
plot.set_xlabel('Position', fontsize=70)
plot.set_ylabel('Log2 Fold Change Sort No AA Matches/Pre-sort', fontsize=70)
plt.axhline(y=0, color='black', linestyle='dashed')

#plot.axhline(y=4, color='black', linestyle='dashed', linewidth=10)
#plot.axhline(y=2, color='black', linestyle='dashed', linewidth=5)
handles=[]
for label, color in palette.items():
    handle=matplotlib.patches.Patch(color=color, label=label)
    handles.append(handle)

plot.legend(handles=handles,fontsize=200)
sns.move_legend(plot, "best")

x=Known_mutant_log2_Analysis_fold_change_category_df['Position']
y=Known_mutant_log2_Analysis_fold_change_category_df['Fold Change Sort No AA Matches']
for i, txt in enumerate(Mutant_List_3):
    plot.annotate(txt, (x[i], (y[i])), fontsize=60)
    
plot.axhline(y=4, color='black', linestyle='dashed', linewidth=10)
plot.axhline(y=2, color='black', linestyle='dashed', linewidth=5)


# In[292]:


plt.figure(figsize=(170,50))
#plt.ylim(1.5, 10.5)
plt.xlim(64.5,178.5)
#plt.xlabel('Position', fontsize=1)
#plt.ylabel('Fold Change Sort AA Matches', fontsize=1)
plt.xticks(rotation=90, fontsize=60)
plt.yticks(fontsize=70)
#sns.set(font_scale=10)
plot=sns.scatterplot(data=log2_Analysis_df_fold_change, x='Position', y='Fold Change Sort No AA Matches', hue=legend2, palette=palette, size=size_legend_2, sizes=(600,1000))
plot.set_xlabel('Position', fontsize=70)
plot.set_ylabel('Log2 Fold Change Sort No AA Matches/Pre-sort', fontsize=70)
plt.axhline(y=0, color='black', linestyle='dashed')

handles=[]
for label, color in palette.items():
    handle=matplotlib.patches.Patch(color=color, label=label)
    handles.append(handle)

plot.legend(handles=handles,fontsize=200)
sns.move_legend(plot, "best")

x=Known_mutant_log2_Analysis_fold_change_category_df['Position']
y=Known_mutant_log2_Analysis_fold_change_category_df['Fold Change Sort No AA Matches']
for i, txt in enumerate(Mutant_List_3):
    plot.annotate(txt, (x[i], (y[i])), fontsize=60)
    
plot.axhline(y=4, color='black', linestyle='dashed', linewidth=10)
plot.axhline(y=2, color='black', linestyle='dashed', linewidth=5)


# In[293]:


plt.figure(figsize=(12,6))
#plt.ylim(1.5, 10.5)
plt.xlim(154.5,162.5)
#plt.xlabel('Position', fontsize=1)
#plt.ylabel('Fold Change Sort AA Matches', fontsize=1)
plt.xticks(rotation=90, fontsize=10)
plt.yticks(fontsize=10)
#sns.set(font_scale=10)
plot=sns.scatterplot(data=log2_Analysis_df_fold_change, x='Position', y='Fold Change Sort No AA Matches', hue=legend2, palette=palette, size=size_legend_2, sizes=(50,100))
plot.set_xlabel('Position', fontsize=10)
plot.set_ylabel('Log2 Fold Change Sort No AA Matches/Pre-sort', fontsize=10)
plt.axhline(y=0, color='black', linestyle='dashed')

handles=[]
for label, color in palette.items():
    handle=matplotlib.patches.Patch(color=color, label=label)
    handles.append(handle)

plot.legend(handles=handles,fontsize=200)
sns.move_legend(plot, "best")

x=Known_mutant_log2_Analysis_fold_change_category_df['Position']
y=Known_mutant_log2_Analysis_fold_change_category_df['Fold Change Sort No AA Matches']
for i, txt in enumerate(Mutant_List_3):
    plot.annotate(txt, (x[i], (y[i]+0.3)), fontsize=10)
    
plot.axhline(y=4, color='black', linestyle='dashed', linewidth=5)
plot.axhline(y=2, color='black', linestyle='dashed', linewidth=2)


# In[294]:


Sort_No_AA_Raw_Totals_P1_df=pd.DataFrame({'Raw Sort Count No AA': S_Total_Matches_sort_No_AA})

# In[295]:


Sort_AA_Raw_Totals_P1_df=pd.DataFrame({'Raw Sort Count AA': S_Total_Matches_sort_AA_P1})

# In[296]:


Sort_No_AA_Raw_Totals_P2_df=pd.DataFrame({'Raw Sort Count No AA':S_Total_Matches_sort_No_AA_3})


# In[297]:


Sort_AA_Raw_Totals_P2_df=pd.DataFrame({'Raw Sort Count AA':S_Total_Matches_sort_AA_3})


# In[298]:


Sort_Raw_Totals_df=pd.concat([Sort_No_AA_Raw_Totals_P1_df,Sort_AA_Raw_Totals_P1_df],axis=1)


# In[299]:


Sort_Raw_Totals_df_2=pd.concat([Sort_No_AA_Raw_Totals_P2_df,Sort_AA_Raw_Totals_P2_df], axis=1)


# In[300]:


Sort_Raw_Totals_df=pd.concat([Sort_Raw_Totals_df, Sort_Raw_Totals_df_2])


# In[301]:


Summary_Plate_Sort_df=Raw_Totals_df.join(Sort_Raw_Totals_df, how='right')
Summary_Plate_Sort_df.replace(0.5, 0, inplace=True)
Summary_Plate_Sort_df=Summary_Plate_Sort_df.join(Analysis_df, how='right')


# In[302]:


Summary_Plate_Sort_df.rename(columns={'Normalized Total_Matches_Pre-Sort': 'Raw Normalized Total_Matches_Pre-Sort', 'Normalized Total_Matches_Sort_No_AA': 'Raw Normalized Total_Matches_Sort_No_AA','Normalized Total_Matches_Sort_AA':'Raw Normalized Total_Matches_Sort_AA',  'Fold Change Sort AA Matches':'Raw Fold Change Sort AA Matches', 'Fold Change Sort No AA Matches' : 'Raw Fold Change Sort No AA Matches' }, inplace=True)


# In[303]:


Summary_Plate_Sort_df.drop(['Position', 'Position Number'], axis=1, inplace=True)


# In[304]:


Summary_Plate_Sort_df=Summary_Plate_Sort_df.join(log2_Analysis_df, how='right')

# In[305]:


Summary_Plate_Sort_df.rename(columns={'Normalized Total_Matches_Pre-Sort':'Log2 Normalized Total_Matches_Pre-Sort', 'Normalized Total_Matches_Sort_No_AA':'Log2 Normalized Total_Matches_Sort_No_AA', 'Normalized Total_Matches_Sort_AA': 'Log2 Normalized Total_Matches_Sort_AA' }, inplace=True)


# In[306]:


Summary_Plate_Sort_df=Summary_Plate_Sort_df.join(log2_Analysis_df_fold_change, how='right')

# In[307]:


Summary_Plate_Sort_df.rename(columns={'Fold Change Sort AA Matches': 'Log2 Fold Change Sort AA Matches', 'Fold Change Sort No AA Matches':'Log2 Fold Change Sort No AA Matches' }, inplace=True)


# In[308]:


Summary_Plate_Sort_df.drop(['Position','Position Number','Mutant Name','Codon', 'Amino', 'Full Mutant Name'], axis=1, inplace=True)


# In[309]:


Summary_Plate_Sort_df=pd.concat([Summary_Plate_Sort_df, Category_S], axis=1)
Summary_Plate_Sort_df.rename(columns={0:'Category'}, inplace=True)


# In[310]:

# In[311]:


Analysis_df_by_AA=Analysis_df.groupby(level=[0,2]).sum()


# In[312]:

# In[313]:


Analysis_df_by_AA['Normalized Total_Matches_Pre-Sort'].value_counts()


# In[314]:


Analysis_df_by_AA['Normalized Total_Matches_Sort_No_AA'].value_counts()


# In[315]:


Analysis_df_by_AA['Normalized Total_Matches_Sort_AA'].value_counts()


# In[316]:


Analysis_df_by_AA.drop(['Fold Change Sort AA Matches', 'Fold Change Sort No AA Matches', 'Position Number'],axis=1, inplace=True)
Analysis_df_by_AA['Fold Change Sort AA Matches'] = Analysis_df_by_AA['Normalized Total_Matches_Sort_AA'] /Analysis_df_by_AA['Normalized Total_Matches_Pre-Sort']
Analysis_df_by_AA['Fold Change Sort No AA Matches'] = Analysis_df_by_AA['Normalized Total_Matches_Sort_No_AA'] /Analysis_df_by_AA['Normalized Total_Matches_Pre-Sort']
Analysis_df_by_AA_fold_change=Analysis_df_by_AA.drop(['Normalized Total_Matches_Sort_AA', 'Normalized Total_Matches_Pre-Sort','Normalized Total_Matches_Sort_No_AA'], axis=1, inplace=False)
Analysis_df_by_AA_fold_change['Position']=Analysis_df_by_AA_fold_change.index.get_level_values(0)
Analysis_df_by_AA_fold_change['Amino Acid']=Analysis_df_by_AA_fold_change.index.get_level_values(1)


# In[317]:


Analysis_df_by_AA_fold_change['Fold Change Sort AA Matches'].value_counts()


# In[318]:


Analysis_df_by_AA_fold_change['Fold Change Sort No AA Matches'].value_counts()


# In[319]:


Analysis_df_by_AA_fold_change.isna().sum()


# In[320]:


Analysis_df_fold_change_by_AA_for_log=Analysis_df_by_AA_fold_change.drop(["Position"], axis=1, inplace=False)
Analysis_df_fold_change_by_AA_for_log=Analysis_df_fold_change_by_AA_for_log.drop(["Amino Acid"], axis=1)


# In[321]:


log2_Analysis_df_fold_change_by_AA=Analysis_df_fold_change_by_AA_for_log.replace(0,0.001)
log2_Analysis_df_fold_change_by_AA.replace(np.inf, 700, inplace=True)
log2_Analysis_df_fold_change_by_AA=log2_Analysis_df_fold_change_by_AA.applymap(np.log2)


# In[322]:


log2_Analysis_df_fold_change_by_AA['Position']=log2_Analysis_df_fold_change_by_AA.index.get_level_values(0)
log2_Analysis_df_fold_change_by_AA['Amino Acid']=log2_Analysis_df_fold_change_by_AA.index.get_level_values(1)


# In[323]:


# In[324]:


# In[325]:


S_Position_6 = pd.Series(log2_Analysis_df_fold_change_by_AA['Position'], dtype="string")
S_Position_sep_6 = S_Position_6.str.split('(\d+)', 1, expand=False)
Positions_from_index_6=S_Position_sep_6.str[1]
PCR_Product_from_index_6=S_Position_sep_6.str[2]


# In[326]:


log2_Analysis_df_fold_change_by_AA['Position Number']=Positions_from_index_6
log2_Analysis_df_fold_change_by_AA['Position Number']=log2_Analysis_df_fold_change_by_AA['Position Number'].astype(int)


# In[327]:


log2_Analysis_df_fold_change_by_AA.sort_values(by=['Position Number'], ascending=True, inplace=True)


# In[328]:

# In[329]:


Analysis_df_index_list_2=list(log2_Analysis_df_fold_change_by_AA.index)


# In[330]:


List_Analysis_4=[]
S_Analysis_df_index_4=pd.Series(Analysis_df_index_list_2)
position_from_index_4=S_Analysis_df_index_4.str[0]
AA_from_index_4=S_Analysis_df_index_4.str[1]

for p2 in range(len(S_Analysis_df_index_4)):
    value_p_4=(position_from_index_4[p2], AA_from_index_4[p2])
    List_Analysis_4.append(value_p_4)


# In[331]:


legend5=[]
for i in List_Analysis_4:
    if i  == ('585P1', 'P'):
        legend5.append('Known mutant')
    elif i == ('603P1', 'I'):
        legend5.append('Known mutant')
    elif i == ('615P1', 'H'):
        legend5.append('Known mutant')
    elif i == ('618P1', 'E'):
        legend5.append('Known mutant')
    elif i == ('618P1', 'D'):
        legend5.append('Known mutant')    
    elif i  == ('972P2', 'A'):
        legend5.append('Known mutant')
    elif i == ('981P2', 'E'):
        legend5.append('Known mutant')
    elif i == ('981P2', 'W'):
        legend5.append('Known mutant')
    elif i == ('981P2', 'R'):
        legend5.append('Known mutant')
    elif i == ('771P1', 'S'):
        legend5.append('Known mutant') 
    elif i == ('771P1', 'G'):
        legend5.append('Known mutant')
    elif i == ('771P2', 'S'):
        legend5.append('Known mutant') 
    elif i == ('771P2', 'G'):
        legend5.append('Known mutant')
    elif i == ('1122P2', 'P'):
        legend5.append('Known mutant') 
    elif i == ('1065P2', 'D'):
        legend5.append('Known mutant')
    elif i[1] == ('*'):
        legend5.append('Stop_Codon')
    elif i in WT_AA_loop:
        legend5.append('WT_AA')
    else:
        legend5.append('unknown')


# In[332]:


size_legend_5=[]
for legend in legend5:
    if legend == ('Known mutant'):
        size_legend_5.append(40)
    if legend == ('Stop_Codon'):
        size_legend_5.append(10)
    if legend == ('WT_AA'):
        size_legend_5.append(10)
    if legend == ('unknown'):
        size_legend_5.append(5)


# In[333]:


Category_S_4=pd.Series(legend5)


# In[334]:


Category_S_4.index = log2_Analysis_df_fold_change_by_AA.index


# In[335]:


log2_Analysis_df_fold_change_by_AA_category_df=pd.concat([log2_Analysis_df_fold_change_by_AA, Category_S_4], axis=1)


# In[336]:


log2_Analysis_df_fold_change_by_AA_category_df.rename(columns = {0:'Category'}, inplace=True)


# In[337]:


unknown_df_2=log2_Analysis_df_fold_change_by_AA_category_df[log2_Analysis_df_fold_change_by_AA_category_df['Category'].str.contains('unknown')]


# In[338]:


known_mutants_df_2=log2_Analysis_df_fold_change_by_AA_category_df[log2_Analysis_df_fold_change_by_AA_category_df['Category'].str.contains('Known mutant')]


# In[339]:


WT_AA_df_2=log2_Analysis_df_fold_change_by_AA_category_df[log2_Analysis_df_fold_change_by_AA_category_df['Category'].str.contains('WT_AA')]


# In[340]:


stop_codons_df_2=log2_Analysis_df_fold_change_by_AA_category_df[log2_Analysis_df_fold_change_by_AA_category_df['Category'].str.contains('Stop_Codon')]


# In[341]:


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


# In[342]:


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
plot.set_xlabel('Log2 Fold Change Sort AA/Pre-sort', fontsize='x-large')
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


# In[343]:


palette6={'WT_AA':'#377eb8' , 'Known mutant':'#984ea3', 'unknown': '#999999'}
A=known_mutants_df_2
B=stop_codons_df_2
C=WT_AA_df_2
D=unknown_df_2
plt.figure(figsize=(12,8))
plot=sns.histplot(data=A, x='Fold Change Sort No AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density', legend=True, color='#984ea3')
#plot=sns.histplot(data=B, x='Fold Change Sort AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density',legend=True, color='#e41a1c')
plot=sns.histplot(data=C, x='Fold Change Sort No AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density',legend=True, color='#377eb8')
plot=sns.histplot(data=D, x='Fold Change Sort No AA Matches' , binwidth=0.5,  multiple='layer',kde=True, stat='density',legend=True, color='#999999')
plot.set_ylabel('Density', fontsize='x-large')
plot.set_xlabel('Log2 Fold Change Sort No AA/Pre-sort', fontsize='x-large')
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


# In[344]:


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


# In[345]:


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


# In[346]:


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


# In[347]:


S_Positions_4 = pd.Series(log2_Analysis_df_fold_change_by_AA_category_df.index.get_level_values(0))
S_Positions_4.reindex(index=None)


# In[348]:


WT_List_4=[]
for i in range(len(S_Positions_4)):
    Positions_4=S_Positions_4[i]
    Amino_Acid_WT_4= WT_dict[Positions_4]
    WT_List_4.append(Amino_Acid_WT_4)

# In[349]:


S_WT_List_4=pd.Series(WT_List_4)


# In[350]:


S_WT_List_4.index=log2_Analysis_df_fold_change_by_AA_category_df.index


# In[351]:


log2_Analysis_df_fold_change_by_AA_category_df=pd.concat([log2_Analysis_df_fold_change_by_AA_category_df, S_WT_List_4], axis=1)
log2_Analysis_df_fold_change_by_AA_category_df.rename(columns = {0:'WT_Amino'}, inplace=True)


# In[352]:


log2_Analysis_df_fold_change_by_AA_category_df['Mutant Name']= (log2_Analysis_df_fold_change_by_AA_category_df['Position Number']/3)+1
log2_Analysis_df_fold_change_by_AA_category_df['Mutant Name']=log2_Analysis_df_fold_change_by_AA_category_df['Mutant Name'].astype(int)
log2_Analysis_df_fold_change_by_AA_category_df['Mutant Name']=log2_Analysis_df_fold_change_by_AA_category_df['Mutant Name'].astype(str)
log2_Analysis_df_fold_change_by_AA_category_df['Full Mutant Name']=log2_Analysis_df_fold_change_by_AA_category_df['Mutant Name'] + log2_Analysis_df_fold_change_by_AA_category_df['Amino Acid'] 


# In[353]:


log2_Analysis_df_fold_change_by_AA_category_df['Full Mutant Name 2']= log2_Analysis_df_fold_change_by_AA_category_df['WT_Amino'] + log2_Analysis_df_fold_change_by_AA_category_df['Full Mutant Name']


# In[354]:


Known_mutant_log2_Analysis_fold_change_by_AA_category_df=log2_Analysis_df_fold_change_by_AA_category_df.loc[log2_Analysis_df_fold_change_by_AA_category_df['Category'] == 'Known mutant']


# In[355]:


Mutant_List_4=list(Known_mutant_log2_Analysis_fold_change_by_AA_category_df['Full Mutant Name 2'])


# In[356]:


plt.figure(figsize=(170,50))
#plt.ylim(0, 30)
plt.xlim(-0.5,413.5)
#plt.xlabel('Position', fontsize=1)
#plt.ylabel('Fold Change Sort AA Matches', fontsize=1)
plt.xticks(rotation=90, fontsize=60)
plt.yticks(fontsize=70)
#sns.set(font_scale=10)
plot=sns.scatterplot(data=log2_Analysis_df_fold_change_by_AA, x='Position', y='Fold Change Sort AA Matches', hue=legend5, palette=palette, size=size_legend_5, sizes=(600,1000))
plot.set_xlabel('Position', fontsize=70)
plot.set_ylabel('Log2 Fold Change Sort AA Matches/Pre-sort', fontsize=70)
plt.axhline(y=0, color='black', linestyle='dashed')

#plot.axhline(y=4, color='black', linestyle='dashed', linewidth=10)
#plot.axhline(y=2, color='black', linestyle='dashed', linewidth=5)
handles=[]
for label, color in palette.items():
    handle=matplotlib.patches.Patch(color=color, label=label)
    handles.append(handle)

plot.legend(handles=handles,fontsize=200)
sns.move_legend(plot, "best")

x=Known_mutant_log2_Analysis_fold_change_by_AA_category_df['Position']
y=Known_mutant_log2_Analysis_fold_change_by_AA_category_df['Fold Change Sort AA Matches']
for i, txt in enumerate(Mutant_List_4):
    plot.annotate(txt, (x[i], (y[i])), fontsize=60)
    
plot.axhline(y=4, color='black', linestyle='dashed', linewidth=10)
plot.axhline(y=2, color='black', linestyle='dashed', linewidth=5)


# In[357]:


plt.figure(figsize=(170,50))
#plt.ylim(0, 30)
plt.xlim(-0.5,413.5)
#plt.xlabel('Position', fontsize=1)
#plt.ylabel('Fold Change Sort AA Matches', fontsize=1)
plt.xticks(rotation=90, fontsize=60)
plt.yticks(fontsize=70)
#sns.set(font_scale=10)
plot=sns.scatterplot(data=log2_Analysis_df_fold_change_by_AA, x='Position', y='Fold Change Sort No AA Matches', hue=legend5, palette=palette, size=size_legend_5, sizes=(600,1000))
plot.set_xlabel('Position', fontsize=70)
plot.set_ylabel('Log2 Fold Change Sort No AA Matches/Pre-sort', fontsize=70)
plt.axhline(y=0, color='black', linestyle='dashed')

#plot.axhline(y=4, color='black', linestyle='dashed', linewidth=10)
#plot.axhline(y=2, color='black', linestyle='dashed', linewidth=5)
handles=[]
for label, color in palette.items():
    handle=matplotlib.patches.Patch(color=color, label=label)
    handles.append(handle)

plot.legend(handles=handles,fontsize=200)
sns.move_legend(plot, "best")

x=Known_mutant_log2_Analysis_fold_change_by_AA_category_df['Position']
y=Known_mutant_log2_Analysis_fold_change_by_AA_category_df['Fold Change Sort No AA Matches']
for i, txt in enumerate(Mutant_List_4):
    plot.annotate(txt, (x[i], (y[i])), fontsize=60)
    
plot.axhline(y=4, color='black', linestyle='dashed', linewidth=10)
plot.axhline(y=2, color='black', linestyle='dashed', linewidth=5)

