#!/usr/bin/env python
# coding: utf-8

# In[1]:


from itertools import product
li = ['A', 'T', 'G', 'C']
keywords = [A+B+C for A,B,C in product(li, repeat=3)]


# In[2]:


from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
reference = SeqIO.read("(fast file name)", "fasta")


# In[3]:


new_seqs_P1_30 = {}
seq = str(reference.seq)
for position in range(438, 824, 3):
    last_chunk_end_P1_30 = min(827, position+33)
    firstchunk_P1_30 = seq[position-30:position]
    lastchunk_P1_30 = seq[position+3:last_chunk_end_P1_30]
    if len(lastchunk_P1_30) < 30:
        add_to_firstchunk_P1 = seq[position-(30-len(lastchunk_P1_30)+30):position-30]
        for codon in keywords:
            new_seq_P1_30 = add_to_firstchunk_P1 + firstchunk_P1_30 + codon + lastchunk_P1_30
            new_seq_name_P1_30 = str(position) + codon
            new_seqs_P1_30[new_seq_name_P1_30] = new_seq_P1_30       
    else:
        pass
        for codon in keywords:
            new_seq_P1_30 = firstchunk_P1_30 + codon + lastchunk_P1_30
            new_seq_name_P1_30 = str(position) + codon
            new_seqs_P1_30[new_seq_name_P1_30] = new_seq_P1_30


# In[4]:


from itertools import islice
dict(islice(new_seqs_P1_30.items(), 0, 1))


# In[5]:

# In[6]:


with open('(name you want for file)', 'w') as ref2:
    ref2.write(str(new_seqs_P1_30))


# In[7]:


new_seqs_P2_30 = {}
seq = str(reference.seq)
for position in range(612, len(seq)-3, 3):
    first_chunk_start = max(612, position-30)
    firstchunk = seq[first_chunk_start:position]
    lastchunk = seq[position+3:position+33]
    if len(firstchunk) < 30:
        add_to_lastchunk = seq[position+33:position+(30-len(firstchunk)+33)]
        for codon in keywords:
            new_seq_P2 = firstchunk + codon + lastchunk + add_to_lastchunk
            new_seq_name_P2 = str(position) + codon
            new_seqs_P2_30[new_seq_name_P2] = new_seq_P2
    elif len(lastchunk) < 30:
        add_to_firstchunk2 = seq[position-(30-len(lastchunk)+30): position-30]
        for codon in keywords:
            new_seq_P2 = add_to_firstchunk2 + firstchunk + codon + lastchunk
            new_seq_name_P2 = str(position) + codon
            new_seqs_P2_30[new_seq_name_P2] = new_seq_P2
    else:
        for codon in keywords:
            new_seq_P2 = firstchunk + codon + lastchunk 
            new_seq_name_P2 = str(position) + codon
            new_seqs_P2_30[new_seq_name_P2] = new_seq_P2
            


# In[8]:

# In[9]:

# In[10]:

# In[11]:

# In[12]:

# In[13]:


from itertools import islice

# In[14]:

# In[15]:

# In[16]:


import re
import gzip


# In[17]:


step = 4
list_FASTQ=[]
with gzip.open('(insert fastq file name)' ,'rt', encoding='utf8', errors='ignore') as handle:
    for lineno, line in enumerate(handle, 1):
        if lineno % step-2 == 0: 
            list_FASTQ.append(line[0:150])


# In[18]:

# In[19]:

# In[20]:


string_list_FASTQ = ',' .join(str(element) for element in list_FASTQ)


# In[21]:

# In[22]:

# In[23]:

# In[24]:

# ## testing out making file into dictionary

# In[25]:


import ast

with open('(file name)') as r:
    reference = r.read()


# In[26]:


ast_reference = ast.literal_eval(reference)


# In[27]:

# In[28]:


from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
reference = SeqIO.read("(FASTA file name)", "fasta")


# In[29]:


new_seqs_P1_30_2 = {}
seq = str(reference.seq)
for position in range(438, 441, 3):
    last_chunk_end_P1_30_2 = min(443, position+33)
    firstchunk_P1_30_2 = seq[position-30:position]
    lastchunk_P1_30_2 = seq[position+3:last_chunk_end_P1_30]
    if len(lastchunk_P1_30_2) < 30:
        add_to_firstchunk_P1_2 = seq[position-(30-len(lastchunk_P1_30)+30):position-30]
        for codon in keywords:
            new_seq_P1_30_2 = add_to_firstchunk_P1 + firstchunk_P1_30 + codon + lastchunk_P1_30
            new_seq_name_P1_30_2 = str(position) + codon
            new_seqs_P1_30_2[new_seq_name_P1_30_2] = new_seq_P1_30_2       
    else:
        pass
        for codon in keywords:
            new_seq_P1_30_2 = firstchunk_P1_30_2 + codon + lastchunk_P1_30_2
            new_seq_name_P1_30_2 = str(position) + codon
            new_seqs_P1_30_2[new_seq_name_P1_30_2] = new_seq_P1_30_2


# In[30]:

# In[31]:


import ast


# In[32]:


with open ('(file name)') as r:
    reference = r.read()
ast_reference = ast.literal_eval(reference)


# In[33]:


from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
reference = SeqIO.read("(file name)", "fasta")


# In[34]:


new_seqs_P1_20 = {}
seq = str(reference.seq)
for position in range(438, 824, 3):
    last_chunk_end_P1_20 = min(827, position+23)
    firstchunk_P1_20 = seq[position-20:position]
    lastchunk_P1_20 = seq[position+3:last_chunk_end_P1_20]
    if len(lastchunk_P1_20) < 20:
        add_to_firstchunk_P1 = seq[position-(20-len(lastchunk_P1_20)+20):position-20]
        for codon in keywords:
            new_seq_P1_20 = add_to_firstchunk_P1 + firstchunk_P1_20 + codon + lastchunk_P1_20
            new_seq_name_P1_20 = str(position) + codon
            new_seqs_P1_20[new_seq_name_P1_20] = new_seq_P1_20       
    else:
        for codon in keywords:
            new_seq_P1_20 = firstchunk_P1_20 + codon + lastchunk_P1_20
            new_seq_name_P1_20 = str(position) + codon
            new_seqs_P1_20[new_seq_name_P1_20] = new_seq_P1_20


# In[35]:

# In[36]:

# In[37]:


new_seqs_P2_20 = {}
seq = str(reference.seq)
for position in range(612, len(seq)-3, 3):
    first_chunk_start_2 = max(612, position-20)
    firstchunk_2 = seq[first_chunk_start_2:position]
    lastchunk_2 = seq[position+3:position+23]
    if len(firstchunk_2) < 20:
        add_to_lastchunk_2 = seq[position+23:position+(20-len(firstchunk_2)+23)]
        for codon in keywords:
            new_seq_P2_2 = firstchunk_2 + codon + lastchunk_2 + add_to_lastchunk_2
            new_seq_name_P2_2 = str(position) + codon
            new_seqs_P2_20[new_seq_name_P2_2] = new_seq_P2_2
    elif len(lastchunk_2) < 20:
        add_to_firstchunk2_2 = seq[position-(20-len(lastchunk_2)+20): position-20]
        for codon in keywords:
            new_seq_P2_2 = add_to_firstchunk2_2 + firstchunk_2 + codon + lastchunk_2
            new_seq_name_P2_2 = str(position) + codon
            new_seqs_P2_20[new_seq_name_P2_2] = new_seq_P2_2
    else:
        for codon in keywords:
            new_seq_P2_2 = firstchunk_2 + codon + lastchunk_2 
            new_seq_name_P2_2 = str(position) + codon
            new_seqs_P2_20[new_seq_name_P2_2] = new_seq_P2_2


# In[38]:

# In[39]:

# In[40]:


from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
reference = SeqIO.read("(insert file name)", "fasta")


# In[41]:


new_seqs_P1_10 = {}
seq = str(reference.seq)
for position in range(438, 824, 3):
    last_chunk_end_P1_10 = min(827, position+13)
    firstchunk_P1_10 = seq[position-10:position]
    lastchunk_P1_10 = seq[position+3:last_chunk_end_P1_10]
    if len(lastchunk_P1_10) < 10:
        add_to_firstchunk_P1_2 = seq[position-(10-len(lastchunk_P1_10)+10):position-10]
        for codon in keywords:
            new_seq_P1_10 = add_to_firstchunk_P1_2 + firstchunk_P1_10 + codon + lastchunk_P1_10
            new_seq_name_P1_10 = str(position) + codon
            new_seqs_P1_10[new_seq_name_P1_10] = new_seq_P1_10       
    else:
        for codon in keywords:
            new_seq_P1_10 = firstchunk_P1_10 + codon + lastchunk_P1_10
            new_seq_name_P1_10 = str(position) + codon
            new_seqs_P1_10[new_seq_name_P1_10] = new_seq_P1_10


# In[45]:


new_seqs_P2_10 = {}
seq = str(reference.seq)
for position in range(612, len(seq)-3, 3):
    first_chunk_start_3 = max(612, position-10)
    firstchunk_3 = seq[first_chunk_start_3:position]
    lastchunk_3 = seq[position+3:position+13]
    if len(firstchunk_3) < 10:
        add_to_lastchunk_3 = seq[position+13:position+(10-len(firstchunk_3)+13)]
        for codon in keywords:
            new_seq_P2_3 = firstchunk_3 + codon + lastchunk_3 + add_to_lastchunk_3
            new_seq_name_P2_3 = str(position) + codon
            new_seqs_P2_10[new_seq_name_P2_3] = new_seq_P2_3
    elif len(lastchunk_3) < 10:
        add_to_firstchunk2_3 = seq[position-(10-len(lastchunk_3)+10): position-10]
        for codon in keywords:
            new_seq_P2_3 = add_to_firstchunk2_3 + firstchunk_3 + codon + lastchunk_3
            new_seq_name_P2_3 = str(position) + codon
            new_seqs_P2_10[new_seq_name_P2_3] = new_seq_P2_3
    else:
        for codon in keywords:
            new_seq_P2_3 = firstchunk_3 + codon + lastchunk_3 
            new_seq_name_P2_3 = str(position) + codon
            new_seqs_P2_10[new_seq_name_P2_3] = new_seq_P2_3



# In[61]:


new_seqs_WT_R206H_10 = {}
seq = str(reference.seq)
for position in range(605, 628, 3):
    firstchunk_4 = seq[605:615]
    lastchunk_4 = seq[618:628]
    WT_codon = seq[615:618]
    new_seq_WT = firstchunk_4 + WT_codon + lastchunk_4
    new_seq_WT_name = '615' + WT_codon
    new_seqs_WT_R206H_10[new_seq_WT_name] = new_seq_WT
    new_seq_R206H = firstchunk_4 + 'CAC' + lastchunk_4
    new_seq_R206H_name = '615' + 'CAC'
    new_seqs_WT_R206H_10[new_seq_R206H_name] = new_seq_R206H
  




