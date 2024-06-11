#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 10 23:16:24 2024

@author: cs
"""

import pandas as pd


# Open list peptides to test

df = pd.read_excel("self-nonself discrimination.xlsx")
df2 = pd.read_excel("self-nonself discrimination.xlsx", sheet_name="self-immunopeptidome")

#%%
viral_peptides = df["Agonist Peptide"]
viral_peptides_9mer = [x for x in viral_peptides if len(x) == 9]

viral_peptides_9mer_H2Db = [x for x, y in zip(df["Agonist Peptide"], df["MHC Class I Restriction"]) if len(x) == 9 and y == "H-2Db"]
viral_peptides_9mer_H2KB = [x for x, y in zip(df["Agonist Peptide"], df["MHC Class I Restriction"]) if len(x) == 9 and y in ["H-2Kb", "H-2kb", "H2-Kb"]]

with open("viral_peptides_H2Db.txt", "w") as F:
    for x in viral_peptides_9mer_H2Db:
        F.write(x)
        F.write("\n")

with open("viral_peptides_H2KB.txt", "w") as F:
    for x in viral_peptides_9mer_H2KB:
        F.write(x)
        F.write("\n")

mouse_peptides_h2db = [x for x, y in zip(df2["Peptide"], df2["MHC"]) if y == "H2Db_B22" and len(x) == 9]
mouse_peptides_h2kb = [x for x, y in zip(df2["Peptide"], df2["MHC"]) if y == "H2Kb_Y3" and len(x) == 9]

with open("mouse_mhc_peptides_h2db.txt", "w") as F:
    for x in mouse_peptides_h2db:
        F.write(x)
        F.write("\n")


with open("mouse_peptides_h2kb.txt", "w") as F:
    for x in mouse_peptides_h2kb:
        F.write(x)
        F.write("\n")
        