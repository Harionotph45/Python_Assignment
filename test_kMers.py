#!/usr/bin/env python3

from kMers import *
import pytest

def test_tidyup():
  ex_tidy = "ACCGTA;:,."
  dna_list = get_dnalst(ex_tidy)
  assert dna_list == "ACCGTA"

def test_create_table():
  import pandas as pd
  ex_dna = "ACCGTA"
  table_seq = kMers(ex_dna,3)
  table = {'k':[1,2,3,'Total'],
         'Obervasi':[4,5,4,13],
         'Possible':[6,5,4,15]},
         'Sequence':[{'A':2,'C':2,'G':1,'T':1},{'AC':1,'CC':1,'CG':1,'GT':1,'TA':1},{'ACC':1,'CCG':1,'CGT':1,'GTA':1},'done']
  assert table_seq == pd.DataFrame(table)

def test_linguistic():
  import pandas as pd
  tab_ling = {'k':[1,2,3,'Total'],
         'Obervasi':[4,5,4,13],
         'Possible':[6,5,4,15]},
         'Sequence':[{'A':2,'C':2,'G':1,'T':1},{'AC':1,'CC':1,'CG':1,'GT':1,'TA':1},{'ACC':1,'CCG':1,'CGT':1,'GTA':1},'done']
  tbl = pd.DataFrame(tab_ling)
  ling_comp = ling_compx(tbl,k)
  assert ling_comp == 13/15
