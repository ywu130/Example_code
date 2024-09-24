
#Modules used to create a fasta file that contains sequences of antibiotic gene of interest from all Staphyloccus species
#This fasta file is compatabile with MEGA and was used for alginnment before phylogenetic analysis
import os
import sys
import pandas as pd
import openpyxl
from pathlib import Path


#This is the list that will contain the species followed by the sequence of the gene of interest
#List will be parsed through a dictionary to remove uneccesary varibales that might disrupt alignment
SEP = []
SDR = []
NORs = []
QUA = []

#This is the list that will contain the species followed by the sequence of the gene of interest
#List will be parsed through a dictionary to remove uneccesary varibales that might disrupt alignment
Sep2Dict = {}
NOR2Dict = {}
QUA2Dict = {}
SDR2Dict = {}

#Enter the name of the excel output file from the '' script
#The program makes the assumption the script and the excel file are in the same location
#The user will manually type this into the terminal
print('enter excel file name that contains sequences(e.g.xxxx.xls)')
excel = str(input())
xlsls = excel
df = pd.read_excel(xlsls)

#From the excel file, you can choose which gene you would like to have sequences developed for
#This requires the user to manually type the gene of interest into the terminal and have a understanding of
#genes detected in the excel sheet
print('enter gene of interest(QUA,SEP,SDR,or NOR):')
name = str(input())

#The pandas module is used to create a sequence dataframe from the user's gene of interest of all species
#This section will add the species and sequence to a list that will be transfomed into a dictionary
#Parsing the list through the dictionary will remove uncessary variables that might disrupt analysis
sep = pd.DataFrame(df, columns= ['Species','sepA']) #Specific specifiic sequence List and Dictionary for sepA gene
for row in sep.itertuples():
    SEP.append(row.Species)
    SEP.append(row.sepA)
for value in SEP:
    if '.txt' in value:
        new = value.replace('.txt', '')
        Final = SEP[SEP.index(value)+1].replace("'",'')
        Sep2Dict[new] = Final

sdr = pd.DataFrame(df, columns= ['Species','sdrM'])#Specific specifiic sequence list and dictionary for sdrM gene
for row in sdr.itertuples():
    SDR.append(row.Species)
    SDR.append(row.sdrM)
for value in SDR:
    if '.txt' in value:
        new = value.replace('.txt', '')
        Final = SDR[SDR.index(value)+1].replace("'",'')
        SDR2Dict[new] = Final

nor = pd.DataFrame(df, columns= ['Species','norC'])#Specific specifiic sequence list and fictionary for norC gene
for row in nor.itertuples():
    NORs.append(row.Species)
    NORs.append(row.norC)
for value in NORs:
    if '.txt' in value:
        new = value.replace('.txt', '')
        Final = NORs[NORs.index(value)+1].replace("'",'')
        NOR2Dict[new] = Final

qua = pd.DataFrame(df, columns= ['Species','qacJ'])#Specific specifiic sequence list and dictionary for qacJ gene
for row in qua.itertuples():
    QUA.append(row.Species)
    QUA.append(row.qacJ)
for value in QUA:
    if '.txt' in value:
        new = value.replace('.txt', '')
        Final = QUA[QUA.index(value)+1].replace("'",'')
        QUA2Dict[new] = Final


#This section of the program creates a fasta sequence file from the species sequence-specific dictionary based on the gene of interest entered
#by the user
FILE = name + '.fasta'
input_file = open(FILE, 'w')

#Individual fasta files for specific antibitoic genes of interest
if 'SDR' == name:
    for key in SDR2Dict:
        sequencez = SDR2Dict.get(key)
        if len(sequencez) >= 3:
            input_file.write(">"+ key + '\n')
            input_file.write(sequencez.replace(']','')+'\n')

if 'NOR' in name:
    for key in SDR2Dict:
        sequencez = NOR2Dict.get(key)
        if len(sequencez) >= 3:
            input_file.write(">"+ key + '\n')
            input_file.write(sequencez.replace(']','')+'\n')

if 'QUA' == name:
    for key in QUA2Dict:
        sequencez = QUA2Dict.get(key)
        if len(sequencez) >= 3:
            input_file.write(">"+ key + '\n')
            input_file.write(sequencez.replace(']','')+'\n')

if 'SEP' == name:
    for key in QUA2Dict:
        sequencez = Sep2Dict.get(key)
        if len(sequencez) >= 3:
            input_file.write(">"+ key + '\n')
            input_file.write(sequencez.replace(']','')+'\n')

input_file.close()
