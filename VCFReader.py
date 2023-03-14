#!/usr/bin/python3#-*- coding : utf-8 -*-
import sys, os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

base_filename = 'VCF1Filtre.txt'
base_filename2 = 'VCF2Filtre.txt'#Output file name
#os.remove(base_filename)  #Delete to refresh file if already existing
#os.remove(base_filename2) 
vcf1=[] #List of every line in the file once filtered(no header, no #CHROM line), used to create the first dataframe (input vcf file 1)
freq1=[] #list with only Allele frequency after extracting it from the column "INFO"
type1=[]#list with only mutation types after extracting it from the column "INFO"
lenght1=[]#list with only mutation lenghts after extracting it from the column "INFO"
DR1 = []
DV1 = []
Coverage1 = []

vcf2=[] #IDEM but for the 2nd dataframe (input vcf file 2)
freq2=[]
type2=[]
lenght2=[]
DR2= []
DV2 = []
Coverage2 = []


FileName = sys.argv[1] #extract 1st file from arguments
FileName2 = sys.argv[2]


fd = open(FileName, "r") 
lines = fd.readlines() #go throught the whole file
fichier = open(FileName, "w") 
for line in lines :
    x = re.findall("^##", line) #Filter the header
    y = re.findall("#CHROM", line)  #Filter the #Chrom line
    if not x and not y: #re-Write the file with only the necessary data
           fichier.write(line)
fd.close()  #Close file after usage
fichier.close()  #Close file after usage

fichierF = open(FileName, "r") #Open the file containing only the needed infos     
reads = fichierF.readlines() #Go throught the whole file
for line in reads :
    tmp = line.split('\t') #Split the file into columns 
    vcf1.append(tmp) #Store the data into a list, dataframe1 !

fd2 = open(FileName2, "r") 
lines2 = fd2.readlines() #go throught the whole file
fichier2 = open(FileName2, "w") 
for line in lines2 :
    x = re.findall("^##", line) #Filter the header
    y = re.findall("#CHROM", line)  #Filter the #Chrom line
    if not x and not y: #re-Write the file with only the necessary data
           fichier2.write(line)
fd2.close()  #Close file after usage
fichier2.close()  #Close file after usage

fichierF2 = open(FileName2, "r") #Open the file containing only the needed infos     
reads2 = fichierF2.readlines() #Go throught the whole file
for line in reads2 :
    tmp2 = line.split('\t') #Split the file into columns 
    vcf2.append(tmp2) #Store the data into a list, dataframe1 !
    
 
    
    
df = pd.DataFrame(vcf1, columns =['ID', 'POS', 'REF', 'ALT', 'ALLELE', 'QUAL', "FILTER", "INFO", "INU2", "INU3"]) #Create the first dataframe
df
df2 = pd.DataFrame(vcf2, columns =['ID', 'POS', 'REF', 'ALT', 'ALLELE', 'QUAL', "FILTER", "INFO", "INU2", "INU3"]) #Create the first dataframe
df2

for tmp in vcf1 :  #Now we extract frenquency, type and lenght from the vcf1 list : already in "info" column but this way its clearer and prettier...
    inf = tmp[7].split(";")
    freq = inf[7].split("=")[-1] # values of frequency obtained !
    freq1.append(freq)#Store into a list
    type = inf[1].split("=")[-1]# Type of mutations obtained !
    type1.append(type)#Store into a list
    lenght=inf[2].split("=")[-1]# values of mutation lenght obtained !
    lenght1.append(lenght)#Store into a list
    DV1a = tmp[9].split(":")[-1]
    DV1.append(DV1a)
    DR1a = tmp[9].split(":")[-2]
    DR1.append(DR1a)
    Cove=inf[5].split("=")[-1]
    Coverage=Cove.split(",")[-3]# values of mutation lenght obtained !
    Coverage1.append(Coverage) #Store into a list
    
fichierF.close()  

for tmp2 in vcf2 :  #Now we extract frenquency, type and lenght from the vcf1 list : already in "info" column but this way its clearer and prettier...
    inf2 = tmp2[7].split(";")
    freqB = inf2[7].split("=")[-1] # values of frequency obtained !
    freq2.append(freqB)#Store into a list
    typeB = inf2[1].split("=")[-1]# Type of mutations obtained !
    type2.append(typeB)#Store into a list
    lenghtB=inf2[2].split("=")[-1]# values of mutation lenght obtained !
    lenght2.append(lenghtB)#Store into a list
    DV2a = tmp2[9].split(":")[-1]
    DV2.append(DV2a)
    DR2a = tmp2[9].split(":")[-2]
    DR2.append(DR2a) 
    Cove2=inf2[5].split("=")[-1]
    CoverageB=Cove2.split(",")[-3]
    Coverage2.append(CoverageB)
fichierF2.close()    


dfreq = pd.DataFrame(freq1, columns=["FREQUENCY"])
dtype = pd.DataFrame(type1, columns=["MUTATION"])
dlenght = pd.DataFrame(lenght1, columns=["MUTLEN"])
dDR1 = pd.DataFrame(DR1, columns=["DR"])
dDV1 = pd.DataFrame(DV1, columns=["DV"])
dCov1 = pd.DataFrame(Coverage1, columns=["COVERAGE"])

dfreqV2 = pd.DataFrame(freq2, columns=["FREQUENCY"])
dtypeV2 = pd.DataFrame(type2, columns=["MUTATION"])
dlenghtV2 = pd.DataFrame(lenght2, columns=["MUTLEN"])
dDR2 = pd.DataFrame(DR2, columns=["DR"])
dDV2 = pd.DataFrame(DV2, columns=["DV"])
dCov2 = pd.DataFrame(Coverage2, columns=["COVERAGE"])

df = df.drop(columns=["INFO","INU2","INU3"]) #Store into a list Drop the last columns, not necessary for the analysis...
df["FREQUENCY"] =dfreq #concatenante every new dataframe into the first one 
df["MUTATION"] =dtype
df["LENGHT"] =dlenght
df["DR"] = dDR1
df["DV"] = dDV1
df["COVERAGE"] = dCov1

df2 = df2.drop(columns=["INFO","INU2","INU3"]) #Store into a list Drop the last columns, not necessary for the analysis...
df2["FREQUENCY"] =dfreqV2 #concatenante every new dataframe into the first one 
df2["MUTATION"] =dtypeV2
df2["LENGHT"] =dlenghtV2
df2["DR"] = dDR2
df2["DV"] = dDV2
df2["COVERAGE"] = dCov2

with open (base_filename, "a+") as outfile : #Create the output file
    df.to_string(outfile) #Transfert the dataframe into a .txt file

outfile.close() #Close the file, end of data processing for the first file


with open (base_filename2, "a+") as outfile : #Create the output file
    df2.to_string(outfile) #Transfert the dataframe into a .txt file

outfile.close() #Close the file, end of data processing for the first file

frames = [df, df2] #Group the two dataframes and concatenate, so we can check row by row if the position of a mutation is present or not, .compare does not work since the two files might not have the same row numbers
result = pd.concat(frames).reset_index(drop=True)   

column_name = 'POS'
Rep = input("Voulez vous trouver des mutations spécifiques aux deux VCF( exemple : entre les deux P30 chauds ou deux P30 froids) ? (Y/N)") # Optional feature, to find duplicate mutations in the two files
if (Rep =="Y") or (Rep =="y") : 
    def find_SameMut(column_name):
        SameMut = result[result.duplicated(subset=column_name, keep=False)]
        if SameMut.empty:
            print(f"There are no duplicates in the column {column_name}.")
        else:
            print(f"Les duplicats dans la colonne {column_name} sont:")
            print(SameMut)
    find_SameMut(column_name)# Call the function and pass a column names
    
else :
    Tempo = input(print("Voulez vous voir si une mutation particulière est présente dans deux fichiers VCF ( exemple entre P30 chauds et l'ensemble des P20 concaténés) ?(Y/N)"))
    
if (Tempo =="Y") or (Tempo =="y") : 
    value = input("Entrez la position de votre mutation !")  # check for every row in the concatenated dataframes if this value exists or not
    def find_duplicate(value):
        duplicate = result[result.isin([value]).any(axis=1)] #I use the isin function, if isin is true it will add a value to duplicates 
        if duplicate.empty: #check if duplicates is empty (no occurences of our value in two or more rows)
            print(f"La mutation a la position {value} est présente dans une autre ligne.") #Print the rows where the duplicates are found
        else:
            print(f"La mutation a la position {value} n'est PAS présente dans une autre ligne.") 
            print(duplicate)

    find_duplicate(value)  # Call the function and pass a value

else : 
        print("Deux fichiers filtrés et triés de vos VCF ont été créés avec succès")  #Else, only output the two files with filtered and sorted VCF datas     



############### PARTIE VISUALISATION ###############


### profondeur des mutants en fonction de la position ###

positions = df["POS"]
mutant_depths = df["DV"].astype(float)


# df = df.sort_values("POS")
# df = df.sort_values("DV")
# plt.plot(positions, mutant_depths)

# plt.gcf().autofmt_xdate()
# plt.xlabel("Position")
# plt.ylabel("Mutant Depth")
# plt.title("Muta
# nt Depth versus Position")
# plt.show()

### répartition du type de mutation pour les deux échantillons donnés ###


mutation_types = df["MUTATION"]
unique_types = mutation_types.unique()
plt.bar(unique_types, mutation_types.value_counts(), width=0.4)
plt.xlabel("Mutation Type")
plt.ylabel("Count")
plt.title("Histogram of Mutation Type Distribution for P50 samples ")
plt.show()