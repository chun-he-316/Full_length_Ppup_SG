import pandas as pd
from collections import defaultdict

adict=defaultdict(list)
newadict=defaultdict(list)
genelistadict=defaultdict(list)
newgenelistadict=defaultdict(list)
trans2geneadict={}
with open("gene.transcript.match.txt","r")as f:
    lines=f.readlines()
for line in lines:
    line=line.strip()
    transid=line.split("\t")[1]
    trans2geneadict[transid]=line.split("\t")[0]

df2=pd.read_excel('P20210100549_PP_SA_3_5.xlsx',header=0,index_col=False)
for index,row in df2.iterrows():
	print(row["Sequence"])
	genelist=[]
	for trans in str(row["Reference"]).split("/")[:]:
		print(trans)
		geneid=trans2geneadict[trans]
		genelist.append(geneid)

	genelist=list(set(genelist))
	for sep in genelist:
		adict[sep].append(row["Sequence"])
		genelistadict[sep].append(genelist)

def number(gene,adict):
	if gene in adict.keys():
		return(len(adict[gene]))
	else:
		return('-')
#uniquemap
def unique(gene,genelistadict):
	if gene in genelistadict.keys():
		uninum=0
		for i in genelistadict[gene]:
			if len(i)==1:
				uninum=uninum+1
		return(uninum)
	else:
		return('-')

df1=pd.read_excel('gene.xlsx',sheet_name="Sheet1",header=0)
df1['peptide_num']=df1['gene'].apply(lambda x : number(x,adict))
df1['unique_num']=df1['gene'].apply(lambda x :unique(x,genelistadict))
df1.to_csv('gene.addpeptide.csv',sep="\t",index=False)
