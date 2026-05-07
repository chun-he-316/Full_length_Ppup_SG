# This code is used to determine the genomic coordinates of peptides.
import pandas as pd 
prosadict={} 
pepname={} 
df=pd.read_excel('P20210100648_PP_SA_55.xlsx',header=0) 
for index,row in df.iterrows(): 
	sep=row["Sequence"] 
	pepname[sep]=row["peptide_id"].strip() 
	prosadict[sep]=row["Reference"] 
def proslist(seq,prosadict): 
	pros=[] 
	for pro in prosadict[seq].split('/')[:]: 
		pros.append(pro) 
	return(pros) 
fastaadict={} 
with open('Ppup.final.pep.fa','r')as f: 
	lines=f.readlines() 
for line in lines: 
	line=line.strip() 
	if line.startswith('>'): 
		sep=line[1:] 
		fastaadict[sep]='' 
	else: 
		fastaadict[sep]+=line 

beddf=pd.read_csv('Ppup.final.bed12',sep='\t',header=None) 
def cdsposi(pro,beddf): 
	global cdsstart 
	beddf1=beddf[beddf[3]==pro] 
	for index,row in beddf1.iterrows():		 
		if row[5]=='+': 
			cdsstart=row[6] 
		else: 
			cdsstart=row[7] 
	return(cdsstart) 
beddf=pd.read_csv('Ppup.final.bed12',sep='\t',header=None) 
def cdsblock(pro,beddf): 
	global i 
	cdsstart=cdsposi(pro,beddf) 
	beddf1=beddf[beddf[3]==pro] 
	for index,row in beddf1.iterrows(): 
		blocksizelist=row[10].split(',')[:] 
		blockstartlist=row[11].split(',')[:] 
		blocknumber=row[9] 
		for i in range(0,blocknumber): 
			if cdsstart<=int(row[1])+int(blocksizelist[i])+int(blockstartlist[i]): 
				break 
	return(i)				 
def txposi(pro,beddf): 
	global utr 
	global l 
	cdsstart=cdsposi(pro,beddf) 
	beddf1=beddf[beddf[3]==pro] 
	number=cdsblock(pro,beddf) 
	for index,row in beddf1.iterrows(): 
		blocksizelist=row[10].split(',')[:] 
		blockstartlist=row[11].split(',')[:] 
		blocknumber=row[9] 
		if row[5]=='+':			 
			if number==0: 
				utr=int(cdsstart)-int(row[1])-int(blockstartlist[number]) 
			if number>0: 
				l=0 
				for a in range(0,number):						 
					l+=int(blocksizelist[a]) 
					utr=int(cdsstart)-int(row[1])-int(blockstartlist[number])+l					 
		if row[5]=='-': 
			if number==int(blocknumber)-1: 
				utr=int(row[1])+int(blockstartlist[number])+int(blocksizelist[number])-cdsstart 
			if number<int(blocknumber)-1: 
				l=0 
				for a in range(number+1,blocknumber): 
					l+=int(blocksizelist[a]) 
					utr=int(row[1])+int(blockstartlist[number])+int(blocksizelist[number])-cdsstart+l 
	return(utr) 

result=[] 
withoutcds=[] 
for seq in df['Sequence']: 
	proslist1=[] 
	proslist1=proslist(seq,prosadict) 
	for pro in proslist1: 
		if pro in beddf[3].tolist(): 
			if pro in fastaadict.keys(): 
				prostart=fastaadict[pro].find(seq) 
				proend=prostart+len(seq) 
				finalutr=txposi(pro,beddf)	 
				txstart=prostart*3+finalutr 
				txend=proend*3+finalutr 
				peptidename=pepname[seq].strip() 
				result.append([pro,txstart,txend,peptidename]) 
		else: 
			withoutcds.append([pepname[seq],pro]) 

dffinal=pd.DataFrame(result) 
dffinal.to_csv('finalallpeptide.tx',sep='\t',index=False,header=False) 
