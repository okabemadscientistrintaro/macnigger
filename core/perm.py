import numpy as np

#generation in lexicographic order - Narayana Pandita
#adapted to fhi-aims input to generate all possibilities of binary metal core
#the example here is for Ni5Ga3

def permute(s):
    a = sorted(s)
    n = len(a) - 1
    while True:
        yield ''.join(a)
        for j in range(n-1, -1, -1):
            if a[j] < a[j + 1]:
                break
        else:
            return
        v = a[j]
        for k in range(n, j, -1):
            if v < a[k]:
                break
        a[j], a[k] = a[k], a[j]
        a[j+1:] = a[j+1:][::-1]
m=0
n=[]
cluster=["Ni","Ni","Ni","Ni","Ni","Ga","Ga","Ga"]   #change elements to be permutated here
for s in permute(cluster): 
	v=[s[i:i+2] for i in range(0, len(s), 2)]
	n.append(v)
with open('geometry1','a+') as f:   #change file here; index guide in line 31  
    lines = f.readlines()
f.close()
lines = [w.replace('\t', '') for w in lines]
lines = [w.replace('\n', '') for w in lines]
for x in range(len(n)):
    y=x+1    #change 1=1,2=57,3=113,4=169,5=225,6=281,7=337,8=393,9=449,10=505... if using several generating geometries (one at a time)
    if y<10:
        yy=str('00'+str(y))
    if y>9:
        yy=str('0'+str(y))
    if y>99:
        yy = str(y)			 
    new= open('Ni5Ga3_'+yy+'.xyz', 'a+')     ##DEFINIR O NOME DOS ARQUIVOS A SEREM GERADOS
    k=np.array(n)[x]
    new.write(str(len(felipe))+'\n'+'\n') #write 8, skip line twice
    for i, l in zip(k,lines):
        new.write(i+'   '+l+'\n') #write element + coordinates
f.close()
#Place elements that will be exchanged above in the script.
#Folder containing file 'geometry.in' WITHOUT ELEMENT INDEXES in the same directory. e.g., line 1 = 0.0034 0.3456 0.2356
#Of course, if n elements (with or without repeating) will be permutated, the geometry.in file must have n lines with their coordinates.
