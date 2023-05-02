import os
import random
import time
wdir=os.path.realpath(__file__).removesuffix('\\'+os.path.basename(__file__))
print(wdir)
os.chdir(wdir)
flist=os.listdir(wdir)
print(*flist)
reads=[]
while True:
    print("Choose first fastq file")
    file=input()
    if file.endswith((".fastq",".fq",".txt")):
        try:
            reads=open(file)
            break
        except: print("Can`t open",file)
    else:
        for ex in [".fastq",".fq",".txt"]:
            try:
                reads=open(file+ex)
                file=file+ex
                print("opened",file)
                break
            except: print("Can`t open",file+ex)
    if reads!=[]:break
n=0
for i in reads:
    print(i)
    n+=1
    if n>7:break
while True:
    print("Choose paired fasta file")
    file2=input()
    if file2.endswith((".fastq",".fq",".txt")):
        try:
            reads2=open(file2)
            break
        except: print("Can`t open",file2)
    else:
        for ex in [".fastq",".fq",".txt"]:
            try:
                reads2=open(file2+ex)
                file2=file2+ex
                print("opened",file2)
                break
            except: print("Can`t open",file2+ex)
    if reads2!=[]:break
n=0
for i in reads2:
    print(i)
    n+=1
    if n>7:break
while True:
    print("Choose percentage of reads to extract (integer 1-100)")
    percentage=input()
    try:
        perc=int(percentage)
        if perc in range(100):break
        else:
            print("Invalid number")
            continue
    except: print("Invalid number")
n=0
newfile1=open(percentage + file, "a")
newfilet2=open(percentage +"t2"+ file2, "a")
reads=open(file)
for i in reads:
    if n==0:
        rn=random.choice(range(100))
        if rn<perc:
            newfile1.write(i)
            newfilet2.write("1")
        else: newfilet2.write("0")
        n+=1
    elif 0<n<3:
        if rn<perc:
            newfile1.write(i)
            newfilet2.write("1")
        else:newfilet2.write("0")
        n+=1
    elif n==3:
        if rn<perc:
            newfile1.write(i)
            newfilet2.write("1")
        else:newfilet2.write("0")
        n=0
newfile1.close()
newfilet2.close()
newfile2=open(percentage + file2, "a")
temp=open(percentage +"t2"+ file2)
for line in temp:
    lineee=line
ttt=0
reads2=open(file2)
for liner in reads2:
    if lineee[ttt]=="1":
        newfile2.write(liner)
    ttt=ttt+1
temp.close()
newfile2.close()
time.sleep(10)