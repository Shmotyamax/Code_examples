import math
import pandas
import os

wdir=os.path.realpath(__file__).removesuffix('\\'+os.path.basename(__file__))
print(wdir)
os.chdir(wdir)
flist=os.listdir(wdir)

KOdict={}
KOMetabolicDict={}
NameSecondaryDict={}
KOSecondaryDict={}
KOHydrolasesDict={}
NameFactorsDict={}
KOAnnotationDict={}
GeneDataDict={}
NameList=[]
PathwaysDict={}
KOIDDict={}

Hand=open("query KO.txt")
for line in Hand:
    if line.startswith("gene"):
        KOline=line[5:].strip()
    else: KOline=line.strip()

    try:
        KOdict[KOline.split("\t")[0]]=KOline.split("\t")[1]
        NameList.append(KOline.split("\t")[0])
    except:
        KOdict[KOline]=None
        NameList.append(KOline)
    #Name   KO
Hand=open("KEGG annotation.txt")
for line in Hand:
    KOAnn=line.strip().split("\t")
    KOAnnotationDict[KOAnn[0]]=KOAnn[1]
    #KO    Annotation

Hand=open("KEGG metabolism annotation.txt")
for line in Hand:
    KOMAnn=line.strip().split("\t")
    KOMetabolicDict[KOMAnn[0]]=KOMAnn[1]
    #KO    Annotation

Hand=open("Secondary metabolism.txt")
for line in Hand:
    NameSec=line.strip().split("\t")
    try:
        NameSecondaryDict[NameSec[0]]=NameSec[1]
    except:
        NameSecondaryDict[NameSec[0]]=False
    #Name    Secondary func

Hand=open("KEGG sec metabolism annotation.txt")
for line in Hand:
    KOSec=line.strip().split("\t")
    KOSecondaryDict[KOSec[0]]=KOSec[1]
    #KO    Secondsry func

Hand=open("Hydrolases.txt")
for line in Hand:
    Hydro=line.strip().split("\t")
    KOHydrolasesDict[Hydro[0]]=Hydro[1]
    #KO    Annotation

Hand=open("Factors and regulation.txt")
for line in Hand:
    NameFac=line.strip().split("\t")
    NameFactorsDict[NameFac[0]]=NameFac[1]
    #Name    Factor

Hand=open("t_data.ctab")
for line in Hand:
    GeneData=line.strip().split("\t")
    if GeneData[0].startswith("t"): continue
    GeneDataDict[int(GeneData[0])]=[GeneData[1],GeneData[3],GeneData[4]]
    #Name    Factor
'''
Hand=open("KOlistIDs")
for line in Hand:
    IDG=line.strip().split("\t")
    KOIDDict[IDG[1]]=IDG[0]
    #KO    Annotation
'''
Hand=open("KOlistPs")
for line in Hand:
    if line.startswith("K")==False:
        Pathway=line.strip()
    else:
        KOPS=line.split()
        for KOP in KOPS:
            PathwaysDict[KOP]=PathwaysDict.get(KOP,[])
            PathwaysDict[KOP].append(Pathway)
for key in PathwaysDict:
    PathwaysDict[key]="    ".join(PathwaysDict[key])
    #KO    Patwhays
Hand.close()



class Gene:
    def __init__(self,ID,Qval,Fcb,Name):
        self.ID=ID
        self.Name=Name
        self.Qval=Qval
        self.Fc=math.log(((math.pow(2,Fcb))-1),2)
        self.KO=None
        self.IsHydrolase=False
        self.IsSecondaryMetabolic=False
        self.IsFactor=False
        self.IsPrimaryMetabolic=False
        self.Description=None
        self.HowRegulated=None
        self.Chromosome=None
        self.ChrPos=None
        self.Pathways=None
        self.KOID=None

for file in flist:
    if file.endswith(".csv"):
        Data=open(file)
        GeneList=[]
        l=0
        for line in Data:
            if line[1]=='"':
                continue
            dataline=line.split(",")
            try:
                GeneList.append(Gene(int(dataline[0].strip('"')),float(dataline[-1]),float(dataline[3]),NameList[l]))
            except:
                GeneList.append(Gene(int(dataline[0].strip('"')),1, 1,NameList[l]))
            l+=1
        for item in GeneList:
            item.Chromosome=GeneDataDict[item.ID][0]
            item.ChrPos=f"{GeneDataDict[item.ID][1]}-{GeneDataDict[item.ID][2]}"
            item.KO=KOdict[item.Name]
            item.IsSecondaryMetabolic = NameSecondaryDict.get(item.Name,False)
            item.IsFactor = NameFactorsDict.get(item.Name,False)
            if item.KO:
                item.Description=KOAnnotationDict.get(item.KO,None)
#                item.KOID=KOIDDict.get(item.KO,None)
                if KOHydrolasesDict.get(item.KO,None):
                    item.IsHydrolase=True
                if item.IsSecondaryMetabolic==False:
                    if KOSecondaryDict.get(item.KO,None):
                        item.IsSecondaryMetabolic=True
                    elif KOMetabolicDict.get(item.KO,None):
                        item.IsPrimaryMetabolic=True
                else:
                    if KOMetabolicDict.get(item.KO,None):
                        item.IsPrimaryMetabolic = True
                item.Pathways=PathwaysDict.get(item.KO,None)

        GeneListSignificant=[]
        for item in GeneList:
            if item.Qval<=0.05 and math.fabs(item.Fc)>=2:
                if item.Fc>=2: item.HowRegulated="Up"
                else: item.HowRegulated="Down"
                GeneListSignificant.append(item)

        up1=down1=0
        for item in GeneListSignificant:
            if item.HowRegulated=="Up": up1+=1
            else: down1+=1

        GeneListProfit=[]
        for item in GeneListSignificant:
            if item.IsHydrolase or item.IsFactor or item.IsSecondaryMetabolic or item.Description:
                GeneListProfit.append(item)

        up2=down2=0
        secup = secdown = hydup = hydown = factup = factdown = keggedup=keggedown = 0
        for item in GeneListProfit:
            if item.HowRegulated=="Up":
                up2+=1
                if item.IsHydrolase: hydup+=1
                if item.IsFactor: factup+=1
                if item.IsSecondaryMetabolic: secup+=1
                if item.KO: keggedup+=1
            else:
                down2+=1
                if item.IsHydrolase: hydown+=1
                if item.IsFactor: factdown+=1
                if item.IsSecondaryMetabolic: secdown+=1
                if item.KO: keggedown+=1


        AllData=pandas.DataFrame({"ID":{item.ID:item.ID for item in GeneList},"Gene name":{item.ID:item.Name for item in GeneList},"Chromosome":{item.ID:item.Chromosome for item in GeneList},"Position":{item.ID:item.ChrPos for item in GeneList},"KO":{item.ID:item.KO for item in GeneList},"KO Description":{item.ID:item.Description for item in GeneList},"KO Pathways":{item.ID:item.Pathways for item in GeneList},"Fold change":{item.ID:item.Fc for item in GeneList},"Q-value":{item.ID:item.Qval for item in GeneList},"How regulated?":{item.ID:item.HowRegulated for item in GeneList},"Secondary metabolism":{item.ID:item.IsSecondaryMetabolic for item in GeneList},"Primary metabolism":{item.ID:item.IsPrimaryMetabolic for item in GeneList},"Hydrolitic enzyme?":{item.ID:item.IsHydrolase for item in GeneList},"Transcription regulator?":{item.ID:item.IsFactor for item in GeneList}})
        SignificantData=pandas.DataFrame({"ID":{item.ID:item.ID for item in GeneListSignificant},"Gene name":{item.ID:item.Name for item in GeneListSignificant},"Chromosome":{item.ID:item.Chromosome for item in GeneListSignificant},"Position":{item.ID:item.ChrPos for item in GeneListSignificant},"KO":{item.ID:item.KO for item in GeneListSignificant},"KO Description":{item.ID:item.Description for item in GeneListSignificant},"KO Pathways":{item.ID:item.Pathways for item in GeneListSignificant},"Fold change":{item.ID:item.Fc for item in GeneListSignificant},"How regulated?":{item.ID:item.HowRegulated for item in GeneListSignificant},"Secondary metabolism":{item.ID:item.IsSecondaryMetabolic for item in GeneListSignificant},"Primary metabolism":{item.ID:item.IsPrimaryMetabolic for item in GeneListSignificant},"Hydrolitic enzyme?":{item.ID:item.IsHydrolase for item in GeneListSignificant},"Transcription regulator?":{item.ID:item.IsFactor for item in GeneListSignificant}})
        ProfitData=pandas.DataFrame({"ID":{item.ID:item.ID for item in GeneListProfit},"Gene name":{item.ID:item.Name for item in GeneListProfit},"Chromosome":{item.ID:item.Chromosome for item in GeneListProfit},"Position":{item.ID:item.ChrPos for item in GeneListProfit},"KO":{item.ID:item.KO for item in GeneListProfit},"KO Description":{item.ID:item.Description for item in GeneListProfit},"KO Pathways":{item.ID:item.Pathways for item in GeneListProfit},"Fold change":{item.ID:item.Fc for item in GeneListProfit},"How regulated?":{item.ID:item.HowRegulated for item in GeneListProfit},"Secondary metabolism":{item.ID:item.IsSecondaryMetabolic for item in GeneListProfit},"Primary metabolism":{item.ID:item.IsPrimaryMetabolic for item in GeneListProfit},"Hydrolitic enzyme?":{item.ID:item.IsHydrolase for item in GeneListProfit},"Transcription regulator?":{item.ID:item.IsFactor for item in GeneListProfit}})
        General=pandas.DataFrame({"N of upregulated features":{"All":up1,"Annotated":up2,"Of them...":"","Have KO":keggedup,"Secondary metabolism":secup,"Hydrolases":hydup,"Factors":factup},"N of downregulated features":{"All":down1,"Annotated":down2,"Of them...":"","Have KO":keggedown,"Secondary metabolism":secdown,"Hydrolases":hydown,"Factors":factdown}})

        def QvalHighlight(qval):
            if qval<=0.05:
                return "background-color:lightgreen"
            else: return None

        def FoldChange(Fold):
            if Fold=="Up":
                return "background-color:#7F6BFF"
            elif Fold=="Down":
                return "background-color:#FF6464"
            else: return None
        def BoolVal(b):
            if b:
                return "background-color:#65FF32"
            else:
                return "background-color:#FF8080"

        def Position(pos):
            if 383<=pos<=536:
                return "background-color:#F5C43D"
            elif 300<=pos<383:
                return "background-color:#FCD568"

        writer=pandas.ExcelWriter(f"{file.removesuffix('csv')}xlsx", engine="openpyxl")
        AllData.style.applymap(QvalHighlight,subset="Q-value").applymap(Position,subset="ID").applymap(FoldChange,subset="How regulated?").applymap(BoolVal, subset=["Secondary metabolism","Primary metabolism","Hydrolitic enzyme?","Transcription regulator?"]).to_excel(writer,sheet_name="All data")
        SignificantData.style.applymap(FoldChange,subset="How regulated?").applymap(Position,subset="ID").applymap(BoolVal, subset=["Secondary metabolism","Primary metabolism","Hydrolitic enzyme?","Transcription regulator?"]).to_excel(writer,sheet_name="Significant change")
        ProfitData.style.applymap(FoldChange,subset="How regulated?").applymap(Position,subset="ID").applymap(BoolVal, subset=["Secondary metabolism","Primary metabolism","Hydrolitic enzyme?","Transcription regulator?"]).to_excel(writer,sheet_name="Significant annotated")
        General.to_excel(writer,sheet_name="General Info")
        writer.close()