data=open("data2.csv")
genelist=[]
class Gene:
    def __init__(self, ID, species, expression_level, sequence):
        self.ID=ID
        self.species=species
        self.expression=expression_level
        self.sequence=sequence
for string in data:
    attributes=string.split(",")
    gen=Gene(attributes[2],attributes[0],attributes[3].strip(),attributes[1])
    genelist.append(gen)
def function0():
    atgc={"a":"t","t":"a","c":"g","g":"c"}
    for gene in genelist:
        comp=[]
        for n in gene.sequence:
            comp.append(atgc[n])
        print(gene.sequence+"\n"+"".join(comp))
def function1():
    for gene in genelist:
        if gene.species == "Drosophila melanogaster" or gene.species == "Drosophila simulans":
            print(gene.ID+" 1")
def function2():
    for gene in genelist:
        if 90<=len(gene.sequence)<=110:
            print(gene.ID+" 2")
def function3():
    for gene in genelist:
        AT=(gene.sequence.count("t")+gene.sequence.count("a"))/len(gene.sequence)
        if AT<0.5 and int(gene.expression)>200:
            print(gene.ID+" 3")
def function4():
    for gene in genelist:
        if str(gene.ID).startswith("h") or str(gene.ID).startswith("k") and gene.species!= "Drosophila melanogaster":
            print(gene.ID+" 4")
def function5():
    for gene in genelist:
        AT =(gene.sequence.count("t")+gene.sequence.count("a")) / len(gene.sequence)
        if AT<0.45: cont=" low "
        elif 0.45<=AT<0.65: cont=" medium "
        else: cont=" high "
        print(str(gene.ID)+cont+str(AT)+" 5")
function0()
function1()
function2()
function3()
function4()
function5()