import copy
import math

import numpy
import matplotlib
import matplotlib.pyplot as plt

capacity = 1000000000
bottleneckprob = 0#0.0001
mutationrate=0.003*0.5*0.75*(10/10) #per genome per generation
possiblemutations=500000
numberofruns=1
#Fill Out Your Mutation Numbers
num_ben_muts=0
ben_mut_f=1.1
weak_del_mut_prob=1#22/25
weak_del_f=1-0.000035#1-0.000028
del_mut_prob=0
del_f=0.99999
strong_del_mut_prob=0
strong_del_f=0.99999
neutral_mut_prob=1-strong_del_mut_prob-del_mut_prob-weak_del_mut_prob #0.225
mut_probs=[[weak_del_mut_prob,del_mut_prob,strong_del_mut_prob,neutral_mut_prob], #1.205
           [del_mut_prob, weak_del_mut_prob, neutral_mut_prob,strong_del_mut_prob],
            [del_mut_prob, weak_del_mut_prob, neutral_mut_prob,strong_del_mut_prob],
            [del_mut_prob,weak_del_mut_prob,neutral_mut_prob,strong_del_mut_prob],
[weak_del_mut_prob,del_mut_prob,strong_del_mut_prob,neutral_mut_prob],
[weak_del_mut_prob,del_mut_prob,strong_del_mut_prob,neutral_mut_prob],]
f_coefs=[ben_mut_f, weak_del_f, del_f, strong_del_f]
bottleneck_f_coefs=[1.25, 1, 1, 1]
transmission_surv=1000
#The class containing the varying "allele"
class allele:
    def __init__(self, size, genotype): #Initialize species class
        self.size=size
        genoholder=genotype.split(":") #Translate from string to array
        self.genotype=[int(numeric_string) for numeric_string in genoholder]

    def update(self, population,availablebens,mutantstoadd,norm): #Increase population size and add mutants
        if population>10:
            self.size=numpy.random.poisson((self.get_nf()/norm)*population) #Number of given offspring. Determined by
        else:
            self.size=int((self.get_nf()/norm)*population)
        #Wright Fisher model and realized as Poisson distribution. Only difference from standard Wright-Fisher Model being
        #this distribution is not conditioned on total population being capacity, though the expected total population
        #is capacity
        if self.size>0:
            numberofmutants = numpy.random.binomial(self.size, mutationrate) #Here mutations are determined by bionomial
        else:
            numberofmutants=0
        #distribution
        mutantclasses=numpy.random.multinomial(numberofmutants, numpy.divide(self.genotype, possiblemutations)) #Specific types
        #of those mutants are then determined by the multinomial distribution
        totalnewmuts=[]
        for i in range(len(mutantclasses)):
            mutantnum=mutantclasses[i]
            newbens=0
            if availablebens>self.genotype[0]: #Beneficial mutations are checked for seperately, generated by binomial
                if mutantnum>0:
                    newbens=numpy.random.binomial(mutantnum, (availablebens-self.genotype[0]) / (possiblemutations))
            else:
                newbens=0
            mutantnum=mutantnum-newbens #Remove beneficial mutations from mutant pool
            newsubmuts=numpy.random.multinomial(mutantnum,mut_probs[0]) #Choose other possible types of mutations
            newmuts=[newbens]+newsubmuts.tolist()+[0] #Get set of mutant types and corresponding numbers
            if i<4: #Count any non neutral mutation to neutral mutation as a reversion
                holder=newmuts[4]
                newmuts[5]=holder
                newmuts[4]=0
            totalnewmuts.append(newmuts) #Count total number of new mutants
        self.size=self.size-numberofmutants #Remove mutants from population
        #Need to update appropriate classes
        for i in range(len(totalnewmuts)):
            for j in range(len(totalnewmuts[i])):
                if i!=j:
                    if totalnewmuts[i][j]>0:
                        added=self.genotype[j]+1
                        lost=self.genotype[i]-1
                        copygene=copy.deepcopy(self.genotype)
                        copygene[j]=added
                        copygene[i]=lost
                        if ((i == 3 and (j == 2 or j == 1)) or (i == 2 and j == 1)):
                            copygene[1]=copygene[1]+1
                            copygene[4]=copygene[4]-1
                        mkey=":".join([str(k) for k in copygene])
                        if mkey in mutantstoadd:
                            mutantstoadd[mkey]=mutantstoadd[mkey]+totalnewmuts[i][j]
                        else:
                            mutantstoadd[mkey]=totalnewmuts[i][j]
        return mutantstoadd
    def get_nf(self): #Returns population multiplied by selective advantage
        toproduct=1
        for i in range(4):
            toproduct=toproduct*(f_coefs[i]**self.genotype[i]) #Calulate selective advantage
        return toproduct*self.size
    def get_f(self):
        toproduct=1
        for i in range(4):
            toproduct=toproduct*(f_coefs[i]**self.genotype[i]) #Calulate selective advantage
        return toproduct
    def increase_size(self,value): #Increasing size (for adding mutants)
        self.size=self.size+value
    def get_size(self): #Return size
        return self.size
    def get_reversions(self): #Return number of reversions multiplied by population (reversion mass?)
        return self.size*self.genotype[5]
    def get_mutations(self):
        tosum=0
        for i in range(1,4):
            tosum=tosum+self.genotype[i] #Calulate selective advantage
        return tosum*self.size
    def get_f_bottleneck(self):
        toproduct = 1
        for i in range(4):
            toproduct = toproduct * (bottleneck_f_coefs[i] ** self.genotype[i])  # Caclulate selective advantage
        return toproduct
    def bottlenecked(self, value):
        self.size=value
        self.genotype[4]=self.genotype[4]+self.genotype[0]
        self.genotype[0]=0
        newgeno = str(self.genotype[0])
        for i in range(len(self.genotype) - 1):
            newgeno = newgeno + ':' + str(self.genotype[i + 1])
        return (newgeno, self.size)
if __name__ == '__main__':
    for jm in range(2,3):
        reversionsz = []
        maxesz = []
        avgfz = []
        maxfz = []
        maxesfz = []
        mutationsz = []
        listopain=[]
        timer=0
        for zed in range(1):
            reversions = [None]*numberofruns
            maxes=[None]*numberofruns
            avgf=[None]*numberofruns
            maxf = [None] * numberofruns
            maxesf=[None]*numberofruns
            mutations=[None]*numberofruns
            for u in range(numberofruns):
                print(numberofruns*zed+u)
                numbergens=5000000
                capacity = 1000000000#10**(6+3*jm)
                bottleneckprob =0#0.0001
                currentpopsize=capacity
                norm=0
                currentclasses={}
                currentclasses["0:0:0:0:" + str(possiblemutations) + ":0"]=allele(capacity, "0:0:0:0:" + str(possiblemutations) + ":0") #Set up initial class
                availablebens=0 #Provide initial number of beneficial mutations
                reversions[u]=numpy.zeros(numbergens)
                transmissionbottleneck=False
                selectivebottleneck=False
                mutations[u]=numpy.zeros(numbergens)
                maxes[u]=[]
                maxesf[u] = []
                avgf[u]=numpy.zeros(numbergens)
                cumulativebeneficialmuts=0
                for k in range(numbergens):
                    currentpopsize=min(math.ceil(currentpopsize+(currentpopsize*(1-(currentpopsize/capacity)))),capacity) #Use logistic growth
                    #model with sharp cut off at capacity
                    norm = 0 #Need to check how one allele will be compared to all other alleles
                    for x in currentclasses:
                        norm=norm+currentclasses[x].get_nf()
                    mutantstoadd={}
                    transmissionmuts=1
                    for x in currentclasses:
                        modifier=1
                        if selectivebottleneck:
                            modifier=(1/1.99)+((1.99*currentclasses[x].genotype[0])/19.9)
                        mutantstoadd=currentclasses[x].update(modifier*currentpopsize,availablebens,mutantstoadd,norm) #With projected
                        #population we can launch our update
                    #Add in mutants
                    for x in mutantstoadd:
                        if x in currentclasses:
                            currentclasses[x].increase_size(mutantstoadd[x])
                        else:
                            currentclasses[x]=allele(mutantstoadd[x], x)
                    marked=[]
                    currentmax=0
                    currentmaxf=0
                    currentmaxfpop=0
                    maxspec=None
                    populationsize=0
                    spec_size_list=[]
                    spec_name_list=[]
                    for x in currentclasses:
                        spec_size=currentclasses[x].get_size()
                        spec_f=currentclasses[x].get_f()
                        populationsize=populationsize+spec_size
                        if spec_size>currentmax:
                            currentmax=spec_size
                            maxspec=x
                        if spec_f>=currentmaxf:
                            if spec_size>currentmaxfpop or spec_f!=currentmaxf:
                                currentmaxf=spec_f
                                currentmaxfpop=spec_size
                        if spec_size<1:
                            marked.append(x)
                        else:
                            spec_size_list.append(spec_size)
                            spec_name_list.append(x)
                        avgf[u][k]=avgf[u][k]+currentclasses[x].get_nf()
                        reversions[u][k]=reversions[u][k]+currentclasses[x].get_reversions()
                        mutations[u][k]=mutations[u][k]+currentclasses[x].get_mutations()
                    for x in marked:
                        del currentclasses[x]
                    print(k)
                    print(len(currentclasses))
                    #print(populationsize)
                    print(maxspec)
                    print(currentmaxf)
                    print(currentmaxfpop/populationsize)
                    reversions[u][k]=reversions[u][k]/populationsize
                    print(reversions[u][k])
                    avgf[u][k]=avgf[u][k]/populationsize
                    maxes[u].append(currentmax)
                    maxesf[u].append(currentmaxfpop)
                    mutations[u][k]=mutations[u][k]/populationsize#+cumulativebeneficialmuts
                    print(mutations[u][k])
                    if transmissionbottleneck:
                        survivors=numpy.random.multinomial(transmission_surv,numpy.divide(spec_size_list,populationsize))
                        recreate=[]
                        newcurrentclasses={}
                        for z in range(len(survivors)):
                            recreate.append(currentclasses[spec_name_list[z]].bottlenecked(survivors[z]))
                        for g in recreate:
                            if g[0] in newcurrentclasses:
                                newcurrentclasses[g[0]].increase_size(g[1])
                            else:
                                newcurrentclasses[g[0]] = allele(g[1], g[0])
                        givensurvivor=('0:0:0:0:400000:0', 1)
                        if givensurvivor[0] in newcurrentclasses:
                            newcurrentclasses[givensurvivor[0]].increase_size(givensurvivor[1])
                        else:
                            newcurrentclasses[givensurvivor[0]] = allele(givensurvivor[1], givensurvivor[0])
                        currentclasses = newcurrentclasses
                        populationsize=transmission_surv
                        transmissionbottleneck=False
                    #if k>0 and k%10000==0:# and k<40001: #Release new beneficial mutations

                    if numpy.random.uniform(0,1)<bottleneckprob:
                        #availablebens=numpy.random.binomial(80,1/16)
                        availablebens=0
                        listopain.append((1-currentmaxfpop/populationsize)**1000)
                        #cumulativebeneficialmuts=cumulativebeneficialmuts+num_ben_muts
                        transmissionbottleneck = True
                        selectivebottleneck=False
                    currentpopsize=populationsize
            reversionsz.append(numpy.mean(reversions,axis=0))
            maxesz.append(numpy.mean(maxes,axis=0))
            avgfz.append(numpy.mean(avgf,axis=0))
          #  maxfz.append(numpy.mean(maxf,axis=0))
            maxesfz.append(numpy.mean(maxesf,axis=0))
            mutationsz.append(numpy.mean(mutations,axis=0))
        #matplotlib.pyplot.plot(mutationsz[0])
        #matplotlib.pyplot.xlabel('Metapopulation Turnovers')
        #matplotlib.pyplot.ylabel('Average Mutational Load')
        #matplotlib.pyplot.show()
        save=[]
        mutationsz = numpy.mean(mutationsz, axis=0)
        reversionsz = numpy.mean(reversionsz, axis=0)
        hitchz = numpy.mean(hitchz, axis=0)
        for i in range(len(mutationsz)):
            if i%100==0:
                save.append(mutationsz[i])
        numpy.savetxt('highermutation.txt',save)
