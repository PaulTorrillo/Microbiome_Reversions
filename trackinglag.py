import copy
import math

import numpy
import matplotlib
import matplotlib.pyplot as plt

capacity=10000000000
bottleneckprob = 1 / 10000
mutationrate=0.003*0.5*0.75 #per genome per generation
possiblemutations=500000
numberofruns=1
#Fill Out Your Mutation Numbers
num_ben_muts=0
ben_f=1.03
del_f=1-0.003
f_coefs=[ben_f, ben_f, del_f]
bottleneck_f_coefs=[1.25, 1, 1, 1]
transmission_surv=10
#The class containing the varying "allele"
class allele:
    def __init__(self, size, genotype): #Initialize species class
        self.size=size
        genoholder=genotype.split(":") #Translate from string to array
        self.genotype=[int(numeric_string) for numeric_string in genoholder]

    def update(self, population,availablebens, availablerevs,mutantstoadd,norm): #Increase population size and add mutants
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
                    #loci likely multiple sites so possibbly easier to find beneificial mutation than single site
                    newbens=numpy.random.binomial(mutantnum, (availablebens-self.genotype[0]) / (possiblemutations/5))
            else:
                newbens=0
            mutantnum = mutantnum - newbens
            newrevs = 0
            if availablerevs > self.genotype[
                1]:  # Beneficial mutations are checked for seperately, generated by binomial
                if mutantnum > 0:
                    # 1/5 penalty as reversions are specfic amino acid subsitutions (from usually 6 possible swtiches but transition transvertion bias then suggests reversion might be easier)
                    newrevs = numpy.random.binomial(mutantnum, (availablerevs - self.genotype[1]) / (5*possiblemutations))
            else:
                newrevs = 0
            mutantnum=mutantnum-newrevs #Remove beneficial mutations from mutant pool
            newsubmuts=numpy.random.multinomial(mutantnum,[0.9,0.1]) #Choose other possible types of mutations
            newmuts=[newbens]+[newrevs]+newsubmuts.tolist()+[0] #Get set of mutant types and corresponding numbers
            if i<3: #Count any non neutral mutation to neutral mutation as a reversion
                holder=newmuts[3]
                newmuts[4]=holder
                newmuts[3]=0
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
                        mkey=":".join([str(k) for k in copygene])
                        if mkey in mutantstoadd:
                            mutantstoadd[mkey]=mutantstoadd[mkey]+totalnewmuts[i][j]
                        else:
                            mutantstoadd[mkey]=totalnewmuts[i][j]
        return mutantstoadd
    def get_nf(self): #Returns population multiplied by selective advantage
        toproduct=1
        for i in range(3):
            toproduct=toproduct*(f_coefs[i]**self.genotype[i]) #Calulate selective advantage
        return toproduct*self.size
    def get_f(self):
        toproduct=1
        for i in range(3):
            toproduct=toproduct*(f_coefs[i]**self.genotype[i]) #Calulate selective advantage
        return toproduct
    def increase_size(self,value): #Increasing size (for adding mutants)
        self.size=self.size+value
    def get_size(self): #Return size
        return self.size
    #Modified uses of how tracking mutations and reversions work
    def get_reversions(self): #Return number of reversions multiplied by population (reversion mass?)
        tosum=self.genotype[1]
        return self.size*tosum
    def get_mutations(self):
        tosum=self.genotype[0]#-self.genotype[1]+self.genotype[2]
        return tosum*self.size
    def get_hitchikers(self):
        tosum = self.genotype[2]
        return tosum * self.size
    def bottlenecked(self, value):
        self.size=value
        #self.genotype[4]=self.genotype[4]+self.genotype[0]
        #self.genotype[0]=0
        newgeno = str(self.genotype[0])
        for i in range(len(self.genotype) - 1):
            newgeno = newgeno + ':' + str(self.genotype[i + 1])
        return (newgeno, self.size)

if __name__ == '__main__':
    reversionsz = []
    hitchz = []
    maxesz = []
    avgfz = []
    maxfz = []
    maxesfz = []
    mutationsz = []
    for zed in range(1):
        reversions = [None]*numberofruns
        hitch = [None] * numberofruns
        maxes=[None]*numberofruns
        avgf=[None]*numberofruns
        maxf = [None] * numberofruns
        maxesf=[None]*numberofruns
        mutations=[None]*numberofruns
        running_sum_forward=0
        running_sum_backward=0
        for u in range(numberofruns):
            print(numberofruns*zed+u)
            numbergens=2000000
            currentpopsize=capacity
            norm=0
            currentclasses={}
            currentclasses["0:0:0:" + str(possiblemutations) + ":0"]=allele(capacity, "0:0:0:" + str(possiblemutations) + ":0") #Set up initial class
            availablebens=0 #Provide initial number of beneficial mutations
            availablerevs=0
            reversions[u]=numpy.zeros(numbergens)
            hitch[u] = numpy.zeros(numbergens)
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
                for x in currentclasses:
                    mutantstoadd=currentclasses[x].update(currentpopsize,availablebens,availablerevs,mutantstoadd,norm) #With projected
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
                    hitch[u][k] = hitch[u][k] + currentclasses[x].get_hitchikers()
                    mutations[u][k]=mutations[u][k]+currentclasses[x].get_mutations()
                for x in marked:
                    del currentclasses[x]
                print(k)
                print(len(currentclasses))
                #print(populationsize)
                print(maxspec)
                print(currentmaxf)
                reversions[u][k]=reversions[u][k]/populationsize
                hitch[u][k] = hitch[u][k] / populationsize
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
                    currentclasses = newcurrentclasses
                    populationsize=transmission_surv
                    transmissionbottleneck=False
                #if k>0 and k%10000==0 and k<40001: #Release new beneficial mutations
                if numpy.random.uniform(0,1)<bottleneckprob:
                    #availablebens=numpy.random.binomial(80,1/16)
                    #cumulativebeneficialmuts=cumulativebeneficialmuts+num_ben_muts
                    transmissionbottleneck = True
                int_array = [int(x) for x in maxspec.split(':')]
                if numpy.random.uniform(0,1)<((1/(840))):
                    int_array = [int(x) for x in maxspec.split(':')]
                    numbertoadd=numpy.random.poisson(1,1)
                    addingrevs=numpy.random.binomial(numbertoadd[0],(availablebens-availablerevs)/55)
                    if addingrevs>(availablebens-availablerevs):
                        addingrevs=availablebens-availablerevs
                    addingbens = numbertoadd[0] - addingrevs
                    availablebens=availablebens+addingbens
                    availablerevs=availablerevs+addingrevs

                print("current lags")
                running_sum_forward=running_sum_forward+(availablebens-mutations[u][k])
                print('lag for forward')
                print(running_sum_forward/(k+1))
                running_sum_backward = running_sum_backward + (availablerevs-reversions[u][k])
                print('lag for reverse')
                print(running_sum_backward / (k + 1))
                currentpopsize=populationsize
        reversionsz.append(numpy.mean(reversions,axis=0))
        hitchz.append(numpy.mean(hitch, axis=0))
        mutationsz.append(numpy.mean(mutations,axis=0))
    matplotlib.pyplot.plot(mutationsz[0])
    matplotlib.pyplot.xlabel('Generations')
    matplotlib.pyplot.ylabel('Average Mutational Load')
    matplotlib.pyplot.xlim([1000,1000000])
    matplotlib.pyplot.xscale('log')
    save1 = []
    save2 = []
    save3=[]
    for i in range(len(mutationsz)):
        if i % 100 == 0:
            save1.append(mutationsz[i])
            save2.append(reversionsz[i])
            save3.append(hitchz[i])
    matplotlib.pyplot.show()