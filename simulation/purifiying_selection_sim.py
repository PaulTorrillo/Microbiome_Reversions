import copy
import math
import numpy

bottleneckprob = 0 #Probability of bottleneck though unused here
mutationrate=0.003*0.5*0.75*(9/10) #Mutation rate per genome per generation,
#The 0.003 is from assuming 1 in a billion mutations per bp per generation multiplied by 3 for to be in terms of codons
#than multiplied to take into account a genome of 1 million codons. But we end up only saying 500,000 codons so we multiply by 0.5
#not all changes are nonsynonymous so multiply by 0.75 and finally we choose to pretake into acount mutations being nonneutral so
#we multily by 9/10
possiblemutations=500000 #The core genome size in terms of codons
numberofruns=1 #If you wanted to do multiple simulations and take the average
ben_mut_f=1.1 #The fitness of those beneficial mutations
weak_del_mut_prob=1 #Only simulating weak deleterious mutations but its possible to do more
weak_del_f=1-0.000035 #The fitness effect of weak deleterious mutation
del_mut_prob=0 #Could add in deleterious mutations
del_f=0.99 #Fitness effect
strong_del_mut_prob=0 #Could add in strongly deleterious mutations
strong_del_f=0.9 #Fitness effect
neutral_mut_prob=1-strong_del_mut_prob-del_mut_prob-weak_del_mut_prob #Whats left over is neutral mutations,
#we dont want back mutations for this simulation so we took them into account in the mutation rate, not here but we dont have to

#this goes unusued ultimately in our purifying simulations to be more efficent, but if there were multiple classes
#of mutations one would want to fill out essentially a mutation transition table below
mut_probs=[[weak_del_mut_prob,del_mut_prob,strong_del_mut_prob,neutral_mut_prob],
           [weak_del_mut_prob,del_mut_prob,strong_del_mut_prob,neutral_mut_prob],
            [weak_del_mut_prob,del_mut_prob,strong_del_mut_prob,neutral_mut_prob],
            [weak_del_mut_prob,del_mut_prob,strong_del_mut_prob,neutral_mut_prob],
            [weak_del_mut_prob,del_mut_prob,strong_del_mut_prob,neutral_mut_prob],
            [weak_del_mut_prob,del_mut_prob,strong_del_mut_prob,neutral_mut_prob]]
f_coefs=[ben_mut_f, weak_del_f, del_f, strong_del_f] #Fitness advantages in list form
transmission_surv=1000 #Number of survivors given a bottleneck occured
#The class containing the varying "allele"
class allele:
    def __init__(self, size, genotype): #Initialize genetic class
        self.size=size
        genoholder=genotype.split(":") #Translate from string to array
        self.genotype=[int(numeric_string) for numeric_string in genoholder]

    def update(self, population,availablebens,mutantstoadd,norm): #Increase population size and add mutants
        if population>10: #We do not want the population to die out due to random fluctuations so population is deterministic
            #if below size 10
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
    def bottlenecked(self, value): #Make sure genetic classes are carried over correctly during bottleneck
        self.size=value
        self.genotype[4]=self.genotype[4]+self.genotype[0]
        self.genotype[0]=0
        newgeno = str(self.genotype[0])
        for i in range(len(self.genotype) - 1):
            newgeno = newgeno + ':' + str(self.genotype[i + 1])
        return (newgeno, self.size)
if __name__ == '__main__':
    for jm in range(0,4): #This is an outer loop that allows us to change capacity in this code

        #Some lists set aside for storing information if needed
        reversionsz = []
        maxesz = []
        avgfz = []
        maxfz = []
        maxesfz = []
        mutationsz = []
        for zed in range(1): #Another loop in case you want to further manipulations (unusued here)
            #Setting aside some lists to fill up for each run
            reversions = [None]*numberofruns
            maxes=[None]*numberofruns
            avgf=[None]*numberofruns
            maxf = [None] * numberofruns
            maxesf=[None]*numberofruns
            mutations=[None]*numberofruns
            for u in range(numberofruns): #If you want to take the average of multiple runs (unneeded in this scenario)
                print(numberofruns*zed+u)
                numbergens=5000000
                capacity = 10**(6+3*jm)
                bottleneckprob =0 #Probability of bottleneck though unused here
                currentpopsize=capacity #Starting population size
                norm=0
                currentclasses={}
                currentclasses["0:0:0:0:" + str(possiblemutations) + ":0"]=allele(capacity, "0:0:0:0:" + str(possiblemutations) + ":0") #Set up initial class
                availablebens=0 #Provide initial number of beneficial mutations
                reversions[u]=numpy.zeros(numbergens)
                #Make sure to reset the bottleneck if it was set to true
                transmissionbottleneck=False
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
                        mutantstoadd=currentclasses[x].update(currentpopsize,availablebens,mutantstoadd,norm) #With projected
                        #population we can launch our update
                    #Add in mutants
                    for x in mutantstoadd:
                        if x in currentclasses:
                            currentclasses[x].increase_size(mutantstoadd[x])
                        else:
                            currentclasses[x]=allele(mutantstoadd[x], x)

                    #This is all mainly for just collecting some statistics
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
                        #We do need to remove genetic classes that are now empty
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

                    #Printing out some summary statistics and adding to lists keeping data
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
                    mutations[u][k]=mutations[u][k]/populationsize
                    print(mutations[u][k])

                    #There are not bottlenecks in the standard purifying selection so the following is never run in
                    # this simulation but may be useful if you want to modify the code for your ends
                    if transmissionbottleneck:
                        #We recreate the survivors from the bottleneck
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
                        #This was created to explore what would happen if somehow during transmissions the least loaded class always made it
                        #regardless.
                        givensurvivor=('0:0:0:0:500000:0', 1)
                        if givensurvivor[0] in newcurrentclasses:
                            newcurrentclasses[givensurvivor[0]].increase_size(givensurvivor[1])
                        else:
                            newcurrentclasses[givensurvivor[0]] = allele(givensurvivor[1], givensurvivor[0])
                        currentclasses = newcurrentclasses
                        populationsize=transmission_surv
                        transmissionbottleneck=False
                    #Checking for random bottlenecks
                    if numpy.random.uniform(0,1)<bottleneckprob:
                        availablebens=0
                        transmissionbottleneck = True
                    currentpopsize=populationsize
            #Saving off averages for the simulation
            reversionsz.append(numpy.mean(reversions,axis=0))
            maxesz.append(numpy.mean(maxes,axis=0))
            avgfz.append(numpy.mean(avgf,axis=0))
            maxesfz.append(numpy.mean(maxesf,axis=0))
            mutationsz.append(numpy.mean(mutations,axis=0))
        save=[]
        mutationsz=mutationsz[0]
        #Saving of the simulations and only keeping every 100th generation to save on space
        for i in range(len(mutationsz)):
            if i%100==0:
                save.append(mutationsz[i])
        numpy.savetxt('purifying_sim_muts' +str(jm)+'.txt',save)