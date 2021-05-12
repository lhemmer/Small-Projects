#load library
library(ggplot2)
library(ggpubr)
library(reshape2)


#initialize population and some conditions
avgNoOffspring <- 2

makePop <- function(numMales, numFemales, numDriver, numMod) {
  startPopMales <- rep(3,numMales)
  startPopMales[sample(1:numMales,numDriver,replace=F)] <- 2
  startPopFemales <- rep(3,numFemales)
  startPopFemales[sample(1:numFemales,numDriver,replace=F)] <- 2
  
  modPopMales <- rep(3,numMales)
  modPopMales[sample(1:numMales,round(numMod*2/3))] <- 2
  modPopMales[sample(1:numMales,round(numMod/3))] <- 1
  modPopFemales <- rep(3,numFemales)
  modPopFemales[sample(1:numFemales,round(numMod*2/3))] <- 2
  modPopFemales[sample(1:numFemales,round(numMod/3))] <- 1
  
  startPopMales <- data.frame(startPopMales,modPopMales,"m")
  colnames(startPopMales) <- c("genotype","mod.genotype","sex")
  startPopFemales <- data.frame(startPopFemales,modPopFemales,"f")
  colnames(startPopFemales) <- c("genotype","mod.genotype","sex")
  startingPop <- rbind.data.frame(startPopMales,startPopFemales)
  
  return(startingPop)
}

## Drive functions

#start pop dataframe, probabilty of getting A when het, probably of male/female, fitness of homozygote, 
#fitness of het, number of generations, turn modifier on/off, trial index
simulateAutosome <- function(startPopulation, probA, probS, numGenerations, AAfitness, Aafitness, modEffect) {
  df.allele <- NULL
  tmpPop <- startPop
  for (gen in 1:numGenerations) {
    tmpMales <- subset(tmpPop, tmpPop$sex == "m")
    tmpFemales <- subset(tmpPop, tmpPop$sex == "f")
    if (sum(tmpPop$sex=="m") == 0 | sum(tmpPop$sex=="f") == 0) {break}
    avgMates <- sum(tmpPop$sex=="m") / sum(tmpPop$sex=="f")
    offSpring <- getOffspring(tmpFemales, tmpMales, avgMates, avgNoOffspring, probA, probS, modEffect)
    if (AAfitness < 1 | Aafitness < 1) {
      offSpring <- offspringSelection(offSpring, AAfitness, Aafitness)
    }
    nextGen <- rbind.data.frame(tmpPop,offSpring)
    if (nrow(nextGen) > nrow(startPopulation)) {
      x <- sample(1:nrow(nextGen),size=(nrow(startPopulation)),replace=F)
      tmpPop <- nextGen[x,]
    } else { tmpPop <- nextGen }
    numA <- calcAllelesInGeneration(tmpPop$genotype)
    numS <- calcSexInGeneration(tmpPop$sex)
    numM <- calcAllelesInGeneration(tmpPop$mod.genotype)
    df.allele <- rbind(df.allele, c(gen, numA, numS, numM))
  } 
  colnames(df.allele) <- c("Generation","Driver","NonDriver","Male","Female","Mod","NonMod")
  df.allele <- as.data.frame(df.allele)
  return(df.allele)
}

## remove offspring based on fitness
offspringSelection <- function(offSpring,wAA,wAa) {
  AAoffspring <- offSpring[sample((which(offSpring$genotype==1)), round(wAA*length(which(offSpring$genotype==1)))),]
  Aaoffspring <- offSpring[sample((which(offSpring$genotype==2)), round(wAa*length(which(offSpring$genotype==2)))),]
  aaoffspring <- offSpring[offSpring$genotype==3,]
  selectedOffspring <- rbind(AAoffspring,Aaoffspring,aaoffspring)
  return(selectedOffspring)
}

#Bookkeeping 
calcAllelesInGeneration <- function(x) {
  AA = sum(x==1)
  AB = sum(x==2)
  BB = sum(x==3)
  #print(unlist(x))
  #print(c(AA, AB, BB))
  num.A <- (2 * AA)   + AB
  num.B <- (2 * BB)   + AB
  return(c(num.A, num.B))
}

#Bookkeeping 
calcSexInGeneration <- function(x) {
  num.M <- sum(x=="m")
  num.F <- sum(x=="f")
  return(c(num.M, num.F))
}

#Given a parent individual, get one of their Alleles in the gamete
getGamete <- function(indiv, prob.big.A) {
  if (indiv == 1) return(1) #AA
  if (indiv == 3) return(0) #aa
  #if Parent is Aa, the gamete is one binomial trial with prob.big.A
  if (indiv == 2) return(rbinom(1, size=1, prob.big.A)) #Aa
}

#Given a parent individual, get one of their Alleles in the gamete
getGameteMom <- function(indiv, prob.big.A=0.5) {
  if (indiv == 1) return(1) #AA
  if (indiv == 3) return(0) #aa
  #if Parent is Aa, the gamete is one binomial trial with prob.big.A
  if (indiv == 2) return(rbinom(1, size=1, prob.big.A)) #Aa
}

#Given a parent individual, get one of their Alleles in the gamete
getGameteDad <- function(indiv, prob.big.A=0.5) {
  if (indiv == 1) return(1) #AA
  if (indiv == 3) return(0) #aa
  #if Parent is Aa, the gamete is one binomial trial with prob.big.A
  if (indiv == 2) return(rbinom(1, size=1, prob.big.A)) #Aa
}

#Two parental Gametes combine to form a zygote
combineGametes <- function(mg,dg) {
  if ((mg == 1) && (dg == 1)) return(1) #AA
  else if ((mg == 0) && (dg == 0)) return(3) #aa
  else return(2) #Aa
}

#Given two individuals, get an offspring for next generation
getOffspringMod <- function(mom, dad,p) {
  momGam <- getGameteMom(mom)
  dadGam <- getGameteDad(dad,p)
  return(combineGametes(momGam, dadGam))  
}

#given two individuals and modifier, get offspring for next generation
getOffspringGeno <- function(mom, dad, dadMod, p, ModEffect) {
  momGam <- getGameteMom(mom)
  if ( (dadMod == 1 | dadMod == 2) && ModEffect==T) {
    dadGam <- getGameteDad(dad,0.5)
  } else { dadGam <- getGameteDad(dad,p) }
  return(combineGametes(momGam, dadGam))  
}

# Get offspring Sex
getOffspringSex <- function(prob) {
  x <- rbinom(1,size=1,prob)
  if (x==1) return("f")
  else return("m")
}

## genotype of offspring
getOffspring <- function(female.df, male.df, avgMates, avgNoOffspring, probA, probS, ModEffect) {
  tmpGeno <- NULL
  tmpMod <- NULL
  tmpSex <- NULL
  for (i in 1:nrow(female.df)) {
    #number of matings
    noMatings <- rpois(1,avgMates)
    if (noMatings > 0) {
      for (j in 1:length(noMatings)) {
        noOffspring <- rpois(1,avgNoOffspring)
        maleGeno <- male.df[sample(1:nrow(male.df),size=1),]
        if (noOffspring > 0) {
          for (k in 1:length(noOffspring)) {
            tmpMod <- c(tmpMod, getOffspringMod(female.df$mod.genotype[i],maleGeno$mod.genotype,0.5))
            tmpGeno <- c(tmpGeno, getOffspringGeno(female.df$genotype[i],maleGeno$genotype,maleGeno$mod.genotype,probA, ModEffect))
            tmpSex <- c(tmpSex, getOffspringSex(probS)) 
          } } } } }
  if (length(tmpGeno) == 0 | length(tmpSex) == 0) { return(NULL) }
  else {
    tmpDF <- cbind.data.frame(tmpGeno,tmpMod,tmpSex)
    colnames(tmpDF) <- c("genotype","mod.genotype","sex")
    return(tmpDF) } 
}


###plotting function
plotAllelesWithTime <- function(df) {  
  colorRange<-colorRampPalette(c(rgb(0,0,1), rgb(1,0.7,0) ))
  p <- ggplot(df, aes(x= Generation, y= value/(4*500), group=variable, color=factor(variable))) + geom_line()
  p <- p + scale_colour_manual(values = colorRange(2), name="Alleles")  
  p <- p + labs(title = "Allele Frequencies Across Generations")
  p <- p + ylab("Frequency of Allele in the Population") + ylim(0, 1)
  return(p)
}






###################
###################

# our initial starting population with number of drive alleles and modifier alleles

startPop <- makePop(numMales = 500, numFemales = 500, numDriver = 10, numMod = 100)

# probA is strength of drive (0.5 fair, > selfish)
# probS is probabilty of male / female (0.5)
# numGenerations is number of Generations in our simulation
# AA fitness is fitness of individual homozygous for the driver (0 dead, 1 no diff)
# AA fitness is fitness of individual heterogyzous for the driver (0 dead, 1 no diff)
# modEffect is whether a modifier makes the driver fair (T) or does nothing (F)

df <- simulateAutosome(startPop, probA=0.95, probS=0.5, numGenerations = 1000, 
                       AAfitness=0, Aafitness=1,modEffect=T)

# illustrate our simulation
df1 <- melt(data = df, id = c("Generation"), measure=c("Driver", "Mod"))
plotAllelesWithTime(df1)










