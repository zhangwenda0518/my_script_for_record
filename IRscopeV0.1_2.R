###Under the license of the University of Helsinki
#Author: Ali Amiryousefi
#email: ali.amiryousefi@helisnki.fi


#loading required packages

library(seqinr)
library(ape)
library(shape)
library(diagram)
library(reutils)
library(snow)
library(snowfall)
library(knitr)
library(shiny)
library(jpeg)

GC.content<- function(fas){
  return((sum(fas=="g")+sum(fas=="c"))/length(fas))
}
FasExtract<- function(gb){
  #' Genome extracter
  #' 
  #' Extracting the fasta format chloroplast genome from the GeneBank File
  #' @param file Name of the GeneBank file of a sequence on the working reposotory
  #' @return The genome of the sequence that is deposited at the end of the GeneBank file in fasta format
  #' @export
  fasta<-gb[(grep("ORIGIN", gb)+1):length(gb)]
  while(fasta[length(fasta)]=="") {
    fasta<- fasta[1:length(fasta)-1]
  }
  while(fasta[length(fasta)]=="//") {
    fasta<- fasta[1:length(fasta)-1]
  }
  fas<-""
  for (i in 1:length(fasta)){
    sort.let<- sort(unique(c(grep("c", strsplit(fasta[i], " ")[[1]]),grep("a", strsplit(fasta[i], " ")[[1]]), grep("t", strsplit(fasta[i], " ")[[1]]), grep("g", strsplit(fasta[i], " ")[[1]]))))
    try(if(length(sort.let)==!6) stop("Check the gb file; the columns of ORIGIN should be 6"))
    fasta[i]<- paste(strsplit(fasta[i], " ")[[1]][sort.let[1]],strsplit(fasta[i], " ")[[1]][sort.let[2]],strsplit(fasta[i], " ")[[1]][sort.let[3]],strsplit(fasta[i], " ")[[1]][sort.let[4]],strsplit(fasta[i], " ")[[1]][sort.let[5]],strsplit(fasta[i], " ")[[1]][sort.let[6]], sep="")
    fas<-paste(fas, fasta[i], sep="")
  }
  #fasta[length(fasta)]<- gsub("NA", "", fasta[length(fasta)])
  fas<-gsub("NA", "", fas)
  strsplit(fas, "")[[1]]
} 
read.gb<- function(file){
  #' R version of the GB file
  #'
  #' Reading the GB file in the working directory producing the RGB format file needed for the functions of the package
  #'@param file The name of the file in the directory in GeneBank format
  #'@return An object of class Rgb
  #'@export
  return(readLines(file))
}
fetch.gb<- function(GI, read=TRUE){#the GI can be also accession number
  p<- efetch(GI, "nucleotide", "gb")
  write(content(p, "text"), file = paste(GI, ".gb", sep=""))
  if (read){
    read.gb(paste(GI, ".gb", sep=""))
  }
}
rdnFixer<- function(gb){
  seq<- FasExtract(gb)
  seq[which(seq=="u")]<-sample(c("t"), length(which(seq=="u")), TRUE)
  seq[which(seq=="r")]<-sample(c("a", "g"), length(which(seq=="r")), TRUE)
  seq[which(seq=="y")]<-sample(c("c", "t"), length(which(seq=="y")), TRUE)
  seq[which(seq=="s")]<-sample(c("c", "g"), length(which(seq=="s")), TRUE)
  seq[which(seq=="w")]<-sample(c("a", "t"), length(which(seq=="w")), TRUE)
  seq[which(seq=="k")]<-sample(c("g", "t"), length(which(seq=="k")), TRUE)
  seq[which(seq=="m")]<-sample(c("c", "a"), length(which(seq=="m")), TRUE)
  seq[which(seq=="b")]<-sample(c("c", "g", "t"), length(which(seq=="b")), TRUE)
  seq[which(seq=="d")]<-sample(c("a", "g", "t"), length(which(seq=="d")), TRUE)
  seq[which(seq=="h")]<-sample(c("c", "a", "t"), length(which(seq=="h")), TRUE)
  seq[which(seq=="v")]<-sample(c("c", "a", "g"), length(which(seq=="v")), TRUE)
  seq[which(seq=="n")]<-sample(c("c", "g", "t", "a"), length(which(seq=="n")), TRUE)
  seq[which(seq=="-")]<-sample(c("c", "g", "t", "a"), length(which(seq=="-")), TRUE)
  return(seq)
}
sp.name<- function(gb){
  paste(strsplit(gb[2], " ")[[1]][3], strsplit(gb[2], " ")[[1]][4], sep=" ")
}
SSCrev<- function(genecord, SSCs, SSCe){
  g<- genecord
  a<- SSCs
  b<- SSCe
  g[, 2]<- as.numeric(g[,2])
  g[, 3]<- as.numeric(g[,3])
  g[, 4]<- as.numeric(g[,4])
  g[, 5]<- as.numeric(g[,5])
  for (i in 1:length(g[,1])){
    
    #if the two values of the second and third column are in the SSC, then reverse them
    if(max(g[i, 2], g[i, 3]) <= b  && min(g[i, 2], g[i, 3]) >= a){
      g[i, 2]<- paste(b-(as.numeric(g[i, 2])-a))
      g[i, 3]<- b-(as.numeric(g[i, 3])-a)
    }
    
    #If the gene is starting from SSC then it should also be fixed
    
    #on the positive strand
    if(g[i, 2] < a & g[i, 3] > a){
      g[i, 2]<- b+(a-as.numeric(g[i,2]))
      g[i, 3]<- b-(as.numeric(g[i,3])-a)
    }
    #reverseve of that
    else if(g[i, 2] > b & g[i, 3] < b){
      g[i, 2]<- a-(as.numeric(g[i,2])-b)
      g[i, 3]<- a+(b-as.numeric(g[i,3]))
    }
    
    #on the negative strand
    else if(g[i, 2] < b & g[i, 3] > b){
      g[i, 2]<- a+(b-as.numeric(g[i,2]))
      g[i, 3]<- a-(as.numeric(g[i,3])-b)
    }
    #reverseve of that
    else if(g[i, 2] > a & g[i, 3] < a){
      g[i, 2]<- b-(as.numeric(g[i,2])-a)
      g[i, 3]<- b+(a-as.numeric(g[i,3]))
    }
    
    
    
    ###for the 4th and 5th columns
    
    #if the two values of the second and third column are in the SSC, then reverse them
    if(max(g[i, 4], g[i, 5]) <= b  && min(g[i, 4], g[i, 5]) >= a){
      g[i, 4]<- paste(b-(as.numeric(g[i, 4])-a))
      g[i, 5]<- b-(as.numeric(g[i, 5])-a)
    }
    
    #If the gene is starting from SSC then it should also be fixed
    
    #on the positive strand
    if(g[i, 4] < a & g[i, 5] > a){
      g[i, 4]<- b+(a-as.numeric(g[i,4]))
      g[i, 5]<- b-(as.numeric(g[i,5])-a)
    }
    #reverseve of that
    else if(g[i, 4] > b & g[i, 5] < b){
      g[i, 4]<- a-(as.numeric(g[i,4])-b)
      g[i, 5]<- a+(b-as.numeric(g[i,5]))
    }
    
    #on the negative strand
    else if(g[i, 4] < b & g[i, 5] > b){
      g[i, 4]<- a+(b-as.numeric(g[i,4]))
      g[i, 4]<- a-(as.numeric(g[i,5])-b)
    }
    #reverseve of that
    else if(g[i, 4] > a & g[i, 5] < a){
      g[i, 4]<- b-(as.numeric(g[i,4])-a)
      g[i, 5]<- b+(a-as.numeric(g[i,5]))
    }
    
    
    
    
    
    
  }
  g[, 2]<- as.character(g[,2])
  g[, 3]<- as.character(g[,3])
  g[, 4]<- as.character(g[,4])
  g[, 5]<- as.character(g[,5])
  return(g)
}

IRinfo<- function(genome, parallel=TRUE){
  
  #' IR information
  #' 
  #' Detecting the starts and the length of the IR regions on the chloroplast genome
  #' 
  #' @param genome The plastid genome of a species as a simple vector of nucleotides 
  #' @return a vector of four elements as the start of the first and second IR region, their length and the total lenght of the genome, respectively
  #' @export 
  
  ###Preliminary functions 
  
  cirtick<<- function(tick, vector){#circular rotative function
    if(tick > length(vector)-1 || tick < 1){
      return(vector)
    }
    else {
      return(c(vector[(tick+1):length(vector)], vector[1:tick]))
    }
  }
  
  genome.comp.rev<<- function(genome){#reverse complement function
    gcr<-genome[length(genome):1]
    gcr<-gsub("a", "T", gcr)
    gcr<-gsub("t", "A", gcr)
    gcr<-gsub("g", "C", gcr)
    gcr<-gsub("c", "G", gcr)
    return(tolower(gcr))
  }
  
  #Checked
  phase.detector<- function(genome){#detecting the phase difference of the two inverted genomes
    shifter=84000 #this shifter is faster mode, to be sure set the value to 80000,
    gcr<- genome.comp.rev(genome)
    genome<-cirtick(shifter, genome) 
    l<-length(genome)
    track<- numeric(l+1)
    track[l+1]<- round(l/4)
    for (i in 1:l){
      track[i]<- sum(cirtick(i, genome)==gcr)
      if ((track[i] - track[l+1])/l > 0.1) {###stable version with 0.07 but changing it for the Guizotia abyssinica, for the parallel p.d fucntion as well
        break
      }
    }
    a<<- which(track==max(track))+shifter
    if ( a > l){
      return(a - l)
    }
    else {
      return(a)
    }
  }
  #the parallel version of the phase.detector function. Set with default 4
  p.d<<- function(genome, nCPU=2){#tested in the GlobEnv, it might fail when embedded in the bigger function
    genome<<-genome
    fun<<- function(shifter){#the function to be passed to the slaves for the parallel computing, the out put is with max 10 second either NA or the phase.detector
      gcr<- genome.comp.rev(genome)
      cir.genome<<-cirtick(shifter, genome)
      l<-length(genome)
      track<- numeric(l+1)
      track[l+1]<- round(l/4)
      s.time<- Sys.time()
      no.value<- FALSE
      for (i in 1:l){
        track[i]<- sum(cirtick(i, cir.genome)==gcr)
        i.time<- Sys.time()
        if ((track[i] - track[l+1])/l > 0.1) {
          break
        }
        if (i.time - s.time > 11){
          no.value<- TRUE
          break
        }
      }
      if (no.value) {
        a<<- NA
        return(a)
      }
      else {
        a<<- which(track==max(track))+shifter
        if ( a > l){
          return(a - l)
        }
        else {
          return(a)
        } 
      }
    }
    ini.forw<<- 84000
    ini.back<<- ini.forw
    mm<<- rep(NA, nCPU)
    sfStop()
    sfInit(parallel=TRUE, cpus=nCPU)
    while(sum(is.na(mm))==length(mm)){
      sfExport("genome", "cirtick", "genome.comp.rev", "fun", "mm", "ini.forw", "ini.back")
      mm<-unlist(sfLapply(c(seq(ini.forw, ini.forw+nCPU/2*1000-1, 1000), seq(ini.back, ini.back- nCPU/2*1000, -1000)[-1]),  fun))
      ini.forw<<- ini.forw+nCPU/2*1000
      ini.back<<-ini.back-nCPU/2*1000
    }
    sfStop()
    return(unique(mm)[which((is.na(unique(mm))==FALSE))])
  }
  
  #Checked
  True.sequence.finder<- function(genome, phase.difference){#finding the cordinate of the IR region
    #phase.difference<-phase.detector(genome)
    true.search<-cirtick(phase.difference, genome)==genome.comp.rev(genome)
    true.arm<- round(length(genome)/100)
    for (i in 1:length(genome)){
      if (sum(true.search[i:(true.arm+i-1)])==true.arm) {
        return(i); break
      }
    }
  }
  
  #Checked
  IR1start<-function(phase.difference,True.sequence.finder){
    return(phase.difference+True.sequence.finder)
  }
  
  #Checked
  IR.length<- function(IR1start, True.sequence.finder, genome){
    s<-IR1start
    t<-True.sequence.finder
    r<- genome.comp.rev(genome)
    T<-cirtick(s, genome)==cirtick(t, r)
    Tl<- list()
    for (i in 1:50){
      Tl[[i]]<- cirtick((s + i), genome)==cirtick(t, r)
    }
    for (i in 1:50){
      Tl[[i+50]]<- cirtick(s , genome)==cirtick((t+i), r)
    }
    count<-1
    while(T[count]){
      count<- count+1
    }
    ###jump from INDEL
    for (i in 1:50){
      if(sum(Tl[[i]][count:(count+9)])==10){
        while (Tl[[i]][count]){
          count<- (count + i)
        }
      }
      if(sum(Tl[[i]][(count+1):(count+10)])==10){
        while (Tl[[i]][(count+1)]){
          count<- ((count+1) + i)
        }
      }
      if(sum(Tl[[i]][(count+2):(count+11)])==10){
        while (Tl[[i]][(count+2)]){
          count<- ((count+2) + i)
        }
      }
      if(sum(Tl[[i]][(count+3):(count+12)])==10){
        while (Tl[[i]][(count+3)]){
          count<- ((count+3) + i)
        }
      }
      if(sum(Tl[[i]][(count+4):(count+13)])==10){
        while (Tl[[i]][(count+4)]){
          count<- ((count+4) + i)
        }
      }
      if(sum(Tl[[i]][(count+5):(count+14)])==10){
        while (Tl[[i]][(count+5)]){
          count<- ((count+5) + i)
        }
      }
    }
    for (i in 1:50){
      if(sum(Tl[[i+50]][count:(count+9)])==10){
        while (Tl[[i+50]][count]){
          count<- (count + i)
        }
      }
      if(sum(Tl[[i+50]][(count+1):(count+10)])==10){
        while (Tl[[i+50]][(count+1)]){
          count<- ((count+1) + i)
        }
      }
      if(sum(Tl[[i+50]][(count+2):(count+11)])==10){
        while (Tl[[i+50]][(count+2)]){
          count<- ((count+2) + i)
        }
      }
      if(sum(Tl[[i+50]][(count+3):(count+12)])==10){
        while (Tl[[i+50]][(count+3)]){
          count<- ((count+3) + i)
        }
      }
      if(sum(Tl[[i+50]][(count+4):(count+13)])==10){
        while (Tl[[i]][(count+4)]){
          count<- ((count+4) + i)
        }
      }
      if(sum(Tl[[i+50]][(count+5):(count+14)])==10){
        while (Tl[[i]][(count+5)]){
          count<- ((count+5) + i)
        }
      }
    }
    ##end of jump indel
    if (sum(T[(count+10):(count+25)]==rep(TRUE, 16))==16){
      count<- count+10
      while(T[count]){
        count<- count+1
      } 
    }
    if (sum(T[(count+10):(count+25)]==rep(TRUE, 16))==16){
      count<- count+10
      while(T[count]){
        count<- count+1
      } 
    }
    ###jump from INDEL
    for (i in 1:50){
      if(sum(Tl[[i]][count:(count+9)])==10){
        while (Tl[[i]][count]){
          count<- (count + i)
        }
      }
      if(sum(Tl[[i]][(count+1):(count+10)])==10){
        while (Tl[[i]][(count+1)]){
          count<- ((count+1) + i)
        }
      }
      if(sum(Tl[[i]][(count+2):(count+11)])==10){
        while (Tl[[i]][(count+2)]){
          count<- ((count+2) + i)
        }
      }
      if(sum(Tl[[i]][(count+3):(count+12)])==10){
        while (Tl[[i]][(count+3)]){
          count<- ((count+3) + i)
        }
      }
      if(sum(Tl[[i]][(count+4):(count+13)])==10){
        while (Tl[[i]][(count+4)]){
          count<- ((count+4) + i)
        }
      }
      if(sum(Tl[[i]][(count+5):(count+14)])==10){
        while (Tl[[i]][(count+5)]){
          count<- ((count+5) + i)
        }
      }
    }
    for (i in 1:50){
      if(sum(Tl[[i+50]][count:(count+9)])==10){
        while (Tl[[i+50]][count]){
          count<- (count + i)
        }
      }
      if(sum(Tl[[i+50]][(count+1):(count+10)])==10){
        while (Tl[[i+50]][(count+1)]){
          count<- ((count+1) + i)
        }
      }
      if(sum(Tl[[i+50]][(count+2):(count+11)])==10){
        while (Tl[[i+50]][(count+2)]){
          count<- ((count+2) + i)
        }
      }
      if(sum(Tl[[i+50]][(count+3):(count+12)])==10){
        while (Tl[[i+50]][(count+3)]){
          count<- ((count+3) + i)
        }
      }
      if(sum(Tl[[i+50]][(count+4):(count+13)])==10){
        while (Tl[[i]][(count+4)]){
          count<- ((count+4) + i)
        }
      }
      if(sum(Tl[[i+50]][(count+5):(count+14)])==10){
        while (Tl[[i]][(count+5)]){
          count<- ((count+5) + i)
        }
      }
    }
    ##end of jump indel
    if (sum(T[(count+10):(count+25)]==rep(TRUE, 16))==16){
      count<- count+10
      while(T[count]){
        count<- count+1
      } 
    }
    if (sum(T[(count+10):(count+25)]==rep(TRUE, 16))==16){
      count<- count+10
      while(T[count]){
        count<- count+1
      } 
    }
    for (i in 1:50){###jump from INDEL
      if(sum(Tl[[i]][count:(count+9)])==10){
        while (Tl[[i]][count]){
          count<- (count + i)
        }
      }
    }
    for (i in 1:50){
      if(sum(Tl[[i+50]][count:(count+9)])==10){
        while (Tl[[i+50]][count]){
          count<- (count + i)
        }
      }
    }##end of jump indel
    if (sum(T[(count+10):(count+25)]==rep(TRUE, 16))==16){
      count<- count+10
      while(T[count]){
        count<- count+1
      } 
    }
    if (sum(T[(count+10):(count+25)]==rep(TRUE, 16))==16){
      count<- count+10
      while(T[count]){
        count<- count+1
      } 
    }
    for (i in 1:50){###jump from INDEL
      if(sum(Tl[[i]][count:(count+9)])==10){
        while (Tl[[i]][count]){
          count<- (count + i)
        }
      }
    }
    for (i in 1:50){
      if(sum(Tl[[i+50]][count:(count+9)])==10){
        while (Tl[[i+50]][count]){
          count<- (count + i)
        }
      }
    }##end of jump indel
    if (sum(T[(count+10):(count+25)]==rep(TRUE, 16))==16){
      count<- count+10
      while(T[count]){
        count<- count+1
      } 
    }
    if (sum(T[(count+10):(count+25)]==rep(TRUE, 16))==16){
      count<- count+10
      while(T[count]){
        count<- count+1
      } 
    }
    for (i in 1:50){###jump from INDEL
      if(sum(Tl[[i]][count:(count+9)])==10){
        while (Tl[[i]][count]){
          count<- (count + i)
        }
      }
    }
    for (i in 1:50){
      if(sum(Tl[[i+50]][count:(count+9)])==10){
        while (Tl[[i+50]][count]){
          count<- (count + i)
        }
      }
    }##end of jump indel
    if (sum(T[(count+10):(count+25)]==rep(TRUE, 16))==16){
      count<- count+10
      while(T[count]){
        count<- count+1
      } 
    }
    if (sum(T[(count+10):(count+25)]==rep(TRUE, 16))==16){
      count<- count+10
      while(T[count]){
        count<- count+1
      } 
    }
    if (sum(T[(count+10):(count+25)]==rep(TRUE, 16))==16){
      count<- count+10
      while(T[count]){
        count<- count+1
      } 
    }
    if (sum(T[(count+10):(count+25)]==rep(TRUE, 16))==16){
      count<- count+10
      while(T[count]){
        count<- count+1
      } 
    }
    if (sum(T[(count+10):(count+25)]==rep(TRUE, 16))==16){
      count<- count+10
      while(T[count]){
        count<- count+1
      } 
    }
    if (sum(T[(count+10):(count+25)]==rep(TRUE, 16))==16){
      count<- count+10
      while(T[count]){
        count<- count+1
      } 
    }
    if (sum(T[(count+10):(count+25)]==rep(TRUE, 16))==16){
      count<- count+10
      while(T[count]){
        count<- count+1
      } 
    }
    if (sum(T[(count+10):(count+25)]==rep(TRUE, 16))==16){
      count<- count+10
      while(T[count]){
        count<- count+1
      } 
    }
    if (sum(T[(count+10):(count+25)]==rep(TRUE, 16))==16){
      count<- count+10
      while(T[count]){
        count<- count+1
      } 
    }
    if (sum(T[(count+10):(count+25)]==rep(TRUE, 16))==16){
      count<- count+10
      while(T[count]){
        count<- count+1
      } 
    }
    if (sum(T[(count+10):(count+75)]==rep(TRUE, 66))==66){
      count<- count+10
      while(T[count]){
        count<- count+1
      } 
    }
    for (i in 1:50){###jump from INDEL
      if(sum(Tl[[i]][count:(count+9)])==10){
        while (Tl[[i]][count]){
          count<- (count + i)
        }
      }
    }
    for (i in 1:50){
      if(sum(Tl[[i+50]][count:(count+9)])==10){
        while (Tl[[i+50]][count]){
          count<- (count + i)
        }
      }
    }##end of jump indel
    if (sum(T[(count+10):(count+75)]==rep(TRUE, 66))==66){
      count<- count+10
      while(T[count]){
        count<- count+1
      } 
    }
    for (i in 1:50){###jump from INDEL
      if(sum(Tl[[i]][count:(count+9)])==10){
        while (Tl[[i]][count]){
          count<- (count + i)
        }
      }
    }
    for (i in 1:50){
      if(sum(Tl[[i+50]][count:(count+9)])==10){
        while (Tl[[i+50]][count]){
          count<- (count + i)
        }
      }
    }##end of jump indel
    if (sum(T[(count+10):(count+75)]==rep(TRUE, 66))==66){
      count<- count+10
      while(T[count]){
        count<- count+1
      } 
    }
    for (i in 1:50){###jump from INDEL
      if(sum(Tl[[i]][count:(count+9)])==10){
        while (Tl[[i]][count]){
          count<- (count + i)
        }
      }
    }
    for (i in 1:50){
      if(sum(Tl[[i+50]][count:(count+9)])==10){
        while (Tl[[i+50]][count]){
          count<- (count + i)
        }
      }
    }##end of jump indel
    if (sum(T[(count+10):(count+175)]==rep(TRUE, 166))==66){
      count<- count+10
      while(T[count]){
        count<- count+1
      } 
    }
    count
  }
  
  IR2start<- function(True.sequence.finder, IR.length, genome){
    return(length(genome)-(True.sequence.finder+IR.length-2))
  }
  
  #declaration
  if (parallel){
    phase.difference<- p.d(genome)
  }
  else{
    phase.difference<- phase.detector(genome)
  }
  Trsf<- True.sequence.finder(genome, phase.difference)
  IR1s<- IR1start(phase.difference, Trsf)
  IR.l<- IR.length(IR1s, Trsf, genome)
  IR2s<- IR2start(Trsf, IR.l, genome)
  
  #calculation
  return(c(IR1s, IR2s, IR.l, length(genome)))#returning the start of IR one and two follow by their lenght and the lenght of genome
  #gives a vector with four elemens as the start of the IRb and IRa and their length and the lenght of the genome sequence
}
IRtab<- function(GBFiles){#making the table of the IR information
  l<-length(GBFiles)
  IRtable<- matrix(0, l, 5)
  for (i in 1:l){
    IRtable[i,1]<-GBFiles[i]
    IRtable[i, 2:5]<- IRinfo(FasExtract(read.gb(GBFiles[i])))
  }
  return(IRtable)
}
intronCordinates<- function(gb){#The Gene bank file as the input an
  t<- gb[grep("  intron ", gb)]
  m<- matrix(0, length(t), 2)
  for (i in 1:length(t)){
    crude<-strsplit(t[i], " ")[[1]]
    cord<-crude[which(crude!=rep("", length(crude)))]
    if(length(cord)!= 2){
      stop("Check the GeneBank file: error with their intron row(s)","\n")
    }
    cord<-cord[2]
    if(length(strsplit(cord, "complement")[[1]])==2){
      cord<-strsplit(cord, "complement")[[1]][2]
      m[i,]<- strsplit(cord, "\\..")[[1]][c(2,1)]
    }
    else if (length(strsplit(cord, "complement")[[1]])==1){
      m[i,]<- strsplit(cord, "\\..")[[1]]
    } 
    else {
      stop("Check the Genebank file: error with cordinates of the intron row(s)")
    }
  }
  m<- gsub(")", "", m)
  m<- gsub("\\(", "", m)
  return(m)#giving the table format of the cordinates of the genes with introns (need imporvment to return the gene name as well and their exoin cordinates too)
}  
gene.name<-function(gb, type){
  if(type=="gene"){
    t <- gb[grep("  gene  ", gb)+1]
    for (i in 1:length(t)){
      crude<-strsplit(t[i], " ")[[1]]
      if(length(grep("=", crude))==0){#if it cannot find the gene name in the following line go to the next
        t[i] <- gb[grep("  gene  ", gb)+2][i]
        crude<- strsplit(t[i], " ")[[1]]
      }
      if(length(grep("=", crude))==0){#if it cannot find the gene name in the following line go to the next
        t[i] <- gb[grep("  gene  ", gb)+3][i]
        crude<- strsplit(t[i], " ")[[1]]
      }
      if(length(grep("=", crude))==0){#if it cannot find the gene name in the following line go to the next
        t[i] <- gb[grep("  gene  ", gb)+4][i]
        crude<- strsplit(t[i], " ")[[1]]
      }
      crude.name<-crude[which(crude!=rep("", length(crude)))]
      t[i]<-gsub("\"", "", strsplit(crude.name, "=")[[1]][2])
    }
    na<-which(is.na(t))
    if(length(na) > 0){
      for (i in 1:length(na)){
        warning(paste(paste(paste("Gene No.", na[i], " "), paste("is not properly named and is deleted from the list."), ""), paste("Check the gb file on line", which(gb==gb[grep("gene ", gb)+1][na]), ""), ""))
      }
    }
  }
  else if(type=="tRNA"){
    t <- gb[grep(" tRNA  ", gb)+1]
    for (i in 1:length(t)){
      crude<-strsplit(t[i], " ")[[1]]
      if(length(grep("=", crude))==0){
        t[i] <- gb[grep("tRNA  ", gb)+2][i]
        crude<- strsplit(t[i], " ")[[1]]
      }
      if(length(grep("=", crude))==0){#if it cannot find the gene name in the following line go to the next
        t[i] <- gb[grep("  tRNA  ", gb)+3][i]
        crude<- strsplit(t[i], " ")[[1]]
      }
      if(length(grep("=", crude))==0){#if it cannot find the gene name in the following line go to the next
        t[i] <- gb[grep("  tRNA  ", gb)+4][i]
        crude<- strsplit(t[i], " ")[[1]]
      }
      crude.name<-crude[which(crude!=rep("", length(crude)))]
      t[i]<-gsub("\"", "", strsplit(crude.name, "=")[[1]][2])
    }
    na<-which(is.na(t))
    if(length(na) > 0){
      for (i in 1:length(na)){
        warning(paste(paste(paste("Gene No.", na[i], " "), paste("is not properly named and is deleted from the list."), ""), paste("Check the gb file on line", which(gb==gb[grep("tRNA ", gb)+1]), ""), ""), "\n")
      }
    }    	
  }
  else if(type=="rRNA"){
    t <- gb[grep("rRNA  ", gb)+1]
    for (i in 1:length(t)){
      crude<-strsplit(t[i], " ")[[1]]
      if(length(grep("=", crude))==0){
        t[i] <- gb[grep("rRNA  ", gb)+2][i]
        crude<- strsplit(t[i], " ")[[1]]
      }
      if(length(grep("=", crude))==0){#if it cannot find the gene name in the following line go to the next
        t[i] <- gb[grep("  rRNA  ", gb)+3][i]
        crude<- strsplit(t[i], " ")[[1]]
      }
      if(length(grep("=", crude))==0){#if it cannot find the gene name in the following line go to the next
        t[i] <- gb[grep("  rRNA  ", gb)+4][i]
        crude<- strsplit(t[i], " ")[[1]]
      }
      crude.name<-crude[which(crude!=rep("", length(crude)))]
      t[i]<-gsub("\"", "", strsplit(crude.name, "=")[[1]][2])
    }
    na<-which(is.na(t))
    if(length(na) > 0){
      for (i in 1:length(na)){
        warning(paste(paste(paste("Gene No.", na[i], " "), paste("is not properly named and is deleted from the list."), ""), paste("Check the gb file on line", which(gb==gb[grep("gene ", gb)+1][na]), ""), ""))
      }
    }
  }
  else if(type=="mRNA"){
    t <- gb[grep("mRNA  ", gb)+1]
    for (i in 1:length(t)){
      crude<-strsplit(t[i], " ")[[1]]
      if(length(grep("=", crude))==0){
        t[i] <- gb[grep("mRNA  ", gb)+2][i]
        crude<- strsplit(t[i], " ")[[1]]
      }
      if(length(grep("=", crude))==0){#if it cannot find the gene name in the following line go to the next
        t[i] <- gb[grep("  mRNA  ", gb)+3][i]
        crude<- strsplit(t[i], " ")[[1]]
      }
      if(length(grep("=", crude))==0){#if it cannot find the gene name in the following line go to the next
        t[i] <- gb[grep("  mRNA  ", gb)+4][i]
        crude<- strsplit(t[i], " ")[[1]]
      }
      crude.name<-crude[which(crude!=rep("", length(crude)))]
      t[i]<-gsub("\"", "", strsplit(crude.name, "=")[[1]][2])
    }
    na<-which(is.na(t))
    if(length(na) > 0){
      for (i in 1:length(na)){
        warning(paste(paste(paste("Gene No.", na[i], " "), paste("is not properly named and is deleted from the list."), ""), paste("Check the gb file on line", which(gb==gb[grep("gene ", gb)+1][na]), ""), ""))
      }
    }
  }
  else {
    stop("The type should be defined as either gene or tRNA")
  }
  t<-t[!is.na(t)]
  return(t)
  #intermediate gene name function to substract the gene names of either gene or tRNA, mRNA or rRNA from their second line information(or third), the input is the gb file.
}
gb.gene.cor<- function(gb){
  gene<- gb[grep("  gene  ", gb)]
  trna<- gb[grep("  tRNA ", gb)]
  rrna<- gb[grep("  rRNA ", gb)]
  t<-c(gene, trna, rrna)
  m<- matrix(0, length(t), 2)
  for (i in 1:length(t)){
    if(strsplit(t[i], "")[[1]][length(strsplit(t[i],"")[[1]])]==","){#check the last characted if its ',' and add the next one from gb file to that
      #previous code--> if(strsplit(t[i], "")[[1]][length(strsplit(t[i],"")[[1]])]=="," && strsplit(t[i], "")[[1]][length(strsplit(t[i],"")[[1]])-1]==")")
      s<-t[i]
      t[i]<-paste(t[i], gsub(" ", "", gb[which(t[i]==gb)+1]), "")
    }
    t[i]<- gsub(", ", ",", t[i])
    t[i]<- gsub(") ", ")", t[i])
    if(strsplit(t[i], "")[[1]][length(strsplit(t[i],"")[[1]])]==","){#check again if the last characted if its ',' and add the next one from gb file to that
      #previous code--> if(strsplit(t[i], "")[[1]][length(strsplit(t[i],"")[[1]])]=="," && strsplit(t[i], "")[[1]][length(strsplit(t[i],"")[[1]])-1]==")")
      t[i]<-paste(t[i], gsub(" ", "", gb[which(s==gb)+2]), "")
    }
    t[i]<- gsub(", ", ",", t[i])
    t[i]<- gsub(") ", ")", t[i])
    m[i,2]<-strsplit(t[i], " ")[[1]][length(strsplit(t[i], " ")[[1]])]
  }
  names<-  c(gene.name(gb, "gene"), gene.name(gb, "tRNA"), gene.name(gb, "rRNA") )
  if(length(t)!=length(names)){
    stop("Error: Check the gene names of the gb file")
  }
  m[,1]<-names
  m #gives the crude table extracted from gb file for genes
}
gene.cordinates<- function(gb){
  gb.gene.cor.out<-gb.gene.cor(gb)
  l<- length(gb.gene.cor.out[,1])
  m<- matrix(0, l, 5)
  m[,1]<- gb.gene.cor.out[,1]
  for (i in 1:l){
    if(length(strsplit(gb.gene.cor.out[i,2], ",")[[1]])==2){
      m[i,2]<-sub("join\\(", "", strsplit(gb.gene.cor.out[i,2], ",")[[1]][1]) 
      m[i,2]<-sub("order\\(", "", m[i,2]) 
      m[i,4]<-sub(")", "", strsplit(gb.gene.cor.out[i,2], ",")[[1]][2])
      cord<-m[i,2]
      if(length(strsplit(cord, "complement")[[1]])==2){
        cord<-strsplit(cord, "complement")[[1]][2]
        m[i,c(2,3)]<- strsplit(cord, "\\..")[[1]][c(2,1)]
      }
      else if (length(strsplit(cord, "complement")[[1]])==1){
        m[i,c(2,3)]<- strsplit(cord, "\\..")[[1]]
      } 
      else {
        stop("Check the Genebank file: error with cordinates of the intron row(s)")
      }
      cord<-m[i,4]
      if(length(strsplit(cord, "complement")[[1]])==2){
        cord<-strsplit(cord, "complement")[[1]][2]
        m[i,c(4,5)]<- strsplit(cord, "\\..")[[1]][c(2,1)]
      }
      else if (length(strsplit(cord, "complement")[[1]])==1){
        m[i,c(4,5)]<- strsplit(cord, "\\..")[[1]]
      } 
      else {
        stop("Check the Genebank file: error with cordinates of the intron row(s)")
      }
    }
    else {
      m[i,2]<-gb.gene.cor.out[i,2]
      cord<-m[i,2]
      if(length(strsplit(cord, "complement")[[1]])==2){
        cord<-strsplit(cord, "complement")[[1]][2]
        m[i,c(2,3)]<- strsplit(cord, "\\..")[[1]][c(2,1)]
      }
      else if (length(strsplit(cord, "complement")[[1]])==1){
        m[i,c(2,3)]<- strsplit(cord, "\\..")[[1]]
      } 
      else {
        stop("Check the Genebank file: error with cordinates of the intron row(s)")
      }
    }
  }		
  m<- gsub(")", "", m)
  m<- gsub("\\(", "", m)
  m<- gsub("<", "", m)
  m<- gsub(">", "", m)
  return(m)#gives a polished table of genes and their cordinates, as the gene names their start and end of first part and the second parts as the first to fifth columns. The orders are reflected with the swapping of values
}   
JunctRadiusFinder<- function(gene.cordinates, IRinfo, J.pos, radius, silence=TRUE){#function to find the genes near the junction cordinate based on the output of the gene.cordinate function
  if(J.pos==1){
    J<- IRinfo[1]
  }
  else if(J.pos==2){
    J<- IRinfo[1]+IRinfo[3]
  }
  else if(J.pos==3){
    J<- IRinfo[2]
  }
  else if(J.pos==4){
    J<- IRinfo[2]+IRinfo[3]
  }
  else {stop("J.pos missing or out bound. It should be either 1,2,3, or 4 for JLB, JSB, JSA, and JLA, recpectively ")}
  g<-gene.cordinates
  l<- ncol(g)
  r<- radius
  gs<- IRinfo[4]
  t1<- subset(g, J-r<as.numeric(g[,2]) & as.numeric(g[,2])<J+r)
  t2<- subset(g, J-r<as.numeric(g[,3]) & as.numeric(g[,3])<J+r)
  t3<- subset(g, J-r<as.numeric(g[,4]) & as.numeric(g[,4])<J+r)
  t4<- subset(g, J-r<as.numeric(g[,5]) & as.numeric(g[,5])<J+r)
  for (i in 1:nrow(g)){
    if(sum(g[i, ]=="0")==2){
      g[i,4:5]<-c(NA, NA)
    }
  }
  t5<- subset(g, J-r<(as.numeric(g[,2])+gs) & (as.numeric(g[,2])+gs)<J+r)
  t6<- subset(g, J-r<(as.numeric(g[,3])+gs) & (as.numeric(g[,3])+gs)<J+r)
  t7<- subset(g, J-r<(as.numeric(g[,4])+gs) & (as.numeric(g[,4])+gs)<J+r)
  t8<- subset(g, J-r<(as.numeric(g[,5])+gs) & (as.numeric(g[,5])+gs)<J+r)
  t<-unique(do.call("rbind", list(t1, t2, t3, t4, t5, t6, t7, t8)))
  if((length(t)+silence)==0){
    warning("There is no gene on the specified junction with the given radius, increasing the radius might help")
  }
  return(t)
}
JunctNoFinder<- function(gene.cordinates, IRinfo, J.pos ,Number){#the same function like it precusor but returns the number of genes found near the junction site
  g<-gene.cordinates
  l<- ncol(g)
  n<-Number
  counter=0
  count=0
  while (nrow(JunctRadiusFinder(g, IRinfo, J.pos ,count)) <= n-1){
    count<-count+1
  }
  t<-JunctRadiusFinder(g, IRinfo, J.pos ,count)
  t[is.na(t)]<- "0"
  tup<- matrix(0, n, 3)
  for (i in 1:n){
    dist<-abs(as.numeric(t[i,][2:5])-J)
    if (which(dist==min(dist))>3){
      tup[i,]<-t[i, c(1,4,5)] 
    }
    else {
      tup[i,]<-t[i, c(1,2,3)]
    }
  }
  return(tup)
}
thFinder<- function(gene.cordinates, IRinfo, J.pos, n){#finding the nth closest gene to the junctions site
  g<-gene.cordinates
  l<- ncol(g)
  if(n==1){
    return(JunctNoFinder(g, IRinfo, J.pos, n))
  }
  else{
    previous<-JunctNoFinder(g, IRinfo, J.pos, n-1)
    current<- JunctNoFinder(g, IRinfo, J.pos, n  )
    return(current[!current[,2] %in% previous[,2],])
  }
}
RadiusNoFinder<- function(gene.cordinates, IRinfo, J.pos, Number){
  g<-gene.cordinates
  l<- ncol(g)
  n<-Number
  counter=0
  count=0
  while (nrow(JunctRadiusFinder(g, IRinfo, J.pos, count)) <= n-1){
    count<-count+1
  }
  return(count)
}
GinRadius<- function(Radius, J.pos, track){#thosse who will be plotted
  if(J.pos==1){
    pc<-15
    J<- IRList[[track]][1]
  }
  else if(J.pos==2){
    pc<-40
    J<- IRList[[track]][1]+IRList[[track]][3]
  }
  else if(J.pos==3){
    pc<-70
    J<- IRList[[track]][2]
  }
  else if(J.pos==4){
    pc<-95
    J<- IRList[[track]][2]+IRList[[track]][3]
  }
  else {stop("J.pos missing or out bound. It should be either 1,2,3, or 4 for JLB, JSB, JSA, and JLA, recpectively ")}
  t<- JunctRadiusFinder(GeneList[[track]], IRList[[track]], J.pos , Radius)
  t[is.na(t)]<- "0"
  n<- length(t[,1])
  tup<- matrix(0, n, 3)
  for (i in 1:n){
    dist<-abs(as.numeric(t[i,][2:5])-J)
    if (which(dist==min(dist))>3){
      tup[i,]<-t[i, c(1,4,5)] 
    }
    else {
      tup[i,]<-t[i, c(1,2,3)]
    }
  }               
  return(tup)
}
JG.plotter<- function(Radius, J.pos, track){#JunctionGene plotter
  #' Junction site gene plotter
  #' 
  #' Plotting the genes in the vicinity of the junction site of the chloroplast. 
  #' @param gene.cordinate The output of the gene.cordinate function with the name of the genes and their cordinates on genome
  #' @param junction.cordinate The cordinate of the junction site for which the genes are to be plotted
  #' @param n The number of the genes that is needed to be plotted
  #' @param J.pos The position of the J site: JLB, JSB, JSA, and JLA for 1,2,3, and 4 respectively
  #' @param margin The margin value to be added to the interval as the proportion of the minimum lenght of the nth gene in vicinity
  #' @param track The track on which the genes are to be plotted, strating from the bottom to up as integers 1,2,...
  #' @export
  if(J.pos==1){
    pc<-15
    J<- IRList[[track]][1]
  }
  else if(J.pos==2){
    pc<-40
    J<- IRList[[track]][1]+IRList[[track]][3]
  }
  else if(J.pos==3){
    pc<-70
    J<- IRList[[track]][2]
  }
  else if(J.pos==4){
    pc<-95
    J<- IRList[[track]][2]+IRList[[track]][3]
  }
  else {stop("J.pos missing or out bound. It should be either 1,2,3, or 4 for JLB, JSB, JSA, and JLA, recpectively ")}
  t<- JunctRadiusFinder(GeneList[[track]], IRList[[track]], J.pos , Radius)    
  t[is.na(t)]<- "0"
  n<- length(t[,1])
  tup<- matrix(0, n, 3)
  for (i in 1:n){
    dist<-abs(as.numeric(t[i,][2:5])-J)
    if (which(dist==min(dist))>3){
      tup[i,]<-t[i, c(1,4,5)] 
    }
    else {
      tup[i,]<-t[i, c(1,2,3)]
    }
  }               
  bw<- 10/Radius
  Rcord<- matrix(0, n, 2)
  Rcord[,1]<-as.numeric(tup[,2])
  Rcord[,2]<-as.numeric(tup[,3])
  Pcord<- (Rcord-J)*bw+pc
  gcol<- function(tup){
    l<-length(tup[,1])
    col<-numeric(l)
    for (i in 1:l){
      if(tup[i,1]=="trnH"){col[i]<- "burlywood4"}
      else if(tup[i,1]=="trnN"){col[i]<- "violet"}
      else if(tup[i,1]=="psbA"){col[i]<- "purple2"}
      else if(tup[i,1]=="rps19"){col[i]<- "firebrick1"}
      else if(tup[i,1]=="ycf1"){col[i]<- "dodgerblue3"}
      else if(tup[i,1]=="ndhF"){col[i]<- "darkred"}
      else if(tup[i,1]=="rpl22"){col[i]<- "navy"}
      else if(tup[i,1]=="rpl2"){col[i]<- "forestgreen"}
      else if(tup[i,1]=="rpl23"){col[i]<- "palevioletred4"}
      else if(tup[i,1]=="ycf2"){col[i]<- "blue3"}
      else if(tup[i,1]=="rrn23"){col[i]<- "deeppink"}
      else if(tup[i,1]=="rrn4"){col[i]<- "deeppink3"}
      else if(tup[i,1]=="rrn5"){col[i]<- "deeppink4"}
      else if(tup[i,1]=="chlL"){col[i]<- "darkorange1"}
      else if(tup[i,1]=="chlN"){col[i]<- "darkorange3"}
      else if(tup[i,1]=="rps12"){col[i]<- "slategray1"}
      else if(tup[i,1]=="rps7"){col[i]<- "slategray3"}
      else if(tup[i,1]=="rps3"){col[i]<- "slategrey"}
      else if(tup[i,1]=="16rrn"){col[i]<- "darkkhaki"}
      else if(tup[i,1]=="trnV"){col[i]<- "slateblue4"}
      else if(tup[i,1]=="trnM"){col[i]<- "steelblue4"}
      else if(tup[i,1]=="trnL"){col[i]<- "steelblue"}
      else if(tup[i,1]=="ccsA"){col[i]<- "goldenrod1"}
      else if(tup[i,1]=="rpl32"){col[i]<- "darkolivegreen4"}
      else if(tup[i,1]=="ndhB"){col[i]<- "lightseagreen"}
      else if(tup[i,1]=="rps15"){col[i]<- "lightsalmon1"}
      else if(tup[i,1]=="ndhH"){col[i]<- "navajowhite"}
      else if(tup[i,1]=="ndhA"){col[i]<- "palevioletred2"}
      else {col[i]<- "black"}
    }
    return(col)
  }
  if (J.pos==4){
    Pcord[Pcord < 0]<- ((Rcord[which(Pcord< 0)] + IRList[[track]][4])-J)*bw+pc
  }
  for (i in 1:n){
    if(Rcord[i,1]>Rcord[i,2]){
      x1<-min(Pcord[i,1], pc+10)
      x2<-max(Pcord[i,2], pc-10)
      for (j in seq(0.10, 0.70, 0.05)){
        segments(x1, track*5+j+5, x2, track*5+j+5, lwd=1, col=paste(gcol(tup)[i]))
      }
    }
    else {
      for (j in seq(1.1, 1.7, 0.05)){
        segments(max(Pcord[i,1], pc-10), track*5-j+5, min(Pcord[i,2], pc+10), track*5-j+5, lwd=1, col=paste(gcol(tup)[i]))
      }
    }
  }
}
GN.plotter<- function(Radius, J.pos, track){#GeneName plotter
  #' Gene Name plotter
  #' 
  #' Plotting the gene names on a given gene which is already plotted on the tracks of the IR plot
  if(J.pos==1){
    pc<-15
    J<- IRList[[track]][1]
  }
  else if(J.pos==2){
    pc<-40
    J<- IRList[[track]][1]+IRList[[track]][3]
  }
  else if(J.pos==3){
    pc<-70
    J<- IRList[[track]][2]
  }
  else if(J.pos==4){
    pc<-95
    J<- IRList[[track]][2]+IRList[[track]][3]
  }
  else {stop("J.pos missing or out bound. It should be either 1,2,3, or 4 for JLB, JSB, JSA, and JLA, recpectively ")}
  txtincol<- "white"
  txtoutcol="black"
  txtfont<- 4
  numfont<- 1
  numcex<- 0.44
  txtcex<- 0.46
  t<- JunctRadiusFinder(GeneList[[track]], IRList[[track]], J.pos , Radius)
  t[is.na(t)]<- "0"
  n<- length(t[,1])
  tup<- matrix(0, n, 3)
  for (i in 1:n){
    dist<-abs(as.numeric(t[i,][2:5])-J)
    if (which(dist==min(dist))>3){
      tup[i,]<-t[i, c(1,4,5)] 
    }
    else {
      tup[i,]<-t[i, c(1,2,3)]
    }
  }               
  bw<- 10/Radius
  Rcord<- matrix(0, n, 2)
  Rcord[,1]<-as.numeric(tup[,2])
  Rcord[,2]<-as.numeric(tup[,3])
  Pcord<- (Rcord-J)*bw+pc
  if (J.pos==4){
    Pcord[Pcord < 0]<- ((Rcord[which(Pcord< 0)] + IRList[[track]][4])-J)*bw+pc
  }
  for (i in 1:n){
    if(Rcord[i,1]>Rcord[i,2]){
      x1<-min(Pcord[i,1], pc+10)
      x2<-max(Pcord[i,2], pc-10)
      if (min(x1, x2) >= 104){
        text(103.5, track*5+1.05+5, tup[i,1], cex=txtcex, col=txtoutcol, font=txtfont)
      }
      else if (abs(x1-x2) >= 9.7){
        text(min(x1, x2)+1.8, track*5+0.35+5, tup[i,1], cex=txtcex, col=txtincol, font=txtfont)
        text(max(x1, x2)+1, track*5+0.3+5, paste(Rcord[i,1]-Rcord[i,2], "bp", " "), cex=numcex, col=txtincol, font=numfont, pos=2)
      }
      else if (abs(x1-x2) >= 3 & abs(x1-x2) < 9.7){
        text(((min(x1, x2)+max(x1, x2))/2), track*5+0.35+5, tup[i,1], cex=txtcex, col=txtincol, font=txtfont)
      }
      else if (abs(x1-x2) < 3){
        text(((min(x1, x2)+max(x1, x2))/2), track*5+1.15+5, tup[i,1], cex=txtcex, col=txtoutcol, font=txtfont)
      }
    }
    else {
      x1<-max(Pcord[i,1], pc-10)
      x2<-min(Pcord[i,2], pc+10)
      if (abs(x1-x2) >= 9.7){
        text(min(x1, x2)+1.8, track*5-1.40+5, tup[i,1], cex=txtcex, col=txtincol, font=txtfont)
        text(max(x1, x2)+1, track*5-1.45+5, paste(Rcord[i,2]-Rcord[i,1], "bp", " "), cex=numcex, col=txtincol, font=numfont, pos=2)
      }
      else if (abs(x1-x2) >= 3 & abs(x1-x2) < 9.7){
        text(((min(x1, x2)+max(x1, x2))/2), track*5-1.40+5, tup[i,1], cex=txtcex, col=txtincol, font=txtfont)
      }
      else if (abs(x1-x2) < 3){
        text(((min(x1, x2)+max(x1, x2))/2), track*5-2.25+5, tup[i,1], cex=txtcex, col=txtoutcol, font=txtfont)
      }
    }
  }
}
OJ.plotter<- function(Radius, J.pos, track){#GeneName plotter
  #' On Junction plotter
  #' 
  #' Plotting the fine tuned narrow lines showing the limits of the genes which are passing through the junction sites with their bp
  if(J.pos==1){
    pc<-15
    J<- IRList[[track]][1]
  }
  else if(J.pos==2){
    pc<-40
    J<- IRList[[track]][1]+IRList[[track]][3]
  }
  else if(J.pos==3){
    pc<-70
    J<- IRList[[track]][2]
  }
  else if(J.pos==4){
    pc<-95
    J<- IRList[[track]][2]+IRList[[track]][3]
  }
  else {stop("J.pos missing or out bound. It should be either 1,2,3, or 4 for JLB, JSB, JSA, and JLA, recpectively ")}
  txtincol<- "white"
  txtoutcol="black"
  txtfont<- 4
  numfont<- 1
  numcex<- 0.44
  txtcex<- 0.46
  t<- JunctRadiusFinder(GeneList[[track]], IRList[[track]], J.pos , Radius)
  t[is.na(t)]<- "0"
  n<- length(t[,1])
  tup<- matrix(0, n, 3)
  for (i in 1:n){
    dist<-abs(as.numeric(t[i,][2:5])-J)
    if (which(dist==min(dist))>3){
      tup[i,]<-t[i, c(1,4,5)] 
    }
    else {
      tup[i,]<-t[i, c(1,2,3)]
    }
  }               
  bw<- 10/Radius
  Rcord<- matrix(0, n, 2)
  Rcord[,1]<-as.numeric(tup[,2])
  Rcord[,2]<-as.numeric(tup[,3])
  Pcord<- (Rcord-J)*bw+pc
  if (J.pos==4){
    Pcord[Pcord < 0]<- ((Rcord[which(Pcord< 0)] + IRList[[track]][4])-J)*bw+pc
  }
  for (i in 1:n){
    if(Rcord[i,1]>Rcord[i,2]){
      x1<-min(Pcord[i,1], pc+10)
      x2<-max(Pcord[i,2], pc-10)
      if (Rcord[i,1] > J & Rcord[i,2] <J ){
        Arrows(min(x1, x2), track*5+1.1+5, pc-0.15, track*5+1.1+5, arr.type = "T", cex=0.5, arr.length = 0.12, lwd=0.5, arr.width = 0.4)
        Arrows(pc-0.15, track*5+1.1+5, min(x1, x2), track*5+1.1+5, arr.type = "T", cex=0.5, arr.length = 0.12, lwd=0.5, arr.width = 0.4)
        Arrows(pc+0.15, track*5+1.1+5, max(x1, x2), track*5+1.1+5, arr.type = "T", cex=0.5, arr.length = 0.12, lwd=0.5, arr.width = 0.4)
        Arrows(max(x1, x2), track*5+1.1+5, pc+0.15, track*5+1.1+5, arr.type = "T", cex=0.5, arr.length = 0.12, lwd=0.5, arr.width = 0.4)
        text(pc-1.4, track*5+1.4+5, paste(Rcord[i, 1]-J+1, "bp", sep=" "), cex=0.4, pos=4)#up-right
        text(pc+1.4, track*5+1.4+5, paste(J-Rcord[i, 2], "bp", sep=" "), cex=0.4, pos=2)#up-left
      }
    }
    else {
      x1<-max(Pcord[i,1], pc-10)
      x2<-min(Pcord[i,2], pc+10)
      if (Rcord[i,2] > J & Rcord[i,1] <J ){
        Arrows(max(x1, x2), track*5-2.1+5, pc+0.15, track*5-2.1+5, arr.type = "T", cex=0.5, arr.length = 0.12, lwd=0.5, arr.width = 0.4)
        Arrows(pc+0.15, track*5-2.1+5, max(x1, x2), track*5-2.1+5, arr.type = "T", cex=0.5, arr.length = 0.12, lwd=0.5, arr.width = 0.4)
        Arrows(pc-0.15, track*5-2.1+5, min(x1, x2), track*5-2.1+5, arr.type = "T", cex=0.5, arr.length = 0.12, lwd=0.5, arr.width = 0.4)
        Arrows(min(x1, x2), track*5-2.1+5, pc-0.15, track*5-2.1+5, arr.type = "T", cex=0.5, arr.length = 0.12, lwd=0.5, arr.width = 0.4)
        text(pc-1.4, track*5-2.6+5, paste(Rcord[i, 2]-J, "bp", sep=" "), cex=0.4, pos=4)#low-right
        text(pc+1.4, track*5-2.6+5, paste(J-Rcord[i, 1], "bp", sep=" "), cex=0.4, pos=2)#low-left
      }
    }
  }
}
JD.plotter<- function(Radius, J.pos, track){#GeneName plotter
  #' Junction Distance plotter
  #' 
  #' plotting the narrow lines of the distance of the genes for the junction sites which are not passing through any gene and their bp
  if(J.pos==1){
    pc<-15
    J<- IRList[[track]][1]
  }
  else if(J.pos==2){
    pc<-40
    J<- IRList[[track]][1]+IRList[[track]][3]
  }
  else if(J.pos==3){
    pc<-70
    J<- IRList[[track]][2]
  }
  else if(J.pos==4){
    pc<-95
    J<- IRList[[track]][2]+IRList[[track]][3]
  }
  else {stop("J.pos missing or out bound. It should be either 1,2,3, or 4 for JLB, JSB, JSA, and JLA, recpectively ")}
  txtincol<- "white"
  txtoutcol="black"
  txtfont<- 4
  numfont<- 1
  numcex<- 0.44
  txtcex<- 0.46
  t<- JunctRadiusFinder(GeneList[[track]], IRList[[track]], J.pos , Radius)
  t[is.na(t)]<- "0"
  n<- length(t[,1])
  tup<- matrix(0, n, 3)
  for (i in 1:n){
    dist<-abs(as.numeric(t[i,][2:5])-J)
    if (which(dist==min(dist))>3){
      tup[i,]<-t[i, c(1,4,5)] 
    }
    else {
      tup[i,]<-t[i, c(1,2,3)]
    }
  }               
  bw<- 10/Radius
  Rcord<- matrix(0, n, 2)
  Rcord[,1]<-as.numeric(tup[,2])
  Rcord[,2]<-as.numeric(tup[,3])
  Pcord<- (Rcord-J)*bw+pc
  if (J.pos==4){
    ind<-Pcord < 0
    Pcord[Pcord < 0]<- ((Rcord[which(Pcord< 0)] + IRList[[track]][4])-J)*bw+pc
    Rcord[ind]<- Rcord[ind] + IRList[[track]][4]
  }
  counter<- 0
  for (i in 1:n){
    if(Rcord[i,1]>Rcord[i,2]){
      x1<-min(Pcord[i,1], pc+10)
      x2<-max(Pcord[i,2], pc-10)
      if (Rcord[i,1] > J & Rcord[i,2] <J ){
        counter<- counter+1
      }
    }
    else {
      x1<-max(Pcord[i,1], pc-10)
      x2<-min(Pcord[i,2], pc+10)
      if (Rcord[i,2] > J & Rcord[i,1] <J ){
        counter<- counter+1
      }
    }
  }
  if (counter==0){#find the closest gene to the junction site
    nearest<- which(abs(Pcord-pc)==min(abs(Pcord-pc)))
    col.cor<- floor((nearest-0.01)/n)+1
    row.cor<- nearest-n*floor((nearest-0.01)/n)###Now we have the row of the nearest gene
    #with this setting the position of the zero (genes tangant to the junction site) will not be plotted. If interested either put "=" for the middle condition of the left or right binary operator or better develop zero only handling if function
    if     (Rcord[row.cor, 1] > Rcord[row.cor, 2] & pc-Pcord[row.cor, col.cor] < 0 & min(Pcord[row.cor, 1], pc+10) - max(Pcord[row.cor, 2], pc-10) > 3 ){#top,right, big
      curvedarrow (from=c(pc+(Pcord[row.cor, col.cor]-pc)/2, track*5+0.3+5), to=c(pc+(Pcord[row.cor, col.cor]-pc)/2 + 3, track*5+1.3+5), curve = -0.21, lwd=0.6, arr.type = "curved", arr.col = "white", arr.length=0.08, arr.lwd=0.4, arr.pos=0.69, endhead=TRUE)
      text(pc+(Pcord[row.cor, col.cor]-pc)/2 + 4, track*5+1.3+0.2+5, paste(Rcord[row.cor, 2]-J, "bp", sep=" "), cex=0.4)
    }
    else if(Rcord[row.cor, 1] > Rcord[row.cor, 2] & pc-Pcord[row.cor, col.cor] < 0 & min(Pcord[row.cor, 1], pc+10) - max(Pcord[row.cor, 2], pc-10) <= 3 ){#top, right, small
      arrows(pc+(Pcord[row.cor, col.cor]-pc)/2, track*5+0.3+5 , pc+(Pcord[row.cor, col.cor]-pc)/2 -2 , track*5+1.3+5, angle = 15, length = 0.05, lwd=0.6)
      text(pc+(Pcord[row.cor, col.cor]-pc)/2 - 4.5 , track*5+1.3+0.2+5, paste(Rcord[row.cor, 2]-J, "bp", sep=" "), cex=0.4)
    }
    else if(Rcord[row.cor, 1] < Rcord[row.cor, 2] & pc-Pcord[row.cor, col.cor] < 0 & abs(max(Pcord[row.cor, 1], pc-10) - min(Pcord[row.cor, 2], pc+10)) > 3 ){#low, right, big
      curvedarrow (from=c(pc+(Pcord[row.cor, col.cor]-pc)/2, track*5-1.3+5), to=c(pc+(Pcord[row.cor, col.cor]-pc)/2 + 3, track*5-2.3+5), curve = 0.21, lwd=0.6, arr.type = "curved", arr.col = "white", arr.length=0.08, arr.lwd=0.4, arr.pos=0.69, endhead=TRUE) 
      text(pc+(Pcord[row.cor, col.cor]-pc)/2 + 4, track*5-2.3-0.2+5, paste(Rcord[row.cor, 1]-J, "bp", sep=" "), cex=0.4)
    }
    else if(Rcord[row.cor, 1] < Rcord[row.cor, 2] & pc-Pcord[row.cor, col.cor] < 0 & abs(max(Pcord[row.cor, 1], pc-10) - min(Pcord[row.cor, 2], pc+10)) <= 3 ){#low, right, small
      arrows(pc+(Pcord[row.cor, col.cor]-pc)/2, track*5-1.3+5 , pc+(Pcord[row.cor, col.cor]-pc)/2 -3 , track*5-2.3+5, angle = 15, length = 0.05, lwd=0.6)
      text(pc+(Pcord[row.cor, col.cor]-pc)/2 -4.5 , track*5-2.3-0.2+5, paste(Rcord[row.cor, 1]-J, "bp", sep=" "), cex=0.4)
    }
    else if(Rcord[row.cor, 1] < Rcord[row.cor, 2] & pc-Pcord[row.cor, col.cor] > 0 & abs(max(Pcord[row.cor, 1], pc-10) - min(Pcord[row.cor, 2], pc+10)) > 3 ){#low, left, big
      curvedarrow (from=c(pc+(Pcord[row.cor, col.cor]-pc)/2, track*5-1.3+5), to=c(pc+(Pcord[row.cor, col.cor]-pc)/2 - 3, track*5-2.3+5), curve = -0.21, lwd=0.6, arr.type = "curved", arr.col = "white", arr.length=0.08, arr.lwd=0.4, arr.pos=0.69, endhead=TRUE) 
      text(pc+(Pcord[row.cor, col.cor]-pc)/2 - 4, track*5-2.3-0.2+5, paste(J-Rcord[row.cor, 2], "bp", sep=" "), cex=0.4)
    }
    else if(Rcord[row.cor, 1] < Rcord[row.cor, 2] & pc-Pcord[row.cor, col.cor] > 0 & abs(max(Pcord[row.cor, 1], pc-10) - min(Pcord[row.cor, 2], pc+10)) <= 3 ){#low, left, small
      arrows(pc+(Pcord[row.cor, col.cor]-pc)/2, track*5-1.3+5 , pc+(Pcord[row.cor, col.cor]-pc)/2 + 3 , track*5-2.3+5, angle = 15, length = 0.05, lwd=0.6)
      text(pc+(Pcord[row.cor, col.cor]-pc)/2 + 4.5 , track*5-2.3-0.2+5, paste(J-Rcord[row.cor, 2], "bp", sep=" "), cex=0.4)
    }
    else if(Rcord[row.cor, 1] > Rcord[row.cor, 2] & pc-Pcord[row.cor, col.cor] > 0 & min(Pcord[row.cor, 1], pc+10) - max(Pcord[row.cor, 2], pc-10) > 3 ){#top, left, big
      curvedarrow (from=c(pc+(Pcord[row.cor, col.cor]-pc)/2, track*5+0.3+5), to=c(pc+(Pcord[row.cor, col.cor]-pc)/2 - 3, track*5+1.3+5), curve = 0.21, lwd=0.6, arr.type = "curved", arr.col = "white", arr.length=0.08, arr.lwd=0.4, arr.pos=0.69, endhead=TRUE)
      text(pc+(Pcord[row.cor, col.cor]-pc)/2 - 4.5 , track*5+1.3+0.2+5, paste(J-Rcord[row.cor, 1], "bp", sep=" "), cex=0.4)
    }
    else if(Rcord[row.cor, 1] > Rcord[row.cor, 2] & pc-Pcord[row.cor, col.cor] > 0 & min(Pcord[row.cor, 1], pc+10) - max(Pcord[row.cor, 2], pc-10) <= 3 ){#top, left, small
      arrows(pc+(Pcord[row.cor, col.cor]-pc)/2, track*5+0.3+5 , pc+(Pcord[row.cor, col.cor]-pc)/2 + 2 , track*5+1.3+5, angle = 15, length = 0.05, lwd=0.6)
      text(pc+(Pcord[row.cor, col.cor]-pc)/2 + 4 , track*5+1.3+0.2+5, paste(J-Rcord[row.cor, 1], "bp", sep=" "), cex=0.4)
    }
  }
}
Max.Radius<-function(J.pos, l, genelist, IRlist){
  if(J.pos==1){
    Radius<-550#680
  }
  if(J.pos==2){
    Radius<-100
  }
  if(J.pos==3){
    Radius<-700#1800
  }
  if(J.pos==4){
    Radius<-700#1000
  }
  R<- numeric(l)
  for (track in 1:l){
    t<- JunctRadiusFinder(GeneList[[track]], IRList[[track]], J.pos , Radius)
    while(nrow(t)==0){
      Radius<- Radius+(0.2)*Radius#0.2
      t<- JunctRadiusFinder(GeneList[[track]], IRList[[track]], J.pos , Radius)
    }
    R[track]<-Radius
  }
  return(round(max(R)+1))
}
trnfixer<- function(Genelist){
  for (i in 1:length(Genelist[,1])){
    
    if (Genelist[i,1] %in% c("trnK-UUU", "tRNA-Lys")){
      Genelist[i,1]<- "trnK"
    }
    
    if (Genelist[i,1] %in% c("trnN-UUA", "tRNA-Asn")){
      Genelist[i,1]<- "trnN"
    }
    
    if (Genelist[i,1] %in% c("trnT-UGA", "tRNA-Thr")){
      Genelist[i,1]<- "trnT"
    }
    
    if (Genelist[i,1] %in% c("trnR-UCU", "tRNA-Arg")){
      Genelist[i,1]<- "trnR"
    }
    
    if (Genelist[i,1] %in% c("trnM-UAC", "tRNA-Met")){
      Genelist[i,1]<- "trnM"
    }
    
    if (Genelist[i,1] %in% c("trnI-UAA", "tRNA-Ile")){
      Genelist[i,1]<- "trnI"
    }
    
    if (Genelist[i,1] %in% c("trnQ-GUU", "tRNA-Gln")){
      Genelist[i,1]<- "trnQ"
    }
    
    if (Genelist[i,1] %in% c("trnH-GUA", "tRNA-His", "trnH-GUG")){
      Genelist[i,1]<- "trnH"
    }
    
    if (Genelist[i,1] %in% c("trnP-GGG", "tRNA-Pro")){
      Genelist[i,1]<- "trnP"
    }
    
    if (Genelist[i,1] %in% c("trnE-CUC", "tRNA-Glu")){
      Genelist[i,1]<- "trnE"
    }
    
    if (Genelist[i,1] %in% c("trnD-CUA", "tRNA-Asp")){
      Genelist[i,1]<- "trnD"
    }
    
    if (Genelist[i,1] %in% c("trnA-CGC", "tRNA-Ala")){
      Genelist[i,1]<- "trnA"
    }
    
    if (Genelist[i,1] %in% c("trnG-CCC", "tRNA-Gly")){
      Genelist[i,1]<- "trnG"
    }
    
    if (Genelist[i,1] %in% c("trnV-CAC", "tRNA-Val")){
      Genelist[i,1]<- "trnV"
    }
    
    if (Genelist[i,1] %in% c("trnY-AUA", "tRNA-Tyr")){
      Genelist[i,1]<- "trnY"
    }
    
    if (Genelist[i,1] %in% c("trnS-AGA", "tRNA-Ser")){
      Genelist[i,1]<- "trnS"
    }
    
    if (Genelist[i,1] %in% c("trnW-ACC", "tRNA-Trp")){
      Genelist[i,1]<- "trnW"
    }
    
    if (Genelist[i,1] %in% c("trnC-ACA", "tRNA-Cys")){
      Genelist[i,1]<- "trnC"
    }
    
    if (Genelist[i,1] %in% c("trnL-AAU", "tRNA-Leu")){
      Genelist[i,1]<- "trnL"
    }
    
    if (Genelist[i,1] %in% c("trnF-AAA", "tRNA-Phe")){
      Genelist[i,1]<- "trnF"
    }
  }
  return(Genelist)
}
trnCut<- function(GL){
  name<- GL[,1]
  for (i in 1:length(name)){
    if (length(strsplit(as.character(name[i]), split="")[[1]])>5){
      name[i]<-paste(strsplit(as.character(name[i]), split="")[[1]][1:4], collapse = "")
    }
  }
  GL[,1]<- as.character(name)
  return(GL)
}#second version for dogma input and fix and output
chr.count<- function(word){ 
  if (length(word)==1){
    return(length(strsplit(word, "")[[1]]))
  }
  else {
    t<- numeric(length(word))
    for (i in 1:length(word)){
      t[i]<- length(strsplit(word[i], "")[[1]])
    }
    return(t)
  }
}
LSC<- function(IRinfo){#a vecotr as a list
  c<-as.vector(IRinfo)
  IR<- c[3]
  SSC<- c[2]-(c[1]+IR)
  return(c[4]-(SSC+2*IR))
}
SSC<- function(IRinfo){
  c<-as.vector(IRinfo)
  return(c[2]-(c[1]+c[3]))
}
gblister<- function(gbFiles){
  gbfiles<- list()
  j<-0
  for (i in 1:length(gbFiles)){
    if (!is.null(gbFiles[[i]])){
      j<-j+1
      gbfiles[[j]]<- gbFiles[[i]]
    }
  }
  return(gbfiles)
}
IRs<- function(gbfiles, Sfiles, file="IR.jpg"){
  l<<- length(gbfiles)#STOP with the proper length condition to be added
  FasList<<-  list()
  IRList<<-   list()
  GeneList<<- list()
  spnames<<-  list()
  nuclw<<- numeric(l)
  for (i in 1:l){
    gb<-gbfiles[[i]]
    rev<- Sfiles[[i]]
    #form one IRListDimp
    FasList[[i]]<<- rdnFixer(gb)
    IRList[[i]]<<-  IRinfo(as.vector(FasList[[i]]))
    GeneList[[i]]<<-gene.cordinates(gb)
    if (rev){GeneList[[i]]<<- SSCrev(GeneList[[i]], (IRList[[i]][1]+IRList[[i]][3]), IRList[[i]][2])}
    spnames[[i]]<<- sp.name(gb)
    nuclw[i]<<- 100/length(FasList[[1]]) 
  }
  for (i in 1:l){
    GeneList[[i]]<<- trnfixer(GeneList[[i]])
  }
  for (i in 1:l){
    GeneList[[i]]<<- trnCut(GeneList[[i]])
  }
}
IRs2<- function(gbfiles, file="IR.jpg"){
  names(FasList)<- unlist(spnames)
  names(IRList)<- unlist(spnames)
  names(GeneList)<- unlist(spnames)
  #copied and paste the following line to make sure produce the graph
  width<-8.3
  height<- l
  dev.new(width, height)
  m<-matrix(rep(c(15, 25, 30, 25, 10), l+2), 5, l+2)
  m[,l+2]<-rep(NA,5)
  m[,1]<-rep(0,5)
  jpeg(file, width=8.3, height=(l+2)*8.3/12, units="in", res=300)#I changed the track to 12 from 11#alternative togglier for height can be l*.83 /// or l ///(l*7.47/10+0.83) ///(l+1)*8.3/11
  par(mai=c(0.5, 0.6+max(chr.count(unlist(spnames)))/20, 0.7, 0.5))# c(bottom, left, top, right) 
  barplot(m, horiz=T, lwd=1,cex.names=0.46, space = 4, border = T, axes = F, col=c("lightblue", "orange1", "lightgreen", "orange1", "lightblue"))
  title(main="Inverted Repeats")
  segments(15, 5*l+8, 15, 7, lty=1, lwd=1, col = "darkgrey")
  segments(40, 5*l+8, 40, 7, lty=1, lwd=1, col = "darkgrey")
  segments(70, 5*l+8, 70, 7, lty=1, lwd=1, col = "darkgrey")
  segments(95, 5*l+8, 95, 7, lty=1, lwd=1, col = "darkgrey")
  text(15, 5*l+9, "JLB", font=4, cex=0.7)
  text(40, 5*l+9, "JSB", font=4, cex=0.7)
  text(70, 5*l+9, "JSA", font=4, cex=0.7)
  text(95, 5*l+9, "JLA", font=4, cex=0.7)
  
  intextcol<- "white"
  innocol<- "blue"
  for (i in seq(9.5, 5*l+9, 5)){
    text(101, i, "LSC", cex=0.57, font=2, col=intextcol)
    text(2.5, i, "LSC", cex=0.57, font=2, col=intextcol)
    text(18, i, "IRb", cex=0.57, font=2, col=intextcol)
    text(43, i, "SSC", cex=0.57, font=2, col=intextcol)
    text(73, i, "IRa", cex=0.57, font=2, col=intextcol)
    points(27.5, i, cex=2.3, pch=18, col="white")
    points(55, i, cex=2.3, pch=18, col="white")
    points(82.5, i, cex=2.3, pch=18, col="white")
    text(27.5, i, "//", cex=0.95, font=1)
    text(55, i, "//", cex=0.95, font=1)
    text(82.5, i, "//", cex=0.95, font=1)
    count<-which(seq(9.5, 5*l+9, 5)==i)
    LSCp<-paste(prettyNum(LSC(IRList[[count]]), big.mark = ","), "bp", "")
    SSCp<-paste(prettyNum(SSC(IRList[[count]]), big.mark = ","), "bp", "")
    IRp<- paste(prettyNum(IRList[[count]][3], big.mark = ","), "bp", "")
    text(10.5, i-0.09, LSCp, font=4, cex=0.52, col=innocol)
    text(64.5, i-0.09, SSCp, font=4, cex=0.52, col=innocol)
    text(35.5, i-0.09, IRp, font=4, cex=0.52, col=innocol)
    text(90, i-0.09, IRp, font=4, cex=0.52, col=innocol)
    axis(2, i-1.5, labels = paste(prettyNum(IRList[[count]][4], big.mark=","), "bp", ""), font = 4, cex.axis=0.7, las=2, tick=F)
  }
  axis(2, seq(9.5, 5*l+9, 5), labels = spnames, font = 4, cex.axis=0.7, las=2, tick=F)
  points(0, 4.5, cex=2.3, pch=18, col="white")
  
  I<-   Max.Radius(1, l, genelist = GeneList, IRlist = IRList)
  II<-  Max.Radius(2, l, genelist = GeneList, IRlist = IRList)
  III<- Max.Radius(3, l, genelist = GeneList, IRlist = IRList)
  IV<-  Max.Radius(4, l, genelist = GeneList, IRlist = IRList)
  
  for (i in 1:l){JG.plotter(I, 1, i)}
  for (i in 1:l){GN.plotter(I, 1, i)}
  for (i in 1:l){JG.plotter(II, 2, i)}#default has been 100
  for (i in 1:l){GN.plotter(II, 2, i)}#
  for (i in 1:l){JG.plotter(III, 3, i)}
  for (i in 1:l){GN.plotter(III, 3, i)}
  for (i in 1:l){JG.plotter(IV, 4, i)}
  for (i in 1:l){GN.plotter(IV, 4, i)}
  for (i in 1:l){OJ.plotter(I, 1, i)}
  for (i in 1:l){OJ.plotter(III, 3, i)}
  for (i in 1:l){OJ.plotter(II, 2, i)}
  for (i in 1:l){OJ.plotter(IV, 4, i)} 
  for (i in 1:l){JD.plotter(I, 1, i)}
  for (i in 1:l){JD.plotter(III, 3, i)}
  for (i in 1:l){JD.plotter(II, 2, i)}
  for (i in 1:l){JD.plotter(IV, 4, i)}
  #copied here
  width<-8.3
  height<- l
  dev.off()
  dev.new(width, height)
  m<-matrix(rep(c(15, 25, 30, 25, 10), l+2), 5, l+2)
  m[,l+2]<-rep(NA,5)
  m[,1]<-rep(0,5)
  jpeg(file, width=8.3, height=(l+2)*8.3/12, units="in", res=300)#I changed the track to 12 from 11#alternative togglier for height can be l*.83 /// or l ///(l*7.47/10+0.83) ///(l+1)*8.3/11
  par(mai=c(0.5, 0.5+max(chr.count(unlist(spnames)))/20, 0.7, 0.5))# c(bottom, left, top, right) 
  barplot(m, horiz=T, lwd=1,cex.names=0.46, space = 4, border = T, axes = F, col=c("lightblue", "orange1", "lightgreen", "orange1", "lightblue"))
  title(main="Inverted Repeats")
  segments(15, 5*l+8, 15, 7, lty=1, lwd=1, col = "darkgrey")
  segments(40, 5*l+8, 40, 7, lty=1, lwd=1, col = "darkgrey")
  segments(70, 5*l+8, 70, 7, lty=1, lwd=1, col = "darkgrey")
  segments(95, 5*l+8, 95, 7, lty=1, lwd=1, col = "darkgrey")
  text(15, 5*l+9, "JLB", font=4, cex=0.7)
  text(40, 5*l+9, "JSB", font=4, cex=0.7)
  text(70, 5*l+9, "JSA", font=4, cex=0.7)
  text(95, 5*l+9, "JLA", font=4, cex=0.7)
  intextcol<- "white"
  innocol<- "blue"
  for (i in seq(9.5, 5*l+9, 5)){
    text(101, i, "LSC", cex=0.57, font=2, col=intextcol)
    text(2.5, i, "LSC", cex=0.57, font=2, col=intextcol)
    text(18, i, "IRb", cex=0.57, font=2, col=intextcol)
    text(43, i, "SSC", cex=0.57, font=2, col=intextcol)
    text(73, i, "IRa", cex=0.57, font=2, col=intextcol)
    points(27.5, i, cex=2.3, pch=18, col="white")
    points(55, i, cex=2.3, pch=18, col="white")
    points(82.5, i, cex=2.3, pch=18, col="white")
    text(27.5, i, "//", cex=0.95, font=1)
    text(55, i, "//", cex=0.95, font=1)
    text(82.5, i, "//", cex=0.95, font=1)
    count<-which(seq(9.5, 5*l+9, 5)==i)
    LSCp<-paste(prettyNum(LSC(IRList[[count]]), big.mark = ","), "bp", "")
    SSCp<-paste(prettyNum(SSC(IRList[[count]]), big.mark = ","), "bp", "")
    IRp<- paste(prettyNum(IRList[[count]][3], big.mark = ","), "bp", "")
    text(10.5, i-0.09, LSCp, font=4, cex=0.52, col=innocol)
    text(64.5, i-0.09, SSCp, font=4, cex=0.52, col=innocol)
    text(35.5, i-0.09, IRp, font=4, cex=0.52, col=innocol)
    text(90, i-0.09, IRp, font=4, cex=0.52, col=innocol)
    axis(2, i-1.5, labels = paste(prettyNum(IRList[[count]][4], big.mark=","), "bp", ""), font = 4, cex.axis=0.7, las=2, tick=F)
  }
  axis(2, seq(9.5, 5*l+9, 5), labels = spnames, font = 4, cex.axis=0.7, las=2, tick=F)
  points(0, 4.5, cex=2.3, pch=18, col="white")
  
  I<-   Max.Radius(1, l, genelist = GeneList, IRlist = IRList)
  II<-  Max.Radius(2, l, genelist = GeneList, IRlist = IRList)
  III<- Max.Radius(3, l, genelist = GeneList, IRlist = IRList)
  IV<-  Max.Radius(4, l, genelist = GeneList, IRlist = IRList)
  
  
  for (i in 1:l){JG.plotter(I, 1, i)}
  for (i in 1:l){GN.plotter(I, 1, i)}
  for (i in 1:l){JG.plotter(II, 2, i)}#default has been 100
  for (i in 1:l){GN.plotter(II, 2, i)}#
  for (i in 1:l){JG.plotter(III, 3, i)}
  for (i in 1:l){GN.plotter(III, 3, i)}
  for (i in 1:l){JG.plotter(IV, 4, i)}
  for (i in 1:l){GN.plotter(IV, 4, i)}
  for (i in 1:l){OJ.plotter(I, 1, i)}
  for (i in 1:l){OJ.plotter(III, 3, i)}
  for (i in 1:l){OJ.plotter(II, 2, i)}
  for (i in 1:l){OJ.plotter(IV, 4, i)} 
  for (i in 1:l){JD.plotter(I, 1, i)}
  for (i in 1:l){JD.plotter(III, 3, i)}
  for (i in 1:l){JD.plotter(II, 2, i)}
  for (i in 1:l){JD.plotter(IV, 4, i)}
  dev.off()
  dev.off()
  dev.off()
}##For the GB files
###THis section contains the functions related to the dogma input
rdnFixerD<- function(FasGenome){
  seq<- FasGenome
  seq[which(seq=="u")]<-sample(c("t"), length(which(seq=="u")), TRUE)
  seq[which(seq=="r")]<-sample(c("a", "g"), length(which(seq=="r")), TRUE)
  seq[which(seq=="y")]<-sample(c("c", "t"), length(which(seq=="y")), TRUE)
  seq[which(seq=="s")]<-sample(c("c", "g"), length(which(seq=="s")), TRUE)
  seq[which(seq=="w")]<-sample(c("a", "t"), length(which(seq=="w")), TRUE)
  seq[which(seq=="k")]<-sample(c("g", "t"), length(which(seq=="k")), TRUE)
  seq[which(seq=="m")]<-sample(c("c", "a"), length(which(seq=="m")), TRUE)
  seq[which(seq=="b")]<-sample(c("c", "g", "t"), length(which(seq=="b")), TRUE)
  seq[which(seq=="d")]<-sample(c("a", "g", "t"), length(which(seq=="d")), TRUE)
  seq[which(seq=="h")]<-sample(c("c", "a", "t"), length(which(seq=="h")), TRUE)
  seq[which(seq=="v")]<-sample(c("c", "a", "g"), length(which(seq=="v")), TRUE)
  seq[which(seq=="n")]<-sample(c("c", "g", "t", "a"), length(which(seq=="n")), TRUE)
  seq[which(seq=="-")]<-sample(c("c", "g", "t", "a"), length(which(seq=="-")), TRUE)
  return(seq)
}##for the Dogma and fasta
sp.nameD<- function(gb){
  paste(strsplit(gb[2], " ")[[1]][3], strsplit(gb[2], " ")[[1]][4], sep=" ")
}##for the dogma input, read the name from the header of fasta files
trnfixerD<- function(Genelist){
  Genelist<- Genelist[,c(3,2,1,4)]
  for (i in 1:length(Genelist[,1])){
    
    if (Genelist[i,1] %in% c("trnK-UUU", "tRNA-Lys")){
      Genelist[i,1]<- "trnK"
    }
    
    if (Genelist[i,1] %in% c("trnN-UUA", "tRNA-Asn")){
      Genelist[i,1]<- "trnN"
    }
    
    if (Genelist[i,1] %in% c("trnT-UGA", "tRNA-Thr")){
      Genelist[i,1]<- "trnT"
    }
    
    if (Genelist[i,1] %in% c("trnR-UCU", "tRNA-Arg")){
      Genelist[i,1]<- "trnR"
    }
    
    if (Genelist[i,1] %in% c("trnM-UAC", "tRNA-Met")){
      Genelist[i,1]<- "trnM"
    }
    
    if (Genelist[i,1] %in% c("trnI-UAA", "tRNA-Ile")){
      Genelist[i,1]<- "trnI"
    }
    
    if (Genelist[i,1] %in% c("trnQ-GUU", "tRNA-Gln")){
      Genelist[i,1]<- "trnQ"
    }
    
    if (Genelist[i,1] %in% c("trnH-GUA", "tRNA-His")){
      Genelist[i,1]<- "trnH"
    }
    
    if (Genelist[i,1] %in% c("trnP-GGG", "tRNA-Pro")){
      Genelist[i,1]<- "trnP"
    }
    
    if (Genelist[i,1] %in% c("trnE-CUC", "tRNA-Glu")){
      Genelist[i,1]<- "trnE"
    }
    
    if (Genelist[i,1] %in% c("trnD-CUA", "tRNA-Asp")){
      Genelist[i,1]<- "trnD"
    }
    
    if (Genelist[i,1] %in% c("trnA-CGC", "tRNA-Ala")){
      Genelist[i,1]<- "trnA"
    }
    
    if (Genelist[i,1] %in% c("trnG-CCC", "tRNA-Gly")){
      Genelist[i,1]<- "trnG"
    }
    
    if (Genelist[i,1] %in% c("trnV-CAC", "tRNA-Val")){
      Genelist[i,1]<- "trnV"
    }
    
    if (Genelist[i,1] %in% c("trnY-AUA", "tRNA-Tyr")){
      Genelist[i,1]<- "trnY"
    }
    
    if (Genelist[i,1] %in% c("trnS-AGA", "tRNA-Ser")){
      Genelist[i,1]<- "trnS"
    }
    
    if (Genelist[i,1] %in% c("trnW-ACC", "tRNA-Trp")){
      Genelist[i,1]<- "trnW"
    }
    
    if (Genelist[i,1] %in% c("trnC-ACA", "tRNA-Cys")){
      Genelist[i,1]<- "trnC"
    }
    
    if (Genelist[i,1] %in% c("trnL-AAU", "tRNA-Leu")){
      Genelist[i,1]<- "trnL"
    }
    
    if (Genelist[i,1] %in% c("trnF-AAA", "tRNA-Phe")){
      Genelist[i,1]<- "trnF"
    }
  }
  return(Genelist)
}##For the dogma and fasta
trnDogma<- function(dogma){##takes in the dogma input and return the same output but with the gene names fixed
  name<- as.character(dogma[,3])
  for (i in 1:length(name)){
    if (length(strsplit(as.character(name[i]), split="")[[1]])>5){
      name[i]<-paste(strsplit(as.character(name[i]), split="")[[1]][1:4], collapse = "")
    }
  }
  dogma[,3]<- as.character(name)
  return(dogma)
}#second version for dogma input and fix and output
GnlBuilder<- function(dogma){##takes the dogma input and return the type of the list suited for the IR and other functions
  l<- length(dogma[,1])
  list<- cbind(as.character(rep(0, l)),as.character(rep(0, l)),as.character(rep(0, l)),as.character(rep(0, l)),as.character(rep(0, l)))
  for (i in 1:l){
    if (as.character(dogma[i,4])=="-"){
      p<-dogma[i, 2]
      dogma[i,2]<- dogma[i,1]
      dogma[i,1]<- p
    }
    list[i,1]<-as.character(dogma[i,3])
    list[i,2]<-as.character(dogma[i,1])
    list[i,3]<-as.character(dogma[i,2])
  }
  return(list)
}
###Suppose that you are given two vectors of IRinp and one value of the genome size.
###IRinp going to be the vector with input cordinates of the JLB, JSB, JSA, JLA, and genome size
#for its first, second, third, fourth, and fifth values respectively 
JunctRadiusFinderDinp<- function(gene.cordinates, IRinp, J.pos, radius, silence=TRUE){#function to find the genes near the junction cordinate based on the output of the gene.cordinate function
  if(J.pos==1){
    J<- IRinp[1]
  }
  else if(J.pos==2){
    J<- IRinp[2]
  }
  else if(J.pos==3){
    J<- IRinp[3]
  }
  else if(J.pos==4){
    J<- IRinp[4]
  }
  else {stop("J.pos missing or out bound. It should be either 1,2,3, or 4 for JLB, JSB, JSA, and JLA, recpectively ")}
  g<-gene.cordinates
  l<- ncol(g)
  r<- radius
  gs<- IRinp[5]
  t1<- subset(g, J-r<as.numeric(g[,2]) & as.numeric(g[,2])<J+r)
  t2<- subset(g, J-r<as.numeric(g[,3]) & as.numeric(g[,3])<J+r)
  t3<- subset(g, J-r<as.numeric(g[,4]) & as.numeric(g[,4])<J+r)
  t4<- subset(g, J-r<as.numeric(g[,5]) & as.numeric(g[,5])<J+r)
  for (i in 1:nrow(g)){
    if(sum(g[i, ]=="0")==2){
      g[i,4:5]<-c(NA, NA)
    }
  }
  t5<- subset(g, J-r<(as.numeric(g[,2])+gs) & (as.numeric(g[,2])+gs)<J+r)
  t6<- subset(g, J-r<(as.numeric(g[,3])+gs) & (as.numeric(g[,3])+gs)<J+r)
  t7<- subset(g, J-r<(as.numeric(g[,4])+gs) & (as.numeric(g[,4])+gs)<J+r)
  t8<- subset(g, J-r<(as.numeric(g[,5])+gs) & (as.numeric(g[,5])+gs)<J+r)
  t<-unique(do.call("rbind", list(t1, t2, t3, t4, t5, t6, t7, t8)))
  if((length(t)+silence)==0){
    warning("There is no gene on the specified junction with the given radius, increasing the radius might help")
  }
  return(t)
}
JunctNoFinderDinp<- function(gene.cordinates, IRinp, J.pos ,Number){#the same function like it precusor but returns the number of genes found near the junction site
  if(J.pos==1){
    J<- IRinp[1]
  }
  else if(J.pos==2){
    J<- IRinp[2]
  }
  else if(J.pos==3){
    J<- IRinp[3]
  }
  else if(J.pos==4){
    J<- IRinp[4]
  }
  else {stop("J.pos missing or out bound. It should be either 1,2,3, or 4 for JLB, JSB, JSA, and JLA, recpectively ")}
  g<-gene.cordinates
  l<- ncol(g)
  n<-Number
  counter=0
  count=0
  while (nrow(JunctRadiusFinderDinp(g, IRinp, J.pos ,count)) <= n-1){
    count<-count+1
  }
  t<-JunctRadiusFinderDinp(g, IRinp, J.pos ,count)
  t[is.na(t)]<- "0"
  tup<- matrix(0, n, 3)
  for (i in 1:n){
    dist<-abs(as.numeric(t[i,][2:5])-J)
    if (which(dist==min(dist))>3){
      tup[i,]<-t[i, c(1,4,5)] 
    }
    else {
      tup[i,]<-t[i, c(1,2,3)]
    }
  }
  return(tup)
}
thFinderDinp<- function(gene.cordinates, IRinp, J.pos, n){#finding the nth closest gene to the junctions site
  g<-gene.cordinates
  l<- ncol(g)
  if(n==1){
    return(JunctNoFinderDinp(g, IRinp, J.pos, n))
  }
  else{
    previous<-JunctNoFinderDinp(g, IRinp, J.pos, n-1)
    current<- JunctNoFinderDinp(g, IRinp, J.pos, n  )
    return(current[!current[,2] %in% previous[,2],])
  }
}
RadiusNoFinderDinp<- function(gene.cordinates, IRinp, J.pos, Number){
  g<-gene.cordinates
  l<- ncol(g)
  n<-Number
  counter=0
  count=0
  while (nrow(JunctRadiusFinderDinp(g, IRinp, J.pos, count)) <= n-1){
    count<-count+1
  }
  return(count)
}
JG.plotterDinp<- function(Radius, J.pos, track){#JunctionGene plotter given an global IRListDinp
  #' Junction site gene plotter
  #' 
  #' Plotting the genes in the vicinity of the junction site of the chloroplast. 
  #' @param gene.cordinate The output of the gene.cordinate function with the name of the genes and their cordinates on genome
  #' @param junction.cordinate The cordinate of the junction site for which the genes are to be plotted
  #' @param n The number of the genes that is needed to be plotted
  #' @param J.pos The position of the J site: JLB, JSB, JSA, and JLA for 1,2,3, and 4 respectively
  #' @param margin The margin value to be added to the interval as the proportion of the minimum lenght of the nth gene in vicinity
  #' @param track The track on which the genes are to be plotted, strating from the bottom to up as integers 1,2,...
  #' @export
  if(J.pos==1){
    pc<-15
    J<- IRListDinp[[track]][1]
  }
  else if(J.pos==2){
    pc<-40
    J<- IRListDinp[[track]][2]
  }
  else if(J.pos==3){
    pc<-70
    J<- IRListDinp[[track]][3]
  }
  else if(J.pos==4){
    pc<-95
    J<- IRListDinp[[track]][4]
  }
  else {stop("J.pos missing or out bound. It should be either 1,2,3, or 4 for JLB, JSB, JSA, and JLA, recpectively ")}
  t<- JunctRadiusFinderDinp(GeneList[[track]], IRListDinp[[track]], J.pos , Radius)    
  t[is.na(t)]<- "0"
  n<- length(t[,1])
  tup<- matrix(0, n, 3)
  for (i in 1:n){
    dist<-abs(as.numeric(t[i,][2:5])-J)
    if (which(dist==min(dist))>3){
      tup[i,]<-t[i, c(1,4,5)] 
    }
    else {
      tup[i,]<-t[i, c(1,2,3)]
    }
  }               
  bw<- 10/Radius
  Rcord<- matrix(0, n, 2)
  Rcord[,1]<-as.numeric(tup[,2])
  Rcord[,2]<-as.numeric(tup[,3])
  Pcord<- (Rcord-J)*bw+pc
  gcol<- function(tup){
    l<-length(tup[,1])
    col<-numeric(l)
    for (i in 1:l){
      if(tup[i,1]=="trnH"){col[i]<- "burlywood4"}
      else if(tup[i,1]=="trnN"){col[i]<- "violet"}
      else if(tup[i,1]=="psbA"){col[i]<- "purple2"}
      else if(tup[i,1]=="rps19"){col[i]<- "firebrick1"}
      else if(tup[i,1]=="ycf1"){col[i]<- "dodgerblue3"}
      else if(tup[i,1]=="ndhF"){col[i]<- "darkred"}
      else if(tup[i,1]=="rpl22"){col[i]<- "navy"}
      else if(tup[i,1]=="rpl2"){col[i]<- "forestgreen"}
      else if(tup[i,1]=="rpl23"){col[i]<- "blue3"}
      else {col[i]<- "black"}
    }
    return(col)
  }
  if (J.pos==4){
    Pcord[Pcord < 0]<- ((Rcord[which(Pcord< 0)] + IRListDinp[[track]][5])-J)*bw+pc
  }
  for (i in 1:n){
    if(Rcord[i,1]>Rcord[i,2]){
      x1<-min(Pcord[i,1], pc+10)
      x2<-max(Pcord[i,2], pc-10)
      for (j in seq(0.10, 0.70, 0.05)){
        segments(x1, track*5+j+5, x2, track*5+j+5, lwd=1, col=paste(gcol(tup)[i]))
      }
    }
    else {
      for (j in seq(1.1, 1.7, 0.05)){
        segments(max(Pcord[i,1], pc-10), track*5-j+5, min(Pcord[i,2], pc+10), track*5-j+5, lwd=1, col=paste(gcol(tup)[i]))
      }
    }
  }
}
GN.plotterDinp<- function(Radius, J.pos, track){#GeneName plotter
  #' Gene Name plotter
  #' 
  #' Plotting the gene names on a given gene which is already plotted on the tracks of the IR plot
  if(J.pos==1){
    pc<-15
    J<- IRListDinp[[track]][1]
  }
  else if(J.pos==2){
    pc<-40
    J<- IRListDinp[[track]][2]
  }
  else if(J.pos==3){
    pc<-70
    J<- IRListDinp[[track]][3]
  }
  else if(J.pos==4){
    pc<-95
    J<- IRListDinp[[track]][4]
  }
  else {stop("J.pos missing or out bound. It should be either 1,2,3, or 4 for JLB, JSB, JSA, and JLA, recpectively ")}
  txtincol<- "white"
  txtoutcol="black"
  txtfont<- 4
  numfont<- 1
  numcex<- 0.44
  txtcex<- 0.46
  t<- JunctRadiusFinderDinp(GeneList[[track]], IRListDinp[[track]], J.pos , Radius)
  t[is.na(t)]<- "0"
  n<- length(t[,1])
  tup<- matrix(0, n, 3)
  for (i in 1:n){
    dist<-abs(as.numeric(t[i,][2:5])-J)
    if (which(dist==min(dist))>3){
      tup[i,]<-t[i, c(1,4,5)] 
    }
    else {
      tup[i,]<-t[i, c(1,2,3)]
    }
  }               
  bw<- 10/Radius
  Rcord<- matrix(0, n, 2)
  Rcord[,1]<-as.numeric(tup[,2])
  Rcord[,2]<-as.numeric(tup[,3])
  Pcord<- (Rcord-J)*bw+pc
  if (J.pos==4){
    Pcord[Pcord < 0]<- ((Rcord[which(Pcord< 0)] + IRListDinp[[track]][5])-J)*bw+pc
  }
  for (i in 1:n){
    if(Rcord[i,1]>Rcord[i,2]){
      x1<-min(Pcord[i,1], pc+10)
      x2<-max(Pcord[i,2], pc-10)
      if (min(x1, x2) >= 104){
        text(103.5, track*5+1.05+5, tup[i,1], cex=txtcex, col=txtoutcol, font=txtfont)
      }
      else if (abs(x1-x2) >= 9.7){
        text(min(x1, x2)+1.8, track*5+0.35+5, tup[i,1], cex=txtcex, col=txtincol, font=txtfont)
        text(max(x1, x2)+1, track*5+0.3+5, paste(Rcord[i,1]-Rcord[i,2], "bp", " "), cex=numcex, col=txtincol, font=numfont, pos=2)
      }
      else if (abs(x1-x2) >= 3 & abs(x1-x2) < 9.7){
        text(((min(x1, x2)+max(x1, x2))/2), track*5+0.35+5, tup[i,1], cex=txtcex, col=txtincol, font=txtfont)
      }
      else if (abs(x1-x2) < 3){
        text(((min(x1, x2)+max(x1, x2))/2), track*5+1.15+5, tup[i,1], cex=txtcex, col=txtoutcol, font=txtfont)
      }
    }
    else {
      x1<-max(Pcord[i,1], pc-10)
      x2<-min(Pcord[i,2], pc+10)
      if (abs(x1-x2) >= 9.7){
        text(min(x1, x2)+1.8, track*5-1.40+5, tup[i,1], cex=txtcex, col=txtincol, font=txtfont)
        text(max(x1, x2)+1, track*5-1.45+5, paste(Rcord[i,2]-Rcord[i,1], "bp", " "), cex=numcex, col=txtincol, font=numfont, pos=2)
      }
      else if (abs(x1-x2) >= 3 & abs(x1-x2) < 9.7){
        text(((min(x1, x2)+max(x1, x2))/2), track*5-1.40+5, tup[i,1], cex=txtcex, col=txtincol, font=txtfont)
      }
      else if (abs(x1-x2) < 3){
        text(((min(x1, x2)+max(x1, x2))/2), track*5-2.25+5, tup[i,1], cex=txtcex, col=txtoutcol, font=txtfont)
      }
    }
  }
}
OJ.plotterDinp<- function(Radius, J.pos, track){#GeneName plotter
  #' On Junction plotter
  #' 
  #' Plotting the fine tuned narrow lines showing the limits of the genes which are passing through the junction sites with their bp
  if(J.pos==1){
    pc<-15
    J<- IRListDinp[[track]][1]
  }
  else if(J.pos==2){
    pc<-40
    J<- IRListDinp[[track]][2]
  }
  else if(J.pos==3){
    pc<-70
    J<- IRListDinp[[track]][3]
  }
  else if(J.pos==4){
    pc<-95
    J<- IRListDinp[[track]][4]
  }
  else {stop("J.pos missing or out bound. It should be either 1,2,3, or 4 for JLB, JSB, JSA, and JLA, recpectively ")}
  txtincol<- "white"
  txtoutcol="black"
  txtfont<- 4
  numfont<- 1
  numcex<- 0.44
  txtcex<- 0.46
  t<- JunctRadiusFinderDinp(GeneList[[track]], IRListDinp[[track]], J.pos , Radius)
  t[is.na(t)]<- "0"
  n<- length(t[,1])
  tup<- matrix(0, n, 3)
  for (i in 1:n){
    dist<-abs(as.numeric(t[i,][2:5])-J)
    if (which(dist==min(dist))>3){
      tup[i,]<-t[i, c(1,4,5)] 
    }
    else {
      tup[i,]<-t[i, c(1,2,3)]
    }
  }               
  bw<- 10/Radius
  Rcord<- matrix(0, n, 2)
  Rcord[,1]<-as.numeric(tup[,2])
  Rcord[,2]<-as.numeric(tup[,3])
  Pcord<- (Rcord-J)*bw+pc
  if (J.pos==4){
    Pcord[Pcord < 0]<- ((Rcord[which(Pcord< 0)] + IRListDinp[[track]][5])-J)*bw+pc
  }
  for (i in 1:n){
    if(Rcord[i,1]>Rcord[i,2]){
      x1<-min(Pcord[i,1], pc+10)
      x2<-max(Pcord[i,2], pc-10)
      if (Rcord[i,1] > J & Rcord[i,2] <J ){
        Arrows(min(x1, x2), track*5+1.1+5, pc-0.15, track*5+1.1+5, arr.type = "T", cex=0.5, arr.length = 0.12, lwd=0.5, arr.width = 0.4)
        Arrows(pc-0.15, track*5+1.1+5, min(x1, x2), track*5+1.1+5, arr.type = "T", cex=0.5, arr.length = 0.12, lwd=0.5, arr.width = 0.4)
        Arrows(pc+0.15, track*5+1.1+5, max(x1, x2), track*5+1.1+5, arr.type = "T", cex=0.5, arr.length = 0.12, lwd=0.5, arr.width = 0.4)
        Arrows(max(x1, x2), track*5+1.1+5, pc+0.15, track*5+1.1+5, arr.type = "T", cex=0.5, arr.length = 0.12, lwd=0.5, arr.width = 0.4)
        text(pc-1.4, track*5+1.4+5, paste(Rcord[i, 1]-J+1, "bp", sep=" "), cex=0.4, pos=4)#up-right
        text(pc+1.4, track*5+1.4+5, paste(J-Rcord[i, 2], "bp", sep=" "), cex=0.4, pos=2)#up-left
      }
    }
    else {
      x1<-max(Pcord[i,1], pc-10)
      x2<-min(Pcord[i,2], pc+10)
      if (Rcord[i,2] > J & Rcord[i,1] <J ){
        Arrows(max(x1, x2), track*5-2.1+5, pc+0.15, track*5-2.1+5, arr.type = "T", cex=0.5, arr.length = 0.12, lwd=0.5, arr.width = 0.4)
        Arrows(pc+0.15, track*5-2.1+5, max(x1, x2), track*5-2.1+5, arr.type = "T", cex=0.5, arr.length = 0.12, lwd=0.5, arr.width = 0.4)
        Arrows(pc-0.15, track*5-2.1+5, min(x1, x2), track*5-2.1+5, arr.type = "T", cex=0.5, arr.length = 0.12, lwd=0.5, arr.width = 0.4)
        Arrows(min(x1, x2), track*5-2.1+5, pc-0.15, track*5-2.1+5, arr.type = "T", cex=0.5, arr.length = 0.12, lwd=0.5, arr.width = 0.4)
        text(pc-1.4, track*5-2.6+5, paste(Rcord[i, 2]-J, "bp", sep=" "), cex=0.4, pos=4)#low-right
        text(pc+1.4, track*5-2.6+5, paste(J-Rcord[i, 1], "bp", sep=" "), cex=0.4, pos=2)#low-left
      }
    }
  }
}
JD.plotterDinp<- function(Radius, J.pos, track){#GeneName plotter
  #' Junction Distance plotter
  #' 
  #' plotting the narrow lines of the distance of the genes for the junction sites which are not passing through any gene and their bp
  if(J.pos==1){
    pc<-15
    J<- IRListDinp[[track]][1]
  }
  else if(J.pos==2){
    pc<-40
    J<- IRListDinp[[track]][2]
  }
  else if(J.pos==3){
    pc<-70
    J<- IRListDinp[[track]][3]
  }
  else if(J.pos==4){
    pc<-95
    J<- IRListDinp[[track]][4]
  }
  else {stop("J.pos missing or out bound. It should be either 1,2,3, or 4 for JLB, JSB, JSA, and JLA, recpectively ")}
  txtincol<- "white"
  txtoutcol="black"
  txtfont<- 4
  numfont<- 1
  numcex<- 0.44
  txtcex<- 0.46
  t<- JunctRadiusFinderDinp(GeneList[[track]], IRListDinp[[track]], J.pos , Radius)
  t[is.na(t)]<- "0"
  n<- length(t[,1])
  tup<- matrix(0, n, 3)
  for (i in 1:n){
    dist<-abs(as.numeric(t[i,][2:5])-J)
    if (which(dist==min(dist))>3){
      tup[i,]<-t[i, c(1,4,5)] 
    }
    else {
      tup[i,]<-t[i, c(1,2,3)]
    }
  }               
  bw<- 10/Radius
  Rcord<- matrix(0, n, 2)
  Rcord[,1]<-as.numeric(tup[,2])
  Rcord[,2]<-as.numeric(tup[,3])
  Pcord<- (Rcord-J)*bw+pc
  if (J.pos==4){
    ind<-Pcord < 0
    Pcord[Pcord < 0]<- ((Rcord[which(Pcord< 0)] + IRListDinp[[track]][5])-J)*bw+pc
    Rcord[ind]<- Rcord[ind] + IRListDinp[[track]][5]
  }
  counter<- 0
  for (i in 1:n){
    if(Rcord[i,1]>Rcord[i,2]){
      x1<-min(Pcord[i,1], pc+10)
      x2<-max(Pcord[i,2], pc-10)
      if (Rcord[i,1] > J & Rcord[i,2] <J ){
        counter<- counter+1
      }
    }
    else {
      x1<-max(Pcord[i,1], pc-10)
      x2<-min(Pcord[i,2], pc+10)
      if (Rcord[i,2] > J & Rcord[i,1] <J ){
        counter<- counter+1
      }
    }
  }
  if (counter==0){#find the closest gene to the junction site
    nearest<- which(abs(Pcord-pc)==min(abs(Pcord-pc)))
    col.cor<- floor((nearest-0.01)/n)+1
    row.cor<- nearest-n*floor((nearest-0.01)/n)###Now we have the row of the nearest gene
    #with this setting the position of the zero (genes tangant to the junction site) will not be plotted. If interested either put "=" for the middle condition of the left or right binary operator or better develop zero only handling if function
    if     (Rcord[row.cor, 1] > Rcord[row.cor, 2] & pc-Pcord[row.cor, col.cor] < 0 & min(Pcord[row.cor, 1], pc+10) - max(Pcord[row.cor, 2], pc-10) > 3 ){#top,right, big
      curvedarrow (from=c(pc+(Pcord[row.cor, col.cor]-pc)/2, track*5+0.3+5), to=c(pc+(Pcord[row.cor, col.cor]-pc)/2 + 3, track*5+1.3+5), curve = -0.21, lwd=0.6, arr.type = "curved", arr.col = "white", arr.length=0.08, arr.lwd=0.4, arr.pos=0.69, endhead=TRUE)
      text(pc+(Pcord[row.cor, col.cor]-pc)/2 + 4, track*5+1.3+0.2+5, paste(Rcord[row.cor, 2]-J, "bp", sep=" "), cex=0.4)
    }
    else if(Rcord[row.cor, 1] > Rcord[row.cor, 2] & pc-Pcord[row.cor, col.cor] < 0 & min(Pcord[row.cor, 1], pc+10) - max(Pcord[row.cor, 2], pc-10) <= 3 ){#top, right, small
      arrows(pc+(Pcord[row.cor, col.cor]-pc)/2, track*5+0.3+5 , pc+(Pcord[row.cor, col.cor]-pc)/2 -2 , track*5+1.3+5, angle = 15, length = 0.05, lwd=0.6)
      text(pc+(Pcord[row.cor, col.cor]-pc)/2 - 4.5 , track*5+1.3+0.2+5, paste(Rcord[row.cor, 2]-J, "bp", sep=" "), cex=0.4)
    }
    else if(Rcord[row.cor, 1] < Rcord[row.cor, 2] & pc-Pcord[row.cor, col.cor] < 0 & abs(max(Pcord[row.cor, 1], pc-10) - min(Pcord[row.cor, 2], pc+10)) > 3 ){#low, right, big
      curvedarrow (from=c(pc+(Pcord[row.cor, col.cor]-pc)/2, track*5-1.3+5), to=c(pc+(Pcord[row.cor, col.cor]-pc)/2 + 3, track*5-2.3+5), curve = 0.21, lwd=0.6, arr.type = "curved", arr.col = "white", arr.length=0.08, arr.lwd=0.4, arr.pos=0.69, endhead=TRUE) 
      text(pc+(Pcord[row.cor, col.cor]-pc)/2 + 4, track*5-2.3-0.2+5, paste(Rcord[row.cor, 1]-J, "bp", sep=" "), cex=0.4)
    }
    else if(Rcord[row.cor, 1] < Rcord[row.cor, 2] & pc-Pcord[row.cor, col.cor] < 0 & abs(max(Pcord[row.cor, 1], pc-10) - min(Pcord[row.cor, 2], pc+10)) <= 3 ){#low, right, small
      arrows(pc+(Pcord[row.cor, col.cor]-pc)/2, track*5-1.3+5 , pc+(Pcord[row.cor, col.cor]-pc)/2 -3 , track*5-2.3+5, angle = 15, length = 0.05, lwd=0.6)
      text(pc+(Pcord[row.cor, col.cor]-pc)/2 -4.5 , track*5-2.3-0.2+5, paste(Rcord[row.cor, 1]-J, "bp", sep=" "), cex=0.4)
    }
    else if(Rcord[row.cor, 1] < Rcord[row.cor, 2] & pc-Pcord[row.cor, col.cor] > 0 & abs(max(Pcord[row.cor, 1], pc-10) - min(Pcord[row.cor, 2], pc+10)) > 3 ){#low, left, big
      curvedarrow (from=c(pc+(Pcord[row.cor, col.cor]-pc)/2, track*5-1.3+5), to=c(pc+(Pcord[row.cor, col.cor]-pc)/2 - 3, track*5-2.3+5), curve = -0.21, lwd=0.6, arr.type = "curved", arr.col = "white", arr.length=0.08, arr.lwd=0.4, arr.pos=0.69, endhead=TRUE) 
      text(pc+(Pcord[row.cor, col.cor]-pc)/2 - 4, track*5-2.3-0.2+5, paste(J-Rcord[row.cor, 2], "bp", sep=" "), cex=0.4)
    }
    else if(Rcord[row.cor, 1] < Rcord[row.cor, 2] & pc-Pcord[row.cor, col.cor] > 0 & abs(max(Pcord[row.cor, 1], pc-10) - min(Pcord[row.cor, 2], pc+10)) <= 3 ){#low, left, small
      arrows(pc+(Pcord[row.cor, col.cor]-pc)/2, track*5-1.3+5 , pc+(Pcord[row.cor, col.cor]-pc)/2 + 3 , track*5-2.3+5, angle = 15, length = 0.05, lwd=0.6)
      text(pc+(Pcord[row.cor, col.cor]-pc)/2 + 4.5 , track*5-2.3-0.2+5, paste(J-Rcord[row.cor, 2], "bp", sep=" "), cex=0.4)
    }
    else if(Rcord[row.cor, 1] > Rcord[row.cor, 2] & pc-Pcord[row.cor, col.cor] > 0 & min(Pcord[row.cor, 1], pc+10) - max(Pcord[row.cor, 2], pc-10) > 3 ){#top, left, big
      curvedarrow (from=c(pc+(Pcord[row.cor, col.cor]-pc)/2, track*5+0.3+5), to=c(pc+(Pcord[row.cor, col.cor]-pc)/2 - 3, track*5+1.3+5), curve = 0.21, lwd=0.6, arr.type = "curved", arr.col = "white", arr.length=0.08, arr.lwd=0.4, arr.pos=0.69, endhead=TRUE)
      text(pc+(Pcord[row.cor, col.cor]-pc)/2 - 4.5 , track*5+1.3+0.2+5, paste(J-Rcord[row.cor, 1], "bp", sep=" "), cex=0.4)
    }
    else if(Rcord[row.cor, 1] > Rcord[row.cor, 2] & pc-Pcord[row.cor, col.cor] > 0 & min(Pcord[row.cor, 1], pc+10) - max(Pcord[row.cor, 2], pc-10) <= 3 ){#top, left, small
      arrows(pc+(Pcord[row.cor, col.cor]-pc)/2, track*5+0.3+5 , pc+(Pcord[row.cor, col.cor]-pc)/2 + 2 , track*5+1.3+5, angle = 15, length = 0.05, lwd=0.6)
      text(pc+(Pcord[row.cor, col.cor]-pc)/2 + 4 , track*5+1.3+0.2+5, paste(J-Rcord[row.cor, 1], "bp", sep=" "), cex=0.4)
    }
  }
}
Max.RadiusDinp<-function(J.pos, l, genelist, IRlistDinp){
  if(J.pos==1){
    Radius<-680#680#100
  }
  if(J.pos==2){
    Radius<-100#100#100
  }
  if(J.pos==3){
    Radius<-1000#1800#100
  }
  if(J.pos==4){
    Radius<-1000#1000#100
  }
  R<- numeric(l)
  for (track in 1:l){
    t<- JunctRadiusFinderDinp(GeneList[[track]], IRListDinp[[track]], J.pos , Radius)
    while(nrow(t)==0){
      Radius<- Radius+(0.2)*Radius
      t<- JunctRadiusFinderDinp(GeneList[[track]], IRListDinp[[track]], J.pos , Radius)
    }
    R[track]<-Radius
  }
  return(round(max(R)+1))
}
LSCDinp<- function(IRinp){#a vecotr as a list
  c<-as.vector(IRinp)
  return(c[5]-((c[2]-c[1])+(c[4]-c[3])+(c[3]-c[2])))
}
SSCDinp<- function(IRinp){
  c<-as.vector(IRinp)
  return(c[3]-c[2])
}
###everything is cleared for Dogma
IRsD<- function(dgfiles, fastafiles, irfiles, nfiles, file="IR.jpg"){
  l<<- length(dgfiles)#STOP with the proper length condition to be added
  FasList<<-  list()
  IRListDinp<<- list()
  GeneList<<- list()
  spnames<<-  list()
  nuclw<<- numeric(l)
  for (i in 1:l){
    dg<-dgfiles[[i]]
    FasList[[i]]<<- fastafiles[[i]]
    GeneList[[i]]<<-dgfiles[[i]]
    spnames[[i]]<<- nfiles[[i]]
    if (irfiles[[i]][1]==999){
      sss<- IRinfo(as.vector(FasList[[i]]))
      IRListDinp[[i]]<<-  c(sss[1],sss[1]+sss[3],sss[2],sss[2]+sss[3],sss[4])
    }
    else {IRListDinp[[i]]<<- irfiles[[i]]}
    nuclw[i]<<- 100/length(FasList[[1]]) 
  }
  for (i in 1:l){
    GeneList[[i]]<<- GeneList[[i]]
  }
}
IRsD2<- function(dgfiles, file="IR.jpg"){
  names(FasList)<- unlist(spnames)
  names(IRListDinp)<- unlist(spnames)
  names(GeneList)<- unlist(spnames)
  #copied and paste the following line to make sure produce the graph
  width<-8.3
  height<- l
  dev.new(width, height)
  m<-matrix(rep(c(15, 25, 30, 25, 10), l+2), 5, l+2)
  m[,l+2]<-rep(NA,5)
  m[,1]<-rep(0,5)
  jpeg(file, width=8.3, height=(l+2)*8.3/12, units="in", res=300)#I changed the track to 12 from 11#alternative togglier for height can be l*.83 /// or l ///(l*7.47/10+0.83) ///(l+1)*8.3/11
  par(mai=c(0.5, 0.6+max(chr.count(unlist(spnames)))/20, 0.7, 0.5))# c(bottom, left, top, right) 
  barplot(m, horiz=T, lwd=1,cex.names=0.46, space = 4, border = T, axes = F, col=c("lightblue", "orange1", "lightgreen", "orange1", "lightblue"))
  title(main="Inverted Repeats")
  segments(15, 5*l+8, 15, 7, lty=1, lwd=1, col = "darkgrey")
  segments(40, 5*l+8, 40, 7, lty=1, lwd=1, col = "darkgrey")
  segments(70, 5*l+8, 70, 7, lty=1, lwd=1, col = "darkgrey")
  segments(95, 5*l+8, 95, 7, lty=1, lwd=1, col = "darkgrey")
  text(15, 5*l+9, "JLB", font=4, cex=0.7)
  text(40, 5*l+9, "JSB", font=4, cex=0.7)
  text(70, 5*l+9, "JSA", font=4, cex=0.7)
  text(95, 5*l+9, "JLA", font=4, cex=0.7)
  
  intextcol<- "white"
  innocol<- "blue"
  for (i in seq(9.5, 5*l+9, 5)){
    text(101, i, "LSC", cex=0.57, font=2, col=intextcol)
    text(2.5, i, "LSC", cex=0.57, font=2, col=intextcol)
    text(18, i, "IRb", cex=0.57, font=2, col=intextcol)
    text(43, i, "SSC", cex=0.57, font=2, col=intextcol)
    text(73, i, "IRa", cex=0.57, font=2, col=intextcol)
    points(27.5, i, cex=2.3, pch=18, col="white")
    points(55, i, cex=2.3, pch=18, col="white")
    points(82.5, i, cex=2.3, pch=18, col="white")
    text(27.5, i, "//", cex=0.95, font=1)
    text(55, i, "//", cex=0.95, font=1)
    text(82.5, i, "//", cex=0.95, font=1)
    count<-which(seq(9.5, 5*l+9, 5)==i)
    LSCp<-paste(prettyNum(LSCDinp(IRListDinp[[count]]), big.mark = ","), "bp", "")
    SSCp<-paste(prettyNum(SSCDinp(IRListDinp[[count]]), big.mark = ","), "bp", "")
    IRpb<- paste(prettyNum((IRListDinp[[count]][2]-IRListDinp[[count]][1]), big.mark = ","), "bp", "")
    IRpa<- paste(prettyNum((IRListDinp[[count]][4]-IRListDinp[[count]][3]), big.mark = ","), "bp", "")
    text(10.5, i-0.09, LSCp, font=4, cex=0.52, col=innocol)
    text(64.5, i-0.09, SSCp, font=4, cex=0.52, col=innocol)
    text(35.5, i-0.09, IRpb, font=4, cex=0.52, col=innocol)
    text(90, i-0.09, IRpa, font=4, cex=0.52, col=innocol)
    axis(2, i-1.5, labels = paste(prettyNum(IRListDinp[[count]][5], big.mark=","), "bp", ""), font = 4, cex.axis=0.7, las=2, tick=F)
  }
  axis(2, seq(9.5, 5*l+9, 5), labels = spnames, font = 4, cex.axis=0.7, las=2, tick=F)
  points(0, 4.5, cex=2.3, pch=18, col="white")
  
  I<-   Max.RadiusDinp(1, l, genelist = GeneList, IRlistDinp = IRListDinp)
  II<-  Max.RadiusDinp(2, l, genelist = GeneList, IRlistDinp = IRListDinp)
  III<- Max.RadiusDinp(3, l, genelist = GeneList, IRlistDinp = IRListDinp)
  IV<-  Max.RadiusDinp(4, l, genelist = GeneList, IRlistDinp = IRListDinp)
  
  
  for (i in 1:l){JG.plotterDinp(I, 1, i)}
  for (i in 1:l){GN.plotterDinp(I, 1, i)}
  for (i in 1:l){JG.plotterDinp(II, 2, i)}#default has been 100
  for (i in 1:l){GN.plotterDinp(II, 2, i)}#
  for (i in 1:l){JG.plotterDinp(III, 3, i)}
  for (i in 1:l){GN.plotterDinp(III, 3, i)}
  for (i in 1:l){JG.plotterDinp(IV, 4, i)}
  for (i in 1:l){GN.plotterDinp(IV, 4, i)}
  for (i in 1:l){OJ.plotterDinp(I, 1, i)}
  for (i in 1:l){OJ.plotterDinp(III, 3, i)}
  for (i in 1:l){OJ.plotterDinp(II, 2, i)}
  for (i in 1:l){OJ.plotterDinp(IV, 4, i)} 
  for (i in 1:l){JD.plotterDinp(I, 1, i)}
  for (i in 1:l){JD.plotterDinp(III, 3, i)}
  for (i in 1:l){JD.plotterDinp(II, 2, i)}
  for (i in 1:l){JD.plotterDinp(IV, 4, i)}
  #copied here
  width<-8.3
  height<- l
  dev.new(width, height)
  m<-matrix(rep(c(15, 25, 30, 25, 10), l+2), 5, l+2)
  m[,l+2]<-rep(NA,5)
  m[,1]<-rep(0,5)
  jpeg(file, width=8.3, height=(l+2)*8.3/12, units="in", res=300)#I changed the track to 12 from 11#alternative togglier for height can be l*.83 /// or l ///(l*7.47/10+0.83) ///(l+1)*8.3/11
  par(mai=c(0.5, 0.6+max(chr.count(unlist(spnames)))/20, 0.7, 0.5))# c(bottom, left, top, right) 
  barplot(m, horiz=T, lwd=1,cex.names=0.46, space = 4, border = T, axes = F, col=c("lightblue", "orange1", "lightgreen", "orange1", "lightblue"))
  title(main="Inverted Repeats")
  segments(15, 5*l+8, 15, 7, lty=1, lwd=1, col = "darkgrey")
  segments(40, 5*l+8, 40, 7, lty=1, lwd=1, col = "darkgrey")
  segments(70, 5*l+8, 70, 7, lty=1, lwd=1, col = "darkgrey")
  segments(95, 5*l+8, 95, 7, lty=1, lwd=1, col = "darkgrey")
  text(15, 5*l+9, "JLB", font=4, cex=0.7)
  text(40, 5*l+9, "JSB", font=4, cex=0.7)
  text(70, 5*l+9, "JSA", font=4, cex=0.7)
  text(95, 5*l+9, "JLA", font=4, cex=0.7)
  
  intextcol<- "white"
  innocol<- "blue"
  for (i in seq(9.5, 5*l+9, 5)){
    text(101, i, "LSC", cex=0.57, font=2, col=intextcol)
    text(2.5, i, "LSC", cex=0.57, font=2, col=intextcol)
    text(18, i, "IRb", cex=0.57, font=2, col=intextcol)
    text(43, i, "SSC", cex=0.57, font=2, col=intextcol)
    text(73, i, "IRa", cex=0.57, font=2, col=intextcol)
    points(27.5, i, cex=2.3, pch=18, col="white")
    points(55, i, cex=2.3, pch=18, col="white")
    points(82.5, i, cex=2.3, pch=18, col="white")
    text(27.5, i, "//", cex=0.95, font=1)
    text(55, i, "//", cex=0.95, font=1)
    text(82.5, i, "//", cex=0.95, font=1)
    count<-which(seq(9.5, 5*l+9, 5)==i)
    LSCp<-paste(prettyNum(LSCDinp(IRListDinp[[count]]), big.mark = ","), "bp", "")
    SSCp<-paste(prettyNum(SSCDinp(IRListDinp[[count]]), big.mark = ","), "bp", "")
    IRpb<- paste(prettyNum((IRListDinp[[count]][2]-IRListDinp[[count]][1]), big.mark = ","), "bp", "")
    IRpa<- paste(prettyNum((IRListDinp[[count]][4]-IRListDinp[[count]][3]), big.mark = ","), "bp", "")
    text(10.5, i-0.09, LSCp, font=4, cex=0.52, col=innocol)
    text(64.5, i-0.09, SSCp, font=4, cex=0.52, col=innocol)
    text(35.5, i-0.09, IRpb, font=4, cex=0.52, col=innocol)
    text(90, i-0.09, IRpa, font=4, cex=0.52, col=innocol)
    axis(2, i-1.5, labels = paste(prettyNum(IRListDinp[[count]][5], big.mark=","), "bp", ""), font = 4, cex.axis=0.7, las=2, tick=F)
  }
  axis(2, seq(9.5, 5*l+9, 5), labels = spnames, font = 4, cex.axis=0.7, las=2, tick=F)
  points(0, 4.5, cex=2.3, pch=18, col="white")
  
  I<-   Max.RadiusDinp(1, l, genelist = GeneList, IRlistDinp = IRListDinp)
  II<-  Max.RadiusDinp(2, l, genelist = GeneList, IRlistDinp = IRListDinp)
  III<- Max.RadiusDinp(3, l, genelist = GeneList, IRlistDinp = IRListDinp)
  IV<-  Max.RadiusDinp(4, l, genelist = GeneList, IRlistDinp = IRListDinp)
  
  
  for (i in 1:l){JG.plotterDinp(I, 1, i)}
  for (i in 1:l){GN.plotterDinp(I, 1, i)}
  for (i in 1:l){JG.plotterDinp(II, 2, i)}#default has been 100
  for (i in 1:l){GN.plotterDinp(II, 2, i)}#
  for (i in 1:l){JG.plotterDinp(III, 3, i)}
  for (i in 1:l){GN.plotterDinp(III, 3, i)}
  for (i in 1:l){JG.plotterDinp(IV, 4, i)}
  for (i in 1:l){GN.plotterDinp(IV, 4, i)}
  for (i in 1:l){OJ.plotterDinp(I, 1, i)}
  for (i in 1:l){OJ.plotterDinp(III, 3, i)}
  for (i in 1:l){OJ.plotterDinp(II, 2, i)}
  for (i in 1:l){OJ.plotterDinp(IV, 4, i)} 
  for (i in 1:l){JD.plotterDinp(I, 1, i)}
  for (i in 1:l){JD.plotterDinp(III, 3, i)}
  for (i in 1:l){JD.plotterDinp(II, 2, i)}
  for (i in 1:l){JD.plotterDinp(IV, 4, i)}
  dev.off()
  dev.off()
  dev.off()
}##For the Fasta and Dogma Files



ui <- fluidPage(
  br(),
  HTML('<!DOCTYPE html>
       <html>
       <head>
       <style type="text/css">
       h1 {color:#267A43;}
       h4 {color:black;}
       h3 {color:navy;}
       p {color:navy;}
       </style>
       </head>
       <body>
       <div style="opacity:0.1;position:absolute;right:16px;left:16px;height:178px;border-radius:300px;background-color:#40B3DF"></div>
       
       <div style="padding:20px;border-radius:300px;border:20px solid #1E8449";">
       
       <div style="opacity:0.3;position:absolute;background-color:#8E44AD"></div>
       
       <h1 style="font-family:Copperplate Gothic Bold;letter-spacing:8px;text-align:center;">IRSCOPE</h1>
       
       <h4 style = "font-family:Charlesworth;text-align:center;color:navy;"><em> Tool for visualizig the junction sites of the chloroplast genome </em> </h4>        
       
       </div>
       </body>
       </html>
       
       '),
  
  # tags$h1(style = "font-family:Copperplate Gothic Bold", "IRSCOPE"),
  #tags$h4(style = "font-family:Charlesworth", tags$em("Tools for analyzing and visualizing the chloroplast genomes.")),
  #tags$hr(),
  tags$br(),
  
  #navlistPanel(widths = c(1, 11), well = FALSE,
              # tabPanel("Home",column(width = 6, h4(tags$em("Welcome to IRscope...")), br(), p(tags$strong("Chloroscope"),  "is a suit for analyzing and visualazing the chloroplast genomes. The tools can be found in the 'Tools' tab. Currently only the IRScope is a functional program of the website and the more information related to that can be found in the related pages."))),
               #tabPanel("Tools",
  fluidRow(column(width = 1, ""), column(width = 11, offset = 1, 
                        tabsetPanel(
                          tabPanel(strong("Home"), br(), 
                                   column( width=6, h4(tags$em("Welcome to IRscope...")), br(), tags$p(strong("IRscope"),"is a tool for  visualizing the genes on the boundaries of the junction sites of the chloroplast genome. The input files can be uploaded using the 'Files upload' tab. After the submission, the program will start searching the", tags$em("Inverted Repeat"), "regions. Note that the species selected will lead to a consensus radius for each junction site for the genes to be plotted. This indicates that if the input files are not representing closely related species, the program may lead to unsatisfactory result. Finally note that the input files can be either manual annotations, GeneBank accession numbers, or related GB files of the chloroplast genomes.")
                                   )),
                          tabPanel("Files upload", br(), fluidRow(column(width = 11, br(), p("You need to provide information related to at least two species for the analysis to proceed. After you have successfully provided your input files, press the 'Submit' button. This will activate its following box. Be patient till you are communicated through that box on how to proceed for downloading your output. Also note that the computations may take up to five minutes. Please be patient. For running the analysis, one can use either of the following cases."), strong("1. GB Files"), p("In this method, you can provide the GB files or the accession number of the related species to obtain the IR plot. This function is active under the 'GB File' tab"),  strong("2. Manual Files"),p("This method is useful for the times when  simple tabular format of the annotations, the genome sequences of the targeted species, and possibly their corresponding inverted repeat regions boundary coordinates are available. This function is active under 'Manual Files' tab."),  p("Further instructions for each case is given specifically under relevant tab."))),
                                   br(), br(),br(),
                                   tabsetPanel(
                                     tabPanel("GB File", br(), br(), p("In this section you can upload either your GeneBank files or provide their corresponding accession numbers if applicable. You can denote if you want to depict the SSC region in the reverse order by ticking the 'SSC', box for a corresponding input file.",br(), br(),"For example, you can test the program with these accession numbers:"), tags$em("NC_007898, NC_008096"), br(), br(),
                                   fluidRow(column(width = 2, textInput(inputId = "acc1", label = "Accession No.", placeholder = "1st Accession ...")),
                                            column(width = 3,  fileInput(inputId = "data1", label = "GeneBank File", placeholder = "or GB file")), column(width = 5, br(), verbatimTextOutput("sp1", placeholder = TRUE)), column(width = 1, strong("SSC"), br(), checkboxInput("S1", label = icon("transfer", lib =  "glyphicon"), value = FALSE)), column(width = 1, br())
                                   ),
                                   fluidRow(column(width = 2, textInput(inputId = "acc2", label = "", placeholder = "2nd Accession ...")),
                                            column(width = 3,  fileInput(inputId = "data2", label = "", placeholder = "or GB file")), column(width = 5, br(), verbatimTextOutput("sp2", placeholder = TRUE)),  column(width = 1, strong(""), br(), checkboxInput("S2", label = icon("transfer", lib =  "glyphicon"), value = FALSE)), column(width = 2, br())
                                   ),
                                   fluidRow(column(width = 2, textInput(inputId = "acc3", label = "", placeholder = "3rd Accession ...")),
                                            column(width = 3,  fileInput(inputId = "data3", label = "", placeholder = "or GB file")), column(width = 5, br(), verbatimTextOutput("sp3", placeholder = TRUE)),  column(width = 1, strong(""), br(), checkboxInput("S3", label = icon("transfer", lib =  "glyphicon"), value = FALSE)), column(width = 2, br())
                                   ),
                                   fluidRow(column(width = 2, textInput(inputId = "acc4", label = "", placeholder = "4th Accession ...")),
                                            column(width = 3,  fileInput(inputId = "data4", label = "", placeholder = "or GB file")), column(width = 5, br(), verbatimTextOutput("sp4", placeholder = TRUE)),  column(width = 1, strong(""), br(), checkboxInput("S4", label = icon("transfer", lib =  "glyphicon"), value = FALSE)), column(width = 2, br())
                                   ),
                                   fluidRow(column(width = 2, textInput(inputId = "acc5", label = "", placeholder = "5th Accession ...")),
                                            column(width = 3,  fileInput(inputId = "data5", label = "", placeholder = "or GB file")), column(width = 5, br(), verbatimTextOutput("sp5", placeholder = TRUE)),  column(width = 1, strong(""), br(), checkboxInput("S5", label = icon("transfer", lib =  "glyphicon"), value = FALSE)), column(width = 2, br())
                                   ),
                                   fluidRow(column(width = 2, textInput(inputId = "acc6", label = "", placeholder = "6th Accession ...")),
                                            column(width = 3,  fileInput(inputId = "data6", label = "", placeholder = "or GB file")), column(width = 5, br(), verbatimTextOutput("sp6", placeholder = TRUE)),  column(width = 1, strong(""), br(), checkboxInput("S6", label = icon("transfer", lib =  "glyphicon"), value = FALSE)), column(width = 2, br())
                                   ),
                                   fluidRow(column(width = 2, textInput(inputId = "acc7", label = "", placeholder = "7th Accession ...")),
                                            column(width = 3,  fileInput(inputId = "data7", label = "", placeholder = "or GB file")), column(width = 5, br(), verbatimTextOutput("sp7", placeholder = TRUE)),  column(width = 1, strong(""), br(), checkboxInput("S7", label = icon("transfer", lib =  "glyphicon"), value = FALSE)), column(width = 2, br())
                                   ),
                                   fluidRow(column(width = 2, textInput(inputId = "acc8", label = "", placeholder = "8th Accession ...")),
                                            column(width = 3,  fileInput(inputId = "data8", label = "", placeholder = "or GB file")), column(width = 5, br(), verbatimTextOutput("sp8", placeholder = TRUE)),  column(width = 1, strong(""), br(), checkboxInput("S8", label = icon("transfer", lib =  "glyphicon"), value = FALSE)), column(width = 2, br())
                                   ),
                                   fluidRow(column(width = 2, textInput(inputId = "acc9", label = "", placeholder = "9th Accession ...")),
                                            column(width = 3,  fileInput(inputId = "data9", label = "", placeholder = "or GB file")), column(width = 5, br(), verbatimTextOutput("sp9", placeholder = TRUE)),  column(width = 1, strong(""), br(), checkboxInput("S9", label = icon("transfer", lib =  "glyphicon"), value = FALSE)), column(width = 2, br())
                                   ),
                                   fluidRow(column(width = 2, textInput(inputId = "acc10", label = "", placeholder = "10th Accession ...")),
                                            column(width = 3,  fileInput(inputId = "data10", label = "", placeholder = "or GB file")), column(width = 5, br(), verbatimTextOutput("sp10", placeholder = TRUE)),  column(width = 1, strong(""), br(), checkboxInput("S10", label = icon("transfer", lib =  "glyphicon"), value = FALSE)), column(width = 2, br())
                                   ),
                                   column(width=10,verbatimTextOutput("stat"), actionButton("Gbsub", "Submit"), br(), br(), p("After submission the section below will turn innactive, meaning that analysis is in process. Be patient until you are prompted to download your plot."), verbatimTextOutput("StatusGB", placeholder = TRUE), br(), verbatimTextOutput("StatusGB2", placeholder = FALSE),  
                                   downloadButton("downloadData", "Download"))),
                                     tabPanel("Manual Files", br(), p("In case your GB files are of low quality or otherwise not available and still you would like to obtain the plot, you can provide your annotations as a simple tabular format (primary", tags$a("DOGMA", href="https://dogma.ccbb.utexas.edu"), "output) beside their genomes in fasta format and proceed to obtain the plot. In case you also have the values for the junction sites coordinates, you can provide them on the 'IR info' section as JLB, JSB, JSA, and JLA respectively. This decreases the processing time dramatically and would be useful in cases where for example the IR regions are not identical.", br(), br(), "A minimal example of a annotation input can be like", downloadLink('downloadDogma', 'this'), "file."),
                                              fluidRow(br(),column(width = 3, fileInput(inputId = "dogma1", label = "Annotations", placeholder = "1st Annotation ...")),
                                                       column(width = 3,  fileInput(inputId = "fas1", label = "Genome Fasta", placeholder = "And Fasta file")), column(width = 6, textInput(inputId = "IR1",label = "IR info (optional)", placeholder = "IRb end, IRb start, IRa end, IRa start")), column(width = 2, br())
                                              ),
                                              fluidRow(column(width = 3, fileInput(inputId = "dogma2", label = "", placeholder = "1st Annotation ...")),
                                                       column(width = 3,  fileInput(inputId = "fas2", label = "", placeholder = "And Fasta file")), column(width = 6, textInput(inputId = "IR2",label = "", placeholder = "IRb end, IRb start, IRa end, IRa start")), column(width = 2, br())
                                              ),
                                              fluidRow(column(width = 3, fileInput(inputId = "dogma3", label = "", placeholder = "1st Annotation ...")),
                                                       column(width = 3,  fileInput(inputId = "fas3", label = "", placeholder = "And Fasta file")), column(width = 6, textInput(inputId = "IR3",label = "", placeholder = "IRb end, IRb start, IRa end, IRa start")), column(width = 2, br())
                                              ),
                                              fluidRow(column(width = 3, fileInput(inputId = "dogma4", label = "", placeholder = "1st Annotation ...")),
                                                       column(width = 3,  fileInput(inputId = "fas4", label = "", placeholder = "And Fasta file")), column(width = 6, textInput(inputId = "IR4",label = "", placeholder = "IRb end, IRb start, IRa end, IRa start")), column(width = 2, br())
                                              ),
                                              fluidRow(column(width = 3, fileInput(inputId = "dogma5", label = "", placeholder = "1st Annotation ...")),
                                                       column(width = 3,  fileInput(inputId = "fas5", label = "", placeholder = "And Fasta file")), column(width = 6, textInput(inputId = "IR5",label = "", placeholder = "IRb end, IRb start, IRa end, IRa start")), column(width = 2, br())
                                              ),
                                              fluidRow(column(width = 3, fileInput(inputId = "dogma6", label = "", placeholder = "1st Annotation ...")),
                                                       column(width = 3,  fileInput(inputId = "fas6", label = "", placeholder = "And Fasta file")), column(width = 6, textInput(inputId = "IR6",label = "", placeholder = "IRb end, IRb start, IRa end, IRa start")), column(width = 2, br())
                                              ),
                                              fluidRow(column(width = 3, fileInput(inputId = "dogma7", label = "", placeholder = "1st Annotation ...")),
                                                       column(width = 3,  fileInput(inputId = "fas7", label = "", placeholder = "And Fasta file")), column(width = 6, textInput(inputId = "IR7",label = "", placeholder = "IRb end, IRb start, IRa end, IRa start")), column(width = 2, br())
                                              ),
                                              fluidRow(column(width = 3, fileInput(inputId = "dogma8", label = "", placeholder = "1st Annotation ...")),
                                                       column(width = 3,  fileInput(inputId = "fas8", label = "", placeholder = "And Fasta file")), column(width = 6, textInput(inputId = "IR8",label = "", placeholder = "IRb end, IRb start, IRa end, IRa start")), column(width = 2, br())
                                              ),
                                              fluidRow(column(width = 3, fileInput(inputId = "dogma9", label = "", placeholder = "1st Annotation ...")),
                                                       column(width = 3,  fileInput(inputId = "fas9", label = "", placeholder = "And Fasta file")), column(width = 6, textInput(inputId = "IR9",label = "", placeholder = "IRb end, IRb start, IRa end, IRa start")), column(width = 2, br())
                                              ),
                                              fluidRow(column(width = 3, fileInput(inputId = "dogma10", label = "", placeholder = "1st Annotation ...")),
                                                       column(width = 3,  fileInput(inputId = "fas10", label = "", placeholder = "And Fasta file")), column(width = 6, textInput(inputId = "IR10",label = "", placeholder = "IRb end, IRb start, IRa end, IRa start")), column(width = 2, br())
                                              ), column(width=10, actionButton("Go", "Submit"), br(), br(), p("After submission the section below will turn innactive, meaning that analysis is in process. Be patient until you are prompted to download your plot."), verbatimTextOutput("Status", placeholder = TRUE), br(), verbatimTextOutput("Status2", placeholder = FALSE), downloadButton("Down","Download file")))
                          )),
                          tabPanel("FAQ", br(), column(width = 1, ""), br(), column(width = 6, h4(strong("Q:"), tags$em(strong("Why the gene names are not displayed properly?"))), p(strong("A:"),"Check your input file(s), IRscope is based on the gene names specified in them."), br(),
                                                                                                     h4(strong("Q:"), tags$em(strong("Why some of the genes are extended to the large portion of the tracks?"))), p(strong("A:"),"This is due to the bad annotation implemented in the input file(s). Often this is related to the genes with introns or trans-spliced genes like rps12. Check and fix your input files."), br(),
                                                                                                     h4(strong("Q:"), tags$em(strong("I have a well annotated genome but not getting the plot?"))), p(strong("A:"),"IRscope can handle the redundant base pairs but if the sequence is of a very poor quality the program may fail in detecting the inverted regions and hence not produce any output."), br(),
                                                                                                     h4(strong("Q:"), tags$em(strong("I see genes overplotting on each other, how can I fix this?"))), p(strong("A:"),"This is mostly because at least one of your species differs significantly from the rest so that it creates too large radius for finding the genes in the vicinity of the junction site and consequently causes the populated genes in these respective areas. Try reducing the species sampling to a smaller group of more closely related species."), br(),
                                                                                                     h4(strong("Q:"), tags$em(strong("Can I modify and run the codes on my own computer?"))), p(strong("A:"),"Yes, please download the file in the 'About' section, run the whole script in your R session. This will set up a pseudo web interface on your computer which resemles and works exactly like the online version. You can further tune the functions as you wish." ), br(),
                                                                                                     h4(strong("Q:"), tags$em(strong("My analysis runs slow. Can I enhance its speed?"))), p(strong("A:"),"Yes, after downloading the file in the 'About' section, alter the 'parallel' arguments of the function 'IRinfo' into 'TRUE'. Then run the whole code on your R console. This will launch an instance with default 4 CPUs. If your machine is equipped with more cores, you can increase the number of CPUs with the 'nCPU' in 'p.d' subfunction inside the 'IRinfo' function."), br(),
                                                                                                     h4(strong("Q:"), tags$em(strong("What should be the format of the annotation files in the 'Manual Files' section?"))), p(strong("A:"),"This should be a plain text format (.txt) of the four tab separated columns, as start of the gene, end of the gene, its name, and the direction of it (+ or -). Note that the file needs to be without header." ), br(),
                                                                                                     h4(strong("Q:"), tags$em(strong("My species names are not plotted correctly, why?"))), p(strong("A:"),"The genome names in the 'GB Files' are read from the 'Organism' line of the file while the first space separated text of the genomes first line in the manual part is a determinant of species name in this section. In the 'Manual Files' section, the names are read from the first line after the '>' of the fasta file uploaded. Please edit them if they do not follow your expectations and rerun." ), br(),
                                                                                                     h4(strong("Q:"), tags$em(strong("My analysis too excessively long time and even after that the plot is of an overall low quality, why is that?"))), p(strong("A:"),"The IRscope is primarily optimized for the angiosperms and it often turns reliable results for other seed plants as well. However, further departure from the Embryophyta will extend this program to its limits. In such cases, you may want to consider use of the manual section with providing the the IR coordinates." ), br(),
                                                                                                     h4(strong("Q:"), tags$em(strong("How can I cite the program"))), p(strong("A:"), "The paper describing this web app in detail with the title", tags$em("IRscope: An online program to visualize the junction sites of chloroplast genomes"), "is submitted to", tags$em("Bioinformatics"), "journal. Please, cite the paper when you use the program."), br()                                                                                          
                                                                                                     
                          )),
                          tabPanel("About", br(), column(width = 6, p("IRscope and its all dependencies are coded in R. The detailed instruction on how to use each subfuction is provided as well. Also the accession numbers that used to test and train the program are listed at the end of the file.", "Click", downloadLink('downloadDat', 'here'), "to download the file after which, you can run all the codes (apart from accession numbers) into your R Console to obtain this app on your PC. Here are two example outputs of the program:" , downloadLink('downloadPIC1', 'example1'), "and", downloadLink('downloadPIC2', 'example2'), "."))),
                                   #uiOutput('markdown', inline = TRUE)
                          #),
                          tabPanel("Contact", br(), column(width = 6,  p("The paper describing this web app in detail with the title", tags$em("IRscope: An online program to visualize the junction sites of chloroplast genomes"), "is submitted to", tags$em("Bioinformatics"), "journal. Please cite that papers in all your correspondence. In case of practical issues please consult the FAQ sections first or email us at", tags$em("ali(dot)amiryousefi@helsinki(dot)fi."))))
                        )
               )
            #tabPanel("Contact")
))

server <- function(input, output) {
 
  
  
  ###GB File section#############################################
  
  ###Status update for the GB file submission
  output$StatusGB<- renderText({
    if(input$Gbsub){
      paste(paste(paste("Thank you for your GB submission at", date(), sep= " "), "!", sep=""), " You may now download your plot via the tab below.", sep="")
    }
  })
  
  
  
  ####Calculating the IR values with IRs in a pseudo text function in the GB file section
  output$StatusGB2<- renderText({
    if(input$Gbsub){
      dist<- isolate(IRs(gb(), SRev(), file=file))
      dist
    }
  })
  
  
  
  ###Downloading the result with the Download botton in the GB files section
  output$downloadData <- downloadHandler(
    filename = "IR.jpg",
    content = function(file) {
      #mboard(datagb(), file = file)
      IRs2(gb(), file=file)
      #jpeg(file,  width=8.3, height=8.9, units="in", res=100)
      #plot(rnorm(100))
      #dev.off()
    }
  )
  
  
  
  ###Reactive list of the objects needed for the calculation of this section i.e. GB files
  gb<- reactive({
    gbFiles<- list()
    if (input$acc1!=""){gbFiles[[1]]<- fetch.gb(input$acc1)}
    if (input$acc2!=""){gbFiles[[2]]<- fetch.gb(input$acc2)}
    if (input$acc3!=""){gbFiles[[3]]<- fetch.gb(input$acc3)}
    if (input$acc4!=""){gbFiles[[4]]<- fetch.gb(input$acc4)}
    if (input$acc5!=""){gbFiles[[5]]<- fetch.gb(input$acc5)}
    if (input$acc6!=""){gbFiles[[6]]<- fetch.gb(input$acc6)}
    if (input$acc7!=""){gbFiles[[7]]<- fetch.gb(input$acc7)}
    if (input$acc8!=""){gbFiles[[8]]<- fetch.gb(input$acc8)}
    if (input$acc9!=""){gbFiles[[9]]<- fetch.gb(input$acc9)}
    if (input$acc10!=""){gbFiles[[10]]<- fetch.gb(input$acc10)}
    inFile<- input$data1
    if (is.data.frame(inFile)){gbFiles[[1]]<-readLines(inFile$datapath)}
    inFile<- input$data2
    if (is.data.frame(inFile)){gbFiles[[2]]<-readLines(inFile$datapath)}
    inFile<- input$data3 
    if (is.data.frame(inFile)){gbFiles[[3]]<-readLines(inFile$datapath)}
    inFile<- input$data4 
    if (is.data.frame(inFile)){gbFiles[[4]]<-readLines(inFile$datapath)}
    inFile<- input$data5 
    if (is.data.frame(inFile)){gbFiles[[5]]<-readLines(inFile$datapath)}
    inFile<- input$data6
    if (is.data.frame(inFile)){gbFiles[[6]]<-readLines(inFile$datapath)}
    inFile<- input$data7
    if (is.data.frame(inFile)){gbFiles[[7]]<-readLines(inFile$datapath)}
    inFile<- input$data8 
    if (is.data.frame(inFile)){gbFiles[[8]]<-readLines(inFile$datapath)}
    inFile<- input$data9 
    if (is.data.frame(inFile)){gbFiles[[9]]<-readLines(inFile$datapath)}
    inFile<- input$data10 
    if (is.data.frame(inFile)){gbFiles[[10]]<-readLines(inFile$datapath)}
    gbFiles[which(gbFiles!="NULL")]
  })
  
  SRev<- reactive({
    SFiles<- list()
    inFile<- input$data1
    if (input$acc1!="" || is.data.frame(inFile)) {SFiles[[1]]<- input$S1}
    inFile<- input$data2
    if (input$acc2!="" || is.data.frame(inFile)) {SFiles[[2]]<- input$S2}
    inFile<- input$data3
    if (input$acc3!="" || is.data.frame(inFile)) {SFiles[[3]]<- input$S3}
    inFile<- input$data4
    if (input$acc4!="" || is.data.frame(inFile)) {SFiles[[4]]<- input$S4}
    inFile<- input$data5
    if (input$acc5!="" || is.data.frame(inFile)) {SFiles[[5]]<- input$S5}
    inFile<- input$data6
    if (input$acc6!="" || is.data.frame(inFile)) {SFiles[[6]]<- input$S6}
    inFile<- input$data7
    if (input$acc7!="" || is.data.frame(inFile)) {SFiles[[7]]<- input$S7}
    inFile<- input$data8
    if (input$acc8!="" || is.data.frame(inFile)) {SFiles[[8]]<- input$S8}
    inFile<- input$data9
    if (input$acc9!="" || is.data.frame(inFile)) {SFiles[[9]]<- input$S9}
    inFile<- input$data10
    if (input$acc10!="" || is.data.frame(inFile)) {SFiles[[10]]<- input$S10}
    SFiles[which(SFiles!="NULL")]
  })
  
  
  
  
  
  ###Dogma File section#############################################
  
  ###Status file for the submission of the files
  output$Status<- renderText({
   if(input$Go){
     paste(paste(paste("Thank you for your Manual submission at", date(), sep= " "), "!", sep=""), "You may now download your plot via the tab below.", sep="")
     }
     })
  
  
  
  ###Calculating the IR with the IRsD in a pseudo text function in the Dogma section
  output$Status2<- renderText({
   if(input$Go){
     dist<- isolate(IRsD(df(), ff(),  ii(), nf(), file=file))
     dist
    }
 })
  
  
  
  
  ###Downloading the result with the Download botton in the DOGMA files section
  output$Down <- downloadHandler(
    filename = "IR.jpg",
    content = function(file) {
      #mboard(datagb(), file = file)
      IRsD2(df(), file=file)
      #jpeg(file,  width=8.3, height=8.9, units="in", res=100)
      #plot(rnorm(100))
      #dev.off()
    }
  )
 
  
  
 
  ###Reactive list of objects for the calculation of this section
  
  ###df for dogma files
  df<- reactive({
    dFiles<- list()
    inFile<- input$dogma1
    if (is.data.frame(inFile)){dFiles[[1]]<-GnlBuilder(trnDogma(read.table(inFile$datapath)))}
    inFile<- input$dogma2
    if (is.data.frame(inFile)){dFiles[[2]]<-GnlBuilder(trnDogma(read.table(inFile$datapath)))}
    inFile<- input$dogma3
    if (is.data.frame(inFile)){dFiles[[3]]<-GnlBuilder(trnDogma(read.table(inFile$datapath)))}
    inFile<- input$dogma4
    if (is.data.frame(inFile)){dFiles[[4]]<-GnlBuilder(trnDogma(read.table(inFile$datapath)))}
    inFile<- input$dogma5
    if (is.data.frame(inFile)){dFiles[[5]]<-GnlBuilder(trnDogma(read.table(inFile$datapath)))}
    inFile<- input$dogma6
    if (is.data.frame(inFile)){dFiles[[6]]<-GnlBuilder(trnDogma(read.table(inFile$datapath)))}
    inFile<- input$dogma7
    if (is.data.frame(inFile)){dFiles[[7]]<-GnlBuilder(trnDogma(read.table(inFile$datapath)))}
    inFile<- input$dogma8
    if (is.data.frame(inFile)){dFiles[[8]]<-GnlBuilder(trnDogma(read.table(inFile$datapath)))}
    inFile<- input$dogma9
    if (is.data.frame(inFile)){dFiles[[9]]<-GnlBuilder(trnDogma(read.table(inFile$datapath)))}
    inFile<- input$dogma10
    if (is.data.frame(inFile)){dFiles[[10]]<-GnlBuilder(trnDogma(read.table(inFile$datapath)))}
    dFiles[which(dFiles!="NULL")]
  })
  
  
  ###ff for fasta files
  ff<- reactive({
    fFiles<- list()
    inFile<- input$fas1
    if (is.data.frame(inFile)){fFiles[[1]]<-rdnFixerD(as.character(read.fasta(inFile$datapath)[[1]]))}
    inFile<- input$fas2
    if (is.data.frame(inFile)){fFiles[[2]]<-rdnFixerD(as.character(read.fasta(inFile$datapath)[[1]]))}
    inFile<- input$fas3
    if (is.data.frame(inFile)){fFiles[[3]]<-rdnFixerD(as.character(read.fasta(inFile$datapath)[[1]]))}
    inFile<- input$fas4
    if (is.data.frame(inFile)){fFiles[[4]]<-rdnFixerD(as.character(read.fasta(inFile$datapath)[[1]]))}
    inFile<- input$fas5
    if (is.data.frame(inFile)){fFiles[[5]]<-rdnFixerD(as.character(read.fasta(inFile$datapath)[[1]]))}
    inFile<- input$fas6
    if (is.data.frame(inFile)){fFiles[[6]]<-rdnFixerD(as.character(read.fasta(inFile$datapath)[[1]]))}
    inFile<- input$fas7
    if (is.data.frame(inFile)){fFiles[[7]]<-rdnFixerD(as.character(read.fasta(inFile$datapath)[[1]]))}
    inFile<- input$fas8
    if (is.data.frame(inFile)){fFiles[[8]]<-rdnFixerD(as.character(read.fasta(inFile$datapath)[[1]]))}
    inFile<- input$fas9
    if (is.data.frame(inFile)){fFiles[[9]]<-rdnFixerD(as.character(read.fasta(inFile$datapath)[[1]]))}
    inFile<- input$fas10
    if (is.data.frame(inFile)){fFiles[[10]]<-rdnFixerD(as.character(read.fasta(inFile$datapath)[[1]]))}
     fFiles[which(fFiles!="NULL")]
  })
  
  
  
  ###ii for IR info
  ii<- reactive({
    iFiles<- list()
    inFile<-input$fas1
    if (input$IR1=="" &&  is.data.frame(inFile)) {iFiles[[1]]<- 999}
    else if(input$IR1!=""){iFiles[[1]]<- c(as.numeric(unlist(strsplit(input$IR1, split = ","))) , length(as.character(read.fasta(inFile$datapath)[[1]])))}
    inFile<-input$fas2
    if (input$IR2=="" &&  is.data.frame(inFile)) {iFiles[[2]]<- 999}
    else if(input$IR2!=""){iFiles[[2]]<- c(as.numeric(unlist(strsplit(input$IR2, split = ","))) , length(as.character(read.fasta(inFile$datapath)[[1]])))}
    inFile<-input$fas3
    if (input$IR3=="" &&  is.data.frame(inFile)) {iFiles[[3]]<- 999}
    else if(input$IR3!=""){iFiles[[3]]<- c(as.numeric(unlist(strsplit(input$IR3, split = ","))) , length(as.character(read.fasta(inFile$datapath)[[1]])))}
    inFile<-input$fas4
    if (input$IR4=="" &&  is.data.frame(inFile)) {iFiles[[4]]<- 999}
    else if(input$IR4!=""){iFiles[[4]]<- c(as.numeric(unlist(strsplit(input$IR4, split = ","))) , length(as.character(read.fasta(inFile$datapath)[[1]])))}
    inFile<-input$fas5
    if (input$IR5=="" &&  is.data.frame(inFile)) {iFiles[[5]]<- 999}
    else if(input$IR5!=""){iFiles[[5]]<- c(as.numeric(unlist(strsplit(input$IR5, split = ","))) , length(as.character(read.fasta(inFile$datapath)[[1]])))}
    inFile<-input$fas6
    if (input$IR6=="" &&  is.data.frame(inFile)) {iFiles[[6]]<- 999}
    else if(input$IR6!=""){iFiles[[6]]<- c(as.numeric(unlist(strsplit(input$IR6, split = ","))) , length(as.character(read.fasta(inFile$datapath)[[1]])))}
    inFile<-input$fas7
    if (input$IR7=="" &&  is.data.frame(inFile)) {iFiles[[7]]<- 999}
    else if(input$IR7!=""){iFiles[[7]]<- c(as.numeric(unlist(strsplit(input$IR7, split = ","))) , length(as.character(read.fasta(inFile$datapath)[[1]])))}
    inFile<-input$fas8
    if (input$IR8=="" &&  is.data.frame(inFile)) {iFiles[[8]]<- 999}
    else if(input$IR8!=""){iFiles[[8]]<- c(as.numeric(unlist(strsplit(input$IR8, split = ","))) , length(as.character(read.fasta(inFile$datapath)[[1]])))}
    inFile<-input$fas9
    if (input$IR9=="" &&  is.data.frame(inFile)) {iFiles[[9]]<- 999}
    else if(input$IR9!=""){iFiles[[9]]<- c(as.numeric(unlist(strsplit(input$IR9, split = ","))) , length(as.character(read.fasta(inFile$datapath)[[1]])))}
    inFile<-input$fas10
    if (input$IR10=="" &&  is.data.frame(inFile)) {iFiles[[10]]<- 999}
    else if(input$IR10!=""){iFiles[[10]]<- c(as.numeric(unlist(strsplit(input$IR10, split = ","))) , length(as.character(read.fasta(inFile$datapath)[[1]])))}
    iFiles[which(iFiles!="NULL")]
  })
  
  
  
  ###nf for the names of the species
  nf<- reactive({
    nFiles<- list()
    inFile<- input$fas1
    if (is.data.frame(inFile)){nFiles[[1]]<-gsub("_", " ", names(read.fasta(inFile$datapath)))}
    inFile<- input$fas2
    if (is.data.frame(inFile)){nFiles[[2]]<-gsub("_", " ", names(read.fasta(inFile$datapath)))}
    inFile<- input$fas3
    if (is.data.frame(inFile)){nFiles[[3]]<-gsub("_", " ", names(read.fasta(inFile$datapath)))}
    inFile<- input$fas4
    if (is.data.frame(inFile)){nFiles[[4]]<-gsub("_", " ", names(read.fasta(inFile$datapath)))}
    inFile<- input$fas5
    if (is.data.frame(inFile)){nFiles[[5]]<-gsub("_", " ", names(read.fasta(inFile$datapath)))}
    inFile<- input$fas6
    if (is.data.frame(inFile)){nFiles[[6]]<-gsub("_", " ", names(read.fasta(inFile$datapath)))}
    inFile<- input$fas7
    if (is.data.frame(inFile)){nFiles[[7]]<-gsub("_", " ", names(read.fasta(inFile$datapath)))}
    inFile<- input$fas8
    if (is.data.frame(inFile)){nFiles[[8]]<-gsub("_", " ", names(read.fasta(inFile$datapath)))}
    inFile<- input$fas9
    if (is.data.frame(inFile)){nFiles[[9]]<-gsub("_", " ", names(read.fasta(inFile$datapath)))}
    inFile<- input$fas10
    if (is.data.frame(inFile)){nFiles[[10]]<-gsub("_", " ", names(read.fasta(inFile$datapath)))}
    nFiles[which(nFiles!="NULL")]
  })
  
  
  
  
  #output$markdown <- renderUI({
  #  HTML(markdown::markdownToHTML(knit('Rmarkdown.Rmd', quiet = TRUE)))
  #})
  
  
  #doog <- readLines("dogma.txt")
  
  output$downloadDogma <- downloadHandler(###Downloading the Dogma file
    filename = function() {
      paste("annotation-", Sys.Date(), ".txt", sep="")
    },
    content = function(file) {
      writeLines(doog, file)
    }
  )
  
  #pic1 <- readJPEG("IR1.jpg")
  #pic2 <- readJPEG("IRprevioustest.jpg")
  
  output$downloadPIC1 <- downloadHandler(###Downloading the first PIC
    filename = function() {
      paste("IR1example-", Sys.Date(), ".jpg", sep="")
    },
    content = function(file) {
      writeJPEG(pic1, file)
    }
  )
  
  output$downloadPIC2 <- downloadHandler(###Downloading the first PIC
    filename = function() {
      paste("IR2example-", Sys.Date(), ".jpg", sep="")
    },
    content = function(file) {
      writeJPEG(pic2, file)
    }
  )
  
  
  
  
 # data <- readLines("IRscope.R")
  
  output$downloadDat <- downloadHandler(###Downloading the R codes
    filename = function() {
      paste("data-", Sys.Date(), ".R", sep="")
    },
    content = function(file) {
      writeLines(data, file)
    }
  )
  
  
  
  
  
  
  ### Updating the section areas of the borowsed files in GB section
  
  
  output$sp1 <- renderText({
    
    if (input$acc1!="") {print(head(fetch.gb(input$acc1))[2])}
    else {
      frame1<- input$data1
      if (is.data.frame(frame1)) {print(readLines(frame1$datapath)[2])}
    }
    #isolate(input$clear1)
    #emp<- isolate(input$sp1)
    #emp
  })
  
  output$sp2 <- renderText({ if (input$acc2!="") {print(head(fetch.gb(input$acc2))[2])}
          else {
            frame2<- input$data2
          if (is.data.frame(frame2)) {print(readLines(frame2$datapath)[2])}
            }
    })
  
  output$sp3 <- renderText({ if (input$acc3!="") {print(head(fetch.gb(input$acc3))[2])}
    else {
      frame3<- input$data3
      if (is.data.frame(frame3)) {print(readLines(frame3$datapath)[2])}
    }
  })
  
  output$sp4 <- renderText({ if (input$acc4!="") {print(head(fetch.gb(input$acc4))[2])}
    else {
      frame4<- input$data4
      if (is.data.frame(frame4)) {print(readLines(frame4$datapath)[2])}
    }
  })
  
  output$sp5 <- renderText({ if (input$acc5!="") {print(head(fetch.gb(input$acc5))[2])}
    else {
      frame5<- input$data5
      if (is.data.frame(frame5)) {print(readLines(frame5$datapath)[2])}
    }
  })
  
  output$sp6 <- renderText({ if (input$acc6!="") {print(head(fetch.gb(input$acc6))[2])}
    else {
      frame6<- input$data6
      if (is.data.frame(frame6)) {print(readLines(frame6$datapath)[2])}
    }
  })
  
  output$sp7 <- renderText({ if (input$acc7!="") {print(head(fetch.gb(input$acc7))[2])}
    else {
      frame7<- input$data7
      if (is.data.frame(frame7)) {print(readLines(frame7$datapath)[2])}
    }
  })
  
  output$sp8 <- renderText({ if (input$acc8!="") {print(head(fetch.gb(input$acc8))[2])}
    else {
      frame8<- input$data8
      if (is.data.frame(frame8)) {print(readLines(frame8$datapath)[2])}
    }
  })
  
  output$sp9 <- renderText({ if (input$acc9!="") {print(head(fetch.gb(input$acc9))[2])}
    else {
      frame9<- input$data9
      if (is.data.frame(frame9)) {print(readLines(frame9$datapath)[2])}
    }
  })
  
  output$sp10 <- renderText({ 
    if (input$acc10!="") {print(head(fetch.gb(input$acc10))[2])}
    else {
      frame10<- input$data10
      if (is.data.frame(frame10)) {print(readLines(frame10$datapath)[2])}
    }
  })
  
}

shinyApp(ui, server)