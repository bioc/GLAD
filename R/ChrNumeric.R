ChrNumeric <- function(Chromosome)
  {
    if (!is.numeric(Chromosome))
      {

        if (is.factor(Chromosome))
          {
            Chromosome <- as.character(Chromosome)
          }
     
        Chromosome <- gsub("([cC][hH][rR][ ]*)","",Chromosome)
        indChrX <- which(Chromosome=="X")
        Chromosome[indChrX] <- 23
        indChrY <- which(Chromosome=="Y")
        Chromosome[indChrY] <- 24
        Chromosome <- as.numeric(Chromosome)
      }

    return(Chromosome)

  }
