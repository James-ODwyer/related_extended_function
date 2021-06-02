
# Edited related functions to allow for more loci to be analysed at once
# Base code for familysim call left the same. The code used here subsets the data into quantities that can be 
# analysed using FORTRAN underlying calls from the package related


#familysim function changed to 
familysimunlimited <- function(input, ninds,nloci){
  
  if (input$nloci < 40) {
    
    
    simdatatotal <- familysim(input$freqs, ninds)
    
  }
  else {
    vector0 <- vector() 
    vector1 <- vector()  
    vector2 <- vector()  
    vector3 <- vector() 
    vector0 <- "ID"
    
    datalist = list()
    datalist2 = list()
    #This can be changed to whatever the ideal number of loci is at a time but is set at roughly 2* the size of
    # each group of subsetted data
    if (nloci < 40) {
      stop("function requires either at least 40 loci to be simulated from the input file or all loci to be simulated")
      
    }
    
    
    if (nloci <  input$nloci) {
      
      loci <- nloci
    }
    else {
      
      loci  <- input$nloci 
    }
    # This could also be fixed to be more clever and not miss out on the final few loci as the code removes the 
    # remainder of loci not divisible by 40 (not divisible by the group size)
    groups <- loci/20
    groups <- floor(groups)
    
    # creating naming vector for final compiled work
    for (j in c(1:(sum(groups)*20)))  {
      vector1[j] <- paste0("locus",j,"a")
      vector2[j] <- paste0("locus",j,"b")
    }
    vector3 <- c(rbind(vector1,vector2))
    vector3 <- c(vector0,vector3)
    
    #Subsetting input data (where the value of how many loci per group can be changed)
    for (j in c(1:sum(groups))) {
      inputX <- input
      inputX$gdata <- input$gdata[,c(1,((j*40)-39):((j*40)+1))]
      inputX$nloci <-((ncol(inputX$gdata)-1)/2)
      inputX$nalleles <- (ncol(inputX$gdata)-1)
      inputX$freqs <- input$freqs[((j*20)-19):(j*20)]
      colnames(inputX$gdata) <- c(1:41)
      
      datalist[[j]] <- inputX
    }
    #The family sim function  from the related package
    for (j in c(1:sum(groups))) {
      simdataX <- familysim(datalist[[j]]$freqs, ninds)
      datalist2[[j]] <- simdataX
    }
    #Removal of extra ID columns
    for (j in c(2:sum(groups))) {
      datalist2[[j]]$ID <- NULL
    }
    #combining together data
    simdatatotal <- do.call(cbind,datalist2)
    colnames(simdatatotal) <- vector3
    return(simdatatotal)
  }
}

#compare estimators changed to
compareestimatorsanyvalue <- function(input, ninds, nloci, trioml = 1, wang = 1, lynchli = 1, lynchrd = 1, ritland = 0, quellergt = 1, dyadml = 0) {
  if (input$nloci < 40) {
    
    
    simdatatotal <- familysim(input$freqs, ninds)
    
  }
  else {
    vector0 <- vector() 
    vector1 <- vector()  
    vector2 <- vector()  
    vector3 <- vector() 
    vector0 <- "ID"
    
    datalist = list()
    datalist2 = list()
    #This can be changed to whatever the ideal number of loci is at a time but is set at roughly 2* the size of
    # each group of subsetted data
    if (nloci < 40) {
      stop("function requires either at least 40 loci to be simulated from the input file or all loci to be simulated")
      
    }
    
    
    if (nloci <  input$nloci) {
      
      loci <- nloci
    }
    else {
      
      loci  <- input$nloci 
    }
    # This could also be fixed to be more clever and not miss out on the final few loci as the code removes the 
    # remainder of loci not divisible by 40 (not divisible by the group size)
    groups <- loci/20
    groups <- floor(groups)
    
    # creating naming vector for final compiled work
    for (j in c(1:(sum(groups)*20)))  {
      vector1[j] <- paste0("locus",j,"a")
      vector2[j] <- paste0("locus",j,"b")
    }
    vector3 <- c(rbind(vector1,vector2))
    vector3 <- c(vector0,vector3)
    
    #Subsetting input data (where the value of how many loci per group can be changed)
    for (j in c(1:sum(groups))) {
      inputX <- input
      inputX$gdata <- input$gdata[,c(1,((j*40)-39):((j*40)+1))]
      inputX$nloci <-((ncol(inputX$gdata)-1)/2)
      inputX$nalleles <- (ncol(inputX$gdata)-1)
      inputX$freqs <- input$freqs[((j*20)-19):(j*20)]
      colnames(inputX$gdata) <- c(1:41)
      
      datalist[[j]] <- inputX
    }
    #The family sim function  from the related package
    for (j in c(1:sum(groups))) {
      simdataX <- familysim(datalist[[j]]$freqs, ninds)
      datalist2[[j]] <- simdataX
    }
    #Removal of extra ID columns
    for (j in c(2:sum(groups))) {
      datalist2[[j]]$ID <- NULL
    }
    #combining together data
    simdatatotal <- do.call(cbind,datalist2)
    colnames(simdatatotal) <- vector3
  }
  #previous code shown just below here, previous code just called family sim function
  #simdata <- familysim(filename$freqs, ninds)
  # function additionally changed here to allow for compare estimators to use all relatedness measures and not just
  # the for moment based estimators
  output <- coancestry(simdatatotal, trioml = trioml, wang = wang, lynchli = lynchli, lynchrd = lynchrd, ritland = ritland, quellergt = quellergt, dyadml = dyadml)
  simrel <- cleanuprvals(output$relatedness, ninds)
  urval <- rep(0, ninds)
  hsval <- rep(0.25, ninds)
  fsval <- rep(0.5, ninds)
  poval <- rep(0.5, ninds)
  relvals <- c(poval, fsval, hsval, urval)
  if (wang == 1) {
    wangcor <- cor(relvals, simrel[, 6])
  }
  if (lynchli == 1) {
    lynchlicor <- cor(relvals, simrel[, 7])
  }
  if (lynchrd == 1) {
    lynchrdcor <- cor(relvals, simrel[, 8])
  }
  if (quellergt == 1) {
    quellergtcor <- cor(relvals, simrel[, 10])
  }
  if (ritland == 1) {
    ritlandcor <- cor(relvals, simrel[, 9])
  }
  if(dyadml == 1) {
    dyadmlcor <- cor(relvals, simrel[, 11])
  }
  if(trioml == 1) {
    triomlcor <- cor(relvals, simrel[, 5])
  }
  cat(sprintf("\nCorrelation Coefficients Between Observed & Expected Values:\n"))
  if (wang == 1) {
    cat(sprintf("wang\t\t%f\n", wangcor))
  }
  if (lynchli == 1) {
    cat(sprintf("lynchli\t\t%f\n", lynchlicor))
  }
  if (lynchrd == 1) {
    cat(sprintf("lynchrd\t\t%f\n", lynchrdcor))
  }
  if (quellergt == 1) {
    cat(sprintf("quellergt\t%f\n", quellergtcor))
  }
  if (ritland == 1) {
    cat(sprintf("ritland\t%f\n", ritlandcor))
  }
  if (dyadml == 1) {
    cat(sprintf("dyadml\t%f\n", dyadmlcor))
  }
  if(trioml == 1) {
    cat(sprintf("trioml\t%f\n", triomlcor))
  }
  
  if (wang == 1) {
    wangpo <- simrel[1:ninds, 6]
    wangfs <- simrel[(ninds + 1):(2 * ninds), 6]
    wanghs <- simrel[((2 * ninds) + 1):(3 * ninds), 6]
    wangur <- simrel[((3 * ninds) + 1):(4 * ninds), 6]
  }
  else {
    wangpo <- NULL
    wangfs <- NULL
    wanghs <- NULL
    wangur <- NULL
  }
  if (lynchli == 1) {
    lynchlipo <- simrel[1:ninds, 7]
    lynchlifs <- simrel[(ninds + 1):(2 * ninds), 7]
    lynchlihs <- simrel[((2 * ninds) + 1):(3 * ninds), 7]
    lynchliur <- simrel[((3 * ninds) + 1):(4 * ninds), 7]
  } else {
    lynchlipo <- NULL
    lynchlifs <- NULL
    lynchlihs <- NULL
    lynchliur <- NULL
  }
  if (lynchrd == 1){
    lynchrdpo <- simrel[1:ninds, 8]
    lynchrdfs <- simrel[(ninds + 1):(2 * ninds), 8]
    lynchrdhs <- simrel[((2 * ninds) + 1):(3 * ninds), 8]
    lynchrdur <- simrel[((3 * ninds) + 1):(4 * ninds), 8]
  } else { 
    lynchrdpo <- NULL
    lynchrdfs <- NULL
    lynchrdhs <- NULL
    lynchrdur <- NULL
  }
  if (quellergt == 1) {
    quellergtpo <- simrel[1:ninds, 10]
    quellergtfs <- simrel[(ninds + 1):(2 * ninds), 10]
    quellergths <- simrel[((2 * ninds) + 1):(3 * ninds), 10]
    quellergtur <- simrel[((3 * ninds) + 1):(4 * ninds), 10]
  } else {
    quellergtpo <- NULL
    quellergtfs <- NULL
    quellergths <- NULL
    quellergtur <- NULL
  }
  if (ritland == 1) {
    ritlandpo <- simrel[1:ninds, 9]
    ritlandfs <- simrel[(ninds + 1):(2 * ninds), 9]
    ritlandhs <- simrel[((2 * ninds) + 1):(3 * ninds), 9]
    ritlandur <- simrel[((3 * ninds) + 1):(4 * ninds), 9]
  }else {
    ritlandpo <- NULL
    ritlandfs <- NULL
    ritlandhs <- NULL
    ritlandur <- NULL
  }
  if (dyadml == 1){
    dyadmlpo <- simrel[1:ninds, 11]
    dyadmlfs <- simrel[(ninds + 1):(2 * ninds), 11]
    dyadmlhs <- simrel[((2 * ninds) + 1):(3 * ninds), 11]
    dyadmlur <- simrel[((3 * ninds) + 1):(4 * ninds), 11]
  }else {
    dyadmlpo <- NULL
    dyadmlfs <- NULL
    dyadmlhs <- NULL
    dyadmlur <- NULL
  }
  if (trioml == 1) {
    triomlpo <- simrel[1:ninds, 5]
    triomlfs <- simrel[(ninds + 1):(2 * ninds), 5]
    triomlhs <- simrel[((2 * ninds) + 1):(3 * ninds), 5]
    triomlur <- simrel[((3 * ninds) + 1):(4 * ninds), 5]
  }else {
    triomlpo <- NULL
    triomlfs <- NULL
    triomlhs <- NULL
    triomlur <- NULL    
  }
  sumstests <- sum(wang, lynchli, lynchrd, quellergt, ritland, dyadml, trioml)
  
  if (wang == 1) {
    wang <- rep("W", ninds)
  } else {
    wang <- NULL
  }
  
  if (lynchli == 1) {
    lynchli <- rep("L & L", ninds)
  } else {
    lynchli <- NULL
  }
  
  if (lynchrd == 1) {
    lynchrd <- rep("L & R", ninds)
  } else {
    lynchrd <- NULL
  }
  
  if (quellergt == 1) {
    quellergt <- rep("Q & G", ninds)
  } else {
    quellergt <- NULL
  }
  if (ritland == 1){
    ritland <- rep("Ritl", ninds)
  } else {
    ritland <- NULL
  }
  
  if (dyadml == 1) {
    dyadml <- rep("Dyadml", ninds)
  } else {
    dyadml <- NULL
  }
  if (trioml == 1) {
    trioml <- rep("Trioml", ninds)
  } else {
    trioml <- NULL
  }
  estimator2 <- c(wang, lynchli, lynchrd, quellergt, ritland, dyadml, trioml)
  Estimator <- rep(estimator2, 4)
  po <- rep("Parent-Offspring", (sumstests * ninds))
  fs <- rep("Full-Sibs", (sumstests * ninds))
  hs <- rep("Half-Sibs", (sumstests * ninds))
  ur <- rep("Unrelated", (sumstests * ninds))
  relationship <- c(po, fs, hs, ur)
  relatednesspo <- c(wangpo, lynchlipo, lynchrdpo, quellergtpo, ritlandpo, dyadmlpo, triomlpo)
  relatednessfs <- c(wangfs, lynchlifs, lynchrdfs, quellergtfs, ritlandfs, dyadmlfs, triomlfs)
  relatednesshs <- c(wanghs, lynchlihs, lynchrdhs, quellergths, ritlandhs, dyadmlhs, triomlhs)
  relatednessur <- c(wangur, lynchliur, lynchrdur, quellergtur, ritlandur, dyadmlur, triomlur)
  Relatedness_Value <- c(relatednesspo, relatednessfs, relatednesshs, 
                         relatednessur)
  
  combineddata <- as.data.frame(cbind(Estimator, relationship, 
                                      Relatedness_Value))
  
  combineddata$Relatedness_Value <- as.numeric(as.character(combineddata$Relatedness_Value))
  ggplot(combineddata, aes(x = Estimator, y = Relatedness_Value), 
         ylim = c(-0.5, 1)) + geom_boxplot() + facet_wrap(~relationship)
}
