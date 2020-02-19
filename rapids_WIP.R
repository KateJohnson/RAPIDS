library(shiny)
library(shinyjs)
library(ggplot2)
library(reshape2)
library(beanplot)
library(htmltools)
library(Rlab)
library(haven)
library(here)
library(plyr)
library(purrr)


rapids <- function(base_biom, base_outcomes, base_x, treatment, betam, horizon) {
  
  if (ncol(treatment)!=horizon+2) {
    stop("treatment not specified over horizon")
  }
  if ((ncol(base_biom)!=2) | (ncol(base_outcomes)!=2)) {
    stop("baseline biomarker levels or outcomes not specified over last two periods")
  }
  
  m <-1000                        ## no. of monte carlo iterations
  n_biom <-9                   ## no. of biomarkers
  n_outcomes <- 13            ## no. of outcomes
  n_x <- 10                   ## No. of demographic factors
  r <- abs(nrow(base_biom)/n_biom)  ## no. of individuals
  s <- 2                      ## Number of parameter draws (added)
  
  ## coeffcients for the power variance of the predicted biomarkers as a function of their respective means
  varcoeff = matrix(c(1.619101, -5.776551, 1.714504, -3.596814, .4956913, 1.600088, 3.749151, -8.312432, 1.181856,-1.514088,1.106989,.511323, 2.294562, -5.855728, 1.407148, .4426071, .2707567, 2.84833 ), 2, 9 )
  
  
  fi_biom = matrix(, 1, r*n_biom*(horizon+2)*m*s) # added s dimension
  dim(fi_biom)=c((horizon+2),n_biom,  r, m, s)
  fi_biom[1:2, 1:n_biom, 1:r, , ] = t(base_biom)

  fi_outcomes = matrix(, 1, r*n_outcomes*(horizon+2)*m*s) # added s dimension
  dim(fi_outcomes)=c((horizon+2),n_outcomes,  r, m, s)
  fi_outcomes[1:2, 1:n_outcomes,1:r,  , ] = t(base_outcomes)
  
  fi_x = matrix(, 1, r*n_x*(horizon+1)*m*s)  # added s dimension
  dim(fi_x)=c((horizon+1), n_x,r, m, s)
  fi_x[1, 1:n_x, 1:r, , ] = base_x
  
  treat = as.vector(unlist(treatment))
  dim(treat) = c(13, horizon+2, r)
  
  ## microsimulation
  for(i in 1:r){
    for(t in 1:horizon) {
      
      ## Update demographics
      fi_x[t+1, 1:n_x , i, , ] =  fi_x[t, 1:n_x, i, , ]     ## Copy demog vector to next time period
      fi_x[t+1, 1, i, , ] = fi_x[t+1, 1, i, , ] + 0.25     ## Increase age by .25 years
      fi_x[t+1, 2, i, , ] = fi_x[t+1, 1, i, , ]^2           ## Canculate Age ^2
      fi_x[t+1, 10, i, , ] = fi_x[t+1, 10, i, , ] + 0.25   ## Increase duration of diabetes by .25 years
      
      ## Death UPDATE
      bc <- 10
      pr <-  betam[149, bc][rep(1,each=2)] +
                              mapply("%*%", replicate(s,t(betam[1:10, bc]),simplify=FALSE), alply(fi_x[t+1,  , i, , ], c(NULL,NULL,3)))
                            + mapply("%*%", replicate(s,t(betam[11:23, bc]),simplify=FALSE), alply(fi_outcomes[t+1, , i, , ], c(NULL,NULL,3)))
                            + mapply("%*%", replicate(s,t(betam[24:32, bc]),simplify=FALSE), alply(fi_biom[t+1, , i, , ], c(NULL,NULL,3)))
                            + mapply("%*%", replicate(s,t(betam[33:45, bc]),simplify=FALSE), replicate(s,treat[ ,t+1, i],simplify=FALSE))
                            + mapply("%*%", replicate(s,t(betam[46:58, bc]),simplify=FALSE), alply(fi_outcomes[t+1, , i, , ], c(NULL,NULL,3))) * 
                                    fi_x[t+1, 1, i, , ]
                            + mapply("%*%", replicate(s,t(betam[59:67, bc]),simplify=FALSE), alply(fi_biom[t+1, , i, , ], c(NULL,NULL,3))) *
                                    fi_x[t+1, 1, i, , ]
                            + mapply("%*%", replicate(s,t(betam[68:80, bc]),simplify=FALSE), replicate(s,treat[ ,t+1, i],simplify=FALSE)) * 
                                    fi_x[t+1, 1, i, , ]
                            + mapply("%*%", replicate(s,t(betam[81:92, bc]),simplify=FALSE), alply(fi_outcomes[t+1, 2:n_outcomes, i, , ], c(NULL,NULL,3)))
                            + mapply("%*%", replicate(s,t(betam[93:101, bc]),simplify=FALSE), alply(fi_biom[t, , i, , ], c(NULL,NULL,3)))
                            + mapply("%*%", replicate(s,t(betam[102:114, bc]),simplify=FALSE), replicate(s,treat[ ,t+1, i],simplify=FALSE))
                            + mapply("%*%", replicate(s,t(betam[115:126, bc]),simplify=FALSE), alply(fi_outcomes[t,2:n_outcomes  , i, , ], c(NULL,NULL,3))) *
                                    fi_x[t+1, 1, i, , ]
                            + mapply("%*%", replicate(s,t(betam[127:135, bc]),simplify=FALSE), alply(fi_biom[t, , i, , ], c(NULL,NULL,3))) *
                                    fi_x[t+1, 1, i, , ]
                            + mapply("%*%", replicate(s,t(betam[136:148, bc]),simplify=FALSE), replicate(s,treat[ ,t+1, i],simplify=FALSE)) *
                                    fi_x[t+1, 1, i, , ]
      
      print(pr)
      pr = ifelse(pr<0, 0, ifelse(pr>1, 1, pr)) 
      pri = lapply(pr, function(j) { rbern(m, j) } )
      outcome.t2 <- pmap(list(alply(fi_outcomes[t+1, 1, i, , ], 2), pri), function(x,y) { ifelse(x==0, y, x) })
      fi_outcomes[t+2, 1, i, , ] <- matrix(unlist(outcome.t2), ncol=length(outcome.t2))
      
      
      ## BIOMARKER UPDATES
      
        for (b in 1:9){
        bc <- b
        mu <- fi_biom[t+1, b, i, , ] + 
                              betam[149, bc][rep(1,each=2)] +
                              mapply("%*%", replicate(s,t(betam[1:10, bc]),simplify=FALSE), alply(fi_x[t+1,  , i, , ], c(NULL,NULL,3)))
                            + mapply("%*%", replicate(s,t(betam[11, bc]),simplify=FALSE), alply(fi_outcomes[t+2, 1, i, , ], c(NULL,2))) 
                            + mapply("%*%", replicate(s,t(betam[12:23, bc]),simplify=FALSE), alply(fi_outcomes[t+1, 2:n_outcomes, i, , ], c(NULL,NULL,3)))
                            + mapply("%*%", replicate(s,t(betam[24:32, bc]),simplify=FALSE), alply(fi_biom[t+1,  , i, , ], c(NULL,NULL,3)))
                            + mapply("%*%", replicate(s,t(betam[33:45, bc]),simplify=FALSE), replicate(s,treat[ ,t+1, i],simplify=FALSE))
                            + mapply("%*%", replicate(s,t(betam[46, bc]),simplify=FALSE), alply(fi_outcomes[t+2, 1, i, , ], c(NULL,2))) *
                                   fi_x[t+1, 1, i, , ]
                            + mapply("%*%", replicate(s,t(betam[47:58, bc]),simplify=FALSE), alply(fi_outcomes[t+1, 2:n_outcomes, i, , ], c(NULL,NULL,3))) *   #this is wrong
                                   fi_x[t+1, 1, i, , ]
                            + mapply("%*%", replicate(s,t(betam[59:67, bc]),simplify=FALSE), alply(fi_biom[t+1,  , i, , ], c(NULL,NULL,3))) *
                                   fi_x[t+1, 1, i, , ]
                            + mapply("%*%", replicate(s,t(betam[68:80, bc]),simplify=FALSE), replicate(s,treat[ ,t+1, i],simplify=FALSE)) * 
                                   fi_x[t+1, 1, i, , ]
                            + mapply("%*%", replicate(s,t(betam[81:92, bc]),simplify=FALSE), alply(fi_outcomes[t, 2:n_outcomes, i, , ], c(NULL,NULL,3)))
                            + mapply("%*%", replicate(s,t(betam[93:101, bc]),simplify=FALSE), alply(fi_biom[t,  , i, , ], c(NULL,NULL,3)))
                            + mapply("%*%", replicate(s,t(betam[102:114, bc]),simplify=FALSE), replicate(s,treat[ ,t+1, i],simplify=FALSE))
                            + mapply("%*%", replicate(s,t(betam[115:126, bc]),simplify=FALSE), alply(fi_outcomes[t, 2:n_outcomes, i, , ], c(NULL,NULL,3))) *
                                  fi_x[t+1, 1, i, , ]
                            + mapply("%*%", replicate(s,t(betam[127:135, bc]),simplify=FALSE), alply(fi_biom[t,  , i, , ], c(NULL,NULL,3))) *
                                  fi_x[t+1, 1, i, , ]
                            + mapply("%*%", replicate(s,t(betam[136:148, bc]),simplify=FALSE), replicate(s,treat[ ,t+1, i],simplify=FALSE)) * 
                                  fi_x[t+1, 1, i, , ]
        
        print(mu)
        mu <- as.list(mu)
        var <- lapply(mu, function(j) {exp(varcoeff[1,b]*log(j) + varcoeff[2,b]) })
        mu.var <- mapply(c, mu, var, SIMPLIFY=FALSE)
        biom.t2 <- lapply(mu.var, function(k) {rnorm(rep(1, m), k[1], sqrt(k[2])) }) # this is generating warnings
        fi_biom[t+2,  b, i, , ] <- matrix(unlist(biom.t2), ncol=length(biom.t2))
        
        }
        ## Respeing range of values for biomarkers
        fi_biom[t+2, 1, i, , ] =ifelse(fi_biom[t+2,  1, i, , ]<0 , 0 , ifelse(fi_biom[t+2,  1, i, , ]>60 , 60 , fi_biom[t+2,  1, i, , ] ) )      ## BMI
        fi_biom[t+2, 2, i, , ] =ifelse(fi_biom[t+2,  2, i, , ]<0 , 0 , ifelse(fi_biom[t+2,  2, i, , ]>20 , 20 , fi_biom[t+2,  2, i, , ] ) )      ## A1C
        fi_biom[t+2, 3, i, , ] =ifelse(fi_biom[t+2,  3, i, , ]<0 , 0 , ifelse(fi_biom[t+2,  3, i, , ]>100 , 100 , fi_biom[t+2,  3, i, , ] ) )    ## HDL
        fi_biom[t+2, 4, i, , ] =ifelse(fi_biom[t+2,  4, i, , ]<0 , 0 , ifelse(fi_biom[t+2,  4, i, , ]>200 , 200 , fi_biom[t+2,  4, i, , ] ) )    ## LDL
        fi_biom[t+2, 5, i, , ] =ifelse(fi_biom[t+2,  5, i, , ]<0 , 0 , ifelse(fi_biom[t+2,  5, i, , ]>400 , 400 , fi_biom[t+2,  5, i, , ] ) )    ## CHOL
        fi_biom[t+2, 6, i, , ] =ifelse(fi_biom[t+2,  6, i, , ]<0 , 0 , ifelse(fi_biom[t+2,  6, i, , ]>600 , 600 , fi_biom[t+2,  6, i, , ] ) )    ## TRIG
        fi_biom[t+2, 7, i, , ] =ifelse(fi_biom[t+2,  7, i, , ]<90 , 90 , ifelse(fi_biom[t+2,  7, i, , ]>250 , 250 , fi_biom[t+2,  7, i, , ] ) )  ## SBP
        fi_biom[t+2, 8, i, , ] =ifelse(fi_biom[t+2,  8, i, , ]<60 , 60 , ifelse(fi_biom[t+2,  8, i, , ]>140 , 140 , fi_biom[t+2,  8, i, , ] ) )  ## DBP
        fi_biom[t+2, 9, i, , ] =ifelse(fi_biom[t+2,  9, i, , ]<0 , 0 , ifelse(fi_biom[t+2,  9, i, , ]>100 , 100 , fi_biom[t+2,  9, i, , ] ) )    ## EGFR
        
        
       ## OTHER OUTCOMES UPDATE
       ## myo angina stroke hypo 
       for (b in 11:14){
          bc <- b
          pr <- betam[149, bc][rep(1,each=2)] +
                          + mapply("%*%", replicate(s,t(betam[1:10, bc]),simplify=FALSE), alply(fi_x[t+1,  , i, , ], c(NULL,NULL,3))) 
                          + mapply("%*%", replicate(s,t(betam[11, bc]),simplify=FALSE), alply(fi_outcomes[t+2, 1, i, , ], c(NULL,2)))  
                          + mapply("%*%", replicate(s,t(betam[12:23, bc]),simplify=FALSE), alply(fi_outcomes[t+1, 2:n_outcomes, i, , ], c(NULL,NULL,3)))
                          + mapply("%*%", replicate(s,t(betam[24:32, bc]),simplify=FALSE), alply(fi_biom[t+1,  , i, , ], c(NULL,NULL,3)))
                          + mapply("%*%", replicate(s,t(betam[33:45, bc]),simplify=FALSE), replicate(s,treat[ ,t+1, i],simplify=FALSE))
                          + mapply("%*%", replicate(s,t(betam[46, bc]),simplify=FALSE), alply(fi_outcomes[t+2, 1, i, , ], c(NULL,2))) *
                                fi_x[t+1, 1, i, , ]
                          + mapply("%*%", replicate(s,t(betam[47:58, bc]),simplify=FALSE), alply(fi_outcomes[t+1, 2:n_outcomes, i, , ], c(NULL,NULL,3))) *
                                fi_x[t+1, 1, i, , ]
                          + mapply("%*%", replicate(s,t(betam[59:67, bc]),simplify=FALSE), alply(fi_biom[t+1,  , i, , ], c(NULL,NULL,3))) *
                                fi_x[t+1, 1, i, , ]
                          + mapply("%*%", replicate(s,t(betam[68:80, bc]),simplify=FALSE), replicate(s,treat[ ,t+1, i],simplify=FALSE)) * 
                                fi_x[t+1, 1, i, , ]
                          + mapply("%*%", replicate(s,t(betam[81:92, bc]),simplify=FALSE), alply(fi_outcomes[t, 2:n_outcomes, i, , ], c(NULL,NULL,3)))
                          + mapply("%*%", replicate(s,t(betam[93:101, bc]),simplify=FALSE), alply(fi_biom[t, , i, , ], c(NULL,NULL,3)))
                          + mapply("%*%", replicate(s,t(betam[102:114, bc]),simplify=FALSE), replicate(s,treat[ ,t+1, i],simplify=FALSE))
                          + mapply("%*%", replicate(s,t(betam[115:126, bc]),simplify=FALSE), alply(fi_outcomes[t, 2:n_outcomes, i, , ], c(NULL,NULL,3))) *
                                fi_x[t+1, 1, i, , ]
                          + mapply("%*%", replicate(s,t(betam[127:135, bc]),simplify=FALSE), alply(fi_biom[t,  , i, , ], c(NULL,NULL,3))) *
                                fi_x[t+1, 1, i, , ]
                          + mapply("%*%", replicate(s,t(betam[136:148, bc]),simplify=FALSE), replicate(s,treat[ ,t+1, i],simplify=FALSE)) * 
                                fi_x[t+1, 1, i, , ]
          
          print(pr)
          pr = ifelse(pr<0, 0, ifelse(pr>1, 1, pr))
          fi_outcomes_b9 <- lapply(pr, function(j) { rbern(rep(1, m), j) } )
          fi_outcomes[t+2,  b-9, i, , ] <- matrix(unlist(fi_outcomes_b9), ncol=length(fi_outcomes_b9))
          fi_outcomes[t+2,  b-5, i, , ] <- ifelse(fi_outcomes[t+1, b-5, i, , ]==0, fi_outcomes[t+2,  b-9, i, , ], fi_outcomes[t+1, b-5, i, , ])  ## history
          
        }
        

        ## chf lea eye esrd
        for (b in 15:18){
          bc <- b
          pr <- betam[149, bc][rep(1,each=2)] + 
                        + mapply("%*%", replicate(s,t(betam[1:10, bc]),simplify=FALSE), alply(fi_x[t+1,  , i, , ], c(NULL,NULL,3))) 
                        + mapply("%*%", replicate(s,t(betam[11, bc]),simplify=FALSE), alply(fi_outcomes[t+2, 1, i, , ], c(NULL,2)))  
                        + mapply("%*%", replicate(s,t(betam[12:23, bc]),simplify=FALSE), alply(fi_outcomes[t+1, 2:n_outcomes, i, , ], c(NULL,NULL,3))) 
                        + mapply("%*%", replicate(s,t(betam[24:32, bc]),simplify=FALSE), alply(fi_biom[t+1,  , i, , ], c(NULL,NULL,3))) 
                        + mapply("%*%", replicate(s,t(betam[33:45, bc]),simplify=FALSE), replicate(s,treat[ ,t+1, i],simplify=FALSE)) 
                        + mapply("%*%", replicate(s,t(betam[46, bc]),simplify=FALSE), alply(fi_outcomes[t+2, 1, i, , ], c(NULL,2))) *
                              fi_x[t+1, 1, i, , ] 
                        + mapply("%*%", replicate(s,t(betam[47:58, bc]),simplify=FALSE), alply(fi_outcomes[t+1, 2:n_outcomes, i, , ], c(NULL,NULL,3))) *
                              fi_x[t+1, 1, i, , ] 
                        + mapply("%*%", replicate(s,t(betam[59:67, bc]),simplify=FALSE), alply(fi_biom[t+1,  , i, , ], c(NULL,NULL,3))) *
                              fi_x[t+1, 1, i, , ] 
                        + mapply("%*%", replicate(s,t(betam[68:80, bc]),simplify=FALSE), replicate(s,treat[ ,t+1, i],simplify=FALSE)) * 
                              fi_x[t+1, 1, i, , ] 
                        + mapply("%*%", replicate(s,t(betam[81:92, bc]),simplify=FALSE), alply(fi_outcomes[t, 2:n_outcomes, i, , ], c(NULL,NULL,3))) 
                        + mapply("%*%", replicate(s,t(betam[93:101, bc]),simplify=FALSE), alply(fi_biom[t, , i, , ], c(NULL,NULL,3))) 
                        + mapply("%*%", replicate(s,t(betam[102:114, bc]),simplify=FALSE), replicate(s,treat[ ,t+1, i],simplify=FALSE)) 
                        + mapply("%*%", replicate(s,t(betam[115:126, bc]),simplify=FALSE), alply(fi_outcomes[t, 2:n_outcomes, i, , ], c(NULL,NULL,3))) *
                              fi_x[t+1, 1, i, , ] 
                        + mapply("%*%", replicate(s,t(betam[127:135, bc]),simplify=FALSE), alply(fi_biom[t, , i, , ], c(NULL,NULL,3))) *
                              fi_x[t+1, 1, i, , ] 
                        + mapply("%*%", replicate(s,t(betam[136:148, bc]),simplify=FALSE), replicate(s,treat[ ,t+1, i],simplify=FALSE)) * 
                              fi_x[t+1, 1, i, , ]
          
          print(pr)
          pr = ifelse(pr<0, 0, ifelse(pr>1, 1, pr))
          pri = lapply(pr, function(j) { rbern(rep(1, m), j) } )
          outcome.t2b <- pmap(list(alply(fi_outcomes[t+1, b-5, i, , ], 2), pri), function(x,y) { ifelse(x==0, y, x) })
          fi_outcomes[t+2,  b-5, i, , ] <- matrix(unlist(outcome.t2b), ncol=length(outcome.t2b))
          
        
        }
        
        
        
      }
    
  }  
  
  print(apply(fi_biom, c(1,2,3),mean))
  print(apply(fi_outcomes, c(1,2,3),mean))

}


##### Run the model #######

betam <- read_dta("Data/finalbetamatrix.dta")
simuldataforR <- read_dta("Data/siumldataforR.dta")

base_x = rbind(simuldataforR[1, 4], c(NA), simuldataforR[2:9, 4] )
base_x = t(base_x)
base_outcomes = simuldataforR[10:22, 3:4]
base_biom = simuldataforR[23:31, 3:4]
horizon = 10

treatment = simuldataforR[32:44, 3:4]
treatment = cbind(treatment, matrix( rep(t((treatment[,2])), horizon), nrow(treatment), horizon))

# model run
rapids(base_biom, base_outcomes, base_x, treatment, betam, horizon) 
  
  