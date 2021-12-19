#k-steps variable selection scheme

library(gamlss)

#Possible distributions: Normal, Gamma, Inverse_Gaussian.
#Normal: (NO)
#Gamma: (GA)
#Inverse_Gaussian: (IG)


k_Steps_2_mod = function(Y, data, distribution, prop_models){
  set.seed(132)
  
  #number_of_steps[1] : verification 
  #number_of_steps[2] : current count 
  number_of_steps_AIC = numeric(2)
  number_of_steps_AICc = numeric(2)
  number_of_steps_BIC = numeric(2)
  number_of_steps_HQC = numeric(2)
  
  #vector_throughs: vector of count of steps for each selected model 
  vector_throughs_AIC = numeric(1)
  vector_throughs_AICc = numeric(1)
  vector_throughs_BIC = numeric(1)
  vector_throughs_HQC = numeric(1)
  
  
  
  samplesize = length(data[,1])
  Logsamplesize = log(samplesize)
  LogLogsamplesize = log(Logsamplesize)
  
  p = ncol(data)
  a = numeric()
  M = matrix()
  
  for(i in 1:p){
    a = c(a,(rep(c(0,1), each = 2^(p-i), times = 2^(i-1))))  
  }
  
  M_1 = matrix(t(a), 2^(p),p)
  
  M = M_1[-1, ]
  
  var_names = c("winnings", "age", "drive", "accuracy", "regulation", "putts", "save", "events")
  
  #Average
  ModMi = character(nrow(M)+1)
  
  ModMi_mod = character(nrow(M)+1)
  
  ModMi[1] = 'data1[,1] ~ 1'
  ModMi_mod[1] = 'winnings ~ 1'
    
  for (i in 1:nrow(M)){
    vetorcarac = character(ncol(M))
    vetorcarac_mod = character(ncol(M))
    
    if (M[i,1] == 1){
      vetorcarac[1] = 'data1[, 2 ]'
      vetorcarac_mod[1] = 'age'
    }
    
    for (j in 2:ncol(M)){
      
      if (M[i,j] == 1 && vetorcarac [j-1] == ''){
        vetorcarac[j] = paste('data1[,',j+1,']')
        vetorcarac_mod[j] = paste(var_names[j+1])
      }
      
      else if (M[i,j] == 1 && vetorcarac [j-1] != ''){
        vetorcarac[j] = paste(vetorcarac[j-1],'+ data1[,',j+1,']')
        vetorcarac_mod[j] = paste(vetorcarac_mod[j-1],'+',var_names[j+1])
      }
      
      else{
        vetorcarac[j] = vetorcarac[j-1]
        vetorcarac_mod[j] = vetorcarac_mod[j-1]
      }
      
    }
    ModMi[i+1] = paste('data1[,1] ~',vetorcarac[ncol(M)])
    ModMi_mod[i+1] = paste('winnings ~',vetorcarac_mod[ncol(M)])
  }
  
  
  np = c(0,apply(M,1,sum)) + 1
  
  npmi = c(rep(np, each = 1, times = 1))
  
  mi = c(rep(ModMi, each = 1, times = 1))
  
  #######################################################################
  
  #variance 
  ModSigma = character(nrow(M)+1)
  ModSigma_mod = character(nrow(M)+1)
  
  ModSigma[1] = '~ 1'
  ModSigma_mod[1] = '~ 1'
  
  for (i in 1:nrow(M)){
    vetorcarac = character(ncol(M))
    vetorcarac_mod = character(ncol(M))
    
    if (M[i,1] == 1){
      vetorcarac[1] = 'data1[, 2 ]'
      vetorcarac_mod[1] = 'age'
    }
    
    for (j in 2:ncol(M)){
      
      if (M[i,j] == 1 && vetorcarac [j-1] == ''){
        vetorcarac[j] = paste('data1[,',j+1,']')
        vetorcarac_mod[j] = paste(var_names[j+1])
      }
      
      else if (M[i,j] == 1 && vetorcarac [j-1] != ''){
        vetorcarac[j] = paste(vetorcarac[j-1],'+ data1[,',j+1,']')
        vetorcarac_mod[j] = paste(vetorcarac_mod[j-1],'+',var_names[j+1])
      }
      
      else{
        vetorcarac[j] = vetorcarac[j-1]
        vetorcarac_mod[j] = vetorcarac_mod[j-1]
      }
      
    }
    ModSigma[i+1] = paste('~',vetorcarac[ncol(M)]) 
    ModSigma_mod[i+1] = paste('~',vetorcarac_mod[ncol(M)])
  }
  
  npsigma = c(rep(np, each = 1, times = 1))
  
  sigma = c(rep(ModSigma, each = 1, times = 1))
  
  #######################################################################
  #Family for the models
  if(distribution  == "Normal"){
    familia = "NO"
  }
  
  if(distribution  == "Gamma"){
    familia = "GA"
  }
  
  if(distribution  == "Inverse_Gaussian"){
    familia = "IG"
  }
  
  #######################################################################
  
  number_of_steps_AIC[1] = 0
  number_of_steps_AICc[1] = 0
  number_of_steps_BIC[1] = 0
  number_of_steps_HQC[1] = 0
  number_of_steps_AIC[2] = 2
  number_of_steps_AICc[2] = 2
  number_of_steps_BIC[2] = 2
  number_of_steps_HQC[2] = 2
  
  
  data1 = cbind(Y, data)
  
  Deviance = numeric(length(mi))
  
  
  for (i in 1:(length(mi))){
    Deviance[i] = gamlss(as.formula(mi[i]), ~1, family=familia, trace = FALSE,  n.cyc = 1000, method = CG())$G.deviance
  }
  
  
  AIC = numeric(length(mi))
  BIC = numeric(length(mi))
  AICc = numeric(length(mi))
  HQC = numeric(length(mi))
  for(i in 1:length(Deviance)){
    
    AIC[i] = Deviance[i] + 2*(npmi[i] + 1)
    
    AICc[i] = Deviance[i] + 2*(npmi[i] + 1)*((samplesize)/(samplesize - (npmi[i] + 1) - 1))
    
    BIC[i] = Deviance[i] + Logsamplesize*(npmi[i] + 1)
    
    HQC[i] = Deviance[i] + 2*(npmi[i] + 1)*LogLogsamplesize
  }
  Pos1 = rank(AIC) 
  Pos2 = rank(AICc)
  Pos3 = rank(BIC)
  Pos4 = rank(HQC)
  
  
  for(i in 1:length(mi)){
    if(Pos1[i] == 1){
      Pos1_final = i
    }
    if(Pos2[i] == 1){
      Pos2_final = i
    }
    if(Pos3[i] == 1){
      Pos3_final = i
    }
    if(Pos4[i] == 1){
      Pos4_final = i
    }
  }
  
  
  #variance
  Deviance_AIC = numeric(length(mi))
  Deviance_AICc = numeric(length(mi))
  Deviance_BIC = numeric(length(mi))
  Deviance_HQC = numeric(length(mi))
  
  
  for (i in 1:(length(mi))){
    Deviance_AIC[i] = gamlss(as.formula(mi[Pos1_final[1]]), as.formula(sigma[i]), family=familia, trace = FALSE,  n.cyc = 1000, method = CG())$G.deviance
    Deviance_AICc[i] = gamlss(as.formula(mi[Pos2_final[1]]), as.formula(sigma[i]), family=familia, trace = FALSE,  n.cyc = 1000, method = CG())$G.deviance
    Deviance_BIC[i] = gamlss(as.formula(mi[Pos3_final[1]]), as.formula(sigma[i]), family=familia, trace = FALSE,  n.cyc = 1000, method = CG())$G.deviance
    Deviance_HQC[i] = gamlss(as.formula(mi[Pos4_final[1]]), as.formula(sigma[i]), family=familia, trace = FALSE,  n.cyc = 1000, method = CG())$G.deviance
  }
  
  
  #Weights for variance
  pesos_AIC2_sigma = rep(1,length(npsigma))
  pesos_AICc2_sigma = rep(1,length(npsigma))
  pesos_BIC2_sigma = rep(1,length(npsigma))
  pesos_HQC2_sigma = rep(1,length(npsigma))
  
  for(j in 1:length(npsigma)){
    for(i in 1:ncol(data)){
      
      #AIC
      if(M_1[Pos1_final[1],i] == 1 && M_1[j,i] == 1){
        pesos_AIC2_sigma[j] = pesos_AIC2_sigma[j] + 1
      }
      
      if(M_1[Pos1_final[1],i] == 0 && M_1[j,i] == 1){
        pesos_AIC2_sigma[j] = pesos_AIC2_sigma[j] + 2
      } 
      
      #AICc
      if(M_1[Pos2_final[1],i] == 1 && M_1[j,i] == 1){
        pesos_AICc2_sigma[j] = pesos_AICc2_sigma[j] + 1
      }
      
      if(M_1[Pos2_final[1],i] == 0 && M_1[j,i] == 1){
        pesos_AICc2_sigma[j] = pesos_AICc2_sigma[j] + 2
      }
      
      #BIC
      if(M_1[Pos3_final[1],i] == 1 && M_1[j,i] == 1){
        pesos_BIC2_sigma[j] = pesos_BIC2_sigma[j] + 1
      }
      
      if(M_1[Pos3_final[1],i] == 0 && M_1[j,i] == 1){
        pesos_BIC2_sigma[j] = pesos_BIC2_sigma[j] + 2
      } 
      
      #HQC
      if(M_1[Pos4_final[1],i] == 1 && M_1[j,i] == 1){
        pesos_HQC2_sigma[j] = pesos_HQC2_sigma[j] + 1
      }
      
      if(M_1[Pos4_final[1],i] == 0 && M_1[j,i] == 1){
        pesos_HQC2_sigma[j] = pesos_HQC2_sigma[j] + 2
      } 
    }
  }
  
  AIC2 = numeric(length(mi))
  BIC2 = numeric(length(mi))
  AICc2 = numeric(length(mi))
  HQC2 = numeric(length(mi))
  
  for(i in 1:length(Deviance)){
    
    AIC2[i] = Deviance_AIC[i] + 2*(npmi[Pos1_final[1]] + pesos_AIC2_sigma[i])
    
    AICc2[i] = Deviance_AICc[i] + 2*(npmi[Pos2_final[1]] + pesos_AICc2_sigma[i])*((samplesize)/(samplesize - (npmi[Pos2_final[1]] + pesos_AICc2_sigma[i]) - 1))
    
    BIC2[i] = Deviance_BIC[i] + Logsamplesize*(npmi[Pos3_final[1]] + pesos_BIC2_sigma[i])
    
    HQC2[i] = Deviance_HQC[i] + 2*(npmi[Pos4_final[1]] + pesos_HQC2_sigma[i])*LogLogsamplesize
  }
  
  Pos12 = rank(AIC2) 
  Pos22 = rank(AICc2)
  Pos32 = rank(BIC2)
  Pos42 = rank(HQC2)
  
  
  for(i in 1:length(sigma)){
    if(Pos12[i] == 1){
      Pos12_final = i
    }
    if(Pos22[i] == 1){
      Pos22_final = i
    }
    if(Pos32[i] == 1){
      Pos32_final = i
    }
    if(Pos42[i] == 1){
      Pos42_final = i
    }
  }
  
  
  #AIC
  while(number_of_steps_AIC[1] != 1){ 
    
    Deviance_AIC_mi = numeric(length(mi))
    
    for (i in 1:(length(mi))){
      Deviance_AIC_mi[i] = gamlss(as.formula(mi[i]), as.formula(sigma[Pos12_final[1]]), family=familia, trace = FALSE,  n.cyc = 1000, method = CG())$G.deviance
    }
    
    AIC_mi = numeric(length(mi)) 
    
    for(i in 1:length(Deviance_AIC_mi)){
      AIC_mi[i] = Deviance_AIC_mi[i] + 2*(npmi[i] + npsigma[Pos12_final[1]])
    }
    
    Pos111 = rank(AIC_mi) 
    
    for(i in 1:length(mi)){
      if(Pos111[i] == 1){
        Pos11_final = i
      }
    }
    
    number_of_steps_AIC[2] = number_of_steps_AIC[2] + 1
    
    #Ticket verification and check if the model is the same as the previous one - mi
    if(Pos1_final == Pos11_final){
      Pos122 = Pos12
      break
    }
    
    
    if((number_of_steps_AIC[2]*(2^ncol(data))) >= ((2^(2*ncol(data))))*prop_models){
      Pos122 = Pos12
      break
    }
    
    #variance
    Deviance_AIC_sigma = numeric(length(mi))
    
    for (i in 1:(length(mi))){
      Deviance_AIC_sigma[i] = gamlss(as.formula(mi[Pos11_final[1]]), as.formula(sigma[i]), family=familia, trace = FALSE,  n.cyc = 1000, method = CG())$G.deviance
    }
    
    
    #Weights for variance
    pesos_AIC22_sigma = rep(1,length(npsigma))
    
    for(j in 1:length(npsigma)){
      for(i in 1:ncol(data)){
        
        #AIC
        if(M_1[Pos11_final[1],i] == 1 && M_1[j,i] == 1){
          pesos_AIC22_sigma[j] = pesos_AIC22_sigma[j] + 1
        }
        
        if(M_1[Pos11_final[1],i] == 0 && M_1[j,i] == 1){
          pesos_AIC22_sigma[j] = pesos_AIC22_sigma[j] + 2
        } 
      }
    }
    
    
    AIC2_sigma = numeric(length(mi))
    
    for(i in 1:length(Deviance)){
      AIC2_sigma[i] = Deviance_AIC_sigma[i] + 2*(npmi[Pos11_final[1]] + pesos_AIC22_sigma[i])
    }
    
    Pos122 = rank(AIC2_sigma) 
    
    for(i in 1:length(sigma)){
      if(Pos122[i] == 1){
        Pos122_final = i
      }
    }
    
    number_of_steps_AIC[2] = number_of_steps_AIC[2] + 1
    
    #Ticket verification and check if the model is the same as the previous one  - sigma
    if(Pos12_final == Pos122_final){
      break
    }
    
    Pos1_final = Pos11_final
    Pos12_final = Pos122_final 
    
    if((number_of_steps_AIC[2]*(2^ncol(data))) >= ((2^(2*ncol(data))))*prop_models){
      break
    }
  }
  
  vector_throughs_AIC = number_of_steps_AIC[2]
  
  
  #AICc
  while(number_of_steps_AICc[1] != 1){ 
    
    Deviance_AICc_mi = numeric(length(mi))
    
    for (i in 1:(length(mi))){
      Deviance_AICc_mi[i] = gamlss(as.formula(mi[i]), as.formula(sigma[Pos22_final[1]]), family=familia, trace = FALSE,  n.cyc = 1000, method = CG())$G.deviance
    }
    
    AICc_mi = numeric(length(mi))
    
    for(i in 1:length(Deviance_AICc_mi)){
      AICc_mi[i] = Deviance_AICc_mi[i] + 2*(npmi[i] + npsigma[Pos22_final[1]])*((samplesize)/(samplesize - (npmi[i] + npsigma[Pos22_final[1]]) - 1))
    }
    
    Pos211 = rank(AICc_mi)
    
    for(i in 1:length(mi)){
      if(Pos211[i] == 1){
        Pos21_final = i
      }
    }
    
    number_of_steps_AICc[2] = number_of_steps_AICc[2] + 1
    
    #Ticket verification and check if the model is the same as the previous one - mi
    if(Pos2_final == Pos21_final){
      Pos222 = Pos22
      break
    }
    
    if((number_of_steps_AICc[2]*(2^ncol(data))) >= ((2^(2*ncol(data))))*prop_models){
      Pos222 = Pos22
      break
    }
    
    #variance
    Deviance_AICc_sigma = numeric(length(mi))  
    
    for (i in 1:(length(mi))){
      Deviance_AICc_sigma[i] = gamlss(as.formula(mi[Pos21_final[1]]), as.formula(sigma[i]), family=familia, trace = FALSE,  n.cyc = 1000, method = CG())$G.deviance
    }
    
    #Weights for variance
    pesos_AICc22_sigma = rep(1,length(npsigma))
    
    for(j in 1:length(npsigma)){
      for(i in 1:ncol(data)){
        
        #AICc
        if(M_1[Pos21_final[1],i] == 1 && M_1[j,i] == 1){
          pesos_AICc22_sigma[j] = pesos_AICc22_sigma[j] + 1
        }
        
        if(M_1[Pos21_final[1],i] == 0 && M_1[j,i] == 1){
          pesos_AICc22_sigma[j] = pesos_AICc22_sigma[j] + 2
        } 
      }
    }
    
    AICc2_sigma = numeric(length(mi))
    
    for(i in 1:length(Deviance)){
      AICc2_sigma[i] = Deviance_AICc_sigma[i] + 2*(npmi[Pos21_final[1]] + pesos_AICc22_sigma[i])*((samplesize)/(samplesize - (npmi[Pos21_final[1]] + pesos_AICc22_sigma[i]) - 1))
    }
    
    Pos222 = rank(AICc2_sigma)
    
    for(i in 1:length(sigma)){
      if(Pos222[i] == 1){
        Pos222_final = i
      }
    }
    
    number_of_steps_AICc[2] = number_of_steps_AICc[2] + 1
    
    #Ticket verification and check if the model is the same as the previous one - sigma
    if(Pos22_final == Pos222_final){
      break
    }
    
    Pos2_final = Pos21_final 
    Pos22_final = Pos222_final 
    
    if((number_of_steps_AICc[2]*(2^ncol(data))) >= ((2^(2*ncol(data))))*prop_models){
      break
    }
  }
  
  vector_throughs_AICc = number_of_steps_AICc[2]
  
  
  #BIC
  while(number_of_steps_BIC[1] != 1){ 
    
    Deviance_BIC_mi = numeric(length(mi))
    
    for (i in 1:(length(mi))){
      Deviance_BIC_mi[i] = gamlss(as.formula(mi[i]), as.formula(sigma[Pos32_final[1]]), family=familia, trace = FALSE,  n.cyc = 1000, method = CG())$G.deviance
    }
    
    BIC_mi = numeric(length(mi))   
    
    for(i in 1:length(Deviance_BIC_mi)){
      BIC_mi[i] = Deviance_BIC_mi[i] + Logsamplesize*(npmi[i] + npsigma[Pos32_final[1]])
    }
    
    Pos311 = rank(BIC_mi)
    
    for(i in 1:length(mi)){
      if(Pos311[i] == 1){
        Pos31_final = i
      }
    }
    
    
    number_of_steps_BIC[2] = number_of_steps_BIC[2] + 1
    
    #Ticket verification and check if the model is the same as the previous one - mi
    if(Pos3_final == Pos31_final){
      Pos322 = Pos32
      break
    }
    
    if((number_of_steps_BIC[2]*(2^ncol(data))) >= ((2^(2*ncol(data))))*prop_models){
      Pos322 = Pos32
      break
    }
    
    #variance    
    Deviance_BIC_sigma = numeric(length(mi))
    
    for (i in 1:(length(mi))){
      Deviance_BIC_sigma[i] = gamlss(as.formula(mi[Pos31_final[1]]), as.formula(sigma[i]), family=familia, trace = FALSE,  n.cyc = 1000, method = CG())$G.deviance
    }
    
    #Weights for variance
    pesos_BIC22_sigma = rep(1,length(npsigma))
    
    for(j in 1:length(npsigma)){
      for(i in 1:ncol(data)){
        
        #BIC
        if(M_1[Pos31_final[1],i] == 1 && M_1[j,i] == 1){
          pesos_BIC22_sigma[j] = pesos_BIC22_sigma[j] + 1
        }
        
        if(M_1[Pos31_final[1],i] == 0 && M_1[j,i] == 1){
          pesos_BIC22_sigma[j] = pesos_BIC22_sigma[j] + 2
        } 
      }
    }
    
    BIC2_sigma = numeric(length(mi))
    
    for(i in 1:length(Deviance)){
      BIC2_sigma[i] = Deviance_BIC_sigma[i] + Logsamplesize*(npmi[Pos31_final[1]] + pesos_BIC22_sigma[i])
    }
    
    Pos322 = rank(BIC2_sigma)    
    
    for(i in 1:length(sigma)){
      if(Pos322[i] == 1){
        Pos322_final = i
      }
    }
    
    
    number_of_steps_BIC[2] = number_of_steps_BIC[2] + 1
    
    #Ticket verification and check if the model is the same as the previous one - sigma
    if(Pos32_final == Pos322_final){
      break
    }
    
    Pos3_final = Pos31_final 
    Pos32_final = Pos322_final
    
    if((number_of_steps_BIC[2]*(2^ncol(data))) >= ((2^(2*ncol(data))))*prop_models){
      break
    }
  }
  
  vector_throughs_BIC = number_of_steps_BIC[2]
  
  
  #HQC
  while(number_of_steps_HQC[1] != 1){ 
    
    Deviance_HQC_mi = numeric(length(mi))
    
    for (i in 1:(length(mi))){
      Deviance_HQC_mi[i] = gamlss(as.formula(mi[i]), as.formula(sigma[Pos42_final[1]]), family=familia, trace = FALSE,  n.cyc = 1000, method = CG())$G.deviance
    }  
    
    HQC_mi = numeric(length(mi))
    
    for(i in 1:length(Deviance_HQC_mi)){      
      HQC_mi[i] = Deviance_HQC_mi[i] + 2*(npmi[i] + npsigma[Pos42_final[1]])*LogLogsamplesize
    }
    
    Pos411 = rank(HQC_mi)
    
    for(i in 1:length(mi)){
      if(Pos411[i] == 1){
        Pos41_final = i
      }
    }
    
    
    number_of_steps_HQC[2] = number_of_steps_HQC[2] + 1
    
    #Ticket verification and check if the model is the same as the previous one - mi
    if(Pos4_final == Pos41_final){
      Pos422 = Pos42
      break
    }
    
    if((number_of_steps_HQC[2]*(2^ncol(data))) >= ((2^(2*ncol(data))))*prop_models){
      Pos422 = Pos42
      break
    }
    
    
    #variance
    Deviance_HQC_sigma = numeric(length(mi))
    
    for (i in 1:(length(mi))){
      Deviance_HQC_sigma[i] = gamlss(as.formula(mi[Pos41_final[1]]), as.formula(sigma[i]), family=familia, trace = FALSE,  n.cyc = 1000, method = CG())$G.deviance
    }
    
    #Weights for variance
    pesos_HQC22_sigma = rep(1,length(npsigma))
    
    for(j in 1:length(npsigma)){
      for(i in 1:ncol(data)){
        
        #HQC
        if(M_1[Pos41_final[1],i] == 1 && M_1[j,i] == 1){
          pesos_HQC22_sigma[j] = pesos_HQC22_sigma[j] + 1
        }
        
        if(M_1[Pos41_final[1],i] == 0 && M_1[j,i] == 1){
          pesos_HQC22_sigma[j] = pesos_HQC22_sigma[j] + 2
        } 
      }
    }
    
    HQC2_sigma = numeric(length(mi))
    
    for(i in 1:length(Deviance)){
      HQC2_sigma[i] = Deviance_HQC_sigma[i] + 2*(npmi[Pos41_final[1]] + pesos_HQC22_sigma[i])*LogLogsamplesize
    }
    
    Pos422 = rank(HQC2_sigma)
    
    for(i in 1:length(sigma)){
      if(Pos422[i] == 1){
        Pos422_final = i
      }
    }
    
    
    number_of_steps_HQC[2] = number_of_steps_HQC[2] + 1
    
    #Ticket verification and check if the model is the same as the previous one - sigma
    if(Pos42_final == Pos422_final){
      break
    }
    
    Pos4_final = Pos41_final 
    Pos42_final = Pos422_final 
    
    if((number_of_steps_HQC[2]*(2^ncol(data))) >= ((2^(2*ncol(data))))*prop_models){
      break
    }
  }
  
  vector_throughs_HQC = number_of_steps_HQC[2]
  
  #################################################################################################
  
  Mod_AIC = c(ModMi_mod[Pos1_final], ModSigma_mod[Pos12_final])
  
  
  Mod_AICc = c(ModMi_mod[Pos2_final], ModSigma_mod[Pos22_final])
  
  
  Mod_BIC = c(ModMi_mod[Pos3_final], ModSigma_mod[Pos32_final])
  
  
  Mod_HQC = c(ModMi_mod[Pos4_final], ModSigma_mod[Pos42_final])
  
  
  return(list("Steps AIC" = vector_throughs_AIC, "Steps AICc" = vector_throughs_AICc, 
              "Steps BIC" = vector_throughs_BIC, "Steps HQC" = vector_throughs_HQC,
              'Mod_AIC' = Mod_AIC, 'Mod_AICc' = Mod_AICc, 'Mod_BIC' = Mod_BIC, 'Mod_HQC' = Mod_HQC))
}






#Application
data = read.table("PGA2004_Ksteps_aplication_variables.txt", header = T)

{
indice = matrix(c(rep(NA, nrow(data)*1)), nrow(data), 1)

set.seed(15)

for(j in 1:1){
  
  indice_aux = sample(1:nrow(data))
  
  for(i in 1:nrow(data)){
    indice[i,j] = indice_aux[i]
  }
}

rep = 1
k = numeric()

EQM_AIC = numeric(rep)
EQM_AICc = numeric(rep)
EQM_BIC = numeric(rep)
EQM_HQC = numeric(rep)

Steps_AIC = numeric(rep)
Steps_AICc = numeric(rep)
Steps_BIC = numeric(rep)
Steps_HQC = numeric(rep)

vars_AIC_mu = matrix(rep(0, 7*rep), 7,rep)
vars_AIC_sigma = matrix(rep(0, 7*rep), 7,rep)

vars_AICc_mu = matrix(rep(0, 7*rep), 7,rep)
vars_AICc_sigma = matrix(rep(0, 7*rep), 7,rep)

vars_BIC_mu = matrix(rep(0, 7*rep), 7,rep)
vars_BIC_sigma = matrix(rep(0, 7*rep), 7,rep)

vars_HQC_mu = matrix(rep(0, 7*rep), 7,rep)
vars_HQC_sigma = matrix(rep(0, 7*rep), 7,rep)
}






#Application study
for(k in 1:rep){

model = k_Steps_2_mod(data[indice[1:(round(nrow(data)*0.7,0)),k],1], data[indice[1:(round(nrow(data)*0.7,0)),k],2:8], "Gamma", 1)

#train=data[indice[1:(round(nrow(data)*0.7,0)),k],]
#test=data[indice[(round(nrow(data)*0.7,0) + 1):(nrow(data)),k],]

#AIC
model_AIC = gamlss(as.formula(model$Mod_AIC[1]),
                     as.formula(model$Mod_AIC[2]), family=GA, method = CG(),
                     trace = FALSE, n.cyc = 1000, data = data[indice[1:(round(nrow(data)*0.7,0)),k],])

EQM_AIC[k] = mean((data[indice[(round(nrow(data)*0.7,0) + 1):(nrow(data)),k],1] - 
                     predict(model_AIC, newdata = data[indice[(round(nrow(data)*0.7,0) + 1):(nrow(data)),k],], type = "response"))^2)

Steps_AIC[k] = model$`Steps AIC`


for(a in 1:length(names(model_AIC$mu.coefficients))-1){
  for(b in 1:7){
    if(names(model_AIC$mu.coefficients)[a+1] == names(data)[b+1]){
      vars_AIC_mu[b,k] = 1
    }
  }
}


for(a in 1:length(names(model_AIC$sigma.coefficients))-1){
  for(b in 1:7){
    if(names(model_AIC$sigma.coefficients)[a+1] == names(data)[b+1]){
      vars_AIC_sigma[b,k] = 1
    }
  }
}



#AICc
model_AICc = gamlss(as.formula(model$Mod_AICc[1]),
                    as.formula(model$Mod_AICc[2]), family=GA, method = CG(),
                    trace = FALSE, n.cyc = 1000, data = data[indice[1:(round(nrow(data)*0.7,0)),k],])

EQM_AICc[k] = mean((data[indice[(round(nrow(data)*0.7,0) + 1):(nrow(data)),k],1] - 
                      predict(model_AICc, newdata = data[indice[(round(nrow(data)*0.7,0) + 1):(nrow(data)),k],], type = "response"))^2)

Steps_AICc[k] = model$`Steps AICc`


for(a in 1:length(names(model_AICc$mu.coefficients))-1){
  for(b in 1:7){
    if(names(model_AICc$mu.coefficients)[a+1] == names(data)[b+1]){
      vars_AICc_mu[b,k] = 1
    }
  }
}


for(a in 1:length(names(model_AICc$sigma.coefficients))-1){
  for(b in 1:7){
    if(names(model_AICc$sigma.coefficients)[a+1] == names(data)[b+1]){
      vars_AICc_sigma[b,k] = 1
    }
  }
}

#BIC
model_BIC = gamlss(as.formula(model$Mod_BIC[1]),
                    as.formula(model$Mod_BIC[2]), family=GA, method = CG(),
                    trace = FALSE, n.cyc = 1000, data = data[indice[1:(round(nrow(data)*0.7,0)),k],])

EQM_BIC[k] = mean((data[indice[(round(nrow(data)*0.7,0) + 1):(nrow(data)),k],1] - 
                     predict(model_BIC, newdata = data[indice[(round(nrow(data)*0.7,0) + 1):(nrow(data)),k],], type = "response"))^2)

Steps_BIC[k] = model$`Steps BIC`


for(a in 1:length(names(model_BIC$mu.coefficients))-1){
  for(b in 1:7){
    if(names(model_BIC$mu.coefficients)[a+1] == names(data)[b+1]){
      vars_BIC_mu[b,k] = 1
    }
  }
}


for(a in 1:length(names(model_BIC$sigma.coefficients))-1){
  for(b in 1:7){
    if(names(model_BIC$sigma.coefficients)[a+1] == names(data)[b+1]){
      vars_BIC_sigma[b,k] = 1
    }
  }
}

#HQC
model_HQC = gamlss(as.formula(model$Mod_HQC[1]),
                    as.formula(model$Mod_HQC[2]), family=GA, method = CG(),
                    trace = FALSE, n.cyc = 1000, data = data[indice[1:(round(nrow(data)*0.7,0)),k],])

EQM_HQC[k] = mean((data[indice[(round(nrow(data)*0.7,0) + 1):(nrow(data)),k],1] - 
                     predict(model_HQC, newdata = data[indice[(round(nrow(data)*0.7,0) + 1):(nrow(data)),k],], type = "response"))^2)

Steps_HQC[k] = model$`Steps HQC`


for(a in 1:length(names(model_HQC$mu.coefficients))-1){
  for(b in 1:7){
    if(names(model_HQC$mu.coefficients)[a+1] == names(data)[b+1]){
      vars_HQC_mu[b,k] = 1
    }
  }
}


for(a in 1:length(names(model_HQC$sigma.coefficients))-1){
  for(b in 1:7){
    if(names(model_HQC$sigma.coefficients)[a+1] == names(data)[b+1]){
      vars_HQC_sigma[b,k] = 1
    }
  }
}

print(k)
}

