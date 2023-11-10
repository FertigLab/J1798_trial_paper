# set up for J1798 analysis
# functions and parameters
# to use the same colors and annotation in the figures

# labels for time points
timeLab  = c("Baseline", "C1D1", "C2D1")

# p-value threshold
pVal = 0.05

# set colors
recistCol = c("NR" = "#F8766D","R" = "#619CFF", "NonResponder" = "#F8766D","Responder" = "#619CFF")
timeCol = setNames(c("black", 'grey45', 'grey75'), timeLab)


#==============================
# Functions
#==============================

geti = function(x,i) {x[i]}

# paired Wilcoxon test
# k is a column in dat
pairedWilcox = function(k, dat)
{
  # create matrix patients vs time points  
  mat <- data.frame(matrix(nrow = length(unique(dat$patientID)),
      ncol = length(unique(dat$timepoint)),
      dimnames = list(unique(dat$patientID), levels(factor(dat$timepoint)))), 
      check.names = F)
  
    for(i in levels(factor(dat$timepoint)))
    {
      d1 = dat %>% filter(timepoint == i)
      mat[d1$patientID, i] = d1[,k]
    }
   #     print(as.numeric(mat[,1]))
      
    # running paired test
    # Baseline vs C1D1
    p1 = wilcox.test(as.numeric(mat[,1]), as.numeric(mat[,2]), paired = T)$p.value
    #  C1D1 vs  C2D1
    p2 = wilcox.test(as.numeric(mat[,2]), as.numeric(mat[,3]), paired = T)$p.value
    #  Baseline vs  C2D1
    p3 = wilcox.test(as.numeric(mat[,1]), as.numeric(mat[,3]), paired = T)$p.value
    
    # mean of delta between time points
    d1 = mean(as.numeric(mat[,2]) -  as.numeric(mat[,1]), na.rm = T)
    d2 = mean(as.numeric(mat[,3]) -  as.numeric(mat[,2]), na.rm = T)
    d3 = mean(as.numeric(mat[,3]) -  as.numeric(mat[,1]), na.rm = T)
   

    return(c(pValue_BvsC1D1 = p1, 'mean(C1D1-B)' = d1, pValue_C1D1vsC2D1 = p2, 
             'mean(C2D1-C1D1)' = d2, pValue_BvsC2D2 = p3,'mean(C2D1-B)' = d3))
}

# print the number of pairs between time points
# the input table should have patientID and timepoint columns
printPairs = function(dat)
{
  library(Matrix)
  mat = sparseMatrix(i = as.numeric(as.factor(dat$patientID)), 
                     j = as.numeric(as.factor(dat$timepoint)), x = 1, 
                     dimnames = list(levels(as.factor(dat$patientID)),
                                     levels(as.factor(dat$timepoint))))
  mat = as.data.frame(as.matrix(mat))
  
  cat("The number of paired samples between Baseline and C1D1:",
      mat %>% filter(Baseline == 1 & C1D1 == 1) %>% nrow,"\n")
  cat("The number of paired samples between Baseline and C2D1:",
      mat %>% filter(Baseline == 1 & C2D1 == 1) %>% nrow,"\n")
  cat("The number of paired samples between C2D1 and C1D1:",
      mat %>% filter(C2D1 == 1 & C1D1 == 1) %>% nrow)
  
}

#### functions to summarize character and numeric variables for 
#### table of patient characteristics
charSummary=function(col,data){
  x=data[,col]
  tb=table(x)
  vTb=as.vector(tb)
  names(vTb)=names(tb)
  na=sum(is.na(x))
  if(na>0) vTb=c(vTb,"NA"=na)
  return(vTb)
}

dateSummary=function(x){
  dateRange=as.character(range(x,na.rm=T))
  na.ct=sum(is.na(x))
  
  return(c(dateRange,na.ct))
}


numSummary=function(col,data){
  x=data[,col]
  med=round(median(x,na.rm=T),3)
  rng=round(range(x,na.rm=T),3)
  na.ct=sum(is.na(x))
  ans=c(col,med,paste0(rng,collapse="-"))
  if(na.ct>0) ans=c(ans,"NA"=na.ct)
  names(ans)[1:3]=c("Characteristic","Median","Range")
  return(ans)
}
### function to format character summary
modChar=function(name,dat=charSums){
  characteristic=c(name,rep("",length(dat[[name]])))
  vals=c("",names(dat[[name]]))
  counts=c("",dat[[name]])
  data=data.frame(characteristic,vals,counts)
  colnames(data)=c("Characteristic","Median","Range")
  rownames(data)=NULL
  return(data)
}
