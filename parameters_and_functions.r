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

# set root folder
rootDir = "C:/Users/Luda/OneDrive - Johns Hopkins/J1798/"


#==============================
# Functions
#==============================

geti = function(x,i) {x[i]}

# paired Wilcoxon test
# k is a column in dat
pairedWilcox = function(k, dat)
{
        mat <- data.frame(matrix(nrow = length(unique(dat$patientID)),
      ncol = 3, dimnames = list(unique(dat$patientID),
                      levels(factor(dat$timepoint)))), check.names = F)
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



