source("surrogate_functions.R")

# Data
n = 2 ^ 12
y = arima.sim(n = n, model = list(ar = c(0.9, -0.9)))

# Spectrogram specs
win.length = 2 ^ 9
overlap = 0.75      
win.type = "tukey"  

S = 1000  # Number of surrogates in test

ptm = proc.time()

# Compute spectrogram
F = spectrogram(y, win.length = win.length,
                overlap = overlap, win.type = win.type,
                type = "ar", plot = FALSE)$spec
F.mean = apply(F, 1, mean)
F.global = as.vector(spec.ar(tukeywindow(length(y), r = 0.1) * y, 
                             n.freq = length(F.mean), 
                             plot = FALSE)$spec)

#####
# Compute distances
#####

# 1. Compute KS/KL/LSD distance for each segment vs. the mean
KS.mean = apply(F, 2, KS, p2 = F.mean)
KL.mean = apply(F, 2, KL, p2 = F.mean)
LSD.mean = apply(F, 2, LSD, p2 = F.mean)

# 2. Compute KS/KL/LSD distance for each segment vs. global PSD (same freqs)
KS.global = apply(F, 2, KS, p2 = F.global)
KL.global = apply(F, 2, KL, p2 = F.global)
LSD.global = apply(F, 2, LSD, p2 = F.global)

#####
# Compute different test statistics
#####

var.KS.m = var(KS.mean)
var.KL.m = var(KL.mean)
var.LSD.m = var(LSD.mean)
var.KS.g = var(KS.global)
var.KL.g = var(KL.global)
var.LSD.g = var(LSD.global)

# Combine test statistics into same object
TS.orig = c(var.KS.m,
            var.KL.m,
            var.LSD.m,
            var.KS.g,
            var.KL.g,
            var.LSD.g)

TS = matrix(NA, nrow = S, ncol = length(TS.orig))

for (i in 1:S) {
  
  print(i)
  
  # Create stationary surrogate
  y.surrogate = surrogate(y)
  
  # Compute spectrogram
  F = spectrogram(y.surrogate, win.length = win.length,
                  overlap = overlap, win.type = win.type,
                  type = "ar", plot = FALSE)$spec
  F.mean = apply(F, 1, mean)
  F.global = as.vector(spec.ar(tukeywindow(length(y.surrogate), r = 0.1) * y.surrogate, 
                               n.freq = length(F.mean), 
                               plot = FALSE)$spec)
  
  #####
  # Compute distances
  #####
  
  # 1. Compute KS/KL distance for each segment vs. the mean
  KS.mean = apply(F, 2, KS, p2 = F.mean)
  KL.mean = apply(F, 2, KL, p2 = F.mean)
  LSD.mean = apply(F, 2, LSD, p2 = F.mean)
  
  # 2. Compute KS/KL distance for each segment vs. global PSD (same freqs)
  KS.global = apply(F, 2, KS, p2 = F.global)
  KL.global = apply(F, 2, KL, p2 = F.global)
  LSD.global = apply(F, 2, LSD, p2 = F.global)
  
  #####
  # Compute different test statistics
  #####
  
  var.KS.m = var(KS.mean)
  var.KL.m = var(KL.mean)
  var.LSD.m = var(LSD.mean)
  var.KS.g = var(KS.global)
  var.KL.g = var(KL.global)
  var.LSD.g = var(LSD.global)
  
  # Combine test statistics into same object
  TS[i, ] = c(var.KS.m,
              var.KL.m,
              var.LSD.m,
              var.KS.g,
              var.KL.g,
              var.LSD.g)
}

TS.95 = apply(TS, 2, quantile, probs = 0.95)  # Compute critival value
TS.orig >= TS.95  # Reject or not?

proc.time() - ptm
