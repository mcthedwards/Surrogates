library(bspec)
library(fields)

surrogate = function(y) {
  
  n = length(y)
  i = complex(real = 0, imaginary = 1)
  
  yf = fft(y)          # Fourier transform
  magnitude = Mod(yf)  # Compute magnitude
  
  # Random phase, iid U[-pi, pi]
  phase.positive = runif(n / 2 + 1, -pi, pi)
  phase.negative = -rev(phase.positive[-c(1, n / 2 + 1)])
  phase.surrogate = c(phase.positive, phase.negative)
  phase.surrogate[c(1, n / 2 + 1)] = 0  # Set first and last positive coef to 0
  
  # Create surrogate in time domain
  y.surrogate = Re(fft(magnitude * exp(i * phase.surrogate), inverse = TRUE)) / n
  
  return(y.surrogate)
  
}

spectrogram = function(data, 
                       win.length, 
                       overlap = 0.75,
                       fs =  2 * pi,
                       win.type = "rectangular",
                       plot = FALSE,
                       type = "ar") {
  
  # fs: Sampling frequency
  
  n = length(data)
  
  increment = win.length * (1 - overlap)
  
  start = seq(1, n - win.length + 1, by = increment)
  end = seq(win.length, n, by = increment)
  
  n.increment = length(start) 
  n.freq = win.length / 2 + 1
  
  # Create window
  if (win.type == "rectangular") {
    ww = rep(1, win.length)
  }
  if (win.type == "tukey") {
    r = (1 - overlap) / 10  # Seems okay?
    ww = tukeywindow(win.length, r = r)
  }
  if (win.type == "hann") {
    ww = hannwindow(win.length)
  }
  
  if (type == "pdgrm") {
    # Compute magnitude of STFT for each segment
    spec = matrix(NA, nrow = n.freq, ncol = n.increment)
    for (j in 1:n.increment) {
      segment = start[j]:end[j]
      data.segment = data[segment]
      spec[, j] = (2 * Mod(fft(ww * data.segment)) ^ 2 / win.length)[1:n.freq]  # Multiply by 2 for 1-sided
    }
  }
  
  if (type == "ar") {
    spec = matrix(NA, nrow = n.freq, ncol = n.increment)
    for (j in 1:n.increment) {
      segment = start[j]:end[j]
      data.segment = data[segment]
      spec[, j] = as.vector(spec.ar(ww * data.segment, n.freq = n.freq, plot = FALSE)$spec)
    }
  }
  
  # Find associated times and frequencies
  deltat = 1 / fs
  time = seq(0, n * deltat, length = ncol(spec))
  freq = seq(0, fs / 2, length = nrow(spec))
  
  # Remove first and last frequency
  bFreq = c(1, length(freq))
  freq = freq[-bFreq]
  spec = spec[-bFreq, ]
  
  # spec = spec / max(spec)  # Normalise
  
  if (plot == TRUE) {
    image.plot(time, freq, t(spec),
               xlab = "Time", ylab = "Frequency", nlevel = 1000)
  }
  
  return(list(time = time,
              freq = freq,
              spec = spec))  
  
}

KS = function(p1, p2) {
  
  # Computes Kolgmorov-Smirnov statistic for two PSDs
  
  # p1 and p2 are PSDs
  
  # Normalise and take cumulative sum to find standardised CDF
  s1 = cumsum(p1 / sum(p1))  
  s2 = cumsum(p2 / sum(p2))
  
  # Compute distance
  D = max(abs(s1 - s2))
  
  return(D)
  
}

KL = function(p1, p2) {
  
  # Computes KL Divergence
  
  # p1 and p2 are PSDs
  
  s1 = p1 / sum(p1)
  s2 = p2 / sum(p2)
  
  # Compute symmetric KL divergence
  K1 = sum(s1 * log(s1 / s2))
  K2 = sum(s2 * log(s2 / s1))
  
  D = (K1 + K2) / 2
  
  return(D)
  
}

LSD = function(p1, p2) {
  
  D = sum(abs(log(p1) - log(p2)))
  
  return(D)
  
}

