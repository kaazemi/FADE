# FADE
An implementation of the FAst DEconvolution (FADE) Algorithm
# Input: 
  Structure called sys with the following fields:
  
  y: observed calcium traces where the first dimension is the number of ROI's and second dimension is time and
  
  p_norm and q_norm: Penalizes the ell_pq norm of spikes, default: p = q = 1
  
  lambda: regularization parameter, default: lambda = 0 (ML Estimation)
  
  Order: Autoregressive model order, defualt: Order = 2
  
  theta: AR parameters, default: estimated using Yule-Walker equations
  
  noise: estimate of the noise standard variation, default:
  
  estimate using power spectral density
  
  min_iters: minimum number of iterations, default: min_iters = 100
  
  num_iters: maximum number of iterations, default: num_iters = 200
  
 #Output: 
 
  A structure containing the following added fields
 
  spikes: the deconvolved spikes, post processing methods should often be used
  
  smoothed_traces: smoothed calcium traces
   
  ds: relative changes in the spikes  
