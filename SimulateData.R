#library(FSCseq)
source('simulationFunctions.R')
#library(FSCseq)

# varying params
nS=3                     # number of studies
K=c(2,4)                 # number of clusters
n_k=c(25,50)             # number of subjects per cluster
LFCg=1                   # log fold change between clusters
pDEg=c(0.025,0.05)       # proportion of differentially expressed genes
beta0=8                  # baseline log2 mean count of genes
phi0=c(0.15,0.35)        # overdispersion of counts across samples for each gene
g=10000                  # number of genes

simAllData = function(nS,K,n_k,LFCg,pDEg,beta0,phi0,B,LFCb,nsims,g){
  # fixed params
  pDEb=0.5
  sigma_g=0.1; sigma_b=0   # Added Gaussian (nonparametric) noise
  n_pred=25
  #save_dir = "../Basic Simulations"
  save_dir = sprintf("Simulations/%f_%f/B%d_LFCb%d",sigma_g,sigma_b,B,LFCb)
  for(a in 1:length(K)){for(b in 1:length(n_k)){for(c in 1:length(LFCg)){for(d in 1:length(pDEg)){for(e in 1:length(beta0)){for(f in 1:length(phi0)){
    # Simulates data into "./Simulations/<sigma_g>_<sigma_b>/<K>_<n>_<LFCg>_<pDEg>_<beta0>_<phi0>_sim<>_data.RData"
    simulateData(nS=nS, K=K[a], B=B, g=g, n=K[a]*n_k[b], LFCg=LFCg[c], pDEg=pDEg[d], sigma_g=sigma_g,
                         LFCb=LFCb, pDEb=pDEb, sigma_b=sigma_b, beta0=beta0[e], phi0=phi0[f],
                         nsims=nsims, n_pred=n_pred, save_dir=NULL,save_file=T,save_pref=NULL)
    #FSCseq::simulateData(K=K[a], B=B, g=g, n=K[a]*n_k[b], LFCg=LFCg[c], pDEg=pDEg[d], sigma_g=sigma_g,
                     #LFCb=LFCb, pDEb=pDEb, sigma_b=sigma_b, beta0=beta0[e], phi0=phi0[f],
                     #nsims=nsims, n_pred=n_pred, save_dir=save_dir,save_file=T)
  }}}}}}
}

# no batch
B=1; LFCb=0
simAllData(nS=nS,K=K, n_k=n_k, LFCg=LFCg, pDEg=pDEg, beta0=beta0, phi0=phi0, B=B, LFCb=LFCb, nsims=100, g=g)