          seed =  -1

       seqfile = ../bpp_input.txt
      Imapfile = ../mapping.txt
       outfile = results/out-IM_Ameivula.txt
      mcmcfile = results/mcmc-IM_Ameivula.txt
       
  speciesdelimitation = 0 * fixed species tree
* speciesdelimitation = 1 0 2    * species delimitation rjMCMC algorithm0 and finetune(e)
* speciesdelimitation = 1 1 2 1 * species delimitation rjMCMC algorithm1 finetune (a m)
*         speciestree = 1 0.1 0.1 0.2    * species tree SPR/SNL

*   speciesmodelprior = 1  * 0: uniform LH; 1:uniform rooted trees; 2: uniformSLH; 3: uniformSRooted

  species&tree = 2  red green
                    222  52 
                 (red,green);

         phase = 0
       usedata = 1  * 0: no data (prior); 1:seq like
         nloci = 3 * number of data sets in seqfile
         model = gtr

     cleandata = 1    * remove sites with ambiguity data (1:yes, 0:no)?

   thetaprior = 3 0.004 e  # Inv-gamma(a, b) for theta (integrated out by default; add E to also sample theta)
     tauprior = 3 0.02    # Inv-gamma(a, b) for root tau
*    thetaprior = gamma 2 2000  # gamma(a, b) for theta
*      tauprior = gamma 2 1000   # gamma(a, b) for root tau
 
      migprior = 10 1
      migration = 1 
	              red green

	locusrate = 1 0 0 5 iid     # (0: No variation, 1: estimate & a_Dirichlet, 2: from file) 
*    heredity = 2 ../heredity.txt
	
       finetune =  1: 5 0.001 0.001  0.001 0.3 0.33 1.0 # auto (0 or 1): finetune for GBtj, GBspr, theta, tau, mix, locusrate, seqerr

         print = 1 0 0 0   * MCMC samples, locusrate, heredityscalars, Genetrees
        burnin = 5000
      sampfreq = 4
       nsample = 10000
       threads = 3 1 1
       scaling = 1 