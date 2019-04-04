# Introduction
Many models of physics beyond the standard model (BSM) predict enhanced or modified couplings of top quarks to other particles. 
There is a considerable interest in the measurement of the tttt cross section because of its sensitivity to BSM physics, including supersymmetry, two Higgs doublet models, top quark compositeness, and models with extra spatial dimensions. The availability of the EFT predictions makes a generic analysis of tttt production possible

This folder contains main analysis ingredients necessary for recasting upper limit on four top production cross section into constraints on four top dim-6 EFT operators.

# Theoretical predictions

## FeynRules model and dim-6 EFT operators
tttt cross section as a function of Wilson coefficients is calculated using MC@NLO generator with FeynRules model http://feynrules.irmp.ucl.ac.be/wiki/4topEFT

The list of dim-6 EFT operators can be found in Eq.(2) of the TOP-17-019 paper. There are four operators

* `O^1_{tt} = (\bar{t}_R \gamma^\mu t_R) (\bar{t}_R \gamma_\mu t_R)`
* `O^1_{QQ} = (\bar{Q}_L \gamma^\mu Q_L) (\bar{Q}_L \gamma_\mu Q_L)`
* `O^1_{Qt} = (\bar{Q}_L \gamma^\mu Q_L) (\bar{t}_R \gamma_\mu t_R)`
* `O^8_{Qt} = (\bar{Q}_L \gamma^\mu T^A Q_L) (\bar{t}_R T^A \gamma_\mu t_R)`

The notation used in the TOP-17-019 paper is different from the original FeynRules in order to be consistent with SMEFT LHC WG note [arXiv:1802.07237](https://arxiv.org/abs/1802.07237). In order to convert one notation into another, following dictionary can be used

* `O_R -> O^1_{tt}`
* `O^{(1)}_L -> O^1_{QQ}`
* `O^{(8)}_L -> O^8_{QQ}`
* `O^{(1)}_B -> O^1_{Qt}`
* `O^{(8)}_B -> O^8_{Qt}`

The obtained predictions are parametrized by a polynomial of second degree in 4 variables, see Eq.(3). Linear and quadratic parametrization coefficients can be found in Tab.4.

## MadGraph5_aMCatNLO card that was used for the parametrization of EFT predictions
see `mg_cards/run_02_tag_1_banner.txt`

## Assumptions

* Acceptance on the tttt production is estimated using MadGraph5_aMCatNLO NLO QCD tttt samples. Following datacards were used
https://github.com/cms-sw/genproductions/tree/master/bin/MadGraph5_aMCatNLO/cards/production/2017/13TeV/TTTT_5f_NLO

* SM tttt cross section, \sigma_{tttt}^{SM} is set to 9.2 fb.

* Since NLO EFT predictions were not available, scales variations as well as PDF uncertainties are assumed to be the same as for SM NLO predictions. Internal study, documented in the AN, showed that SM and EFT uncertainties agree very well at LO accuracy.

* No explicit assumption on EFT cut-off. It was shown, for example, in [arXiv:1901.05965](https://arxiv.org/abs/1901.05965), see the discussion of m_{tttt} cut in the last paragraph of Sec.5 on p.60.
