# EFT constraints

## Requirements


`SL_OS_SS_combined_unblind_v3.0.json` contains the upper limit on tttt cross section.
Edit this file if you want to get different constrains.

## independent and marginal limits tables
result .tex tables will stored in `build` folder. To compile latex definitions from `resources/common/misc/definitions.tex` have to be included
```bash
mkdir build
python EFT.py -c conf_4op_unblind_paper_v1_cff.py --dir=build
```

## API for EFT cross sections
1. Select 13 TeV of 14 TeV predictions in **mg_calcuations.py**
```python
sig_SM=sig_SM_13TeV
sig_SM=sig_SM_14TeV
# and
MG_SM=MG_SM_13TeV
MG_SM=MG_SM_14TeV
```

2. To calculate tttt cross sections as function of EFT parameters

```python
from eft_coefficients import EftPredictions
from mg_calculations import wilson_coefficients, MG_SM, sig_SM

eft = EftPredictions(wilson_coefficients, MG_SM, sig_SM)

tttt_xs = eft.gen_eft_xs([C_OR,C_OL1,C_OL8,C_B1,CB8])
# or
tttt_xs = eft.vgen_eft_xs(C_OR,C_OL1,C_OL8,C_B1,CB8) # optimzed for numpy
```
alternatively one can use 13 and 14 TeV constants directly
```python
from eft_coefficients import EftPredictions
from mg_calculations import wilson_coefficients, MG_SM_13TeV, sig_SM_13TeV, MG_SM_14TeV, sig_SM_14TeV

eft13 = EftPredictions(wilson_coefficients, MG_SM_13TeV, sig_SM_13TeV)
eft14 = EftPredictions(wilson_coefficients, MG_SM_14TeV, sig_SM_14TeV)
tttt_xs_14_to_13_ratio = eft13.gen_eft_xs([C_OR,C_OL1,C_OL8,C_B1,CB8])/eft14.gen_eft_xs([C_OR,C_OL1,C_OL8,C_B1,CB8])

# test
eft13.gen_eft_xs([0.,0.,0.,0.,0.]) # should give 9.201
eft14.gen_eft_xs([0.,0.,0.,0.,0.]) # should give 11.31723
```

