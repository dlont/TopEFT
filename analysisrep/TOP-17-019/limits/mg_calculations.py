# SM tttt cross section value from MG
sig_SM_13TeV = 0.009201 * 1000. # fb

kFactor13TeV = 0.009201 / 0.009031

# (INFO) RT gg 14/13 TeV parton lumi ratio (https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageAt1314TeV2014#Parton_luminosity_ratio)
#            1.227607773
sig_SM_14TeV = 0.01192 * 1000. # fb
sig_SM_14TeV = sig_SM_14TeV * kFactor13TeV

sig_SM_27TeV = 0.1036 * 1000.
sig_SM_27TeV = sig_SM_27TeV * kFactor13TeV

# Predefined values of Wilson coefs. for which the EFT tttt cross section was calculated.
# One has to make sure that resulting matrix for sigma_i and sigma_ij is not degenerate
wilson_coefficients = [
    [1, 0, 0, 0, 0],
    [0, 1, 0, 0, 0],
    [0, 0, 1, 0, 0],
    [0, 0, 0, 1, 0],
    [0, 0, 0, 0, 1],
    [1, 1, 1, 1, 1],
    [-1, -1, 1, 1, 1],
    [-1, -1, 1, 0, 1],
    [0, 1, 0, 0, -1],
    [0, 1, 1, 1, 0],
    [1, 0, -1, 1, 0],
    [-1, 0, 0, 1, -1],
    [-1, 0, 0, -1, 1],
    [0, 1, -1, 1, -1],
    [0, 1, 0, -1, 0],
    [0, 0, -1, -1, 1],
    [1, -1, 0, -1, 0],
    [1, 1, 0, -1, 1],
    [0, 1, 0, -1, 1],
    [1, -1, -1, 0, 1]
]

# EFT cross section values for different values of Ci in the vectors c1, c2, c3, ... c19, c20
MG_SM_13TeV = [0.01557, 0.01564, 0.0102, 0.01116, 0.01022, 0.02873, 0.0203, 0.01704, 0.01527, 0.02066, 0.0168, 0.01741, 0.01567, 0.01342, 0.01839, 0.01208, 0.02113, 0.0283, 0.01983, 0.02386]  # pb
MG_SM_14TeV = [0.02087, 0.02087, 0.0133, 0.01464, 0.01323, 0.0389, 0.02744, 0.02276, 0.0206, 0.02768, 0.0225, 0.02341, 0.02088, 0.01782, 0.0247, 0.01604, 0.02874, 0.03798, 0.02641, 0.03234]  # pb
MG_SM_27TeV = [0.2251, 0.2247, 0.1189, 0.144, 0.1172, 0.48, 0.332, 0.2688, 0.2237, 0.3327, 0.2628, 0.2745, 0.2457, 0.1939, 0.2768, 0.1607, 0.3564, 0.4429, 0.2944, 0.4082] # pb

##############################################
sig_SM=sig_SM_13TeV
MG_SM=MG_SM_13TeV
sig_SM=sig_SM_13TeV
