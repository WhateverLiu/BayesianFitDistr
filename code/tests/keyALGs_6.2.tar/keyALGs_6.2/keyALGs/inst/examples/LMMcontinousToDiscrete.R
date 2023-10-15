# =====================================================================================
# Discretize gamma distribution with shape = 2, rate = 1 / 2.
# Theoretical mean = 2 / (1 / 2) = 4.
# Theoretical variance = 2 / (1 / 2) ^ 2 = 8.
# =====================================================================================
support = seq(0, 15, len = 64)
disctre1 = keyALGs::LMMcontinousToDiscrete(dgamma, 2, 1 / 2, gd = support)
keyALGs::Mean(disctre1); keyALGs::Var(disctre1)
# 3.937484
# 7.183387


support = seq(0, 30, len = 64)
disctre2 = keyALGs::LMMcontinousToDiscrete(dgamma, 2, 1 / 2, gd = support)
keyALGs::Mean(disctre2); keyALGs::Var(disctre2)
# 3.999862
# 7.996146


support = seq(0, 60, len = 64)
disctre3 = keyALGs::LMMcontinousToDiscrete(dgamma, 2, 1 / 2, gd = support)
keyALGs::Mean(disctre3); keyALGs::Var(disctre3)
# 4
# 8


support = seq(0, qgamma(1 - 1e-10, 2, 1 / 2), len = 64)
disctre4 = keyALGs::LMMcontinousToDiscrete(dgamma, 2, 1 / 2, gd = support)
keyALGs::Mean(disctre4); keyALGs::Var(disctre4)
# 4
# 8
