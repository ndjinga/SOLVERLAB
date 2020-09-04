import numpy as np

# # TSat(P) takes a constant value

# # compute hSat_k; k= liq, vap
# hSat_k(P) = q_k + cp_k* Tsat(P)

# # compute rhoSat_k; k= liq, vap
# rhoSat_k(P) = (P + pi_k) / ((cp_k - cv_k)* Tsat(P))

class ThermoStiffenedGas_CstTsat:
    def __init__(self,
                 cvLiq, gammaLiq, piLiq, qLiq, qPrimeLiq,
                 cvVap, gammaVap, piVap, qVap, qPrimeVap
                 ):
        self.gammaLiq = gammaLiq
        self.gammaVap = gammaVap
        self.cvLiq = cvLiq
        self.cvVap = cvVap
        self.cpLiq = self.cvLiq * self.gammaLiq
        self.cpVap = self.cvVap * self.gammaVap
        self.piLiq = piLiq
        self.piVap = piVap
        self.qLiq = qLiq
        self.qVap = qVap
        self.qPrimeLiq = qPrimeLiq
        self.qPrimeVap = qPrimeVap
        self.A = self.cpVap - self.cvVap
        self.B = self.cpLiq - self.cvLiq
        self.C = self.cpVap - self.cpLiq
        self.D = self.qVap - self.qLiq

        return

# this function P-> Tsat(P) is a "crude" approximation
# that replaces the resolution the fixed point
# A*ln( p + piVap) - B*ln(p + piLiq) - C*(ln(T) - 1) + D/T + qPrimeLiq = 0
# A = cpVap - cvVap
# B = cpLiq - cvLiq
# C = cpVap - cpLiq
# D = qVap - qLiq
    def P_to_Tsat(self, P):
        return 656.0
    
    def T_to_hLiq(self, T):
        return self.qLiq + self.cpLiq * T

    def T_to_hVap(self, T):
        return self.qVap + self.cpVap * T
    
    def P_to_hSatLiq(self, P):
        return self.T_to_hLiq(self.P_to_Tsat(P))
    
    def P_to_hSatVap(self, P):
        return self.T_to_hVap(self.P_to_Tsat(P))
    
    def P_to_rhoSatLiq(self, P):
        return (P + self.piLiq) / ((self.cpLiq - self.cvLiq) * self.P_to_Tsat(P))
    
    def P_to_rhoSatVap(self, P):
        return (P + self.piVap) / ((self.cpVap - self.cvVap) * self.P_to_Tsat(P))
    
    def h_P_to_rhoLiq(self, h, P):
        return self.gammaLiq * (P + self.piLiq) / ((self.gammaLiq - 1) * (h - self.qLiq))

    def h_P_to_rhoVap(self, h, P):
        return self.gammaVap * (P + self.piVap) / ((self.gammaVap - 1) * (h - self.qVap))
