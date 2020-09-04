#Testing equation of state with constant saturation temperature
import testEnv
import ThermoStiffenedGas_CstTsat as eos

thermo = eos.ThermoStiffenedGas_CstTsat(
                                cvLiq=1816.2,
                                gammaLiq=2.35,
                                piLiq=1.0e9,
                                qLiq=-1167.056e3,
                                qPrimeLiq=0.0,
                                cvVap=1040.14,
                                gammaVap=1.43,
                                piVap=0.0,
                                qVap=2030.255e3,
                                qPrimeVap=-23310.0,
                                )
P = 155e5
print "Tsat(P=%10.20e) = %10.20e" % (P, thermo.P_to_Tsat(P))
print "hSatLiq(P=%10.20e) = %10.20e" % (P, thermo.P_to_hSatLiq(P))
print "hSatVap(P=%10.20e) = %10.20e" % (P, thermo.P_to_hSatVap(P))
print "rhoSatLiq(P=%10.20e) = %10.20e" % (P, thermo.P_to_rhoSatLiq(P))
print "rhoSatVap(P=%10.20e) = %10.20e" % (P, thermo.P_to_rhoSatVap(P))

P = 145e5
print "Tsat(P=%10.20e) = %10.20e" % (P, thermo.P_to_Tsat(P))
print "hSatLiq(P=%10.20e) = %10.20e" % (P, thermo.P_to_hSatLiq(P))
print "hSatVap(P=%10.20e) = %10.20e" % (P, thermo.P_to_hSatVap(P))
print "rhoSatLiq(P=%10.20e) = %10.20e" % (P, thermo.P_to_rhoSatLiq(P))
print "rhoSatVap(P=%10.20e) = %10.20e" % (P, thermo.P_to_rhoSatVap(P))


print "done"
