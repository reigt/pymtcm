import mtcm

from math import pi as pi

# Input
sigma_sr = 250 # MPa
Lc = 5000 # mm
hc = 300 # mm
wc = 1000 # mm
cover = 50 # mm
phi_s = 25 # mm
cc_phi_s = 150 # mm
Es = 200000 # MPa
Ecm = 30000 # MPa
fctm = 3.2 # MPa

# Crack width parameters
eps_sr = sigma_sr/Es
n_phi_s = wc/cc_phi_s
As_tot = n_phi_s*pi*phi_s**2/4
hc_ef = min([2.5*(cover+phi_s/2), hc/2])
rho_s = As_tot/(wc*hc_ef)

# MTCM stress
test = mtcm.mtcm(Lc,As_tot,n_phi_s,phi_s,rho_s,Es,Ecm,fctm)
test.stress(eps_sr)

# MTCM strain
test2 = mtcm.mtcm(Lc,As_tot,n_phi_s,phi_s,rho_s,Es,Ecm,fctm)
test2.strain(test.eps_sm)

# Plotting
# plot_test = mtcm.PlotMTCM(test)
# plot_test.chord_distribution()
plot_test2 = mtcm.PlotMTCM(test2)
plot_test2.chord_distribution()
plot_test2.stress_from_strain()

# Print results
results_test = {
    'concept': test.concept,
    'xcr0': test.xcr0,
    'eps_sm': test.eps_sm,
    'Lt': test.Lt/2,
    'wcr': test.wcr,
}
results_test2 = {
    'concept': test2.concept,
    'xcr0': test2.xcr0,
    'eps_sm': test2.eps_sm,
    'Lt': test2.Lt/2,
    'wcr': test2.wcr,
}
print(results_test)
print(results_test2)

