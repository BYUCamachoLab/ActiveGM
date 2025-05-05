from load_sparams import get_GM
from utils import *
import numpy as np
import matplotlib.pyplot as plt

randomseed = 5
np.random.seed(randomseed)

print("Finished imports")
print("Beginning Circuit Construction...")

# Get the general structure of the circuit from load_sparams without including any specific fabrication parameters
circuit = get_GM()

# Coupling imbalance: The transmission of the couplers on the MZIs (.5 for ideal couplers)
coupling_imbalance = 0.4

# Extra coupling error: Standard deviation of errors on the couplers
extra_coupling_error = 0.01

coupling = np.random.normal(loc=coupling_imbalance, scale=extra_coupling_error, size=(8, 2, 4))

# These are the default phase shifts in the MZIs to create a green machine transformation
phi = np.array([[np.pi/2, np.pi/2, np.pi/2, np.pi/2], [0, 0, 0, 0], [np.pi, 0, 0, np.pi], [0, 0, 0, 0], [np.pi/2, np.pi/2, np.pi/2, np.pi/2], [0, np.pi, 0, np.pi], [np.pi/2, np.pi/2, np.pi/2, np.pi/2], [np.pi, np.pi, np.pi, np.pi]])

# Any phase shift errors (radians) for all 8 modes after each of the 8 layers of the Green Machine
internal_theta = np.random.rand(8, 8) * 0.5

# Any loss errors (dB) for all 8 modes at all 8 layers of the Green Machine
internal_loss = np.ones((8, 8)) * 0.2

# The loss in one modulator inside the MZI (maybe poorly named)
closs = 0.3

# The modulators inside of the chip. Phi_corrections are inside the MZI's, theta_corrections are on the outside arms on half of the modes
phi_corrections = np.zeros((8, 4))
# There are only phase shifters on every other arm (0, 2, 4, 6) not all 8
theta_corrections = np.zeros((8, 8))

# These phase shifts work really well for coupling = np.random.normal(loc=0.4, scale=0.01, size=(8, 2, 4)), internal_theta = np.random.rand(8, 8) * 0.5, internal_loss = np.ones((8, 8)) * 0.2, closs = 0.3 with numpy random seed of 5
#phi_corrections = np.array([[-0.04487989,-0.07479983,-0.02991993,-0.04487988],[-0.29919931,0.11967972,0.16455961,-0.94247782],[-0.28423934,-0.2094395,0.14959965,0.26927938],[0.11967973,0.22439948,-0.32911924,-0.94247782],[-0.14959966,-0.17951959,-0.32911925,-0.02991992],[-0.02991993,-0.34407919,-0.02991993,-0.94247782],[-0.20943951,0.2094395,-0.433839,0.],[-0.35903917,-0.67319844,0.19447955,-0.94247782]])
#theta_corrections = np.array([[-7.47998300e-02,0.00000000e+00,3.59039170e-01,0.00000000e+00,-7.47998300e-02,0.00000000e+00,4.93678855e-01,0.00000000e+00],[1.34639695e-01,0.00000000e+00,1.19679730e-01,0.00000000e+00,3.73999125e-01,0.00000000e+00,5.98398601e-02,0.00000000e+00],[2.09439520e-01,0.00000000e+00,2.39359450e-01,0.00000000e+00,-2.99199351e-02,0.00000000e+00,-7.18078340e-01,0.00000000e+00],[2.84239355e-01,0.00000000e+00,-7.47998265e-01,0.00000000e+00,6.88158395e-01,0.00000000e+00,-4.03919065e-01,0.00000000e+00],[2.99199315e-01,0.00000000e+00,-5.98398650e-02,0.00000000e+00,4.93678855e-01,0.00000000e+00,-5.68478670e-01,0.00000000e+00],[3.88959100e-01,0.00000000e+00,5.98398650e-02,0.00000000e+00,-5.98398650e-02,0.00000000e+00,-4.18879030e-01,0.00000000e+00],[2.39359455e-01,0.00000000e+00,7.92878155e-01,0.00000000e+00,5.05128275e-09,0.00000000e+00,-3.29119225e-01,0.00000000e+00],[-3.59039175e-01,0.00000000e+00,-4.03919065e-01,0.00000000e+00,1.34639700e-01,0.00000000e+00,9.42477815e-01,0.00000000e+00]])

# Get ideal sparams (ideal couplers, no internal phase errors, but still including loss) and the first attempt sparams with all of our errors
ideal = get_sparams(circuit, np.ones((8, 2, 4)) * 0.5, phi, closs, internal_loss, np.zeros((8, 8)))
og = get_sparams(circuit, coupling, phi + phi_corrections, closs, internal_loss, internal_theta + theta_corrections)

print("Current sparams:")
print_sparams(og, bukon=0, n=8)
print(quality_func(og))
print("Ideal sparams:")
print_sparams(ideal, bukon=0, n=8)
print(quality_func(ideal))

# Simulation / Optimization parameters (wl doesn't change much of anything, most internal sparams are independent of this value but simphony still requires it)
iterations = 6
step_size = np.pi/210
wl=1.550

# An iterative method for finding phase shifts that lead to optimal coupling / phase relationships. Currently, it does better with the coupling than phase errors, this could be changed by updating the quality function in utils.py
for i in range(iterations):
    sparams = get_sparams(circuit, coupling, phi + phi_corrections, closs, internal_loss, internal_theta + theta_corrections)
    defaultquality = quality_func(sparams)
    print(defaultquality)
    print("")

    delta_phi_corrections = np.zeros((8, 4))
    for p1 in range(8):
        for p2 in range(4):
            print(f'Starting phi modulator ({p1}, {p2})')
            max = -1
            bestdiff = -1
            for d in range(7):
                diff = step_size * (d-3)
                delta_phi_corrections[p1][p2] += diff
                new_sparams = get_sparams(circuit, coupling, phi + phi_corrections + delta_phi_corrections, closs, internal_loss, internal_theta + theta_corrections)
                newquality = quality_func(new_sparams)
                delta_phi_corrections[p1][p2] -= diff
                if newquality > max:
                    max = newquality
                    bestdiff = diff

            delta_phi_corrections[p1][p2] = bestdiff
            defaultquality = max
            print(f'Best shift: {bestdiff}')
            print(f'Best q factor: {defaultquality}')

    print("")
    print("Done with phi corrections")
    print(delta_phi_corrections)
    print("")

    phi_corrections += delta_phi_corrections

    delta_theta_corrections = np.zeros((8, 8))
    for p1 in range(8):
        for p2 in range(4):
            print(f'Starting theta modulator ({p1}, {p2})')
            max = defaultquality
            bestdiff = 0
            for d in range(7):
                diff = step_size * (d-3)
                delta_theta_corrections[p1][2*p2] += diff
                new_sparams = get_sparams(circuit, coupling, phi + phi_corrections, closs, internal_loss, internal_theta + theta_corrections + delta_theta_corrections)
                newquality = quality_func(new_sparams)
                delta_theta_corrections[p1][2*p2] -= diff
                if newquality > max:
                    max = newquality
                    bestdiff = diff

            delta_theta_corrections[p1][2*p2] = bestdiff
            defaultquality = max
            print(f'Best shift: {bestdiff}')
            print(f'Best q factor: {defaultquality}')

    theta_corrections += delta_theta_corrections

    print("")
    print("Done with theta corrections")
    print(delta_theta_corrections)
    print("")

print("")
print("")
print("We started with: ")
print_sparams(og)
print(quality_func(og))
print("")

final = get_sparams(circuit, coupling, phi + phi_corrections, closs, internal_loss, internal_theta + theta_corrections)
print("And were able to correct to: ")
print_sparams(final)
print(quality_func(final))
print("")

print("By using phase shifts of: ")
print(phi_corrections)
print(theta_corrections)
print("")

"""
# This is a gradient descent algorithm for optimizing the phase shifts. It doesn't work as well as the iterative method above.
for i in range(iterations):
    new_sparams = get_sparams_bukon(circuit, coupling, phi + phi_corrections, bukon_coupling, bukon_phi, closs, internal_loss, internal_theta + theta_corrections)
    defaultquality = quality_func(new_sparams)
    print(defaultquality)
    print("")

    print("Iteration: " + str(i))

    delta_phi_corrections = np.random.normal(0, step_size, (samples, 8, 4))

    delta_theta_corrections = np.zeros((samples, 8, 14))
    for n in range(4):
        delta_theta_corrections[:, :, 2*n] = np.random.normal(0, step_size, (samples, 8))

    delta_prime = np.zeros(samples)

    phi_correction = np.zeros((8, 4))
    theta_correction = np.zeros((8, 14))

    valid = False

    for j in range(samples):
        if j % 10 == 0:
            print(f'Testing descent #{j}')
        sparams = get_sparams_bukon(circuit, coupling, phi + phi_corrections + delta_phi_corrections[j], bukon_coupling, bukon_phi, closs, internal_loss, internal_theta + theta_corrections + delta_theta_corrections[j])
        
        newquality = quality_func(sparams)
        delta_prime[j] = newquality - defaultquality

        if delta_prime[j] > 0:
            phi_corrections += 5 * delta_prime[j] * delta_phi_corrections[j]
            theta_corrections += 5 * delta_prime[j] * delta_theta_corrections[j]

    print("Total Correction:")
    print(phi_corrections)
    print(theta_corrections)
"""