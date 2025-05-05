import numpy as np
import matplotlib.pyplot as plt

# Two functions to generate a GM (with / without Bukon ports)

# This function allows you to set parameters for the bokun MZIs if the bokun model is being used.
def get_sparams_bukon(circuit, coupling, phi, bukon_coupling, bukon_phi, closs, internal_loss, internal_theta):
    current_bukon = 0
    mask = np.array([[0,0,0,0],[0,1,1,0],[1,1,0,0],[1,1,1,1],[1,1,0,0],[0,1,1,0],[0,0,0,0],[0,0,0,0]])

    layer_descriptor = []
    for i in range(8):
        total_couplers = 6 + (i%2)
        central_couplers = 4 - (i%2)
        bukon_couplers = total_couplers - central_couplers

        coupling1 = np.zeros(total_couplers)
        coupling2 = np.zeros(total_couplers)
        new_phi = np.zeros(total_couplers)
        new_closs = np.ones(total_couplers) * closs

        coupling1[:central_couplers] = coupling[i][0][:central_couplers]
        coupling2[:central_couplers] = coupling[i][1][:central_couplers]
        new_phi[:central_couplers] = phi[i][:central_couplers]

        for b in range(bukon_couplers):
            if mask[i][b]:
                coupling1[central_couplers + b] = bukon_coupling[current_bukon][0]
                coupling2[central_couplers + b] = bukon_coupling[current_bukon][1]

                new_phi[central_couplers + b] = bukon_phi[current_bukon]

                current_bukon += 1

            else: # Not coupled, no loss
                coupling1[central_couplers + b] = .5
                coupling2[central_couplers + b] = .5
                new_phi[central_couplers + b] = np.pi 
                new_closs[central_couplers + b] = 0

        layer_descriptor.append({"coupling1": coupling1, "coupling2": coupling2, "loss": internal_loss[i], "theta": internal_theta[i], "phi": new_phi, "closs": new_closs})

    sparams = circuit(wl=1.549,
                  layer0=layer_descriptor[0], layer1=layer_descriptor[1], layer2=layer_descriptor[2], layer3=layer_descriptor[3],
                  layer4=layer_descriptor[4], layer5=layer_descriptor[5], layer6=layer_descriptor[6], layer7=layer_descriptor[7])

    return sparams

# This function gets a GM without the extra bokun ports (under the hood it uses the same model but doesn't couple the bokun ports into the GM)
def get_sparams(circuit, coupling, phi, closs, internal_loss, internal_theta):
    length = 0
    loss = 0

    bukon_coupling = np.ones((12, 2)) * 0.5
    bukon_phi = np.ones(12) * np.pi
    new_internal_loss = np.zeros((8, 14))
    new_internal_loss[:, :8] = internal_loss
    new_internal_theta = np.zeros((8, 14))
    new_internal_theta[:, :8] = internal_theta

    return get_sparams_bukon(circuit, coupling, phi, bukon_coupling, bukon_phi, closs, new_internal_loss, new_internal_theta)


# Different functions to print out sparams (set bukon to 0 if printing out a normal GM)
def sparams_phase(sparams, bukon = 0, n = 8):
    result = np.zeros((n + bukon, n + bukon))
    for i in range(n + bukon):
        for o in range(n + bukon):
            inseed = 'in'
            outseed = 'out'
            if i >= n:
                inseed += 'p'
            if o >= n:
                outseed += 'p'

            result[i][o] = np.angle(sparams[(inseed + str(i%n), outseed + str(o%n))])

    return result

def sparams_power(sparams, bukon = 0, n = 8):
    result = np.zeros((n + bukon, n + bukon))
    for i in range(n + bukon):
        for o in range(n + bukon):
            inseed = 'in'
            outseed = 'out'
            if i >= n:
                inseed += 'p'
            if o >= n:
                outseed += 'p'

            result[i][o] = np.abs(sparams[(inseed + str(i%n), outseed + str(o%n))])**2
            if result[i][o] < 0.00001:
                result[i][o] = 0

    return result

def sparams_total(sparams, bukon = 0, n = 8):
    result = np.zeros((n + bukon, n + bukon), dtype=complex)
    for i in range(n + bukon):
        for o in range(n + bukon):
            inseed = 'in'
            outseed = 'out'
            if i >= n:
                inseed += 'p'
            if o >= n:
                outseed += 'p'

            result[i][o] = sparams[(inseed + str(i%n), outseed + str(o%n))]
            if result[i][o] < 0.00001:
                result[i][o] = 0

    return result

def print_sparams(sparams, bukon = 0, n = 8):
    powr = sparams_power(sparams, bukon, n).T
    phase = sparams_phase(sparams, bukon, n).T
    print(powr)
    print(phase)
    plt.imshow(powr, vmin=.015, vmax=.16)
    plt.colorbar()
    plt.xlabel("Input Port")
    plt.ylabel("Output Port")
    plt.show()




# Different functions to evaluate the quality of a GM transformation
def phase_quality_func(phases):
    for i in range(8)[::-1]:
        for o in range(8):
            phases[o][i] -= phases[o][0]

    for o in range(8)[::-1]:
        for i in range(8):
            phases[o][i] -= phases[0][i]

    phases = (phases + 2 * np.pi) % (2*np.pi)

    distances = np.zeros((7, 7))

    for o in range(7):
        for i in range(7):
            distances[o][i] = np.min(np.abs(np.array([0, np.pi, 2*np.pi]) - phases[o+1][i+1]))

    return distances

def quality_func(sparams):
    s1_dB = sparams_power(sparams)
    s1_phases = sparams_phase(sparams)

    phase_error = phase_quality_func(s1_phases)

    dB_std = np.std(s1_dB)
    phase_std = np.mean(phase_error)

    return np.exp(-20*dB_std) - (.5*phase_std)