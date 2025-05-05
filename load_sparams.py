import numpy as np
from pathlib import Path

import warnings
import sax

from simphony.libraries import siepic, ideal

from typing import List

def pic_sublayer_type1_bukon(coupling1: List[float] = None, coupling2: List[float] = None, loss: List[float] = None, theta: List[float] = None, phi: List[float] = None, closs : List[float] = None) -> sax.SDict:
    if coupling1 is None:
         coupling1 = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5]

    if coupling2 is None:
         coupling2 = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5]

    if loss is None:
        loss = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

    if theta is None:
        theta = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

    if phi is None:
        phi = [0, 0, 0, 0, 0, 0]

    if closs is None:
        closs = [0, 0, 0, 0, 0, 0]

    coupling1 = np.array(coupling1)
    coupling2 = np.array(coupling2)
    loss = np.array(loss)
    theta = np.array(theta)
    phi = np.array(phi)
    closs = 10**(-np.array(closs)/10)
    coupling_theta = - np.pi/2

    kappa0 = (coupling1**0.5) * (1-coupling2)**0.5 * np.exp(1j*coupling_theta) + (1-coupling1)**0.5 * (coupling2**0.5) * (closs**0.5) * np.exp(1j*(phi + coupling_theta))
    tau0 = (1- coupling1)**0.5 * (1 - coupling2)**0.5 * (closs**0.5) * np.exp(1j*phi) + (coupling1**0.5) * (coupling2**0.5) * np.exp(2j*coupling_theta)
    kappa1 = (coupling1**0.5) * (1-coupling2)**0.5 * (closs**0.5) * np.exp(1j*(coupling_theta+phi)) + (1-coupling1)**0.5 * (coupling2**0.5) * np.exp(1j*coupling_theta)
    tau1 = (1- coupling1)**0.5 * (1 - coupling2)**0.5 + (coupling1**0.5) * (coupling2**0.5) * (closs**0.5) * np.exp(1j*(2*coupling_theta+phi))
    gamma = 10**(-loss/20) * np.exp(1j*theta)

    sdict = sax.reciprocal(
        {
            ("in0", "out0"): gamma[0] * tau0[0],
            ("in0", "out1"): gamma[0] * kappa0[0],
            ("in1", "out0"): gamma[1] * kappa1[0],
            ("in1", "out1"): gamma[1] * tau1[0],
            ("in2", "out2"): gamma[2] * tau0[1],
            ("in2", "out3"): gamma[2] * kappa0[1],
            ("in3", "out2"): gamma[3] * kappa1[1],
            ("in3", "out3"): gamma[3] * tau1[1],
            ("in4", "out4"): gamma[4] * tau0[2],
            ("in4", "out5"): gamma[4] * kappa0[2],
            ("in5", "out4"): gamma[5] * kappa1[2],
            ("in5", "out5"): gamma[5] * tau1[2],
            ("in6", "out6"): gamma[6] * tau0[3],
            ("in6", "out7"): gamma[6] * kappa0[3],
            ("in7", "out6"): gamma[7] * kappa1[3],
            ("in7", "out7"): gamma[7] * tau1[3],

            ("inp0", "outp0"): gamma[8],

            ("inp1", "outp1"): gamma[9] * tau0[4],
            ("inp1", "outp2"): gamma[9] * kappa0[4],
            ("inp2", "outp1"): gamma[10] * kappa1[4],
            ("inp2", "outp2"): gamma[10] * tau1[4],

            ("inp3", "outp3"): gamma[11] * tau0[5],
            ("inp3", "outp4"): gamma[11] * kappa0[5],
            ("inp4", "outp3"): gamma[12] * kappa1[5],
            ("inp4", "outp4"): gamma[12] * tau1[5],

            ("inp5", "outp5"): gamma[13],
        }
    )

    return sdict

def pic_sublayer_type2_bukon(coupling1: List[float] = None, coupling2: List[float] = None, loss: List[float] = None, theta: List[float] = None, phi: List[float] = None, closs: List[float] = None) -> sax.SDict:
    if coupling1 is None:
         coupling1 = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]

    if coupling2 is None:
         coupling2 = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]

    if loss is None:
        loss = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

    if theta is None:
        theta = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

    if phi is None:
        phi = [0, 0, 0, 0, 0, 0, 0]

    if closs is None:
        closs = [0, 0, 0, 0, 0, 0, 0]

    coupling1 = np.array(coupling1)
    coupling2 = np.array(coupling2)
    loss = np.array(loss)
    theta = np.array(theta)
    phi = np.array(phi)
    closs = 10**(-np.array(closs)/10)
    coupling_theta = - np.pi/2

    kappa0 = (coupling1**0.5) * (1-coupling2)**0.5 * np.exp(1j*coupling_theta) + (1-coupling1)**0.5 * (coupling2**0.5) * (closs**0.5) * np.exp(1j*(phi + coupling_theta))
    tau0 = (1- coupling1)**0.5 * (1 - coupling2)**0.5 * (closs**0.5) * np.exp(1j*phi) + (coupling1**0.5) * (coupling2**0.5) * np.exp(2j*coupling_theta)
    kappa1 = (coupling1**0.5) * (1-coupling2)**0.5 * (closs**0.5) * np.exp(1j*(coupling_theta+phi)) + (1-coupling1)**0.5 * (coupling2**0.5) * np.exp(1j*coupling_theta)
    tau1 = (1- coupling1)**0.5 * (1 - coupling2)**0.5 + (coupling1**0.5) * (coupling2**0.5) * (closs**0.5) * np.exp(1j*(2*coupling_theta+phi))
    gamma = 10**(-loss/20) * np.exp(1j*theta)

    sdict = sax.reciprocal(
        {
            ("in1", "out1"): gamma[1] * tau0[0],
            ("in1", "out2"): gamma[1] * kappa0[0],
            ("in2", "out1"): gamma[2] * kappa1[0],
            ("in2", "out2"): gamma[2] * tau1[0],
            ("in3", "out3"): gamma[3] * tau0[1],
            ("in3", "out4"): gamma[3] * kappa0[1],
            ("in4", "out3"): gamma[4] * kappa1[1],
            ("in4", "out4"): gamma[4] * tau1[1],
            ("in5", "out5"): gamma[5] * tau0[2],
            ("in5", "out6"): gamma[5] * kappa0[2],
            ("in6", "out5"): gamma[6] * kappa1[2],
            ("in6", "out6"): gamma[6] * tau1[2],

            ("inp0", "outp0"): gamma[8] * tau0[3],
            ("inp0", "outp1"): gamma[8] * kappa0[3],
            ("inp1", "outp0"): gamma[9] * kappa1[3],
            ("inp1", "outp1"): gamma[9] * tau1[3],

            ("inp2", "outp2"): gamma[10] * tau0[4],
            ("inp2", "out0"): gamma[10] * kappa0[4],
            ("in0", "outp2"): gamma[0] * kappa1[4],
            ("in0", "out0"): gamma[0] * tau1[4],

            ("in7", "out7"): gamma[7] * tau0[5],
            ("in7", "outp3"): gamma[7] * kappa0[5],
            ("inp3", "out7"): gamma[11] * kappa1[5],
            ("inp3", "outp3"): gamma[11] * tau1[5],

            ("inp4", "outp4"): gamma[12] * tau0[6],
            ("inp4", "outp5"): gamma[12] * kappa0[6],
            ("inp5", "outp4"): gamma[13] * kappa1[6],
            ("inp5", "outp5"): gamma[13] * tau1[6],
        }
    )

    return sdict

def get_GM():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        instance_dict = {}
        connection_dict = {}
        port_dict = {}

        for i in range(8):
            instance_dict[f'layer{i}'] = f'layer_type{i%2+1}'

            #tst = f'layer_type{i%2+1}'
            #print(f'Layer {i} assigned: {tst}')

            if i != 7:
                for j in range(8):
                    connection_dict[f'layer{i}, out{j}'] = f'layer{i+1}, in{j}'

                for j in range(6):
                    connection_dict[f'layer{i}, outp{j}'] = f'layer{i+1}, inp{j}'

            port_dict[f'in{i}'] = f'layer0, in{i}'
            port_dict[f'out{i}'] = f'layer7, out{i}'

        for i in range(6):
            port_dict[f'inp{i}'] = f'layer0, inp{i}'
            port_dict[f'outp{i}'] = f'layer7, outp{i}'

        GM_bukon = {"instances": instance_dict, "connections": connection_dict, "ports": port_dict}

        eve, info = sax.circuit(
            netlist=GM_bukon,
            models={
                "layer_type1": pic_sublayer_type1_bukon,
                "layer_type2": pic_sublayer_type2_bukon
            }
        )

        return eve
