###
#   This module is a collection of cirq quantum circuits for the sequential
#   emission/preparation of some graph states of interest
###

### loading some moduels
import numpy as np

import cirq

import qutrit_utils

###
#   Utilities to build the circuit to prepare different states
###
#---------------------------------------------
#           1D-Cluster state
#---------------------------------------------
##-- prepare the ideal cluster state in 1D
def ideal_cluster_state_1D(Nqubits):
    ##-- defining the qubits
    register = [cirq.LineQid(0, dimension=3),
            *cirq.LineQubit.range(1, Nqubits+1)]

    ##-- initializing the circuit
    ics_circuit = cirq.Circuit()

    for ii in range(Nqubits, 1, -1):
        #- hadamard
        ry = qutrit_utils.Q3_H(np.pi).on(register[0])
        ics_circuit.append(ry)

        #- cnot
        pi_ef = qutrit_utils.Q3_PI_ef(np.pi).on(register[0])
        ics_circuit.append(pi_ef)
        cn = qutrit_utils.Q3_CNOT(0.0, 0.0).on(register[0], register[ii])
        ics_circuit.append(cn)
    #
    #-hadamard
    ha = qutrit_utils.Q3_H(np.pi).on(register[0])
    ics_circuit.append(ha)

    #-swap
    sw = qutrit_utils.Q3_SWAP().on(register[0], register[1])
    ics_circuit.append(sw)
    return ics_circuit
#

##-- prepare the noisy cluster state in 1D
def noisy_cluster_state_1D(Nqubits, wait_ts, ctimes_s1, noise_params):
    ##-- parameters
    t1, t2, t3, t4, t5 = wait_ts
    l1_cnot, sq_gamma = noise_params

    ##-- defining the qubits
    register = [cirq.LineQid(0, dimension=3),
            *cirq.LineQubit.range(1, Nqubits+1)]

    ##-- initializing the circuit
    cs_circuit = cirq.Circuit()

    for ii in range(Nqubits, 1, -1):
        #- hadamard
        qutrit_utils.Q3_idle_time(t1/2, ctimes_s1, register[0], cs_circuit)
        ry = qutrit_utils.Q3_H(np.pi - sq_gamma).on(register[0])
        cs_circuit.append(ry)
        qutrit_utils.Q3_idle_time(t1/2, ctimes_s1, register[0], cs_circuit)

        #- cnot
        #pi pulse
        qutrit_utils.Q3_idle_time(t2/2, ctimes_s1, register[0], cs_circuit)
        pi_ef = qutrit_utils.Q3_PI_ef(np.pi - sq_gamma).on(register[0])
        cs_circuit.append(pi_ef)
        qutrit_utils.Q3_idle_time(t2/2, ctimes_s1, register[0], cs_circuit)

        # excitation transfer
        qutrit_utils.Q3_idle_time(t3/2, ctimes_s1, register[0], cs_circuit)
        cn = qutrit_utils.Q3_CNOT(l1_cnot, 0.0).on(register[0], register[ii])
        cs_circuit.append(cn)
        qutrit_utils.Q3_idle_time(t3/2 + t4, ctimes_s1, register[0], cs_circuit)
    #
    ## The last emission
    #-hadamard
    qutrit_utils.Q3_idle_time(t1/2, ctimes_s1, register[0], cs_circuit)
    ha = qutrit_utils.Q3_H(np.pi - sq_gamma).on(register[0])
    cs_circuit.append(ha)
    qutrit_utils.Q3_idle_time(t1/2, ctimes_s1, register[0], cs_circuit)

    #-swap
    qutrit_utils.Q3_idle_time(t5/2, ctimes_s1, register[0], cs_circuit)
    sw = qutrit_utils.Q3_SWAP().on(register[0], register[1])
    cs_circuit.append(sw)
    qutrit_utils.Q3_idle_time(t5/2, ctimes_s1, register[0], cs_circuit)

    return cs_circuit
#

#------------------------------------------------------------------------------
#
#                     2D-Cluster state / ladder graph states
#
#
#                     graphs of the form 2 by n
#
#------------------------------------------------------------------------------

##-- prepare ideal cluster state in 2D
def leak_cluster_state_2D(Nqubits, gamma, l1_cz, l1_cnot):

    ##-- defining the qubits
    half = int(Nqubits/2)
    storages = cirq.LineQid.range(0, 2, dimension=3)

    qubits = cirq.LineQubit.range(2, Nqubits+2)
    row1 = qubits[::2]
    row2 = qubits[1::2]

    ##-- initializing the circuit
    cs_circuit = cirq.Circuit()

    for ii in range(half-1, 0, -1):
        #-- the hadamards
        ry = qutrit_utils.Q3_H(np.pi).on_each(storages[0], storages[1])
        cs_circuit.append(ry)

        cz = qutrit_utils.Q3_CZ(np.pi - gamma, l1_cz, 0.0).on(storages[0], storages[1])
        cs_circuit.append(cz)

        ##-- add the CNOT
        pi_ef1 = qutrit_utils.Q3_PI_ef(np.pi).on(storages[0])
        cs_circuit.append(pi_ef1)
        cn1 = qutrit_utils.Q3_CNOT(l1_cnot, 0.0).on(storages[0], row1[ii])
        cs_circuit.append(cn1)

        pi_ef2 = qutrit_utils.Q3_PI_ef(np.pi).on(storages[1])
        cs_circuit.append(pi_ef2)
        cn2 = qutrit_utils.Q3_CNOT(l1_cnot, 0.0).on(storages[1], row2[ii])
        cs_circuit.append(cn2)

    #
    #-- the hadamards
    ry = qutrit_utils.Q3_H(np.pi).on_each(storages[0], storages[1])
    cs_circuit.append(ry)

    cz = qutrit_utils.Q3_CZ(np.pi - gamma, l1_cz, 0.0).on(storages[0], storages[1])
    cs_circuit.append(cz)

    sw1 = qutrit_utils.Q3_SWAP().on(storages[0], row1[0])
    cs_circuit.append(sw1)

    sw2 = qutrit_utils.Q3_SWAP().on(storages[1], row2[0])
    cs_circuit.append(sw2)

    return cs_circuit
#

##-- prepare noisy cluster state in 2D
def noisy_cluster_state_2D(Nqubits, wait_ts, ctimes_s1, ctimes_s2, noise_params):

    ##-- parameters
    t1, t2, t3, t4, t5, tsw1, tsw2 = wait_ts
    gamma, l1_cz, l1_cnot, sq_gamma = noise_params

    ##-- defining the qubits
    half = int(Nqubits/2)
    storages = cirq.LineQid.range(0, 2, dimension=3)

    qubits = cirq.LineQubit.range(2, Nqubits+2)
    row1 = qubits[::2]
    row2 = qubits[1::2]

    ##-- initializing the circuit
    cs_circuit = cirq.Circuit()

    for ii in range(half-1, 0, -1):
        #-- the hadamards
        qutrit_utils.Q3_idle_time(t1/2, ctimes_s1, storages[0], cs_circuit)
        qutrit_utils.Q3_idle_time(t1/2, ctimes_s2, storages[1], cs_circuit)

        ry = qutrit_utils.Q3_H(np.pi - sq_gamma).on_each(storages[0], storages[1])
        cs_circuit.append(ry)

        qutrit_utils.Q3_idle_time(t1/2, ctimes_s1, storages[0], cs_circuit)
        qutrit_utils.Q3_idle_time(t1/2, ctimes_s2, storages[1], cs_circuit)

        #-- the CZ
        qutrit_utils.Q3_idle_time(t2/2, ctimes_s1, storages[0], cs_circuit)
        qutrit_utils.Q3_idle_time(t2/2, ctimes_s2, storages[1], cs_circuit)

        cz = qutrit_utils.Q3_CZ(np.pi - gamma, l1_cz, 0.0).on(storages[0], storages[1])
        cs_circuit.append(cz)

        qutrit_utils.Q3_idle_time(t2/2, ctimes_s1, storages[0], cs_circuit)
        qutrit_utils.Q3_idle_time(t2/2, ctimes_s2, storages[1], cs_circuit)

        ##-- add the CNOT
        #- first storage
        qutrit_utils.Q3_idle_time(t3/2, ctimes_s1, storages[0], cs_circuit)
        pi_ef1 = qutrit_utils.Q3_PI_ef(np.pi - sq_gamma).on(storages[0])
        cs_circuit.append(pi_ef1)
        qutrit_utils.Q3_idle_time(t3/2, ctimes_s1, storages[0], cs_circuit)

        qutrit_utils.Q3_idle_time(t4/2, ctimes_s1, storages[0], cs_circuit)
        cn1 = qutrit_utils.Q3_CNOT(l1_cnot, 0.0).on(storages[0], row1[ii])
        cs_circuit.append(cn1)
        qutrit_utils.Q3_idle_time(t4/2 + t5, ctimes_s1, storages[0], cs_circuit)

        qutrit_utils.Q3_idle_time(t3/2, ctimes_s2, storages[1], cs_circuit)
        pi_ef2 = qutrit_utils.Q3_PI_ef(np.pi - sq_gamma).on(storages[1])
        cs_circuit.append(pi_ef2)
        qutrit_utils.Q3_idle_time(t3/2, ctimes_s2, storages[1], cs_circuit)

        qutrit_utils.Q3_idle_time(t4/2, ctimes_s2, storages[1], cs_circuit)
        cn2 = qutrit_utils.Q3_CNOT(l1_cnot, 0.0).on(storages[1], row2[ii])
        cs_circuit.append(cn2)
        qutrit_utils.Q3_idle_time(t4/2 + t5, ctimes_s2, storages[1], cs_circuit)
    #
    #-- the hadamards
    qutrit_utils.Q3_idle_time(t1/2, ctimes_s1, storages[0], cs_circuit)
    qutrit_utils.Q3_idle_time(t1/2, ctimes_s2, storages[1], cs_circuit)

    ry = qutrit_utils.Q3_H(np.pi - sq_gamma).on_each(storages[0], storages[1])
    cs_circuit.append(ry)

    qutrit_utils.Q3_idle_time(t1/2, ctimes_s1, storages[0], cs_circuit)
    qutrit_utils.Q3_idle_time(t1/2, ctimes_s2, storages[1], cs_circuit)

    #-- the cz
    qutrit_utils.Q3_idle_time(t2/2, ctimes_s1, storages[0], cs_circuit)
    qutrit_utils.Q3_idle_time(t2/2, ctimes_s2, storages[1], cs_circuit)

    cz = qutrit_utils.Q3_CZ(np.pi - gamma, l1_cz, 0.0).on(storages[0], storages[1])
    cs_circuit.append(cz)

    qutrit_utils.Q3_idle_time(t2/2, ctimes_s1, storages[0], cs_circuit)
    qutrit_utils.Q3_idle_time(t2/2, ctimes_s2, storages[1], cs_circuit)

    #-- the swaps
    qutrit_utils.Q3_idle_time(tsw1/2, ctimes_s1, storages[0], cs_circuit)
    sw1 = qutrit_utils.Q3_SWAP().on(storages[0], row1[0])
    cs_circuit.append(sw1)
    qutrit_utils.Q3_idle_time(tsw1/2, ctimes_s1, storages[0], cs_circuit)

    qutrit_utils.Q3_idle_time(tsw2/2, ctimes_s2, storages[1], cs_circuit)
    sw2 = qutrit_utils.Q3_SWAP().on(storages[1], row2[0])
    cs_circuit.append(sw2)
    qutrit_utils.Q3_idle_time(tsw2/2, ctimes_s2, storages[1], cs_circuit)

    return cs_circuit
#

#------------------------------------------------------------------------------
#
#                     2D-Cluster state / grid graph states
#
#
#                     graphs of the form 3 by n
#
#------------------------------------------------------------------------------

##-- prepare ideal cluster state in 2D
def leak_cluster_state_2D_3S(Nqubits, gamma, l1_cz, l1_cnot):

    ##-- defining the qubits
    half = int(Nqubits/3)
    storages = cirq.LineQid.range(0, 3, dimension=3)

    qubits = cirq.LineQubit.range(3, Nqubits+3)
    row1 = qubits[::3]
    row2 = qubits[1::3]
    row3 = qubits[2::3]

    ##-- initializing the circuit
    cs_circuit = cirq.Circuit()

    for ii in range(half-1, 0, -1):
        #-- the hadamards
        ry = qutrit_utils.Q3_H(np.pi).on_each(storages[0], storages[1], storages[2])
        cs_circuit.append(ry)

        #- CZ between sources 1 and 2
        cz12 = qutrit_utils.Q3_CZ(np.pi - gamma, l1_cz, 0.0).on(storages[0], storages[1])
        cs_circuit.append(cz12)

        #- cz between sources 2 and 3
        cz23 = qutrit_utils.Q3_CZ(np.pi - gamma, l1_cz, 0.0).on(storages[1], storages[2])
        cs_circuit.append(cz23)

        ##-- add the CNOT
        #- first source-emitter pair
        pi_ef1 = qutrit_utils.Q3_PI_ef(np.pi).on(storages[0])
        cs_circuit.append(pi_ef1)
        cn1 = qutrit_utils.Q3_CNOT(l1_cnot, 0.0).on(storages[0], row1[ii])
        cs_circuit.append(cn1)

        #- second source-emitter pair
        pi_ef2 = qutrit_utils.Q3_PI_ef(np.pi).on(storages[1])
        cs_circuit.append(pi_ef2)
        cn2 = qutrit_utils.Q3_CNOT(l1_cnot, 0.0).on(storages[1], row2[ii])
        cs_circuit.append(cn2)

        #- third source-emitter pair
        pi_ef3 = qutrit_utils.Q3_PI_ef(np.pi).on(storages[2])
        cs_circuit.append(pi_ef3)
        cn3 = qutrit_utils.Q3_CNOT(l1_cnot, 0.0).on(storages[2], row3[ii])
        cs_circuit.append(cn3)
    #
    #-- the hadamards
    ry = qutrit_utils.Q3_H(np.pi).on_each(storages[0], storages[1], storages[2])
    cs_circuit.append(ry)

    #- CZ between sources 1 and 2
    cz12 = qutrit_utils.Q3_CZ(np.pi - gamma, l1_cz, 0.0).on(storages[0], storages[1])
    cs_circuit.append(cz12)

    #- cz between sources 2 and 3
    cz23 = qutrit_utils.Q3_CZ(np.pi - gamma, l1_cz, 0.0).on(storages[1], storages[2])
    cs_circuit.append(cz23)

    sw1 = qutrit_utils.Q3_SWAP().on(storages[0], row1[0])
    cs_circuit.append(sw1)

    sw2 = qutrit_utils.Q3_SWAP().on(storages[1], row2[0])
    cs_circuit.append(sw2)

    sw3 = qutrit_utils.Q3_SWAP().on(storages[2], row3[0])
    cs_circuit.append(sw3)

    return cs_circuit
#

##-- prepare noisy cluster state in 2D
def noisy_cluster_state_2D_3S(Nqubits, wait_ts, ctimes_s1, ctimes_s2, ctimes_s3, noise_params):

    ##-- parameters
    t1, t2, t3, t4, t5, tsw1, tsw2, tsw3 = wait_ts
    gamma, l1_cz, l1_cnot, sq_gamma = noise_params

    ##-- defining the qubits
    half = int(Nqubits/3)
    storages = cirq.LineQid.range(0, 3, dimension=3)

    qubits = cirq.LineQubit.range(3, Nqubits+3)
    row1 = qubits[::3]
    row2 = qubits[1::3]
    row3 = qubits[2::3]

    ##-- initializing the circuit
    cs_circuit = cirq.Circuit()

    for ii in range(half-1, 0, -1):
        #-- the hadamards
        qutrit_utils.Q3_idle_time(t1/2, ctimes_s1, storages[0], cs_circuit)
        qutrit_utils.Q3_idle_time(t1/2, ctimes_s2, storages[1], cs_circuit)
        qutrit_utils.Q3_idle_time(t1/2, ctimes_s3, storages[2], cs_circuit)
        ry = qutrit_utils.Q3_H(np.pi - sq_gamma).on_each(storages[0], storages[1], storages[2])
        cs_circuit.append(ry)
        qutrit_utils.Q3_idle_time(t1/2, ctimes_s1, storages[0], cs_circuit)
        qutrit_utils.Q3_idle_time(t1/2, ctimes_s2, storages[1], cs_circuit)
        qutrit_utils.Q3_idle_time(t1/2, ctimes_s3, storages[2], cs_circuit)

        #-- the CZ
        #- betwewn sources 1 and 2, the 3 is idling
        qutrit_utils.Q3_idle_time(t2/2, ctimes_s1, storages[0], cs_circuit)
        qutrit_utils.Q3_idle_time(t2/2, ctimes_s2, storages[1], cs_circuit)
        qutrit_utils.Q3_idle_time(t2/2, ctimes_s3, storages[2], cs_circuit)
        cz12 = qutrit_utils.Q3_CZ(np.pi - gamma, l1_cz, 0.0).on(storages[0], storages[1])
        cs_circuit.append(cz12)
        qutrit_utils.Q3_idle_time(t2/2, ctimes_s1, storages[0], cs_circuit)
        qutrit_utils.Q3_idle_time(t2/2, ctimes_s2, storages[1], cs_circuit)
        qutrit_utils.Q3_idle_time(t2/2, ctimes_s3, storages[2], cs_circuit)

        #- betwewn sources 2 and 3, the 1 is idling
        qutrit_utils.Q3_idle_time(t2/2, ctimes_s1, storages[0], cs_circuit)
        qutrit_utils.Q3_idle_time(t2/2, ctimes_s2, storages[1], cs_circuit)
        qutrit_utils.Q3_idle_time(t2/2, ctimes_s3, storages[2], cs_circuit)
        cz23 = qutrit_utils.Q3_CZ(np.pi - gamma, l1_cz, 0.0).on(storages[1], storages[2])
        cs_circuit.append(cz23)
        qutrit_utils.Q3_idle_time(t2/2, ctimes_s1, storages[0], cs_circuit)
        qutrit_utils.Q3_idle_time(t2/2, ctimes_s2, storages[1], cs_circuit)
        qutrit_utils.Q3_idle_time(t2/2, ctimes_s3, storages[2], cs_circuit)

        ##-- add the CNOT
        #- first storage
        qutrit_utils.Q3_idle_time(t3/2, ctimes_s1, storages[0], cs_circuit)
        pi_ef1 = qutrit_utils.Q3_PI_ef(np.pi - sq_gamma).on(storages[0])
        cs_circuit.append(pi_ef1)
        qutrit_utils.Q3_idle_time(t3/2, ctimes_s1, storages[0], cs_circuit)

        qutrit_utils.Q3_idle_time(t4/2, ctimes_s1, storages[0], cs_circuit)
        cn1 = qutrit_utils.Q3_CNOT(l1_cnot, 0.0).on(storages[0], row1[ii])
        cs_circuit.append(cn1)
        qutrit_utils.Q3_idle_time(t4/2 + t5, ctimes_s1, storages[0], cs_circuit)

        #- the second storage
        qutrit_utils.Q3_idle_time(t3/2, ctimes_s2, storages[1], cs_circuit)
        pi_ef2 = qutrit_utils.Q3_PI_ef(np.pi - sq_gamma).on(storages[1])
        cs_circuit.append(pi_ef2)
        qutrit_utils.Q3_idle_time(t3/2, ctimes_s2, storages[1], cs_circuit)

        qutrit_utils.Q3_idle_time(t4/2, ctimes_s2, storages[1], cs_circuit)
        cn2 = qutrit_utils.Q3_CNOT(l1_cnot, 0.0).on(storages[1], row2[ii])
        cs_circuit.append(cn2)
        qutrit_utils.Q3_idle_time(t4/2 + t5, ctimes_s2, storages[1], cs_circuit)

        #- the third storage
        qutrit_utils.Q3_idle_time(t3/2, ctimes_s3, storages[2], cs_circuit)
        pi_ef3 = qutrit_utils.Q3_PI_ef(np.pi - sq_gamma).on(storages[2])
        cs_circuit.append(pi_ef3)
        qutrit_utils.Q3_idle_time(t3/2, ctimes_s3, storages[2], cs_circuit)

        qutrit_utils.Q3_idle_time(t4/2, ctimes_s3, storages[2], cs_circuit)
        cn3 = qutrit_utils.Q3_CNOT(l1_cnot, 0.0).on(storages[2], row3[ii])
        cs_circuit.append(cn3)
        qutrit_utils.Q3_idle_time(t4/2 + t5, ctimes_s3, storages[2], cs_circuit)
    #
    #-- the hadamards
    qutrit_utils.Q3_idle_time(t1/2, ctimes_s1, storages[0], cs_circuit)
    qutrit_utils.Q3_idle_time(t1/2, ctimes_s2, storages[1], cs_circuit)
    qutrit_utils.Q3_idle_time(t1/2, ctimes_s3, storages[2], cs_circuit)
    ry = qutrit_utils.Q3_H(np.pi - sq_gamma).on_each(storages[0], storages[1], storages[2])
    cs_circuit.append(ry)
    qutrit_utils.Q3_idle_time(t1/2, ctimes_s1, storages[0], cs_circuit)
    qutrit_utils.Q3_idle_time(t1/2, ctimes_s2, storages[1], cs_circuit)
    qutrit_utils.Q3_idle_time(t1/2, ctimes_s3, storages[2], cs_circuit)

    #-- the cz
    #- betwewn sources 1 and 2, the 3 is idling
    qutrit_utils.Q3_idle_time(t2/2, ctimes_s1, storages[0], cs_circuit)
    qutrit_utils.Q3_idle_time(t2/2, ctimes_s2, storages[1], cs_circuit)
    qutrit_utils.Q3_idle_time(t2/2, ctimes_s3, storages[2], cs_circuit)
    cz12 = qutrit_utils.Q3_CZ(np.pi - gamma, l1_cz, 0.0).on(storages[0], storages[1])
    cs_circuit.append(cz12)
    qutrit_utils.Q3_idle_time(t2/2, ctimes_s1, storages[0], cs_circuit)
    qutrit_utils.Q3_idle_time(t2/2, ctimes_s2, storages[1], cs_circuit)
    qutrit_utils.Q3_idle_time(t2/2, ctimes_s3, storages[2], cs_circuit)

    #- betwewn sources 2 and 3, the 1 is idling
    qutrit_utils.Q3_idle_time(t2/2, ctimes_s1, storages[0], cs_circuit)
    qutrit_utils.Q3_idle_time(t2/2, ctimes_s2, storages[1], cs_circuit)
    qutrit_utils.Q3_idle_time(t2/2, ctimes_s3, storages[2], cs_circuit)
    cz23 = qutrit_utils.Q3_CZ(np.pi - gamma, l1_cz, 0.0).on(storages[1], storages[2])
    cs_circuit.append(cz23)
    qutrit_utils.Q3_idle_time(t2/2, ctimes_s1, storages[0], cs_circuit)
    qutrit_utils.Q3_idle_time(t2/2, ctimes_s2, storages[1], cs_circuit)
    qutrit_utils.Q3_idle_time(t2/2, ctimes_s3, storages[2], cs_circuit)

    #-- the swaps
    qutrit_utils.Q3_idle_time(tsw1/2, ctimes_s1, storages[0], cs_circuit)
    sw1 = qutrit_utils.Q3_SWAP().on(storages[0], row1[0])
    cs_circuit.append(sw1)
    qutrit_utils.Q3_idle_time(tsw1/2, ctimes_s1, storages[0], cs_circuit)

    qutrit_utils.Q3_idle_time(tsw2/2, ctimes_s2, storages[1], cs_circuit)
    sw2 = qutrit_utils.Q3_SWAP().on(storages[1], row2[0])
    cs_circuit.append(sw2)
    qutrit_utils.Q3_idle_time(tsw2/2, ctimes_s2, storages[1], cs_circuit)

    qutrit_utils.Q3_idle_time(tsw3/2, ctimes_s3, storages[2], cs_circuit)
    sw3 = qutrit_utils.Q3_SWAP().on(storages[2], row3[0])
    cs_circuit.append(sw3)
    qutrit_utils.Q3_idle_time(tsw3/2, ctimes_s3, storages[2], cs_circuit)

    return cs_circuit
#

#------------------------------------------------------------------------------
#
#                           Ring cluster state
#
#------------------------------------------------------------------------------

##-- prepare ideal ring cluster state
def leak_ring_cluster_state(Nqubits, gamma, l1_cz, l1_cnot):

    ##-- defining the qubits
    half = int(Nqubits/2)
    storages = cirq.LineQid.range(0, 2, dimension=3)

    qubits = cirq.LineQubit.range(2, Nqubits+2)
    row1 = qubits[::2]
    row2 = qubits[1::2]

    ##-- find the number of steps in the loop
    nsteps = len(row2)

    ##-- initializing the circuit
    cs_circuit = cirq.Circuit()

    ##-- create the transversal link on the right end
    #-- the hadamards
    ry = qutrit_utils.Q3_H(np.pi).on_each(storages[0], storages[1])
    cs_circuit.append(ry)

    #-- the CZ
    cz = qutrit_utils.Q3_CZ(np.pi - gamma, l1_cz, 0.0).on(storages[0], storages[1])
    cs_circuit.append(cz)

    ##-- emit the first photon on the longer side (only when sides are not same length)
    if len(row1) != len(row2):
        ##-- add the CNOT
        pi_ef1 = qutrit_utils.Q3_PI_ef(np.pi).on(storages[0])
        cs_circuit.append(pi_ef1)
        cn1 = qutrit_utils.Q3_CNOT(l1_cnot, 0.0).on(storages[0], row1[-1])
        cs_circuit.append(cn1)

        #-- the hadamards
        ry = qutrit_utils.Q3_H(np.pi).on(storages[0])
        cs_circuit.append(ry)
    #

    for ii in range(nsteps-1, 0, -1):
        #print(ii)
        ##-- add the CNOT
        pi_ef1 = qutrit_utils.Q3_PI_ef(np.pi).on(storages[0])
        cs_circuit.append(pi_ef1)
        cn1 = qutrit_utils.Q3_CNOT(l1_cnot, 0.0).on(storages[0], row1[ii])
        cs_circuit.append(cn1)

        pi_ef2 = qutrit_utils.Q3_PI_ef(np.pi).on(storages[1])
        cs_circuit.append(pi_ef2)
        cn2 = qutrit_utils.Q3_CNOT(l1_cnot, 0.0).on(storages[1], row2[ii])
        cs_circuit.append(cn2)

        #-- the hadamards
        ry = qutrit_utils.Q3_H(np.pi).on_each(storages[0], storages[1])
        cs_circuit.append(ry)
    #
    cz = qutrit_utils.Q3_CZ(np.pi - gamma, l1_cz, 0.0).on(storages[0], storages[1])
    cs_circuit.append(cz)

    sw1 = qutrit_utils.Q3_SWAP().on(storages[0], row1[0])
    cs_circuit.append(sw1)

    sw2 = qutrit_utils.Q3_SWAP().on(storages[1], row2[0])
    cs_circuit.append(sw2)

    return cs_circuit
#

##-- prepare noisy ring cluster state
def noisy_ring_cluster_state(Nqubits, wait_ts, ctimes_s1, ctimes_s2, noise_params):

    ##-- parameters
    t1, t2, t3, t4, t5, tsw1, tsw2 = wait_ts
    gamma, l1_cz, l1_cnot, sq_gamma = noise_params

    ##-- defining the qubits
    half = int(Nqubits/2)
    storages = cirq.LineQid.range(0, 2, dimension=3)

    qubits = cirq.LineQubit.range(2, Nqubits+2)
    row1 = qubits[::2]
    row2 = qubits[1::2]

    ##-- find the number of steps in the loop
    nsteps = len(row2)

    ##-- initializing the circuit
    cs_circuit = cirq.Circuit()

    #-- the hadamards
    qutrit_utils.Q3_idle_time(t1/2, ctimes_s1, storages[0], cs_circuit)
    qutrit_utils.Q3_idle_time(t1/2, ctimes_s2, storages[1], cs_circuit)
    ry = qutrit_utils.Q3_H(np.pi - sq_gamma).on_each(storages[0], storages[1])
    cs_circuit.append(ry)
    qutrit_utils.Q3_idle_time(t1/2, ctimes_s1, storages[0], cs_circuit)
    qutrit_utils.Q3_idle_time(t1/2, ctimes_s2, storages[1], cs_circuit)

    #-- the CZ
    qutrit_utils.Q3_idle_time(t2/2, ctimes_s1, storages[0], cs_circuit)
    qutrit_utils.Q3_idle_time(t2/2, ctimes_s2, storages[1], cs_circuit)
    cz = qutrit_utils.Q3_CZ(np.pi - gamma, l1_cz, 0.0).on(storages[0], storages[1])
    cs_circuit.append(cz)
    qutrit_utils.Q3_idle_time(t2/2, ctimes_s1, storages[0], cs_circuit)
    qutrit_utils.Q3_idle_time(t2/2, ctimes_s2, storages[1], cs_circuit)

    ##-- emit the first photon on the longer side (only when sides are not same length)
    if len(row1) != len(row2):
        ##-- add the CNOT
        qutrit_utils.Q3_idle_time(t3/2, ctimes_s1, storages[0], cs_circuit)
        pi_ef1 = qutrit_utils.Q3_PI_ef(np.pi).on(storages[0])
        cs_circuit.append(pi_ef1)
        qutrit_utils.Q3_idle_time(t3/2, ctimes_s1, storages[0], cs_circuit)

        qutrit_utils.Q3_idle_time(t4/2, ctimes_s1, storages[0], cs_circuit)
        cn1 = qutrit_utils.Q3_CNOT(l1_cnot, 0.0).on(storages[0], row1[-1])
        cs_circuit.append(cn1)
        qutrit_utils.Q3_idle_time(t4/2 + t5, ctimes_s1, storages[0], cs_circuit)

        ## the second storage idles while we emit the first photon out of storage 1
        qutrit_utils.Q3_idle_time(t3 + t4 + t5, ctimes_s1, storages[1], cs_circuit)

        #-- the hadamard in storage 1 after emitting the first photon
        qutrit_utils.Q3_idle_time(t1/2, ctimes_s1, storages[0], cs_circuit)
        ry = qutrit_utils.Q3_H(np.pi).on(storages[0])
        cs_circuit.append(ry)
        qutrit_utils.Q3_idle_time(t1/2, ctimes_s1, storages[0], cs_circuit)

        ## the second storage idles while we apply hadamard on storage 1
        qutrit_utils.Q3_idle_time(t1, ctimes_s1, storages[1], cs_circuit)
    #

    for ii in range(nsteps-1, 0, -1):
        ##-- add the CNOT
        #- first storage
        qutrit_utils.Q3_idle_time(t3/2, ctimes_s1, storages[0], cs_circuit)
        pi_ef1 = qutrit_utils.Q3_PI_ef(np.pi - sq_gamma).on(storages[0])
        cs_circuit.append(pi_ef1)
        qutrit_utils.Q3_idle_time(t3/2, ctimes_s1, storages[0], cs_circuit)

        qutrit_utils.Q3_idle_time(t4/2, ctimes_s1, storages[0], cs_circuit)
        cn1 = qutrit_utils.Q3_CNOT(l1_cnot, 0.0).on(storages[0], row1[ii])
        cs_circuit.append(cn1)
        qutrit_utils.Q3_idle_time(t4/2 + t5, ctimes_s1, storages[0], cs_circuit)

        #- second storage
        qutrit_utils.Q3_idle_time(t3/2, ctimes_s2, storages[1], cs_circuit)
        pi_ef2 = qutrit_utils.Q3_PI_ef(np.pi - sq_gamma).on(storages[1])
        cs_circuit.append(pi_ef2)
        qutrit_utils.Q3_idle_time(t3/2, ctimes_s2, storages[1], cs_circuit)

        qutrit_utils.Q3_idle_time(t4/2, ctimes_s2, storages[1], cs_circuit)
        cn2 = qutrit_utils.Q3_CNOT(l1_cnot, 0.0).on(storages[1], row2[ii])
        cs_circuit.append(cn2)
        qutrit_utils.Q3_idle_time(t4/2 + t5, ctimes_s2, storages[1], cs_circuit)

        #-- the hadamards
        qutrit_utils.Q3_idle_time(t1/2, ctimes_s1, storages[0], cs_circuit)
        qutrit_utils.Q3_idle_time(t1/2, ctimes_s2, storages[1], cs_circuit)
        ry = qutrit_utils.Q3_H(np.pi - sq_gamma).on_each(storages[0], storages[1])
        cs_circuit.append(ry)
        qutrit_utils.Q3_idle_time(t1/2, ctimes_s1, storages[0], cs_circuit)
        qutrit_utils.Q3_idle_time(t1/2, ctimes_s2, storages[1], cs_circuit)
    #

    #-- the cz
    qutrit_utils.Q3_idle_time(t2/2, ctimes_s1, storages[0], cs_circuit)
    qutrit_utils.Q3_idle_time(t2/2, ctimes_s2, storages[1], cs_circuit)
    cz = qutrit_utils.Q3_CZ(np.pi - gamma, l1_cz, 0.0).on(storages[0], storages[1])
    cs_circuit.append(cz)
    qutrit_utils.Q3_idle_time(t2/2, ctimes_s1, storages[0], cs_circuit)
    qutrit_utils.Q3_idle_time(t2/2, ctimes_s2, storages[1], cs_circuit)

    #-- the swaps
    qutrit_utils.Q3_idle_time(tsw1/2, ctimes_s1, storages[0], cs_circuit)
    sw1 = qutrit_utils.Q3_SWAP().on(storages[0], row1[0])
    cs_circuit.append(sw1)
    qutrit_utils.Q3_idle_time(tsw1/2, ctimes_s1, storages[0], cs_circuit)

    qutrit_utils.Q3_idle_time(tsw2/2, ctimes_s2, storages[1], cs_circuit)
    sw2 = qutrit_utils.Q3_SWAP().on(storages[1], row2[0])
    cs_circuit.append(sw2)
    qutrit_utils.Q3_idle_time(tsw2/2, ctimes_s2, storages[1], cs_circuit)

    return cs_circuit
#


#------------------------------------------------------------------------------
#
#                           Tree graph state
#
#------------------------------------------------------------------------------

##-- prepare ideal tree graph state
def leak_tree_graph_state(Nbranches, gamma, l1_cz, l1_cnot):

    ##-- defining the qubits
    Nqubits = 3*Nbranches + 1
    storages = cirq.LineQid.range(0, 2, dimension=3)

    qubits = cirq.LineQubit.range(2, Nqubits+2)
    row1 = qubits[0]
    row2 = qubits[1:]

    ##-- initializing the circuit
    cs_circuit = cirq.Circuit()

    ##-- create the transversal link on the right end
    #-- the hadamards
    ry = qutrit_utils.Q3_H(np.pi).on(storages[0])
    cs_circuit.append(ry)


    for ii in range(Nbranches, 0, -1):

        ###-- first photon in the branch
        #- hadamard
        ry = qutrit_utils.Q3_H(np.pi).on(storages[1])
        cs_circuit.append(ry)

        #- controlled emission
        pi_ef2 = qutrit_utils.Q3_PI_ef(np.pi).on(storages[1])
        cs_circuit.append(pi_ef2)
        cn2 = qutrit_utils.Q3_CNOT(l1_cnot, 0.0).on(storages[1], row2[3*ii - 1])
        cs_circuit.append(cn2)

        ##-- second photon in the branch
        #- hadamard
        ry = qutrit_utils.Q3_H(np.pi).on(storages[1])
        cs_circuit.append(ry)

        #- CZ
        cz = qutrit_utils.Q3_CZ(np.pi - gamma, l1_cz, 0.0).on(storages[0], storages[1])
        cs_circuit.append(cz)

        #- controlled emission
        pi_ef2 = qutrit_utils.Q3_PI_ef(np.pi).on(storages[1])
        cs_circuit.append(pi_ef2)
        cn2 = qutrit_utils.Q3_CNOT(l1_cnot, 0.0).on(storages[1], row2[3*ii - 2])
        cs_circuit.append(cn2)

        ##-- third photon in the branch
        #- hadamard
        ry = qutrit_utils.Q3_H(np.pi).on(storages[1])
        cs_circuit.append(ry)

        #- swap
        sw2 = qutrit_utils.Q3_SWAP().on(storages[1], row2[3*ii - 3])
        cs_circuit.append(sw2)
    #

    ##-- emit the root photon
    sw1 = qutrit_utils.Q3_SWAP().on(storages[0], row1)
    cs_circuit.append(sw1)

    return cs_circuit
#

##-- prepare noisy tree graph state
def noisy_tree_graph_state(Nbranches, wait_ts, ctimes_s1, ctimes_s2, noise_params):

    ##-- parameters
    t1, t2, t3, t4, t5, tsw1, tsw2 = wait_ts
    gamma, l1_cz, l1_cnot, sq_gamma = noise_params

    ##-- defining the qubits
    Nqubits = 3*Nbranches + 1
    storages = cirq.LineQid.range(0, 2, dimension=3)

    qubits = cirq.LineQubit.range(2, Nqubits+2)
    row1 = qubits[0]
    row2 = qubits[1:]

    ##-- initializing the circuit
    cs_circuit = cirq.Circuit()

    ##-- create the transversal link on the right end
    #-- the hadamards
    qutrit_utils.Q3_idle_time(t1/2, ctimes_s1, storages[0], cs_circuit)
    ry = qutrit_utils.Q3_H(np.pi).on(storages[0])
    cs_circuit.append(ry)
    qutrit_utils.Q3_idle_time(t1/2, ctimes_s1, storages[0], cs_circuit)

    for ii in range(Nbranches, 0, -1):

        ###-- first photon in the branch
        #- hadamard
        qutrit_utils.Q3_idle_time(t1/2, ctimes_s1, storages[1], cs_circuit)
        ry = qutrit_utils.Q3_H(np.pi).on(storages[1])
        cs_circuit.append(ry)
        qutrit_utils.Q3_idle_time(t1/2, ctimes_s1, storages[1], cs_circuit)

        #- storage 1 idles while creating the branches
        qutrit_utils.Q3_idle_time(t1, ctimes_s1, storages[0], cs_circuit)

        #- controlled emission
        qutrit_utils.Q3_idle_time(t3/2, ctimes_s1, storages[1], cs_circuit)
        pi_ef2 = qutrit_utils.Q3_PI_ef(np.pi).on(storages[1])
        cs_circuit.append(pi_ef2)
        qutrit_utils.Q3_idle_time(t3/2, ctimes_s1, storages[1], cs_circuit)

        qutrit_utils.Q3_idle_time(t4/2, ctimes_s1, storages[1], cs_circuit)
        cn2 = qutrit_utils.Q3_CNOT(l1_cnot, 0.0).on(storages[1], row2[3*ii - 1])
        cs_circuit.append(cn2)
        qutrit_utils.Q3_idle_time(t4/2 + t5, ctimes_s1, storages[1], cs_circuit)

        #- storage 1 idles while creating the branches
        qutrit_utils.Q3_idle_time(t3 + t4 + t5, ctimes_s1, storages[0], cs_circuit)

        ##-- second photon in the branch
        #- hadamard
        qutrit_utils.Q3_idle_time(t1/2, ctimes_s1, storages[1], cs_circuit)
        ry = qutrit_utils.Q3_H(np.pi).on(storages[1])
        cs_circuit.append(ry)
        qutrit_utils.Q3_idle_time(t1/2, ctimes_s1, storages[1], cs_circuit)

        #- storage 1 idles while creating the branches
        qutrit_utils.Q3_idle_time(t1, ctimes_s1, storages[0], cs_circuit)

        #- CZ
        qutrit_utils.Q3_idle_time(t2/2, ctimes_s2, storages[1], cs_circuit)
        cz = qutrit_utils.Q3_CZ(np.pi - gamma, l1_cz, 0.0).on(storages[0], storages[1])
        cs_circuit.append(cz)
        qutrit_utils.Q3_idle_time(t2/2, ctimes_s2, storages[1], cs_circuit)

        #- storage 1 idles while creating the branches
        qutrit_utils.Q3_idle_time(t2, ctimes_s1, storages[0], cs_circuit)

        #- controlled emission
        qutrit_utils.Q3_idle_time(t3/2, ctimes_s1, storages[1], cs_circuit)
        pi_ef2 = qutrit_utils.Q3_PI_ef(np.pi).on(storages[1])
        cs_circuit.append(pi_ef2)
        qutrit_utils.Q3_idle_time(t3/2, ctimes_s1, storages[1], cs_circuit)

        qutrit_utils.Q3_idle_time(t4/2 + t5, ctimes_s1, storages[1], cs_circuit)
        cn2 = qutrit_utils.Q3_CNOT(l1_cnot, 0.0).on(storages[1], row2[3*ii - 2])
        cs_circuit.append(cn2)
        qutrit_utils.Q3_idle_time(t4/2 + t5, ctimes_s1, storages[1], cs_circuit)

        #- storage 1 idles while creating the branches
        qutrit_utils.Q3_idle_time(t3 + t4 + t5, ctimes_s1, storages[0], cs_circuit)

        ##-- third photon in the branch
        #- hadamard
        qutrit_utils.Q3_idle_time(t1/2, ctimes_s1, storages[1], cs_circuit)
        ry = qutrit_utils.Q3_H(np.pi).on(storages[1])
        cs_circuit.append(ry)
        qutrit_utils.Q3_idle_time(t1/2, ctimes_s1, storages[1], cs_circuit)

        #- storage 1 idles while creating the branches
        qutrit_utils.Q3_idle_time(t1, ctimes_s1, storages[0], cs_circuit)

        #- swap
        qutrit_utils.Q3_idle_time(tsw2/2, ctimes_s2, storages[1], cs_circuit)
        sw2 = qutrit_utils.Q3_SWAP().on(storages[1], row2[3*ii - 3])
        cs_circuit.append(sw2)
        qutrit_utils.Q3_idle_time(tsw2/2, ctimes_s2, storages[1], cs_circuit)
    #

    ##-- emit the root photon
    qutrit_utils.Q3_idle_time(tsw1/2, ctimes_s1, storages[0], cs_circuit)
    sw1 = qutrit_utils.Q3_SWAP().on(storages[0], row1)
    cs_circuit.append(sw1)
    qutrit_utils.Q3_idle_time(tsw1/2, ctimes_s1, storages[0], cs_circuit)

    return cs_circuit
#
