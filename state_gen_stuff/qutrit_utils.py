###
#   This is a module with all the functionalities to simulate leakage, and
#   amplitude and phase damping of qutrits using Cirq.
###

### loading some moduels
import numpy as np
from scipy.linalg import expm
from math import sqrt

import itertools
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple, Union, TYPE_CHECKING

import cirq
from cirq import protocols

###
#   Qutrit functionalities
###
##-- A hadamard for the first two levels of a qutrit
class Q3_H(cirq.Gate):
    """
    A gate that implements a Hadamard between the first two levels of a qutrit.

    """
    def __init__(self, angle):
        super(Q3_H, self)
        self.angle = angle

    def _qid_shape_(self):
        # By implementing this method this gate implements the
        # cirq.qid_shape protocol and will return the tuple (3,)
        # when cirq.qid_shape acts on an instance of this class.
        # This indicates that the gate acts on a single qutrit.
        return (3,)

    def _unitary_(self):
        # Since the gate acts on three level systems it has a unitary
        # effect which is a three by three unitary matrix.
        return np.array([[1j*np.cos(0.5*self.angle) + np.sin(0.5*self.angle)/np.sqrt(2), np.sin(0.5*self.angle)/np.sqrt(2), 0],
                         [np.sin(0.5*self.angle)/np.sqrt(2), 1j*np.cos(0.5*self.angle) - np.sin(0.5*self.angle)/np.sqrt(2), 0],
                         [0, 0, 1]])

    def _circuit_diagram_info_(self, args):
        return '[H]'
#

##-- Pi pulse for the ef transition of the storages
class Q3_PI_ef(cirq.Gate):
    """
    A gate that implements a population transfer between the |e> and |f> levels
    of a qutrit.

    """

    def __init__(self, angle):
        super(Q3_PI_ef, self)
        self.angle = angle

    def _qid_shape_(self):
        # By implementing this method this gate implements the
        # cirq.qid_shape protocol and will return the tuple (3,)
        # when cirq.qid_shape acts on an instance of this class.
        # This indicates that the gate acts on a single qutrit.
        return (3,)

    def _unitary_(self):
        # Since the gate acts on three level systems it has a unitary
        # effect which is a three by three unitary matrix.
        unitary =  np.array([[1, 0, 0],
                             [0, -1j*np.cos(0.5*self.angle), np.sin(0.5*self.angle)],
                             [0, np.sin(0.5*self.angle), -1j*np.cos(0.5*self.angle)]])
        return unitary

    def _circuit_diagram_info_(self, args):
        return r"[pi_ef]"
#

##-- A CNOT between a qutrit and a qubit. Use the qubit subspace of the qutrit as control
class Q3_CNOT(cirq.Gate):
    """
    A gate that implements a population swap between the levels |f0> and |e1> of
    a qutrit-qubit system.
    We include leakage!

    """

    def __init__(self, leakage_rate, leakage_phase):
        super(Q3_CNOT, self)
        self.leakage_rate = leakage_rate
        self.leakage_phase = leakage_phase

    def _qid_shape_(self):
        # By implementing this method this gate implements the
        # cirq.qid_shape protocol and will return the tuple (3,)
        # when cirq.qid_shape acts on an instance of this class.
        # This indicates that the gate acts on a single qutrit.
        return (3, 2)

    def _unitary_(self):
        # Since the gate acts on three level systems it has a unitary
        # effect which is a three by three unitary matrix.
        cnot = np.diag(np.ones(6, dtype = complex))
        cnot[3][3] = 0
        cnot[4][4] = 0
        cnot[3][4] = 1
        cnot[4][3] = 1

        leakage = 4*self.leakage_rate

        generator = np.zeros((6,6), dtype = complex)
        generator[3][4] = -1j * np.arcsin(sqrt(leakage)) * np.exp(1j * self.leakage_phase)
        generator[4][3] = 1j * np.arcsin(sqrt(leakage)) * np.exp(-1j * self.leakage_phase)

        noisy_unitary = expm(1j * generator)

        cnot_unitary = cnot @ noisy_unitary

        return cnot_unitary

    def _circuit_diagram_info_(self, args: 'cirq.CircuitDiagramInfoArgs') -> 'cirq.CircuitDiagramInfo':
        return protocols.CircuitDiagramInfo(wire_symbols=('[CX_Q3_Q2]', '@'))
#

##-- A SWAP between the qubit subspace of a qutrit and a qubit
class Q3_SWAP(cirq.Gate):
    """
    A gate that implements a swap between the qubit subspace of a qutrit and a qubit.

    """

    def _qid_shape_(self):
        # By implementing this method this gate implements the
        # cirq.qid_shape protocol and will return the tuple (3,)
        # when cirq.qid_shape acts on an instance of this class.
        # This indicates that the gate acts on a single qutrit.
        return (3, 2)

    def _unitary_(self):
        # Since the gate acts on three level systems it has a unitary
        # effect which is a three by three unitary matrix.
        return np.array([[1, 0, 0, 0, 0, 0],
                         [0, 0, 1, 0, 0, 0],
                         [0, 1, 0, 0, 0, 0],
                         [0, 0, 0, 1, 0, 0],
                         [0, 0, 0, 0, 1, 0],
                         [0, 0, 0, 0, 0, 1]])

    def _circuit_diagram_info_(self, args: 'cirq.CircuitDiagramInfoArgs') -> 'cirq.CircuitDiagramInfo':
        return protocols.CircuitDiagramInfo(wire_symbols=('X', 'X'))
#

##-- A CZ between two qutrits
class Q3_CZ(cirq.Gate):
    """
    A gate that implements a CPHASE between two qutrits exploiting the
    |ee> to |f0> transition

    """
    def __init__(self, angle, leakage_rate, leakage_phase):
        super(Q3_CZ, self)
        self.angle = angle
        self.leakage_rate = leakage_rate
        self.leakage_phase = leakage_phase

    def _qid_shape_(self):
        # By implementing this method this gate implements the
        # cirq.qid_shape protocol and will return the tuple (3,)
        # when cirq.qid_shape acts on an instance of this class.
        # This indicates that the gate acts on a single qutrit.
        return (3, 3)

    def _unitary_(self):
        # Since the gate acts on three level systems it has a unitary
        # effect which is a three by three unitary matrix.
        phases = np.zeros(9, dtype = complex)
        phases[2] = -self.angle
        phases[4] = self.angle
        ideal_unitary = expm(1j * np.diag(phases))


        leakage = 4*self.leakage_rate

        generator = np.zeros((9,9), dtype = complex)
        generator[2][4] = 1j * np.arcsin(sqrt(leakage)) * np.exp(1j * self.leakage_phase)
        generator[4][2] = -1j * np.arcsin(sqrt(leakage)) * np.exp(-1j * self.leakage_phase)

        noisy_unitary = expm(1j * generator)

        cz_unitary = ideal_unitary @ noisy_unitary

        return cz_unitary

    def _circuit_diagram_info_(self, args: 'cirq.CircuitDiagramInfoArgs') -> 'cirq.CircuitDiagramInfo':
            return protocols.CircuitDiagramInfo(wire_symbols=('[Q3_CZ]', '@'))
#

##-- An amplitude damping channel for qutrits
class Q3_AmplitudeDampingChannel(cirq.Gate):

    def __init__(self, pro1: float, pro2: float) -> None:
        """Construct a channel that dampens qubit phase.

        Args:
            gamma: The damping constant.
        Raises:
            ValueError: if gamma is not a valid probability.
        """
        self.pro1 = pro1
        self.pro2 = pro2
    #

    def _qid_shape_(self):
        return (3,)
    #

    def _num_qubits_(self) -> int:
        return 1

    def _kraus_(self) -> Iterable[np.ndarray]:
        return (
            np.array([[0, np.sqrt(self.pro1), 0], [0, 0, 0], [0, 0, 0]]),    # decay |e) -> |g)
        np.array([[0, 0, 0], [0, 0, np.sqrt(self.pro2)], [0, 0, 0]]),    # decay |f) -> |e)
        np.array([[1, 0, 0], [0, np.sqrt(1-self.pro1), 0], [0, 0, np.sqrt(1-self.pro2)]])
        )

    def _has_kraus_(self) -> bool:
        return True
    def _circuit_diagram_info_(self, args):
        return '[Q3_AD]'
#

##-- A phase damping channel for qutrits
class Q3_PhaseDampingChannel(cirq.Gate):

    def __init__(self, pro1: float, pro2: float) -> None:
        """Construct a channel that dampens qubit phase.

        Args:
            gamma: The damping constant.
        Raises:
            ValueError: if gamma is not a valid probability.
        """
        self.pro1 = pro1
        self.pro2 = pro2
    #

    def _qid_shape_(self):
        return (3,)
    #

    def _num_qubits_(self) -> int:
        return 1

    def _kraus_(self) -> Iterable[np.ndarray]:
        return (
            np.array([[0, 0, 0], [0, np.sqrt(self.pro1), 0], [0, 0, 0]]),
            np.array([[0, 0, 0], [0, 0, 0], [0, 0, np.sqrt(self.pro2)]]),
            np.array([[1, 0, 0], [0, np.sqrt(1-self.pro1), 0], [0, 0, np.sqrt(1-self.pro2)]])
        )

    def _has_kraus_(self) -> bool:
        return True
    def _circuit_diagram_info_(self, args):
        return '[Q3_PD]'
#

###
#    Some other utilities
###
##-- amplitude and phase dmaping probs
def pad(t, T):
    return 1-np.exp(-t/T)
def ppd(t, T1, T2):
    return 1-np.exp(t/T1)*np.exp(-2*t/T2)

##-- defining a damping during idle time
def Q3_idle_time(time, coherence_times, qubit, circuit: cirq.Circuit):
    T1_e, T2_e, T1_f, T2_f = coherence_times

    pad_e = pad(time, T1_e)
    ppd_e = ppd(time, T1_e, T2_e)
    pad_f = pad(time, T1_f)
    ppd_f = ppd(time, T1_f, T2_f)

    circuit.append(Q3_AmplitudeDampingChannel(pad_e, pad_f).on(qubit))
    circuit.append(Q3_PhaseDampingChannel(ppd_e, ppd_f).on(qubit))
#
