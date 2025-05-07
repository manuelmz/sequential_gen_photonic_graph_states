# sequential_gen_photonic_graph_states
This repository contains the tools for the noisy quantum system simulations used in **Detecting genuine multipartite entanglement in multi-qubit devices with restricted measurements** [arXiv:2504.21076](https://arxiv.org/abs/2504.21076). 

The simulations were implented using the open source library Cirq. The module qutrit_utils.py contains useful functions to simulate amplitue and phase damping on qutrits, and single-qutrit, two-qutrit, single-qubit, and qubit-qutrit gates, under the effect of coherent errors, i.e., leakage and under/over rotations. The module graph_state_gen_circuits.py contains functions to build sequential generation circuits for several graph states of interest, namely path, ring, tree, and 2D graph states. Here, a certain number of source qutrits is use to sequentially prepare the desire graph state on a register of $N$ qubits. 

The different notebooks illustrate the usage of these two modules for the generation of these graph states with different number of qubits. 

