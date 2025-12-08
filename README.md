# sequential_gen_photonic_graph_states
This repository contains the tools for the noisy quantum system simulations used in **Detecting genuine multipartite entanglement in multi-qubit devices with restricted measurements** [arXiv:2504.21076](https://arxiv.org/abs/2504.21076). 

The main functionalities are found inside the folder state_gen_stuff in the form of two modules, qutrit_utils.py and graph_state_gen_circuits.py, which are implemented using the open source librayr Cirq.

-  qutrit_utils.py contains useful functions to simulate amplitue and phase damping on qutrits, and single-qutrit, two-qutrit, single-qubit, and qubit-qutrit gates, under the effect of coherent errors, i.e., leakage and under/over rotations.
-  graph_state_gen_circuits.py contains functions to build sequential generation circuits for several graph states of interest, namely path, ring, tree, and 2D graph states. Here, a certain number of source qutrits is use to sequentially prepare the desire graph state on a register of $N$ qubits. 

The different notebooks inside the folder state_gen_stuff illustrate the usage of qutrit_utils.py and graph_state_gen_circuits.py for the generation of the simulated graph states with different number of qubits.

The folder expectation_values contains notebooks used to evaluate expectation values of stabilizers on the simulated noisy graph states.

- 

