# sequential_gen_photonic_graph_states
This repository contains the tools for the noisy quantum system simulations used in **Detecting genuine multipartite entanglement in multi-qubit devices with restricted measurements** [arXiv:2504.21076](https://arxiv.org/abs/2504.21076). 

The main functionalities are found inside the folder state_gen_stuff in the form of two modules, qutrit_utils.py and graph_state_gen_circuits.py, which are implemented using the open source librayr Cirq.

-  qutrit_utils.py contains useful functions to simulate amplitue and phase damping on qutrits, and single-qutrit, two-qutrit, single-qubit, and qubit-qutrit gates, under the effect of coherent errors, i.e., leakage and under/over rotations.
-  graph_state_gen_circuits.py contains functions to build sequential generation circuits for several graph states of interest, namely path, ring, tree, and 2D graph states. Here, a certain number of source qutrits is use to sequentially prepare the desire graph state on a register of $N$ qubits. 

The different notebooks inside the folder state_gen_stuff illustrate the usage of qutrit_utils.py and graph_state_gen_circuits.py for the generation of the simulated graph states with different number of qubits.

The folder expectation_values contains notebooks used to evaluate expectation values of stabilizers on the simulated noisy graph states.

- test_witness_on_dms.ipynb is used to compute the expectation value of the stabilizer generators on nodes and their products on edges, which are neecssary to evaluate the entanglement witness, on the noisy density matrices generated with the notbeooks in state_gen_stuff

- the notbeooks build_all_stabilizers_graph_state.ipynb uses the julia pakcage PauliStrings.jl to generate all the operators in the stabilizer group for a given target graph state and a given number of qubits. The Pauli strings representing these operators are then exported.
  
- the notebook evaluate_all_stabilizers.ipynb uses the pauli strings of the stabilizer group generated with build_all_stabilizers_graph_state.ipynb and evaluates the expectation value of all the stabilizers in the stabilizer group of the target graph states on the simulated noisy states.

- Both build_all_stabilizers_graph_state.ipynb and evaluate_all_stabilizers.ipynb should be used together.

The overall workflow is:
1. use one of the notebooks in state_gen_stuff to generate the density matrix for a given target graph state and number of qubits
2. use test_witness_on_dms.ipynb to evaluate the expectation values of the stabilizers involved in the witness evaluation
3. evaluate the witness from these expectation values
