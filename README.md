#====2019-07-16====

basic input data are:
(1) device hamiltonian calculated by cp2k
(2) device coordinate
(3) electrode hamiltonian calculated by cp2k
(4) electrode coordinate

system specific input data are:
(1) basis set number the of deleted part of device (in device hamiltonian)
(2) basis set number of the lead principle layer (in device hamiltonian)
(3) total basis set number of the electrode (in electrode hamiltonian)

kwant-transmission.py is used for calculating the transmission data

kwant-current-density.py is used for calculating the current density data

plot-3d-quiver.py is used for generating a 3d current density figure with the output of kwant-current-density.py

#====end 2019-07-16====
