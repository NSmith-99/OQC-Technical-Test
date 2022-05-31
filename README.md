# OQC-Technical-Test
Technical Test code for Oxford Quantum Computing job application

## About
This program was made for the technical test for Oxford Quantum Computings application process. It is an object oriented framework within python used to quickly set up a quantum computing processor hamiltonian from a data file containing information about such a system. Multiple usage options are given, as well as capabilities for Coaxamon and Rectanglemon type Qubits.

## Instructions

In order to make use of this small application, define a data file similar to that in `example_setup.dat`. Some things to note:
    - numerical values on one line should be comma seperated.
    - comments and notations on the file should be denoted using C-style `//` at the end of a line but before the comment itself.
    - whitespace is cleaned away
    - keyword/flag setting lines should be *started* with `#`. Style is these keywords should be all upper case, but this isn't enforced.

Comments detailing how the specific sections should be formatted and what values goes where is given in the comments of `example_setup.dat`.

Once the data file is defined, run the following line in the terminal,making not of your working directory,

```bash
python {$PATH_TO_DIR}/Qubit_system.py -F {$PATH_TO_DIR}/data_file.dat
```

This will run the script in data file mode, and return an expression of the Hamiltonian of the system.