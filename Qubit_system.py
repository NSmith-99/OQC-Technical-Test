## Python script defining objects for use in OQC Technical tests

import numpy as np
import sys
import os


class qubit:
    ## Base class for qubits, inherited by more special qubit types
    frequency = 0
    def __init__(self, omega):
        # Constructor, defines the object
        self.frequency = omega

    def omega(self):
        # symbolic way to retreive qubit frequency
        return self.frequency

class coaxmon(qubit):
    # inherited class specialised for coaxamon type qubits.
    def __init__(self, r_c, r_in, r_out, omega = None):
        # constructor takes qubit parameters and makes internal definitions
        self.center_radius = r_c
        self.inner_radius = r_in
        self.outer_radius = r_out

        # Calculates areas
        A_i = np.pi*r_c**2
        A_r = np.pi*(r_out**2 - r_in**2)

        # If no overide, assigns coaxamon frequency via area ratio
        if omega == None:
            super().__init__(A_i/A_r)
        else:
            super().__init__(omega)       

class rectanglemon(qubit):
    # inherited class specialised for rectanglemon type qubits.
    def __init__(self, H1, L1, H2, L2, omega = None):
        # constructor takes parameters and makes internal definitions 
        self.height_1 = H1
        self.length_1 = L1
        self.height_2 = H2
        self.length_2 = L2

        # calculates areas 
        A1 = H1*L1
        A2 = H2*L2

        # if no override, assigns rectanglemon frequency via area sum.
        if omega == None:
            super().__init__(A1 + A2)
        else:
            super().__init__(omega)

class Q_processor:
    # Class that holds the overall system/quantum processor, containing all information about the qubits and their interactions

    # data declarations
    qubits = []
    matrix = np.array([[]])

    def __init__(self):
        # Default constructor
        return

    def __init__(self, system_matrix, types = None):
        # constructor, takes in matrix of system parameters and qubit types

        # asserts matrix is correct
        assert isinstance(system_matrix, np.ndarray) and len(system_matrix.shape) == 2 and system_matrix.shape[0] == system_matrix.shape[1], "System matrix must be a square matrix of complex numbers (NxN complex numpy array)!"

        # if types is given as an array
        if isinstance(types, np.ndarray):
            # asserts this array is correct
            assert types.shape == (system_matrix.shape[0], 5) and np.all(np.logical_and(np.greater_equal(types[:,0], 0), np.less_equal(types[:,0], 2))), "Types array must be an list made up of 0's, 1's or 2's!"

            # for each type given, checks for an override and assigns frequency accordingly.
            for i, ty in enumerate(types):
                # Checks if matrix has frequency overide, assigns None if not (or normal frequency if standard qubit), sets override otherwise
                if system_matrix[i,i] == 0:
                    freq = None
                    if ty[0] == 0:
                        freq = ty[1]
                else:
                    freq = system_matrix[i,i]

                # Appends correct qubit to the list using given values.
                if ty[0] == 0:
                    self.qubits.append(qubit(freq))
                elif ty[0] == 1:
                    self.qubits.append(coaxmon(ty[1], ty[2], ty[3], omega = freq))
                elif ty[0] == 2:
                    self.qubits.append(rectanglemon(ty[1], ty[2], ty[3], ty[4], omega = freq))
    
        else:
            # if no types array is given, construct system purely via matrix
            for i in range(system_matrix.shape[0]):
                self.qubits.append(qubit(system_matrix[i,i]))

        # takes in matrix either as upper-triangular or as self-transpose, then assigns to class value such that it is always self-transpose.
        if np.any(np.not_equal(np.tril(system_matrix, -1), np.zeros_like(system_matrix))):
            assert np.all(np.transpose(system_matrix) == system_matrix), "System matrix must be either upper-triangular (elements below diagonal all zero) or self transpose (i.e. element J_n,m = J_m,n)!"
            self.matrix = system_matrix
        else:
            self.matrix = system_matrix + np.transpose(np.triu(system_matrix, 1))

        # modifies diagonal elements to final values of qubit frequency
        for i, q in enumerate(self.qubits):
            self.matrix[i,i] = q.frequency

    def add_qubit(self, frequency, couplings, type = 0):
        # function to add an arbitrary qubit to the system. only takes in a frequency value, for all qubit types. defaults to standard.
        if type == 0:
            self.qubits.append(qubit(frequency))
        elif type == 1:
            self.qubits.append(coaxmon(0,0,0, omega = frequency))
        elif type == 2:
            self.qubits.append(rectanglemon(0,0,0,0, omega = frequency))
        else:
            print("Given qubit type value does not exist!")
            return
        
        # inserts couplings in to system matrix
        np.append(self.matrix, couplings, axis = 0)
        np.append(self.matrix, np.append(couplings, frequency), axis = 1)

    def add_coaxmon(self, r_in, r_out, r_c, couplings):
        # adds a coaxamon with correct parameters to the system, alongside couplings
        self.qubits.append(coaxmon(r_c, r_in,r_out))
        np.append(self.matrix, couplings, axis = 0)
        np.append(self.matrix, np.append(couplings, self.qubits[-1].frequency), axis = 1)
        return

    def add_rectanglemon(self, H1, L1, H2, L2, couplings):
        # adds a rectanglemon with correct parameters to the system, alongside couplings
        self.qubits.append(rectanglemon(H1, L1, H2, L2))
        np.append(self.matrix, couplings, axis = 0)
        np.append(self.matrix, np.append(couplings, self.qubits[-1].frequency), axis = 1)
        return
    
    def hamiltonian(self):

        return

    def hamiltonian_string(self):
        # prints out formatted string stating overall Hamiltonian operator of system.
        output = "H = "

        # For each qubit, append sigma_z term to operator
        for i, qubit in enumerate(self.qubits):
            output = output + "{}Z{} + ".format(qubit.frequency/2, i)

        # Takes just one copy of each coupling term, where non-zero adds term to expression
        couplings = np.triu(self.matrix, 1)
        args = np.argwhere(couplings != 0)
        for arg in args[:-1]:
            output = output + "{}X{}X{} + ".format(self.matrix[arg[0], arg[1]], arg[0], arg[1])

        # adds final coupling term without trailing addition symbol.
        output = output + "{}X{}X{}".format(self.matrix[args[-1, 0], args[-1, 1]], args[-1, 0], args[-1, 1])

        return output

# when run from command line.
if __name__ == "__main__":

    # Checks if using file is specified
    if sys.argv[1] == "-F" or sys.argv[1] == "--file":
        # finds file name and checls validity, then opens
        file_name = sys.argv[2]
        assert os.path.exists(file_name), "File given does not exist!"
        file = open(file_name)

        # value initialisation
        INPUT_TYPE_READ = False
        INPUT_TYPE = 0
        QUBITS_READ = False
        MATRIX_READ = False

        type_array = np.zeros((0,5))
        matrix = []

        # for each line
        for line in file:
            formatted_line = line.split("//")[0].strip() # Removes c-style comments and leading/trailing whitespace in the data file from lines
            
            if not formatted_line: continue # If blank line

            if formatted_line[0] == '#':
                # sets flags if "#" is used
                if formatted_line[1:] == "INPUT_TYPE":
                    INPUT_TYPE_READ = True
                    QUBITS_READ = False
                    MATRIX_READ = False
                    continue
                elif formatted_line[1:] == "QUBITS":
                    INPUT_TYPE_READ = False
                    QUBITS_READ = True
                    MATRIX_READ = False
                    continue
                elif formatted_line[1:] == "MATRIX":
                    INPUT_TYPE_READ = False
                    QUBITS_READ = False
                    MATRIX_READ = True
                    continue
                else:
                    print("Keyword in file not valid!")
                    input("Press Enter key to exit...")
                    sys.exit()

            # if an input read line, assigns int flag
            if INPUT_TYPE_READ:
                INPUT_TYPE = int(formatted_line)

            # if a qubit read line, splits the array, casts to floats, appends trailing zeros for same size arrays. then appends result to type list
            elif QUBITS_READ:
                array = formatted_line.split(",")
                array = [float(x) for x in array]
                while len(array) < 5:
                    array.append(0)
                type_array = np.append(type_array, np.array([array]), axis = 0)
            
            # if a matrix read line, splits the line, casts to floats, then appends this to the system matrix
            elif MATRIX_READ:
                array = formatted_line.split(",")
                array = [float(x) for x in array]
                matrix.append(array)
        
        matrix = np.array(matrix)
        
        # constructs system class based on input type value
        system = 0
        if INPUT_TYPE == 0:
            system = Q_processor(matrix, type_array)
        elif INPUT_TYPE == 1:
            system = Q_processor(matrix)
        else:
            print("Given input type is not valid!!")
            input("Press Enter key to exit...")
            sys.exit()

        # Nicely prints out required expression for the Hamiltonian.
        print("The Hamiltonian for the given system is:")
        print(system.hamiltonian_string())
        input("Press Enter key to exit...")


###### Code from here was attempt to create a manual entry system for Qubits. This is not particularly required and was taking too much time to make nicely. Usage of the data file method should be more than enough. ######

    # If manual flag given
    elif sys.argv[1] in ["-M", "--manual"]:
        # set up to begin taking in commands
        system = Q_processor(np.zeros((0,0))) # default constructor for Q_processor isn't working for some reason?????

        # input loop
        while True:
            # If previously entered qubits, first defines couplings for the next one
            couplings = np.array([])
            if len(system.qubits) > 0:
                not_entered = True
                coupling_list = []
                #input loop
                while not_entered:
                    coupling_list = input("Please enter couplings to all previous Qubits, seperating with a space.").split()
                    if len(coupling_list) == len(system.qubits):
                        print("number of coupling values not correct! Please try again.")
                    else:
                        not_entered = False
                couplings = np.array([float(x) for x in coupling_list])

            # asks user which qubit to add
            print("Please enter which Qubit type to add.")
            print("0: Standard Qubit.")
            print("1: Coaxamon Qubit.")
            print("2: Rectanglemon Qubit.")
            answer = int(input("Please enter the number of the Qubit"))

            # asks for required parameters, then adds new qubit with couplings.
            if answer == 0:
                frequency = float(input("Please enter the qubit frequency."))
                system.add_qubit(frequency, couplings)
            elif answer == 1:
                r_c = float(input("Please enter the Coaxamon centre radius."))
                r_in = float(input("Please enter the Coaxamon inner radius."))
                r_out = float(input("Please enter the Coaxamon outer radius."))
                system.add_coaxmon(r_in, r_out, r_c, couplings)
            elif answer == 2:
                h1 = float(input("Please enter rectanglemon height 1."))
                l1 = float(input("Please enter rectanglemon length 1."))
                h2 = float(input("Please enter rectanglemon height 2."))
                l2 = float(input("Please enter rectanglemon length 1."))
                system.add_rectanglemon(h1, l1, h2, l2, couplings)
            else:
                print("This is not a correct response, please try again.")
                continue
            
            # asks if another loop should be performed, breaks if not
            if input("Would you like to add another qubit? [Y/N]") in ["Y", "y"]:
                continue
            else:
                break
        
        # Nicely Prints the resultant Hamiltonian expression
        print("The Hamiltonian for the given system is:")
        print(system.hamiltonian_string())
        input("Press Enter key to exit...")


            
