# -*- coding: utf-8 -*-
"""
Created on Thu Sep 10 19:15:16 2020

@author: castromi
"""

from qiskit import QuantumCircuit, QuantumRegister,quantum_info, Aer,BasicAer,execute
import numpy as np

from scipy.optimize import minimize

import pandas as pd
from matplotlib import pyplot as plt


def quantum_simulator(
    theta,
    random_state,
    layers,
    backend,
    qr,
    process_stats):

    """
        Simulates the circuit and for a number of layers and given a set of thetas,
        and a random state vector computes epsilon.

        Paramaters:
            theta: 8x1 array of thetas
            random_state: random statevector needed to compute epsilon_data
            backend: Quantum simulator object from Qiskit
            qr: Quantum register object required by Qiskit
            process_stats: Python object where the stats from the function call will be save during the call
                            for scipy.optimize.minimize()
        Returns:
            epsilon
    """
    # Creating Quantum Circuit object
    circ=QuantumCircuit(qr)
    # bind random state to local variable
    random_state = random_state
    # Reshape theta into 8 by 2
    theta=theta.reshape((layers,8))
    # Build up the system appending the blocks to the circuits
    for i in range(0,layers):
        create_even_block(circ,theta)
        create_odd_block(circ,theta)
    # execute the circuit and get the resulting statevector
    circ_state = execute(circ, backend).result().get_statevector(circ)
    # take the difference of the wave modulus square of the wavefunction with the random state
    difference_vector=(circ_state*np.conjugate(circ_state)-random_state*np.conjugate(random_state))
    # compute epsilon by taking the sum of the square components
    epsilon=(difference_vector*np.conjugate(difference_vector)).sum()

    theta=theta.reshape((8*layers))

    #add value of epsilon to our record
    process_stats.add_data(epsilon.real,theta)


    # print descriptive headers in the first iteration
    if process_stats.func_eval==0:
        print(circ,'\n')
        print('\n\n\t','      Minimizing Circuit Parameters','\t\n')
        print('\t Might take some time if L is larger than 10')
        print("ObjectiveF Evaluations",'\t\t\t','Epsilon value')

    # print value every 10 evaluations
    if (process_stats.func_eval%10 == 0 and process_stats.func_eval !=0) :
        print('\t',process_stats.func_eval,'\t\t\t\t',epsilon.real)


    process_stats.increase_func_eval()

    # return the real part of epsilon the imaginary part is guarantee to be 0
    return epsilon.real



def create_even_block(circ,theta):

    """ This functions appends the even blocks to the quantum circuit
            it takes as a parameters the quantum circuit object (from the Qiskit
            package) and an 8x2 arrays of the corresponding theta angles

        Parameters:
                circ: Quantum Circuit object to which we will append the gates
                theta: 8by2 array of angles for the rotation gates
    """
    for i in range(0,len(circ.qubits)):
        circ.rz(theta[0][i],circ.qubits[i])

       # Append the cz gates
    for i in range(0,len(circ.qubits)-1):
        target_qubit_index=i+1

        while (target_qubit_index) < 4:
             circ.cz(i,i+1)
             target_qubit_index=target_qubit_index+1

       #create barrier for more tidy printing of the circuit

    circ.barrier()
    pass

   #  Append the rz gates to the circuit


def create_odd_block(circ,theta):
    """ This functions appends the odd blocks to the quantum circuit
                it takes as a parameters the quantum circuit object (from the Qiskit
                package) and an 8x2 arrays of the corresponding theta angles

        Parameters:
                circ: Quantum Circuit object to which we will append the gates
                theta: 8by2 array of angles for the rotation gates
    """

    for i in range(0,len(circ.qubits)):
        circ.rx(theta[0][i+4],circ.qubits[i])
    circ.barrier()





def init_theta(layers):
    """ this functions initialize random intial array of thetas
    """
    theta= np.pi*2 * np.random.rand(8*layers)
    return theta

def plot_process(optimization_data):
    """ This functions takes our optimization_stats object and plot the info it
        recorded during the minimization process
        Parameters:
            optimization_data: object of class optimization_stats(see below) where we recorded
            metrics from the optimzation process
    """

    iteration_number_array = optimization_data.iloc[:,0]
    epsilon_array=optimization_data.iloc[:,1]
    plt.title("Evolution of minimization of epsilon per iteration")
    plt.xlabel("Iteration number")
    plt.ylabel("Sum of modulus square of the difference")
    plt.xticks(np.linspace(0,len(iteration_number_array),10,dtype=int))
    plt.ylabel("Sum of modulus square of the difference")
    plt.plot(iteration_number_array,epsilon_array,'-o')
    plt.grid()
    plt.show()

class optimization_stats:
    """ Python object that we will use to keep track of the number of function evals of
        the minimization process and the value of epsilon
    """

    def __init__(self):
        self.func_eval=0
        self.data=[]
        pass

    def increase_func_eval(self):
        self.func_eval+=1
        pass

    def add_data(self,epsilon,theta):
        new_entry= tuple([self.func_eval])+tuple([epsilon.real])+ tuple(theta)
        self.data.append(new_entry)

    def create_data_frame(self,layers):
        column_names=["# of evaluations","epsilon"]
        column_names=column_names + [ 'theta_'+ str(i) for i in range(0,8*layers)]
        optimization_data = pd.DataFrame(self.data,columns=column_names)
        return optimization_data
    pass

def main():

    print("Loading Qiskit Packages...")

    # getting usser input
    display_report = input("Display the full report at the end y/n? ")
    output_file= input("Name your output file (.csv): ")

    layers=int(input("input an integer number of layers: "))

    # Initializing required parameters for call of minimize
    backend =Aer.get_backend('statevector_simulator')
    random_state=quantum_info.random_state(16,seed=101)
    process_stats=optimization_stats()
    angle_bounds=[(0,2*np.pi)]*layers*8
    qr=QuantumRegister(4)
    theta=init_theta(layers)

    print("Optimization Starts (method: L-BFGS-B): ")

    # call for minimize on quantum_simulator() see above for specifications
    result=minimize(quantum_simulator,theta,method='L-BFGS-B',bounds=angle_bounds,args=(random_state,layers,backend,qr,process_stats),options={'ftol':2.22e-08,'gtol':1e-8,'maxiter':100000})
    optimization_data=process_stats.create_data_frame(layers)

    # print report based on user input
    if display_report.upper() != 'Y':
        print("optimization success: ", bool(result.success))
    else:
        print('\n\n',result)
        print(result.x)
        print("epsilon :", result.fun)
        print(optimization_data)

    plot_process(optimization_data)
    optimization_data.to_csv(output_file)






main()
