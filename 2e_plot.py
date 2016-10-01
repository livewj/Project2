#plots the results produced in "2e.cpp"

import numpy as np
import matplotlib.pyplot as plt

#function that reads the data produced in "2e.cpp"

def read(filename):
	infile = open(filename, 'r')
	rho = []
	u1 = []
	u2 = []
	u3 = []
	relevant_lines = infile.readlines()[2:] #skips the irrelevant lines
	for line in relevant_lines:
		data = line.split()
		rho.append(float(data[0]))
		u1.append(float(data[1]))
		u2.append(float(data[2]))
		u3.append(float(data[3]))
	infile.close()
	rho = np.array(rho)
	u1 = np.array(u1)
	u2 = np.array(u2)
	u3 = np.array(u3)
	return rho, u1, u2, u3

#extracting data
rho, u1, u2, u3 = read('wavefunc_omega_5_ON_not.txt')


#plotting
plt.plot(rho, u1, rho, u2, rho, u3)
plt.xlabel('$\\rho$')
plt.ylabel('$|u(\\rho)|^2$')
plt.legend(['$|u_1(\\rho)|^2$', '$|u_2(\\rho)|^2$', '$|u_3(\\rho)|^2$'])
plt.title('Radial wavefunction for omega = 5 with Coulomb interaction')
#plt.savefig('2e_5_ON.png')
plt.show()
