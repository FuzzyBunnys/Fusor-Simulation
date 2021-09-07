import numpy
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import seaborn as sns
import csv #added this to write data for analysis in excel while troubleshooting

#ELECTRIC FIELD PORTION OF SIMULATION

#Potential is a 2d grid of ~10000 points - expanding to 1001x1001 takes a very long time to converge
Nx=301
Ny=201
potential = np.zeros((Ny,Nx))

#setup initial voltagew
Vsource= -13800
VsourceP1 = 180000
VsourceP2 = 90000
VsourceP3 = 82200


#perform a finite difference approximation
res = 0.1

#set the resolution in x and y coordinates
hx= res
hy= res

#precalculate kx and ky
kx= (hx**2)/(2*(hx**2+hy**2))
ky= (hy**2)/(2*(hx**2+hy**2))

#create our voltage matrix
voltmat = np.zeros((Ny,Nx))
#specify the voltage sources
#negative - electron reflector
voltmat[20, 50] = Vsource
voltmat[21, 50] = Vsource
voltmat[60, 50] = Vsource
voltmat[61, 50] = Vsource
voltmat[100, 50] = Vsource
voltmat[101, 50] = Vsource
voltmat[140, 50] = Vsource
voltmat[141, 50] = Vsource
voltmat[180, 50] = Vsource
voltmat[181, 50] = Vsource

#positive - first venetian blind
voltmat[23,147] = VsourceP2
voltmat[22,148] = VsourceP2
voltmat[21,149] = VsourceP2
voltmat[20,150] = VsourceP2
voltmat[19,151] = VsourceP2
voltmat[18,152] = VsourceP2
voltmat[17,153] = VsourceP2

voltmat[63,147] = VsourceP2
voltmat[62,148] = VsourceP2
voltmat[61,149] = VsourceP2
voltmat[60,150] = VsourceP2
voltmat[59,151] = VsourceP2
voltmat[58,152] = VsourceP2
voltmat[57,153] = VsourceP2

voltmat[103,147] = VsourceP2
voltmat[102,148] = VsourceP2
voltmat[101,149] = VsourceP2
voltmat[100,150] = VsourceP2
voltmat[99,151] = VsourceP2
voltmat[98,152] = VsourceP2
voltmat[97,153] = VsourceP2

voltmat[143,147] = VsourceP2
voltmat[142,148] = VsourceP2
voltmat[141,149] = VsourceP2
voltmat[140,150] = VsourceP2
voltmat[139,151] = VsourceP2
voltmat[138,152] = VsourceP2
voltmat[137,153] = VsourceP2

voltmat[183,147] = VsourceP2
voltmat[182,148] = VsourceP2
voltmat[181,149] = VsourceP2
voltmat[180,150] = VsourceP2
voltmat[179,151] = VsourceP2
voltmat[178,152] = VsourceP2
voltmat[177,153] = VsourceP2
#positive - second positive grid

voltmat[22,175] = VsourceP3
voltmat[62,175] = VsourceP3
voltmat[102,175] = VsourceP3
voltmat[142,175] = VsourceP3
voltmat[182,175] = VsourceP3

#positive - second venetian blind
voltmat[17,247] = VsourceP1
voltmat[18,248] = VsourceP1
voltmat[19,249] = VsourceP1
voltmat[20,250] = VsourceP1
voltmat[21,251] = VsourceP1
voltmat[22,252] = VsourceP1
voltmat[23,253] = VsourceP1

voltmat[57,247] = VsourceP1
voltmat[58,248] = VsourceP1
voltmat[59,249] = VsourceP1
voltmat[60,250] = VsourceP1
voltmat[61,251] = VsourceP1
voltmat[62,252] = VsourceP1
voltmat[63,253] = VsourceP1

voltmat[97,247] = VsourceP1
voltmat[98,248] = VsourceP1
voltmat[99,249] = VsourceP1
voltmat[100,250] = VsourceP1
voltmat[101,251] = VsourceP1
voltmat[102,252] = VsourceP1
voltmat[103,253] = VsourceP1

voltmat[137,247] = VsourceP1
voltmat[138,248] = VsourceP1
voltmat[139,249] = VsourceP1
voltmat[140,250] = VsourceP1
voltmat[140,251] = VsourceP1
voltmat[141,252] = VsourceP1
voltmat[142,253] = VsourceP1

voltmat[177,247] = VsourceP1
voltmat[178,248] = VsourceP1
voltmat[179,249] = VsourceP1
voltmat[180,250] = VsourceP1
voltmat[181,251] = VsourceP1
voltmat[182,252] = VsourceP1
voltmat[183,253] = VsourceP1


#correctionfactor
overcorrection_fact = 1.94 #picking this because I was told too start with it

#initial conditions for our while loop
tol = 0.1 #set to low resolution atm. want it to be 0.001 for actual
dsum=1
iter=0

while dsum > tol:

	#grab the squared value of our last try
	square1 = np.square(potential)
	sum1 = square1.sum()
	#iterate through rows
	for ny in range(0,200):
		#iterate through columns
		for nx in range(0,300):
			if nx != 0 and nx != 300 and ny != 0 and ny != 200: #leave the edges alone! - these are our dirichlet boundary conditions
			#calculate potential
				oldpotential = potential[ny,nx] # so we grab the old value
				#calculate the new value
				newpotential = ky * (potential[ny, nx + 1] + potential[ny, nx - 1]) + kx * (potential[ny + 1, nx] + potential[ny - 1, nx]) + kx*voltmat[ny,nx]
				delta = newpotential - oldpotential
				potential[ny,nx] += delta*overcorrection_fact
				#end
			#end
		#end
	#end

	#sum our potential grid again
	square2=np.square(potential)
	sum2=square2.sum()
	#print(sum2)

	#check the absolute value difference between our before and after grid
	dsum = sum2-sum1
	dsum = abs(dsum)
	print(dsum)
	with open('programtimingGS.csv', mode = 'a') as timingresults:
			timing_writer = csv.writer(timingresults, delimiter = ',')
			timing_writer.writerow([dsum,iter])#for performance investigations
	#count iterations
	#count iterations
	iter = iter+1

print("iterations for electrical field=",iter)
# lets save this data
np.savetxt("potential.csv",potential,delimiter=",")
# then we can monkey with it in a plot later on
plt.figure(1)

sns.heatmap(graphdata)
plt.show()


# Improvements  
# - timing macro/function of some kind - CHECK
# - boundary conditions - CHECK
# - plot print in sublime - CHECK
# - Gauss Seidel optimization - CHECK
# - Overcorrection (successive over-relaxation) - CHECK
# - Poisson solution - CHECK
# - Simulate final configuration for venetian blind - CHECK
# - For visualisation values have been shifted by -13800

# make it appear like the 2d slice of the direct energy conversion doohickey 