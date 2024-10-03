import RiemannSolver
from RiemannSolver import *
from numpy import empty
from pylab import plot,  figure, suptitle
from six.moves import range



print("Testing for different number of time steps ...")
noOfTestSamples = 10
max_mx = 101 
max_timeSteps = 101 

timesStepsValues = empty((noOfTestSamples ))
timeResultsVectorized =empty((noOfTestSamples ))
timeResultsPointwise =empty((noOfTestSamples ))
# error when max_mx is divisible by noOfTestSamples. size of arrays should be noOfTestSamples-1, use j to decide what to plot






rs = RiemannSolver(timeSteps = max_timeSteps/noOfTestSamples, mwaves = 2, mx=max_mx , meqn = 2, maux = 2)  


rs.q = random.random((rs.meqn,rs.mx))
rs.aux = random.random((rs.maux,rs.mx))

j = 0
for i in range(max_timeSteps/noOfTestSamples, max_timeSteps, max_timeSteps/noOfTestSamples):
	timesStepsValues[j] = i
	rs.timeSteps = i
	print()
	print("Iteration",j+1,"of", noOfTestSamples)
	print("Number of time steps is", i)
	
	waves1, s1  = rs.solveVectorized(timer = True)
	timeResultsVectorized[j] = rs.elapsedTime
	
	
	
	waves2, s2  = rs.solvePointwize(timer = True)
	timeResultsPointwise[j] = rs.elapsedTime
	
	
	j= j+1


figure(1)	
suptitle("Figure 1 shows the the execution time for different timeStep values\n for the vectorized solver in green and pointwize solver in blue")	

plot(timesStepsValues, timeResultsVectorized, "go-")

#figure(2)
#suptitle( "Figure 2 shows the the execution time for different timeStep values\n for the pointwize solver, mx = {0}".format(max_mx))	
plot(timesStepsValues, timeResultsPointwise, "bo-")
	
	




print("Testing for different sizes of mx ...")
noOfTestSamples = 10
max_mx = 101 
max_timeSteps = 101 

mxValues = empty((noOfTestSamples ))
timeResultsVectorized =empty((noOfTestSamples ))
timeResultsPointwise =empty((noOfTestSamples ))
# error when max_mx is divisible by noOfTestSamples. size of arrays should be noOfTestSamples-1, use j to decide what to plot






rs = RiemannSolver(timeSteps = max_timeSteps, mwaves = 2, mx=max_mx/noOfTestSamples , meqn = 2, maux = 2)  


rs.q = random.random((rs.meqn,rs.mx))
rs.aux = random.random((rs.maux,rs.mx))





 




j = 0
for i in range(max_mx/noOfTestSamples, max_mx,max_mx/noOfTestSamples):
	mxValues[j] = i
	rs.mx = i
	print()
	print("Iteration",j+1,"of", noOfTestSamples)
	print("mx is", i)
	
	rs.q = random.random((rs.meqn,rs.mx))
	rs.aux = random.random((rs.maux,rs.mx))
	
	waves1, s1  = rs.solveVectorized(timer = True)
	timeResultsVectorized[j] = rs.elapsedTime
	
	
	
	waves2, s2  = rs.solvePointwize(timer = True)
	timeResultsPointwise[j] = rs.elapsedTime
	
	
	j= j+1
	
figure(2)
suptitle("Figure 2 shows the the execution time for different mx values\n for the vectorized solver in green and pointwize solver in blue")	
plot(mxValues, timeResultsVectorized, "go-")

#figure(4)
#suptitle("Figure 4 shows the the execution time for different mx values\n for the pointwize solver, timeSteps = {0}".format(max_timeSteps))	
plot(mxValues, timeResultsPointwise, "bo-")
