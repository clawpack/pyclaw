import subprocess 
import os
from pyclaw.solution import Solution

def acoustics2D(use_PETSc=True,kernel_language='Fortran',iplot=False,petscPlot=False,useController=True,htmlplot=False,soltype='classic', np =4):
    
    # Pathes needed to specify the output directory location
    pyclawPath = os.environ['PYCLAW'] 
    examplePath = pyclawPath + "/test/acoustics/2d/homogeneous"
    outputPath = examplePath + "/ParallelOutput"

    # clean the previous output
    run_command = "rm", "-r", outputPath
    p = subprocess.Popen(run_command, stdout=subprocess.PIPE ,stderr=subprocess.STDOUT)
    (stdout_data, ignore) = p.communicate()

    # run subprocess with the required number of processes
    run_command = "mpiexec", "-n", str(np) ,"python","-c", "import sys; sys.path.append('"+str(examplePath)+"'); import acoustics; acoustics.acoustics2D(use_PETSc=True, useController=True, outdir='"+str(outputPath)+"')"
    p = subprocess.Popen(run_command, stdout=subprocess.PIPE ,stderr=subprocess.STDOUT) 
    (stdout_data, ignore) = p.communicate()
    
    sol = Solution(10, format='petsc', path=outputPath )
    return sol.grid.q[0,:,:]

if __name__ == "__main__":
    acoustics2D()
