##################################################################################################
# This python script contains regression and unit tests for PYCLAW and PETCLAW
# These test are executed with nose (http://somethingaboutorange.com/mrl/projects/nose/1.0.0/)
# Attributes are use to tag each regression and unit test. This way nose can run different type 
# of subclasses regression tests.
#
#
# Use: nosetests -a attribute-1 = string-1,attribute-2 = string-2,attribute-3 = string-3
#      for logic AND
#
# Use: nosetests -a attribute-1 = string-1 -a attribute-2 = string-2 -a attribute-3 = string-3
#      for logic OR
#
##################################################################################################


from nose.plugins.attrib import attr
import util                      

if __name__=="__main__":
    import nose
    nose.main()


#=======================================================
#    1D tests
#=======================================================

# Regression test: 1D acoustic in homegeneous material
#@attr(testType ='regression')
@attr(sd='sharpclaw')
@attr(kernel_language='fortran')
@attr(petsc=True)
@attr(ts='implicit')
@attr(speed='fast')
def test_1D_acoustic_homogeneous_1a():
    path           = './test/acoustics/1d/homogeneous'
    target_name    = 'all'
    module_name    = 'acoustics_implicit'
    problem_name    = 'acoustics'
    method_options = {'sclaw': 1, 'petscts': 1, '-ts_type': 'theta'}
    verifier       = lambda error: abs(error-0.00220809553637)<1.e-5
    yield(util.build_run_verify, path, target_name, module_name, problem_name, verifier, method_options)


# Regression test: 1D acoustic in homegeneous material
#@attr(testType ='regression')
@attr(sd='classic')
@attr(kernel_language='python')
@attr(petsc=False)
@attr(ts='explicit')
@attr(speed='fast')
def test_1D_acoustic_homogeneous_1b():
    path           = './test/acoustics/1d/homogeneous'
    module_name    = 'acoustics'
    problem_name   = 'acoustics'
    target_name    = 'classic1.so'
    method_options = {'kernel_language' : 'Fortran', 'use_PETSc' : False, 'soltype' : 'classic'}
    verifier       = lambda error: abs(error-0.00104856594174)<1.e-5
    yield(util.build_run_verify, path, target_name, module_name, problem_name, verifier, method_options)
    

# Regression test: 1D acoustic in homegeneous material
#@attr(testType ='regression')
@attr(sd='classic')
@attr(kernel_language='python')
@attr(petsc=False)
@attr(ts='explicit')
@attr(speed='fast')
def test_1D_acoustic_homogeneous_1c():
    path           = './test/acoustics/1d/homogeneous'
    module_name    = 'acoustics'
    problem_name   = 'acoustics'
    method_options = {'kernel_language' : 'Python', 'use_PETSc' : False, 'soltype' : 'classic'}
    verifier       = lambda error: abs(error-0.00104856594174)<1.e-5
    yield(util.run_verify, path, module_name, problem_name, verifier, method_options)
    

# Regression test: 1D acoustic in homegeneous material
#@attr(testType ='regression')
@attr(sd='classic')
@attr(kernel_language='fortran')
@attr(petsc=True)
@attr(ts='explicit')
@attr(speed='fast')
def test_1D_acoustic_homogeneous_1d():
    path           = './test/acoustics/1d/homogeneous'
    target_name    = 'classic1.so'
    module_name    = 'acoustics'
    problem_name    = 'acoustics'
    method_options = {'kernel_language' : 'Fortran', 'use_PETSc' : True, 'soltype' : 'classic'}
    verifier       = lambda error: abs(error-0.00104856594174)<1.e-5
    yield(util.build_run_verify, path, target_name, module_name, problem_name, verifier, method_options)
    

# Regression test: 1D acoustic in homegeneous material
#@attr(testType ='regression')
@attr(sd='classic')
@attr(kernel_language='python')
@attr(petsc=True)
@attr(ts='explicit')
@attr(speed='fast')
def test_1D_acoustic_homogeneous_1e():
    path           = './test/acoustics/1d/homogeneous'
    module_name    = 'acoustics'
    problem_name   = 'acoustics'
    method_options = {'kernel_language' : 'Python', 'use_PETSc' : True, 'soltype' : 'classic'}
    verifier       = lambda error: abs(error-0.00104856594174)<1.e-5
    yield(util.run_verify, path, module_name, problem_name, verifier, method_options)

# Regression test: 1D acoustic in homegeneous material
#@attr(testType ='regression')
@attr(sd='sharpclaw')
@attr(kernel_language='python')
@attr(petsc=False)
@attr(ts='explicit')
@attr(speed='fast')
def test_1D_acoustic_homogeneous_1f():  
    path           = './test/acoustics/1d/homogeneous'
    module_name    = 'acoustics'
    problem_name   = 'acoustics'
    target_name = 'sharpclaw1.so'
    method_options = {'kernel_language' : 'Python', 'use_PETSc' : False, 'soltype' : 'sharpclaw'}
    verifier       = lambda error: abs(error-0.000818286913339)<1.e-5
    yield(util.run_verify, path, module_name, problem_name, verifier, method_options)


# Regression test: 1D acoustic in homegeneous material
#@attr(testType ='regression')
@attr(sd='sharpclaw')
@attr(kernel_language='fortran')
@attr(petsc=False)
@attr(ts='explicit')
@attr(speed='fast')
def test_1D_acoustic_homogeneous_1g():  
    path           = './test/acoustics/1d/homogeneous'
    target_name    = 'sharpclaw1.so'
    module_name    = 'acoustics'
    problem_name   = 'acoustics'
    method_options = {'kernel_language' : 'Fortran', 'use_PETSc' : False, 'soltype' : 'sharpclaw'}
    verifier       = lambda error: abs(error-0.000818286913339)<1.e-5
    yield(util.build_run_verify, path, target_name, module_name, problem_name, verifier, method_options)
    

# Regression test: 1D acoustic in homegeneous material
#@attr(testType ='regression')
@attr(sd='sharpclaw')
@attr(kernel_language='python')
@attr(petsc=True)
@attr(ts='explicit')
@attr(speed='fast')
def test_1D_acoustic_homogeneous_1h():  
    path           = './test/acoustics/1d/homogeneous'
    module_name    = 'acoustics'
    problem_name   = 'acoustics'
    method_options = {'kernel_language' : 'Python', 'use_PETSc' : True, 'soltype' : 'sharpclaw'}
    verifier       = lambda error: abs(error-0.000818286913339)<1.e-5
    yield(util.run_verify, path, module_name, problem_name, verifier, method_options)


# Regression test: 1D acoustic in homegeneous material
#@attr(testType ='regression')
@attr(sd='sharpclaw')
@attr(kernel_language='fortran')
@attr(petsc=True)
@attr(ts='explicit')
@attr(speed='fast')
def test_1D_acoustic_homogeneous_1i(): 
    path           = './test/acoustics/1d/homogeneous'
    target_name    = 'sharpclaw1.so'
    module_name    = 'acoustics'
    problem_name   = 'acoustics'
    method_options = {'kernel_language' : 'Fortran', 'use_PETSc' : True, 'soltype' : 'sharpclaw'}
    verifier       = lambda error: abs(error-0.000818286913339)<1.e-5
    yield(util.build_run_verify, path, target_name, module_name, module_name, verifier, method_options)    





#=======================================================
#    2D tests
#=======================================================
# Regression test: 2D acoustic in homegeneous material
#@attr(testType ='regression')
#@attr(sd='sharpclaw')
@attr(kernel_language='fortran')
@attr(petsc=True)
@attr(ts='explicit')
@attr(speed='fast')
def test_2D_acoustic_homogeneous_1a(test_x): 
    path           = './test/acoustics/2d/homogeneous'
    target_name = 'classic2.so'
    module_name    = 'acoustics'
    problem_name   = 'acoustics2D'

    def verify_acoustics2D_classic(test_x):
        import numpy
        verify_x=numpy.loadtxt('test/acoustics2D_solution')
        return (numpy.linalg.norm(test_x-verify_x)<1.e-14)

    method_options = {'use_PETSc' : True, 'soltype' : 'classic' }
    yield(util.build_run_verify, path, target_name, module_name, method_name, verify_acoustics2D_classic, method_options)

    method_options = {'use_PETSc' : False, 'soltype' : 'classic' }
    yield(util.build_run_verify, path, target_name, module_name, method_name, verify_acoustics2D_classic, method_options)    


    target_name = 'sharpclaw2.so'

    def verify_acoustics2D_sharpclaw(test_x):
        import numpy
        verify_x=numpy.loadtxt('test/ac_sc_solution')
        print numpy.linalg.norm(test_x-verify_x)
        return (numpy.linalg.norm(test_x-verify_x)<1.e-10)


    method_options = {'use_PETSc' : True, 'soltype' : 'sharpclaw' }
    yield(util.build_run_verify, path, target_name, module_name, method_name, verify_acoustics2D_sharpclaw, method_options)

    method_options = {'use_PETSc' : False, 'soltype' : 'sharpclaw' }
    yield(util.build_run_verify, path, target_name, module_name, method_name, verify_acoustics2D_sharpclaw, method_options)    
   


# Regression test: 2D acoustic in homegeneous material
#@attr(testType ='regression')
#@attr(sd='sharpclaw')
#@attr(kernel_language='fortran')
#@attr(petsc=False)
#@attr(ts='explicit')
#@attr(speed='fast')
def test_2D_acoustic_homogeneous_1b(test_x): 
    path           = './test/acoustics/2d/homogeneous'
    target_name    = 'classic2.so'
    module_name    = 'acoustics'
    problem_name   = 'acoustics2D'

    def verify_acoustics2D(test_x):
        import numpy
        verify_x=numpy.loadtxt('test/acoustics2D_solution')
        return (numpy.linalg.norm(test_x-verify_x)<1.e-14)

    method_options = {'use_PETSc' : False}
    yield(util.build_run_verify, path, target_name, module_name, problem_name, verify_acoustics2D, method_options)   

# Regression test: 2D acoustic in homegeneous material
#@attr(testType ='regression')
#@attr(sd='sharpclaw')
#@attr(kernel_language='fortran')
#@attr(petsc=False)
#@attr(ts='explicit')
#@attr(speed='fast')
def test_2D_shockbubble_1a(test_x): 
    path           = './test/euler/2d'
    target_name = 'classic2.so'
    module_name    = 'shockbubble'
    problem_name   = 'shockbubble'

    def verify_shockbubble(test_x):
        import numpy
        verify_x=numpy.loadtxt('test/sb_density')
        print numpy.max(abs(test_x-verify_x))
        return numpy.max(abs(test_x-verify_x))<1.e-14

    method_options = {'use_PETSc' : False}
    yield(util.build_run_verify, path, target_name, module_name, problem_name, verify_shockbubble, method_options)


# Regression test: 2D acoustic in homegeneous material
#@attr(testType ='regression')
#@attr(sd='sharpclaw')
#@attr(kernel_language='fortran')
#@attr(petsc=False)
#@attr(ts='explicit')
#@attr(speed='fast')
def test_2D_shockbubble_1b(test_x): 
    path           = './test/euler/2d'
    target_name    = 'classic2.so'
    module_name    = 'shockbubble'
    problem_name   = 'shockbubble'

    def verify_shockbubble(test_x):
        import numpy
        verify_x=numpy.loadtxt('test/sb_density')
        print numpy.max(abs(test_x-verify_x))
        return numpy.max(abs(test_x-verify_x))<1.e-14

    method_options = {'use_PETSc' : True}
    yield(util.build_run_verify, path, target_name, module_name, problem_name, verify_shockbubble, method_options)

