###############################################################################################
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
################################################################################################


from nose.plugins.attrib import attr
import util                      

if __name__=="__main__":
    import nose
    nose.main()


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
    method_options = {'sclaw': 1, 'petscts': 1, '-ts_type': theta}
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
    target_name    = 'step1.so'    
    method_options = {'kernel_language' : 'Fortran', 'use_PETSc' : False, 'soltype' : 'classic'}
    verifier       = lambda error: abs(error-0.00104856594174)<1.e-5
    yield(util.build_run_verify, path, target_name, name, name, verifier, method_options)
    

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
    target_name    = 'all'
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
    target_name    = 'flux1.so'
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
    target_name    = 'flux1.so'
    module_name    = 'acoustics'
    problem_name   = 'acoustics'
    method_options = {'kernel_language' : 'Fortran', 'use_PETSc' : True, 'soltype' : 'sharpclaw'}
    verifier       = lambda error: abs(error-0.000818286913339)<1.e-5
    yield(util.build_run_verify, path, target_name, module_name, module_name, verifier, method_options)    



