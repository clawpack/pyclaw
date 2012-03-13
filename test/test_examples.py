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

from nose.plugins.skip import Skip, SkipTest

if __name__=="__main__":
    import nose
    nose.main()


#=======================================================
#    1D tests
#=======================================================

# Regression test: 1D acoustics in homogeneous material
#@attr(testType ='regression')
@attr(solver_type='sharpclaw')
@attr(kernel_language='fortran')
@attr(petsc=True)
@attr(time_stepping_mode='implicit')
@attr(time_stepping_method='theta')
@attr(speed='fast')
def test_1D_acoustics_1a():
    r"""This test is skipped for now because we're no longer supporting
    the implicit approach used by Jed Brown in implementing it."""
    raise SkipTest
    module_path           = './test/acoustics/1d/homogeneous'
    build_target    = 'all'
    module_name    = 'acoustics_implicit'
    problem_name    = 'acoustics'
    method_options = {'sclaw': 1, 'petscts': 1, '-ts_type': 'theta'}
    verifier       = lambda error: abs(error-0.00220809553637)<1.e-5
    yield(util.build_run_verify, build_path, build_target, module_path, module_name, problem_name, verifier, method_options)
   

# Regression test: 1D acoustics in homogeneous material
#@attr(testType ='regression')
@attr(solver_type='classic')
@attr(kernel_language='fortran')
@attr(petsc=False)
@attr(time_stepping_mode='explicit')
@attr(speed='fast')
def test_1D_acoustics_1b():
    module_path           = './test/acoustics/1d/homogeneous'
    module_name    = 'acoustics'
    problem_name   = 'acoustics'
    build_path           = './src/pyclaw/clawpack'
    build_target    = 'classic1.so'
    method_options = {'kernel_language' : 'Fortran', 'use_petsc' : False, 'solver_type' : 'classic'}
    verifier       = lambda error: abs(error-0.00104856594174)<1.e-5
    yield(util.build_run_verify, build_path, build_target, module_path, module_name, problem_name, verifier, method_options)
    

# Regression test: 1D acoustics in homogeneous material
#@attr(testType ='regression')
@attr(solver_type='classic')
@attr(kernel_language='python')
@attr(petsc=False)
@attr(time_stepping_mode='explicit')
@attr(speed='fast')
def test_1D_acoustics_1c():
    path           = './test/acoustics/1d/homogeneous'
    module_name    = 'acoustics'
    problem_name   = 'acoustics'
    method_options = {'kernel_language' : 'Python', 'use_petsc' : False, 'solver_type' : 'classic'}
    verifier       = lambda error: abs(error-0.00104856594174)<1.e-5
    yield(util.run_verify, path, module_name, problem_name, verifier, method_options)
    

# Regression test: 1D acoustics in homogeneous material
#@attr(testType ='regression')
@attr(solver_type='classic')
@attr(kernel_language='fortran')
@attr(petsc=True)
@attr(time_stepping_mode='explicit')
@attr(speed='fast')
def test_1D_acoustics_1d():
    module_path    = './test/acoustics/1d/homogeneous'
    module_name    = 'acoustics'
    problem_name   = 'acoustics'
    build_path     = './src/pyclaw/clawpack'
    build_target   = 'classic1.so'
    method_options = {'kernel_language' : 'Fortran', 'use_petsc' : True, 'solver_type' : 'classic'}
    verifier       = lambda error: abs(error-0.00104856594174)<1.e-5
    yield(util.build_run_verify, build_path, build_target, module_path, module_name, problem_name, verifier, method_options)
    

# Regression test: 1D acoustics in homogeneous material
#@attr(testType ='regression')
@attr(solver_type='classic')
@attr(kernel_language='python')
@attr(petsc=True)
@attr(time_stepping_mode='explicit')
@attr(speed='fast')
def test_1D_acoustics_1e():
    path           = './test/acoustics/1d/homogeneous'
    module_name    = 'acoustics'
    problem_name   = 'acoustics'
    method_options = {'kernel_language' : 'Python', 'use_petsc' : True, 'solver_type' : 'classic'}
    verifier       = lambda error: abs(error-0.00104856594174)<1.e-5
    yield(util.run_verify, path, module_name, problem_name, verifier, method_options)

# Regression test: 1D acoustics in homogeneous material
#@attr(testType ='regression')
@attr(solver_type='sharpclaw')
@attr(kernel_language='python')
@attr(petsc=False)
@attr(time_stepping_mode='explicit')
@attr(time_stepping_method='SSP33')
@attr(speed='fast')
def test_1D_acoustics_1f():  
    module_path           = './test/acoustics/1d/homogeneous'
    module_name    = 'acoustics'
    problem_name   = 'acoustics'
    method_options = {'kernel_language' : 'Python', 'use_petsc' : False, 'solver_type' : 'sharpclaw'}
    verifier       = lambda error: abs(error-0.000298935748775)<1.e-5
    yield(util.run_verify, module_path, module_name, problem_name, verifier, method_options)


# Regression test: 1D acoustics in homogeneous material
#@attr(testType ='regression')
@attr(solver_type='sharpclaw')
@attr(kernel_language='fortran')
@attr(petsc=False)
@attr(time_stepping_mode='explicit')
@attr(time_stepping_method='SSP33')
@attr(speed='fast')
def test_1D_acoustics_1g():  
    module_path    = './test/acoustics/1d/homogeneous'
    module_name    = 'acoustics'
    problem_name   = 'acoustics'
    method_options = {'kernel_language' : 'Fortran', 'use_petsc' : False, 'solver_type' : 'sharpclaw'}
    build_target   = 'sharpclaw1.so'
    build_path     = './src/pyclaw/sharpclaw'
    verifier       = lambda error: abs(error-0.000298935748775)<1.e-5
    yield(util.build_run_verify, build_path, build_target, module_path, module_name, problem_name, verifier, method_options)


# Regression test: 1D acoustics in homogeneous material
# This one tests high order WENO kernels from PyWENO
#@attr(testType ='regression')
@attr(solver_type='sharpclaw')
@attr(kernel_language='fortran')
@attr(petsc=False)
@attr(time_stepping_mode='explicit')
@attr(time_stepping_method='SSP104')
@attr(speed='fast')
def test_1D_acoustics_with_weno17():  
    module_path    = './test/acoustics/1d/homogeneous'
    module_name    = 'acoustics'
    problem_name   = 'acoustics'
    build_target   = 'sharpclaw1.so'
    build_path     = './src/pyclaw/sharpclaw'
    method_options = {'kernel_language' : 'Fortran', 'use_petsc' : False, 'solver_type' : 'sharpclaw', 'weno_order' : 17}
    verifier       = lambda error: abs(error-0.000163221216565)<1.e-5
    yield(util.build_run_verify, build_path, build_target, module_path, module_name, problem_name, verifier, method_options)
    

# Regression test: 1D acoustics in homogeneous material
# This one tests high order WENO kernels from PyWENO
#@attr(testType ='regression')
@attr(solver_type='sharpclaw')
@attr(kernel_language='fortran')
@attr(petsc=True)
@attr(time_stepping_mode='explicit')
@attr(time_stepping_method='SSP104')
@attr(speed='fast')
def test_1D_acoustics_with_weno17_and_PETSc():  
    module_path    = './test/acoustics/1d/homogeneous'
    build_path     = './src/pyclaw/sharpclaw'
    build_target   = 'sharpclaw1.so'
    module_name    = 'acoustics'
    problem_name   = 'acoustics'
    method_options = {'kernel_language' : 'Fortran', 'use_petsc' : True, 'solver_type' : 'sharpclaw', 'weno_order' : 17}
    verifier       = lambda error: abs(error-0.000163221216565)<1.e-5
    yield(util.build_run_verify, build_path, build_target, module_path, module_name, problem_name, verifier, method_options)
    

# Regression test: 1D acoustics in homogeneous material
#@attr(testType ='regression')
@attr(solver_type='sharpclaw')
@attr(kernel_language='python')
@attr(petsc=True)
@attr(time_stepping_mode='explicit')
@attr(time_stepping_method='SSP33')
@attr(speed='fast')
def test_1D_acoustics_1h():  
    module_path           = './test/acoustics/1d/homogeneous'
    module_name    = 'acoustics'
    problem_name   = 'acoustics'
    method_options = {'kernel_language' : 'Python', 'use_petsc' : True, 'solver_type' : 'sharpclaw'}
    verifier       = lambda error: abs(error-0.000298935748775)<1.e-5
    yield(util.run_verify, module_path, module_name, problem_name, verifier, method_options)


# Regression test: 1D acoustics in homogeneous material
#@attr(testType ='regression')
@attr(solver_type='sharpclaw')
@attr(kernel_language='fortran')
@attr(petsc=True)
@attr(time_stepping_mode='explicit')
@attr(time_stepping_method='SSP33')
@attr(speed='fast')
def test_1D_acoustics_1i(): 
    module_path           = './test/acoustics/1d/homogeneous'
    build_target   = 'sharpclaw1.so'
    build_path     = './src/pyclaw/sharpclaw'
    module_name    = 'acoustics'
    problem_name   = 'acoustics'
    method_options = {'kernel_language' : 'Fortran', 'use_petsc' : True, 'solver_type' : 'sharpclaw'}
    verifier       = lambda error: abs(error-0.000298935748775)<1.e-5
    yield(util.build_run_verify, build_path, build_target, module_path, module_name, problem_name, verifier, method_options)



#=======================================================
#    2D tests
#=======================================================
# Regression test: 2D acoustics in homogeneous material
#@attr(testType ='regression')
@attr(solver_type='classic')
@attr(kernel_language='fortran')
@attr(petsc=True)
@attr(time_stepping_mode='explicit')
@attr(speed='fast')
def test_2D_acoustics_1a(): 
    module_path    = './test/acoustics/2d/homogeneous'
    build_target   = 'classic2.so'
    build_path     = './src/pyclaw/clawpack'
    module_name    = 'acoustics'
    problem_name   = 'acoustics2D'

    def verify_acoustics2D_classic(test_x):
        import numpy
        verify_x=numpy.loadtxt('test/acoustics2D_solution')
        diff = numpy.linalg.norm(test_x-verify_x)
        if diff>1.e-12:
            raise Exception('Difference between expected and computed solutions: %s' % diff)
        else: return True

    method_options = {'use_petsc' : True, 'solver_type' : 'classic' }
    yield(util.build_run_verify, build_path, build_target, module_path, module_name, problem_name, verify_acoustics2D_classic, method_options)


# Regression test: Parallel 2D acoustics in homogeneous material
#@attr(testType ='regression')
@attr(solver_type='classic')
@attr(kernel_language='fortran')
@attr(petsc=True)
@attr(time_stepping_mode='explicit')
@attr(speed='fast')
def test_2D_acoustics_1a_parallel():
    module_path           = './test/acoustics/2d/homogeneous'
    build_target    = 'classic2.so'
    build_path           = './src/pyclaw/clawpack'
    module_name    = 'acoustics'
    problem_name   = 'acoustics2D'

    def verify_acoustics2D_classic(sol):
        import numpy
        test_x = sol.q[0,:,:]
        verify_x=numpy.loadtxt('test/acoustics2D_solution')
        return (numpy.linalg.norm(test_x-verify_x)<1.e-5)

    method_options = {'use_petsc' : True, 'solver_type' : 'classic', 'np':6, 'num_output_times':10}
    yield(util.build_run_verify, build_path, build_target, module_path, module_name, problem_name, verify_acoustics2D_classic, method_options)


# Regression test: Parallel 2D shock bubble
#@attr(testType ='regression')
@attr(solver_type='classic')
@attr(kernel_language='fortran')
@attr(petsc=True)
@attr(time_stepping_mode='explicit')
@attr(speed='fast')
def test_2D_shockbubble_classic_parallel():
    raise SkipTest
    module_path           = './test/euler/2d'
    build_target    = 'classic2.so'
    build_path           = './src/pyclaw/clawpack'
    module_name    = 'shockbubble'
    problem_name   = 'shockbubble'

    def verify_shockbubble(test_x):
        import numpy
        verify_x=numpy.loadtxt('test/sb_density')
        return numpy.max(abs(test_x-verify_x))<1.e-14

    method_options = {'use_petsc' : True, 'solver_type' : 'classic', 'np':6, 'num_output_times':10}
    yield(util.build_run_verify, build_path, build_target, module_path, module_name, problem_name, verify_shockbubble, method_options)


# Regression test: 2D acoustics in homogeneous material
#@attr(testType ='regression')
@attr(solver_type='classic')
@attr(kernel_language='fortran')
@attr(petsc=False)
@attr(time_stepping_mode='explicit')
@attr(speed='fast')
def test_2D_acoustics_1b(): 
    module_path           = './test/acoustics/2d/homogeneous'
    build_target    = 'classic2.so'
    build_path           = './src/pyclaw/clawpack'
    module_name    = 'acoustics'
    problem_name   = 'acoustics2D'

    def verify_acoustics2D(test_x):
        import numpy
        verify_x=numpy.loadtxt('test/acoustics2D_solution')
        return (numpy.linalg.norm(test_x-verify_x)<1.e-12)

    method_options = {'use_petsc' : False, 'solver_type' : 'classic' }
    yield(util.build_run_verify, build_path, build_target, module_path, module_name, problem_name, verify_acoustics2D, method_options)


# Regression test: 2D acoustics in homogeneous material
#@attr(testType ='regression')
@attr(solver_type='sharpclaw')
@attr(kernel_language='fortran')
@attr(petsc=True)
@attr(time_stepping_mode='explicit')
@attr(time_stepping_method='SSP33')
@attr(speed='fast')
def test_2D_acoustics_1c(): 
    module_path           = './test/acoustics/2d/homogeneous'
    build_target    = 'sharpclaw2.so'
    build_path     = './src/pyclaw/sharpclaw'
    module_name    = 'acoustics'
    problem_name   = 'acoustics2D'

    def verify_acoustics2D_sharpclaw(test_x):
        import numpy
        verify_x=numpy.loadtxt('test/ac_sc_solution')
        diff = numpy.linalg.norm(test_x-verify_x)
        if diff>1.e-4:
            raise Exception("""In test_2D_acoustics_1c: 
                        Difference between expected and computed solutions: %s""" % diff)
        else: return True

    method_options = {'use_petsc' : True, 'solver_type' : 'sharpclaw' }
    yield(util.build_run_verify, build_path, build_target, module_path, module_name, problem_name, verify_acoustics2D_sharpclaw, method_options)


# Regression test: 2D acoustics in homogeneous material
#@attr(testType ='regression')
@attr(solver_type='sharpclaw')
@attr(kernel_language='fortran')
@attr(petsc=False)
@attr(time_stepping_mode='explicit')
@attr(time_stepping_method='SSP33')
@attr(speed='fast')
def test_2D_acoustics_1d(): 
    module_path           = './test/acoustics/2d/homogeneous'
    build_target    = 'sharpclaw2.so'
    build_path     = './src/pyclaw/sharpclaw'
    module_name    = 'acoustics'
    problem_name   = 'acoustics2D'

    def verify_acoustics2D_sharpclaw(test_x):
        import numpy
        verify_x=numpy.loadtxt('test/ac_sc_solution')
        diff = numpy.linalg.norm(test_x-verify_x)
        if diff>1.e-4:
            raise Exception("""In test_2D_acoustics_1c: 
                        Difference between expected and computed solutions: %s""" % diff)
        else: return True

    method_options = {'use_petsc' : False, 'solver_type' : 'sharpclaw' }
    yield(util.build_run_verify, build_path, build_target, module_path, module_name, problem_name, verify_acoustics2D_sharpclaw, method_options)
   

# Regression test: 2D Euler shock-bubble interaction
@attr(testType ='regression')
@attr(solver_type='classic')
@attr(kernel_language='fortran')
@attr(petsc=False)
@attr(time_stepping_mode='explicit')
def test_2D_shockbubble_1a(): 
    module_path           = './test/euler/2d'
    build_target    = 'classic2.so'
    build_path           = './src/pyclaw/clawpack'
    module_name    = 'shockbubble'
    problem_name   = 'shockbubble'

    def verify_shockbubble(test_x):
        import numpy
        verify_x=numpy.loadtxt('test/sb_density')
        return numpy.max(abs(test_x-verify_x))<1.e-12

    method_options = {'use_petsc' : False}
    yield(util.build_run_verify, build_path, build_target, module_path, module_name, problem_name, verify_shockbubble, method_options)


# Regression test: 2D Euler shock-bubble interaction
@attr(testType ='regression')
@attr(solver_type='classic')
@attr(kernel_language='fortran')
@attr(petsc=True)
@attr(time_stepping_mode='explicit')
def test_2D_shockbubble_petclaw_classic(): 
    module_path           = './test/euler/2d'
    build_target    = 'classic2.so'
    build_path           = './src/pyclaw/clawpack'
    module_name    = 'shockbubble'
    problem_name   = 'shockbubble'

    def verify_shockbubble(test_x):
        import numpy
        verify_x=numpy.loadtxt('test/sb_density')
        return numpy.max(abs(test_x-verify_x))<1.e-12

    method_options = {'use_petsc' : True}
    yield(util.build_run_verify, build_path, build_target, module_path, module_name, problem_name, verify_shockbubble, method_options)


# Regression test: 2D p-system
# This tests 
#       advanced output options
#       aux BCs
# and several other things not covered by the other tests.
@attr(testType ='regression')
@attr(solver_type='classic')
@attr(kernel_language='fortran')
@attr(petsc=True)
@attr(time_stepping_mode='explicit')
def test_psystem_petclaw_classic(): 
    module_path           = './test/psystem'
    build_target    = 'classic2.so'
    build_path           = './src/pyclaw/clawpack'
    module_name    = 'psystem'
    problem_name   = 'psystem2D'

    def verify_psystem(test_x):
        import numpy
        #verify_x=numpy.loadtxt(path+'/_output/F.txt')
        #return numpy.abs(verify_x[-1]-0.)<1.e-12
        return True

    method_options = {'use_petsc' : True}
    yield(util.build_run_verify, build_path, build_target, module_path, module_name, problem_name, verify_psystem, method_options)


# Regression test: 2D shallow water equations on a sphere
#@attr(testType ='regression')
@attr(testType ='regression')
@attr(solver_type='classic')
@attr(kernel_language='fortran')
@attr(petsc=False)
@attr(time_stepping_mode='explicit')
@attr(speed='slow')
@attr(which='this')
def test_2D_shallowwatersphere(): 
    module_path           = './test/shallow_sphere'
    build_target    = 'classic2.so problem.so'
    build_path           = './test/shallow_sphere'
    module_name    = 'shallow_4_Rossby_Haurwitz_wave'
    problem_name   = 'shallow_4_Rossby_Haurwitz'

    def verify_shallowwatersphere(test_x):
        import numpy
        verify_x=numpy.loadtxt('test/swsphere_height')
        diff = numpy.linalg.norm(test_x-verify_x)
        if diff>1.e-4:
            raise Exception("""test_2D_shallowwatersphere: 
                        Difference between expected and computed solutions: %s""" % diff)
        else: return True

    method_options = {}
    yield(util.build_run_verify, build_path, build_target, module_path, module_name, problem_name, verify_shallowwatersphere, method_options)

# Regression test: 3D acoustics in homogeneous material
#@attr(testType ='regression')
@attr(solver_type='classic')
@attr(kernel_language='fortran')
@attr(petsc=False)
@attr(num_dim=3)
@attr(speed='fast')
def test_3D_acoustics_homogeneous():  
    module_path           = './test/acoustics/3d/'
    build_target    = 'classic3.so'
    build_path           = './src/pyclaw/clawpack'
    module_name    = 'acoustics'
    problem_name   = 'acoustics3D'
    method_options = {'use_petsc' : False, 'test' : 'hom'}
    verifier       = lambda error: (abs(error-0.00286)<1.e-4)# and (abs(error[1]-3.2)<1.e-4))
    yield(util.build_run_verify, build_path, build_target, module_path, module_name, problem_name, verifier, method_options)

# Regression test: 3D acoustics in heterogeneous material
#@attr(testType ='regression')
@attr(solver_type='classic')
@attr(kernel_language='fortran')
@attr(petsc=False)
@attr(num_dim=3)
@attr(speed='fast')
def test_3D_acoustics_heterogeneous():
    module_path           = './test/acoustics/3d/'
    build_target    = 'classic3.so'
    build_path           = './src/pyclaw/clawpack'
    module_name    = 'acoustics'
    problem_name   = 'acoustics3D'
    method_options = {'use_petsc' : False, 'test' : 'het'}

    def verify_3D_het_acoustics(test_x):
        import numpy
        verify_x=numpy.loadtxt('test/pressure_3D.txt')
        diff = numpy.linalg.norm(test_x-verify_x)
        if diff>1.e-4:
            raise Exception("""test_3D_acoustics_heterogeneous(): 
                        Difference between expected and computed solutions: %s""" % diff)
        else: return True


    yield(util.build_run_verify, build_path, build_target, module_path, module_name, problem_name, verify_3D_het_acoustics, method_options)


