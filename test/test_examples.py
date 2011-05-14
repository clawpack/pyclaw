from nose.plugins.attrib import attr
import util                      

if __name__=="__main__":
    import nose
    nose.main()

@attr(examples='implicit')
def test_implicit():
    path           = './test/acoustics/1d/homogeneous'
    target_name    = 'all'
    module_name    = 'acoustics_implicit'
    method_name    = 'acoustics'
    method_options = {'sclaw': 1, 'petscts': 1}
    verifier       = lambda error: abs(error-0.00220809553637)<1.e-5
    yield(util.build_run_verify, path, target_name, module_name, method_name, verifier, method_options)

def test_examples():
    path        = './test/acoustics/1d/homogeneous'
    name        = 'acoustics'
    target_name = 'step1.so'
    method_options = {'kernel_language' : 'Fortran', 'use_PETSc' : False, 'soltype' : 'classic'}
    verifier    = lambda error: abs(error-0.00104856594174)<1.e-5
    yield(util.build_run_verify, path, target_name, name, name, verifier, method_options)

    method_options = {'kernel_language' : 'Python', 'use_PETSc' : False, 'soltype' : 'classic'}
    yield(util.run_verify, path, name, name, verifier, method_options)

    method_options = {'kernel_language' : 'Fortran', 'use_PETSc' : True, 'soltype' : 'classic'}
    yield(util.build_run_verify, path, target_name, name, name, verifier, method_options)

    method_options = {'kernel_language' : 'Python', 'use_PETSc' : True, 'soltype' : 'classic'}
    yield(util.run_verify, path, name, name, verifier, method_options)

    verifier    = lambda error: abs(error-0.000818286913339)<1.e-5
    target_name = 'flux1.so'
    method_options = {'kernel_language' : 'Python', 'use_PETSc' : False, 'soltype' : 'sharpclaw'}
    yield(util.run_verify, path, name, name, verifier, method_options)

    method_options = {'kernel_language' : 'Fortran', 'use_PETSc' : False, 'soltype' : 'sharpclaw'}
    yield(util.build_run_verify, path, target_name, name, name, verifier, method_options)

    method_options = {'kernel_language' : 'Python', 'use_PETSc' : True, 'soltype' : 'sharpclaw'}
    yield(util.run_verify, path, name, name, verifier, method_options)

    method_options = {'kernel_language' : 'Fortran', 'use_PETSc' : True, 'soltype' : 'sharpclaw'}
    yield(util.build_run_verify, path, target_name, name, name, verifier, method_options)

    
    path        = './test/acoustics/2d/homogeneous'
    module_name = 'acoustics'
    target_name = ''
    method_name = 'acoustics2D'

    def verify_acoustics2D(test_x):
        import numpy
        verify_x=numpy.loadtxt('test/acoustics2D_solution')
        return (numpy.linalg.norm(test_x-verify_x)<1.e-12)

    yield(util.build_run_verify, path, target_name, module_name,
           method_name, verify_acoustics2D, {})

    method_options = {'use_PETSc' : False}
    yield(util.build_run_verify, path, target_name, module_name,
           method_name, verify_acoustics2D, {})
