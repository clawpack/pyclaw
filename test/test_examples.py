import util                      

if __name__=="__main__":
   nose.main()

def test_examples():
   path        = './test/acoustics/1d/homogeneous'
   name        = 'acoustics'
   target_name = ''
   verifier    = lambda error: abs(error-0.00104856594174)<1.e-5
   yield(util.build_run_verify, path, target_name, name, name, verifier)

   path        = './test/acoustics/2d/homogeneous'
   module_name = 'acoustics'
   target_name = ''
   method_name = 'acoustics2D'
   def verify_acoustics2D(test_x):
      import numpy
      verify_x=numpy.loadtxt('test/acoustics2D_solution')
      return (numpy.linalg.norm(test_x-verify_x)<1.e-12)
   yield(util.build_run_verify, path, target_name, module_name,
         method_name, verify_acoustics2D)
   
