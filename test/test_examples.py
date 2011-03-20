import util                      

if __name__=="__main__":
   nose.main()

def test_examples():
   path        = './apps/acoustics/1d/homogeneous'
   name        = 'acoustics'
   target_name = ''
   verifier    = lambda error: abs(error-0.0044043)<1.e-5
   yield util.build_run_verify, path, target_name, name, name, verifier

