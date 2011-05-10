import sys

def build_run_verify(path, target_name, module_name, method_name, verifier, options):
    sys.path.append(path)
    try:    
        build_command = "make -B -C %s %s" % (path, target_name)
        build_info = "[path: %s, name: %s]" % (path, target_name)
        build_it(build_info, build_command)
        module = __import__(module_name)
        run_method = getattr(module,method_name)
        verify_it(run_method, build_info, verifier, options)
    finally:
        sys.path.remove(path)
        if module_name in sys.modules:
            del(sys.modules[module_name])
            del module

def run_verify(path, module_name, method_name, verifier, options):
    """For tests that need no build (pure Python)"""
    sys.path.append(path)
    try:    
        module = __import__(module_name)
        run_method = getattr(module,method_name)
        verify_it(run_method, method_name, verifier, options)
    finally:
        sys.path.remove(path)
        if module_name in sys.modules:
            del(sys.modules[module_name])
            del module


def build_it(build_info, build_command):
    import subprocess, sys
    p = subprocess.Popen(build_command, executable='/bin/bash',shell=True, stdout=subprocess.PIPE ,stderr=subprocess.STDOUT)
    (stdout_data, ignore) = p.communicate()
    rc = p.returncode
    
    if rc is not None and rc % 256:
        err = "Error building %s\n" % (build_info)
        sys.stderr.write(err)
        sys.stderr.write(build_command)
        sys.stderr.write(stdout_data)
        raise BuildError(err)
 
def verify_it(run_method, method_name, verifier, options):
    output = run_method(**options)

    if not verifier(output):
        err = "Error verifying %s\n" % method_name
        sys.stderr.write("output: %s\n" % output)
        raise VerifyError(err)
    return

class Error(Exception):
    def __init__(self, expr):
        self.expr = expr
    def __str__(self):
        return repr(self.expr)

class BuildError(Exception):
    pass

class VerifyError(Exception):
    pass
