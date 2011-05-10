"""
Utility functions for parsing command line arguments and passing them to Python functions.
"""
import json
import sys

def _method_info_from_argv(argv=None):
    """Command-line -> method call arg processing.
    
    - positional args:
            a b -> method('a', 'b')
    - intifying args:
            a 123 -> method('a', 123)
    - json loading args:
            a '["pi", 3.14, null]' -> method('a', ['pi', 3.14, None])
    - keyword args:
            a foo=bar -> method('a', foo='bar')
    - using more of the above
            1234 'extras=["r2"]'  -> method(1234, extras=["r2"])
    
    @param argv {list} Command line arg list. Defaults to `sys.argv`.
    @returns (<method-name>, <args>, <kwargs>)
    """
    if argv is None:
        argv = sys.argv

    method_name, arg_strs = argv[1], argv[2:]
    args = []
    kwargs = {}
    for s in arg_strs:
        if s.count('=') == 1:
            key, value = s.split('=', 1)
        else:
            key, value = None, s
        try:
            value = json.loads(value) 
        except ValueError:
            pass
        if value=='True': value=True
        if value=='False': value=False
        if key:
            kwargs[key] = value
        else:
            args.append(value)
    return method_name, args, kwargs

def _info_from_argv(argv=None):
    """Command-line -> method call arg processing.
    
    - positional args:
            a b -> method('a', 'b')
    - intifying args:
            a 123 -> method('a', 123)
    - json loading args:
            a '["pi", 3.14, null]' -> method('a', ['pi', 3.14, None])
    - keyword args:
            a foo=bar -> method('a', foo='bar')
    - using more of the above
            1234 'extras=["r2"]'  -> method(1234, extras=["r2"])
    
    @param argv {list} Command line arg list. Defaults to `sys.argv`.
    @returns (<method-name>, <args>, <kwargs>)
    """
    if argv is None:
        argv = sys.argv

    arg_strs = argv[1:]
    args = []
    kwargs = {}
    for s in arg_strs:
        if s.count('=') == 1:
            key, value = s.split('=', 1)
        else:
            key, value = None, s
        try:
            value = json.loads(value) 
        except ValueError:
            pass
        if value=='True': value=True
        if value=='False': value=False
        if key:
            kwargs[key] = value
        else:
            args.append(value)
    return args, kwargs

