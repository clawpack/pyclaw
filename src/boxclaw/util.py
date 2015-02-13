"""
BoxClaw related utilities.
"""

def boxlib_build_variant_arg_dicts(kernel_languages=('Fortran',)):
    import itertools

    # test fboxlib only if it is available
    try:
        import fboxlib
    except ImportError:
        return []

    opt_names = 'use_boxlib','kernel_language'
    opt_product = itertools.product((True,),kernel_languages)
    arg_dicts = [dict(zip(opt_names,argset)) for argset in opt_product]

    return arg_dicts
