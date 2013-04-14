def configuration(parent_package='',top_path=None):
    import os

    from numpy.distutils.misc_util import Configuration

    config = Configuration('next', parent_package, top_path)

    # these are currently work-in-progress
    # iso_src_dir = os.path.join(os.path.dirname(__file__),'iso_c')

    # config.add_extension('advection_iso_c',
    #                      [os.path.join(iso_src_dir, 'rp1_advection.f90')])

    # config.add_extension('classic1_iso_c',
    #                      [os.path.join(iso_src_dir, 'limiter.f90'),
    #                       os.path.join(iso_src_dir, 'philim.f90'),
    #                       os.path.join(iso_src_dir, 'step1.f90')])

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration().todict())
