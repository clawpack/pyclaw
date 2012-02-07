"""
Module for making galleries of thumbnails allowing easy browsing of
applications directories.

These tools assume that the examples have already been run and the plots
produced using "make .plots".  

You should use the script python/run_examples.py to do this first.
"""

import os


# Main root for html links:
#claw_html_root='http://localhost:50005'     
claw_html_root='http://numerics.kaust.edu.sa/pyclaw'


# Determine Clawpack directory:
clawdir_default = os.environ.get('PYCLAW',None)
if clawdir_default is None:
    print "*** Error: set environment variable PYCLAW"

# Location for gallery files:
gallery_dir_default = os.path.join(clawdir_default,'doc/gallery')  

remake = False   # True ==> remake all thumbnails even if they exist.

class GalleryItem(object):
    
    def __init__(self, appdir, plotdir, description, images):
        self.appdir = appdir
        self.plotdir = plotdir
        self.description = description
        self.images = images

class GallerySection(object):
    
    def __init__(self,title,description=""):
        self.title = title
        self.description = description
        self.items = []

    def new_item(self, appdir, plotdir, description, thumbs):
        gitem = GalleryItem(appdir, plotdir, description, thumbs)
        self.items.append(gitem)
        return gitem


class Gallery(object):
    
    def __init__(self, title, clawdir=clawdir_default):
        self.title = title
        self.clawdir = clawdir
        self.claw_html_root = claw_html_root
        self.sections = []

    def new_section(self, title, description=""):
        gsec = GallerySection(title, description)
        gsec.items = []
        self.sections.append(gsec)
        return gsec
    
    def create(self,fname,gallery_dir=None):

        # Directory for gallery files:
        if gallery_dir is None:
            gallery_dir = gallery_dir_default

        print "Gallery files will be created in directory "
        print "   ", gallery_dir

        try:
            if not os.path.isdir(gallery_dir):
                os.system('mkdir %s' % gallery_dir)
                print "Created directory ",gallery_dir
            start_dir = os.getcwd()
            os.chdir(gallery_dir)
        except:
            print "*** Error moving to directory ",gallery_dir
            print "*** Gallery not created"
            raise

        gfile = open(fname, 'w')
        gfile.write(".. _gallery:\n\n")
        gfile.write("==========================\n")
        gfile.write("Application gallery\n")
        gfile.write("==========================\n\n")

        gfile.write("%s\n" % self.title)
        gfile.write("="*len(self.title)+"\n\n")
        for gsec in self.sections:
            gfile.write("%s\n" % gsec.title)
            gfile.write("="*len(gsec.title)+"\n")
            gfile.write("%s\n" % gsec.description)

            for gitem in gsec.items:
                gfile.write('\n')
                desc = gitem.description.lstrip().replace('\n',' ')
                gfile.write('%s\n\n' % desc)
                code = os.path.join(claw_html_root, gitem.appdir)
                plotindex = os.path.join(claw_html_root, gitem.appdir, \
                               gitem.plotdir, '_PlotIndex.html')
                gfile.write('$PYCLAW/%s ... \n`README <%sREADME.rst>`__ ... \n`Plots <%s>`__\n' \
                      % (gitem.appdir,code,plotindex))

                if not os.path.exists('./'+gitem.appdir):
                    os.makedirs ('./'+gitem.appdir)
                os.system('cp %sREADME.rst %s' % ('$PYCLAW/'+gitem.appdir, './'+gitem.appdir))
                gfile.write('\n\n')

                for image in gitem.images:

                    src_name = os.path.join(gitem.appdir, gitem.plotdir, image)
                    thumb_name = src_name.replace('/','_')
                    src_html = os.path.join(claw_html_root,src_name) + '.html'
                    src_name = os.path.join(self.clawdir,src_name)
                    src_png = src_name + '.png'
                    if not os.path.isdir('thumbnails'):
                        print "Creating directory thumbnails"
                        os.system('mkdir thumbnails')
                    thumb_file = os.path.join('thumbnails',thumb_name + '.png')
                    if os.path.isfile(thumb_file) and (not remake):
                        print "Thumbnail exists: ",thumb_file
                    else:
                        scale = 0.3
                        make_thumb(src_png ,thumb_file, scale)
                    gfile.write('.. image:: %s\n   :width: 5cm\n' % thumb_file)
                    gfile.write('   :target: %s\n' % src_html)
                gfile.write('\n\n')

        print "Created ",fname, " in directory ", os.getcwd()
        os.chdir(start_dir)
                    

def make_thumb(src_file, thumb_file, scale):
    
    from numpy import floor 
    if not os.path.exists(src_file):
        print '*** Error in make_thumb: cannot find file ',src_file
    else:
        # convert scale to percent:
        scale = int(floor(scale*100))
        os.system('convert -resize %s' % scale + '% ' + \
            '%s %s' % (src_file, thumb_file))
        print "Converted ",src_file
        print "   to     ",thumb_file

   

def test():
    gallery = Gallery(title="Test Gallery")
    plotdir = '_plots'
    gsec = gallery.new_section('1-dimensional advection')

    appdir = 'apps/advection/1d'
    description = """
        Advecting Gaussian with periodic boundary."""
    images = ('frame0000fig1', 'frame0004fig1')
    gsec.new_item(appdir, plotdir, description, images)
       
    gallery.create('test.rst')



def make_1d():
    gallery = Gallery(title="Gallery of 1d PyClaw applications")
    plotdir = '_plots'


    #----------------------------------------------
    gsec = gallery.new_section('1-dimensional advection')
    appdir = 'apps/advection_1d/'
    description = """
         Advecting Gaussian with periodic boundary."""
    images = ('frame0000fig1', 'frame0004fig1', 'frame0010fig1')
    gsec.new_item(appdir, plotdir, description, images)
    #----------------------------------------------
    gsec = gallery.new_section('1-dimensional variable-velocity advection')
    appdir = 'apps/advection_1d_variable'
    description = """
         Advecting Gaussian and square wave with periodic boundary."""
    images = ('frame0000fig1', 'frame0004fig1', 'frame0008fig1')
    gsec.new_item(appdir, plotdir, description, images)
     #----------------------------------------------
    gsec = gallery.new_section('1-dimensional acoustics')
    appdir = 'apps/acoustics_1d_homogeneous'
    description = """
         Acoustics equations with wall boundary at left and extrap at
         right."""
    images = ('frame0000fig1', 'frame0005fig1', 'frame0010fig1')
    gsec.new_item(appdir, plotdir, description, images)
    #----------------------------------------------
    gsec = gallery.new_section("1-dimensional Burgers' equation")
    appdir = 'apps/burgers_1d/'
    description = """
        Burgers' equation with sinusoidal initial data, steepening to
        N-wave.  """
    images = ('frame0000fig0', 'frame0003fig0', 'frame0006fig0')
    gsec.new_item(appdir, plotdir, description, images)
    #----------------------------------------------
    gsec = gallery.new_section("1-dimensional shallow water equation")
    appdir = 'apps/shallow_1d/'
    description = """Shallow water shock tube."""
    images = ('frame0000fig0', 'frame0003fig0', 'frame0006fig0')
    gsec.new_item(appdir, plotdir, description, images)
    #----------------------------------------------
    gsec = gallery.new_section("1-dimensional nonlinear elasticity")
    appdir = 'apps/stegoton_1d'
    description = """
        Evolution of two trains of solitary waves from an initial gaussian.
        """
    images = ('frame0000fig1', 'frame0003fig1', 'frame0005fig1')
    gsec.new_item(appdir, plotdir, description, images)
    #----------------------------------------------
    gsec = gallery.new_section("1-dimensional Euler equations")
    appdir = 'apps/euler_1d'
    description = """
        Woodward-Colella blast-wave interaction problem.
        """
    images = ('frame0000fig0', 'frame0003fig0', 'frame0010fig0')
    gsec.new_item(appdir, plotdir, description, images)
    #----------------------------------------------
       
    gallery.create('gallery_1d.rst')
    return gallery


def make_2d():
    gallery = Gallery("Gallery of 2d PyClaw applications")
    plotdir = '_plots'

    #----------------------------------------------
    gsec = gallery.new_section('2-dimensional advection')
    #----------------------------------------------
    appdir = 'apps/advection_2d'
    description = """
        Advecting square with periodic boundary conditions."""
    images = ('frame0000fig0', 'frame0002fig0', 'frame0004fig0')
    gsec.new_item(appdir, plotdir, description, images)
    #----------------------------------------------

    #----------------------------------------------
    gsec = gallery.new_section('2-dimensional variable-coefficient advection')
    #----------------------------------------------
    appdir = 'apps/advection_2d_annulus'
    description = """
        Advection in an annular region."""
    images = ('frame0000fig0', 'frame0004fig0', 'frame0008fig0')
    gsec.new_item(appdir, plotdir, description, images)
    #----------------------------------------------

    #----------------------------------------------
    gsec = gallery.new_section('2-dimensional acoustics')
    #----------------------------------------------
    appdir = 'apps/acoustics_2d_homogeneous'
    description = """
        Expanding radial acoustic wave in a homogeneous medium."""
    images = ('frame0000fig0', 'frame0002fig0', 'frame0004fig0')
    gsec.new_item(appdir, plotdir, description, images)
    #----------------------------------------------

    #----------------------------------------------
    gsec = gallery.new_section('2-dimensional variable-coefficient acoustics')
    #----------------------------------------------
    appdir = 'apps/acoustics_2d_variable'
    description = """
        Expanding radial acoustic wave in a two-material medium with an interface."""
    images = ('frame0000fig0', 'frame0010fig0', 'frame0020fig0')
    gsec.new_item(appdir, plotdir, description, images)
    #----------------------------------------------

    #----------------------------------------------
    gsec = gallery.new_section('2-dimensional shallow water equations')
    #----------------------------------------------
    appdir = 'apps/shallow_2d/'
    description = """Radial dam-break."""
    images = ('frame0000fig0', 'frame0004fig0', 'frame0010fig0')
    gsec.new_item(appdir, plotdir, description, images)
    #----------------------------------------------

    #----------------------------------------------
    gsec = gallery.new_section('2-dimensional shallow water on the sphere')
    #----------------------------------------------
    appdir = 'apps/shallow-sphere/'
    description = """Wavenumber 4 Rossby-Haurwitz wave on a rotating sphere."""
    images = ('frame0000fig0', 'frame0004fig0', 'frame0010fig0')
    gsec.new_item(appdir, plotdir, description, images)
    #----------------------------------------------

    #----------------------------------------------
    gsec = gallery.new_section('2-dimensional Euler equations')
    #----------------------------------------------
    appdir = 'apps/euler_2d'
    description = """
        Shock-bubble interaction."""
    images = ('frame0000fig0', 'frame0004fig0', 'frame0010fig0')
    gsec.new_item(appdir, plotdir, description, images)
    #----------------------------------------------

    #----------------------------------------------
    gsec = gallery.new_section('2-dimensional KPP equation')
    #----------------------------------------------
    appdir = 'apps/kpp/'
    description = """
        Non-convex flux example."""
    images = ('frame0000fig1', 'frame0004fig1', 'frame0010fig1')
    gsec.new_item(appdir, plotdir, description, images)
    #----------------------------------------------

    #----------------------------------------------
    gsec = gallery.new_section('2-dimensional p-system')
    #----------------------------------------------
    appdir = 'apps/psystem_2d/'
    description = """
        Radial wave in a checkerboard-like medium."""
    images = ('frame0000fig0', 'frame0004fig0', 'frame0010fig0')
    gsec.new_item(appdir, plotdir, description, images)
    #----------------------------------------------
        
    gallery.create('gallery_2d.rst')
    return gallery

def make_all():
    gallery_1d = make_1d()
    gallery_2d = make_2d()

    # make gallery of everything:
    gallery_all = Gallery(title="Gallery of all PyClaw applications")
    gallery_all.sections = gallery_1d.sections + gallery_2d.sections 
    gallery_all.create('gallery_all.rst')

if __name__ == "__main__":
    make_all()
