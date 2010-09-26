"""
Module for making galleries of thumbnails allowing easy browsing of
applications directories.

These tools assume that the examples have already been run and the plots
produced using "make .plots".  

You should use the script python/run_examples.py to do this first.
"""

import os


# Main root for html links:
claw_html_root='http://localhost:50005'     


# Determine Clawpack directory:
clawdir_default = os.environ.get('CLAW',None)
if clawdir_default is None:
    print "*** Error: set environment variable CLAW"

# Location for gallery files:
gallery_dir_default = os.path.join(clawdir_default,'myclaw','gallery')  

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
        gfile.write("""<html>
              <BODY BGCOLOR="#FFFFE8" LINK="#7F0000" VLINK="#7F0000">
              <font FACE="TREBUCHET MS,HELVETICA,ARIAL">
              <a href="http://www.clawpack.org">
              <IMG SRC="%s/doc/images/clawlogo.jpg" WIDTH="200" HEIGHT="70"
                VSPACE="0" HSPACE="0" ALT="CLAWPACK" BORDER="0" LOOP="0"> </a>
              """ % claw_html_root)

        gfile.write("<h1>%s</h1>" % self.title)
        for gsec in self.sections:
            gfile.write("<p><h2>%s</h2>\n %s\n<p>\n" \
                  % (gsec.title, gsec.description))

            for gitem in gsec.items:
                gfile.write('<p>')
                readme = os.path.join(claw_html_root, gitem.appdir, \
                               'README.html')
                plotindex = os.path.join(claw_html_root, gitem.appdir, \
                               gitem.plotdir, '_PlotIndex.html')
                gfile.write('<p><b>$CLAW/%s</b> ... <a href="%s">README</a> ... <a href="%s">Plot Index</a><p>' \
                      % (gitem.appdir,readme,plotindex))
                gfile.write('<p>\n%s\n<p>\n' % gitem.description)
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
                    gfile.write('&nbsp;&nbsp;<a href="%s"><img src="%s"></a>' \
                        % (src_html, thumb_file))
                gfile.write('<p><hr><p>')

        gfile.write("\n</html>\n")
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

    appdir = 'apps/advection/1d/example1'
    description = """
        Advecting Gaussian with outflow boundary."""
    images = ('frame0000fig1', 'frame0004fig1')
    gsec.new_item(appdir, plotdir, description, images)
       
    gallery.create('test.html')



def make_1d():
    gallery = Gallery(title="Gallery of 1d applications")
    plotdir = '_plots'


    #----------------------------------------------
    gsec = gallery.new_section('1-dimensional advection')
    appdir = 'apps/advection/1d/example1'
    description = """
         Advecting Gaussian with outflow boundary."""
    images = ('frame0000fig1', 'frame0004fig1', 'frame0008fig1')
    gsec.new_item(appdir, plotdir, description, images)
    #----------------------------------------------
    gsec = gallery.new_section('1-dimensional acoustics')
    appdir = 'apps/acoustics/1d/example2'
    description = """
         Acoustics equations with reflecting boundary at left and outflow at
         right."""
    images = ('frame0000fig1', 'frame0005fig1', 'frame0010fig1')
    gsec.new_item(appdir, plotdir, description, images)
    #----------------------------------------------
    gsec = gallery.new_section("1-dimensional Burgers' equation")
    appdir = 'apps/burgers/1d/sine2n'
    description = """
        Burgers' equation with sinusoidal initial data, steepening to
        N-wave.  """
    images = ('frame0000fig0', 'frame0003fig0', 'frame0006fig0')
    gsec.new_item(appdir, plotdir, description, images)
    #----------------------------------------------
    gsec = gallery.new_section('1-dimensional Euler equations')
    appdir = 'apps/euler/1d/hump'
    description = """
        Initial hump in density and pressure propagating outwards and
        forming shocks.  """
    images = ('frame0000fig0', 'frame0002fig0', 'frame0004fig0')
    gsec.new_item(appdir, plotdir, description, images)
    #----------------------------------------------
       
    gallery.create('gallery_1d.html')
    return gallery


def make_2d():
    gallery = Gallery("Gallery of 2d applications")
    plotdir = '_plots'

    #----------------------------------------------
    gsec = gallery.new_section('2-dimensional advection')
    #----------------------------------------------
    appdir = 'apps/advection/2d/example1'
    description = """
        Advecting square with periodic boundary conditions."""
    images = ('frame0000fig0', 'frame0002fig0', 'frame0004fig0')
    gsec.new_item(appdir, plotdir, description, images)
    #----------------------------------------------
    appdir = 'apps/advection/2d/example1/amr'
    description = """
        Advecting square with periodic boundary conditions. Using AMR"""
    images = ('frame0000fig0', 'frame0002fig0', 'frame0004fig0', 'frame0004fig2')
    gsec.new_item(appdir, plotdir, description, images)
    #----------------------------------------------
    appdir = 'apps/advection/2d/example1_gauges/amr'
    description = """
        Advecting square with periodic boundary conditions. Using AMR and
        plotting solution at several gauges"""
    images = ('frame0000fig1', 'frame0006fig1', 'gauge0001fig300', 'gauge0002fig300')
    gsec.new_item(appdir, plotdir, description, images)
    #----------------------------------------------
    appdir = 'apps/advection/2d/swirl'
    description = """
        Advection equation with a time-dependent velocity field that swirls
        one way and then the other."""
    images = ('frame0000fig0', 'frame0005fig0', 'frame0010fig0')
    gsec.new_item(appdir, plotdir, description, images)
    #----------------------------------------------
    appdir = 'apps/advection/2d/swirl/amr'
    description = """
        Advection equation with a time-dependent velocity field that swirls
        one way and then the other.  Using AMR"""
    images = ('frame0000fig0', 'frame0005fig0', 'frame0010fig0')
    gsec.new_item(appdir, plotdir, description, images)
    #----------------------------------------------
    appdir = 'apps/advection/2d/annulus'
    description = """
        Advection equation with solid body rotation in an annular region,
        solved using polar coordinates."""
    images = ('frame0000fig0', 'frame0007fig0', 'frame0007fig1')
    gsec.new_item(appdir, plotdir, description, images)
    #----------------------------------------------
    appdir = 'apps/advection/2d/annulus/amr'
    description = """
        Advection equation with solid body rotation in an annular region,
        solved using polar coordinates.  Using AMR"""
    images = ('frame0000fig0', 'frame0007fig0', 'frame0007fig1', 'frame0007fig2')
    gsec.new_item(appdir, plotdir, description, images)
    #----------------------------------------------
    gsec = gallery.new_section("2-dimensional Burgers'")
    #----------------------------------------------
    appdir = 'apps/burgers/2d/pwconst'
    description = """
        2d Burgers' equation with piecewise constant initial data and
        periodic boundary conditions."""
    images = ('frame0000fig0', 'frame0005fig0', 'frame0010fig0')
    gsec.new_item(appdir, plotdir, description, images)
    #----------------------------------------------
       
    gallery.create('gallery_2d.html')
    return gallery

def make_geoclaw():
    gallery = Gallery("Gallery of GeoClaw applications")
    plotdir = '_plots'

    #----------------------------------------------
    gsec = gallery.new_section('Tsunami models using GeoClaw')
    #----------------------------------------------
    appdir = 'apps/tsunami/bowl-radial'
    description = """
        Radially symmetric solution in a bowl with Gaussian hump initial
        data.  Solved in Cartesian coordinates."""
    images = ('frame0003fig0', 'frame0009fig0', 'frame0014fig0')
    gsec.new_item(appdir, plotdir, description, images)
    #----------------------------------------------
    appdir = 'apps/tsunami/bowl-slosh'
    description = """
        Sloshing water in a parabolic bowl, initialized so the solution 
        should agree with a known exact solution in which the surface is
        always linear."""
    images = ('frame0000fig0', 'frame0002fig0', 'frame0004fig0', 'frame0002fig1')
    gsec.new_item(appdir, plotdir, description, images)
    #----------------------------------------------
    appdir = 'apps/tsunami/chile2010'
    description = """
        Tsunami of 27 February 2010 off the coast of Chile, with a
        comparison to DART buoy data.
        """
    images = ('frame0006fig0', 'frame0012fig0', 'frame0018fig0', 'gauge32412fig300')
    gsec.new_item(appdir, plotdir, description, images)
    #----------------------------------------------
    gallery.create('gallery_geoclaw.html')
    return gallery

def make_book():

    fvmhp_link = """<a href="http://www.amath.washington.edu/~claw/book.html"> 
		  FVMHP</a>"""
    gallery = Gallery(title="Gallery of examples from the book " + fvmhp_link)

    plotdir = '_plots'


    #----------------------------------------------
    gsec = gallery.new_section('Chapter 3: Characteristics and Riemann' + \
        'Problems for Linear Hyperbolic Equations.')
    appdir = 'book/chap3/acousimple'
    description = """Acoustics example with data breaking into left and
        right-going waves, with solid wall at left and outflow at right.
        """
    images = ('frame0000fig1', 'frame0005fig1', 'frame0010fig1')
    gsec.new_item(appdir, plotdir, description, images)
    #----------------------------------------------
    gsec = gallery.new_section('Chapter 6: High-resolution methods')
    appdir = 'book/chap6/compareadv'
    description = """
         Advection equation with first order upwind method.
         Periodic boundary conditions, looping around through unit period.
         """
    images = ('frame0000fig1', 'frame0005fig1', 'frame0010fig1')
    gsec.new_item(appdir, plotdir, description, images)
    #----------------------------------------------
    appdir = 'book/chap6/wavepacket'
    description = """
         Advection equation with first order upwind method.
         Periodic boundary conditions, looping around through unit period.
         """
    images = ('frame0000fig1', 'frame0005fig1', 'frame0010fig1')
    gsec.new_item(appdir, plotdir, description, images)
    #----------------------------------------------
    gsec = gallery.new_section('Chapter 7: Boundary Conditions and Ghost Cells')
    appdir = 'book/chap7/advinflow'
    description = """
         Advection equation with inflow boundary conditions.
         """
    images = ('frame0000fig0', 'frame0003fig0', 'frame0006fig0')
    gsec.new_item(appdir, plotdir, description, images)
    #----------------------------------------------
    appdir = 'book/chap7/acouinflow'
    description = """
         Acoustics equations with inflow boundary condition.
         """
    images = ('frame0000fig0', 'frame0005fig0', 'frame0010fig0')
    gsec.new_item(appdir, plotdir, description, images)
    #----------------------------------------------
    appdir = 'book/chap7/standing'
    description = """
         Acoustics eqautions with solid wall boundary conditions and a
         standing wave solution.
         """
    images = ('frame0000fig1', 'frame0005fig1', 'frame0010fig1', 'frame0015fig1')
    gsec.new_item(appdir, plotdir, description, images)
    #----------------------------------------------
    gsec = gallery.new_section('Chapter 9: Variable-Coefficient Linear Equations')
    appdir = 'book/chap9/acoustics/interface'
    description = """
         Acoustics equations with variable coefficients and a single
         interface  with no impedance jump.
         """
    images = ('frame0000fig1', 'frame0006fig1', 'frame0012fig1')
    gsec.new_item(appdir, plotdir, description, images)
    #----------------------------------------------
    appdir = 'book/chap9/acoustics/layered'
    description = """

         Acoustics equations with variable coefficients in periodic
         layered medium.
         """
    images = ('frame0000fig1', 'frame0006fig1', 'frame0012fig1')
    gsec.new_item(appdir, plotdir, description, images)
    #----------------------------------------------
    gsec = gallery.new_section('Chapter 10: Other Approaches to High Resolution.')
    appdir = 'book/chap10/tvb'
    description = """
         TVB method on the advection equation 
         with wave-packet data and periodic boundary conditions.
         """
    images = ('frame0000fig0', 'frame0005fig0', 'frame0010fig0')
    gsec.new_item(appdir, plotdir, description, images)

    #----------------------------------------------
    gsec = gallery.new_section('Chapter 11: Nonlinear Scalar Conservation Laws.')
    appdir = 'book/chap11/burgers'
    description = """
         Burgers' equation with oscillatory initial data to demonstrate
         decay to an N wave. 
         """
    images = ('frame0000fig1', 'frame0002fig1', 'frame0004fig1')
    #----------------------------------------------
    gsec.new_item(appdir, plotdir, description, images)
    appdir = 'book/chap11/greenlight'
    description = """
         Traffic flow equations.
         Riemann problem for cars starting at a green light. 
         """
    images = ('frame0000fig1', 'frame0005fig1', 'frame0010fig1')
    gsec.new_item(appdir, plotdir, description, images)
    #----------------------------------------------
    appdir = 'book/chap11/redlight'
    description = """
         Traffic flow equations.
         Riemann problem for cars stopping at a red light. 
         """
    images = ('frame0000fig1', 'frame0005fig1', 'frame0010fig1')
    gsec.new_item(appdir, plotdir, description, images)
    #----------------------------------------------
    gsec = gallery.new_section('Chapter 12: Finite Volume Methods for' + \
    'Nonlinear Scalar Conservation Laws.')
    appdir = 'book/chap12/efix'
    description = """
         Burgers' equation with a transonic rarefaction wave. Comparison of
         results obtained with or without the entropy fix.
         """
    images = ('frame0000fig0', 'frame0002fig0', 'frame0004fig0')
    gsec.new_item(appdir, plotdir, description, images)
    #----------------------------------------------
    appdir = 'book/chap12/llf'
    description = """
         Burgers' equation with a transonic rarefaction wave.
         Using the Local Lax-Friedrichs method.
         """
    images = ('frame0000fig0', 'frame0002fig0', 'frame0004fig0')
    gsec.new_item(appdir, plotdir, description, images)
    #----------------------------------------------
    appdir = 'book/chap12/nonconservative'
    description = """
         Burgers' equation solved with a nonconservative upwind method, to
         demonstrate that this does not approximate the weak solution
         properly. 
         """
    images = ('frame0000fig0', 'frame0010fig0', 'frame0020fig0')
    gsec.new_item(appdir, plotdir, description, images)
    #----------------------------------------------
    gsec = gallery.new_section('Chapter 13: Nonlinear Systems of Conservation Laws.')
    appdir = 'book/chap13/collide'
    description = """
         1d Shallow water equations with initial data consisting of two
         2-shocks, which collide and produce a 1-rarefaction and 2-shock. 
         """
    images = ('frame0000fig0', 'frame0003fig0', 'frame0009fig0')
    gsec.new_item(appdir, plotdir, description, images)
    #----------------------------------------------
    gsec = gallery.new_section('Chapter 17: Source Terms and Balance Laws.')
    appdir = 'book/chap17/advdiff'
    description = """
         Advection-diffusion equation solved with a fractional step method. 
         """
    images = ('frame0000fig1', 'frame0004fig1', 'frame0008fig1')
    gsec.new_item(appdir, plotdir, description, images)
    #----------------------------------------------
    gsec = gallery.new_section('Chapter 20: Multidimensional Scalar Equations.')
    appdir = 'book/chap20/burgers'
    description = """
         2D Burgers' equation at angle theta to the x-axis. 
         """
    images = ('frame0000fig0', 'frame0002fig0', 'frame0004fig0')
    #----------------------------------------------
    gsec.new_item(appdir, plotdir, description, images)
    appdir = 'book/chap20/rotate'
    description = """
         Advection equation for solid body rotation 
         """
    images = ('frame0000fig0', 'frame0001fig0', 'frame0002fig0', 'frame0003fig0')
    gsec.new_item(appdir, plotdir, description, images)
    #----------------------------------------------

    gsec = gallery.new_section('Chapter 21: Multidimensional Systems.')
    appdir = 'book/chap21/radialdam'
    description = """
         Radial dam-break problem for 2D shallow water equations. 
         """
    images = ('frame0000fig0', 'frame0001fig0', 'frame0002fig0', 'frame0003fig0')
    gsec.new_item(appdir, plotdir, description, images)

    #----------------------------------------------
    gallery.create('gallery_book.html')
    return gallery

def make_all():
    gallery_1d = make_1d()
    gallery_2d = make_2d()
    gallery_geoclaw = make_geoclaw()
    gallery_book = make_book()

    # make gallery of everything:
    gallery_all = Gallery(title="Gallery of all applications")
    gallery_all.sections = gallery_1d.sections + gallery_2d.sections \
                         + gallery_geoclaw.sections + gallery_book.sections
    gallery_all.create('gallery_all.html')

if __name__ == "__main__":
    make_all()
