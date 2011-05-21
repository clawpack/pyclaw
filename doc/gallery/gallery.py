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
clawdir_default = os.environ.get('PYCLAW',None)
if clawdir_default is None:
    print "*** Error: set environment variable PYCLAW"

# Location for gallery files:
gallery_dir_default = os.path.join(clawdir_default,'doc/gallery')  

remake = True   # True ==> remake all thumbnails even if they exist.

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
              <IMG SRC="%s/doc/_static/clawlogo.jpg" WIDTH="200" HEIGHT="70"
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
                gfile.write('<p><b>$PYCLAW/%s</b> ... <a href="%s">README</a> ... <a href="%s">Plot Index</a><p>' \
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

    appdir = 'apps/advection/1d'
    description = """
        Advecting Gaussian with periodic boundary."""
    images = ('frame0000fig1', 'frame0004fig1')
    gsec.new_item(appdir, plotdir, description, images)
       
    gallery.create('test.html')



def make_1d():
    gallery = Gallery(title="Gallery of 1d applications")
    plotdir = '_plots'


    #----------------------------------------------
    gsec = gallery.new_section('1-dimensional advection')
    appdir = 'apps/advection/1d/constant'
    description = """
         Advecting Gaussian with periodic boundary."""
    images = ('frame0000fig1', 'frame0004fig1', 'frame0010fig1')
    gsec.new_item(appdir, plotdir, description, images)
    #----------------------------------------------
    gsec = gallery.new_section('1-dimensional variable-velocity advection')
    appdir = 'apps/advection/1d/variable'
    description = """
         Advecting Gaussian and square wave with periodic boundary."""
    images = ('frame0000fig1', 'frame0004fig1', 'frame0008fig1')
    gsec.new_item(appdir, plotdir, description, images)
     #----------------------------------------------
    gsec = gallery.new_section('1-dimensional acoustics')
    appdir = 'apps/acoustics/1d/homogeneous'
    description = """
         Acoustics equations with reflecting boundary at left and outflow at
         right."""
    images = ('frame0000fig1', 'frame0005fig1', 'frame0010fig1')
    gsec.new_item(appdir, plotdir, description, images)
    #----------------------------------------------
    gsec = gallery.new_section("1-dimensional Burgers' equation")
    appdir = 'apps/burgers/1d/'
    description = """
        Burgers' equation with sinusoidal initial data, steepening to
        N-wave.  """
    images = ('frame0000fig0', 'frame0003fig0', 'frame0006fig0')
    gsec.new_item(appdir, plotdir, description, images)
    #----------------------------------------------
    gsec = gallery.new_section("1-dimensional nonlinear elasticity")
    appdir = 'apps/elasticity/1d/stegoton'
    description = """
        Evolution of two trains of solitary waves from an initial gaussian.
        """
    images = ('frame0000fig1', 'frame0003fig1', 'frame0005fig1')
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

    #----------------------------------------------
    gsec = gallery.new_section('2-dimensional acoustics')
    #----------------------------------------------
    appdir = 'apps/acoustics/2d/homogeneous'
    description = """
        Expanding radial acoustic wave."""
    images = ('frame0000fig0', 'frame0002fig0', 'frame0004fig0')
    gsec.new_item(appdir, plotdir, description, images)
    #----------------------------------------------

    #----------------------------------------------
    gsec = gallery.new_section('2-dimensional Euler equations')
    #----------------------------------------------
    appdir = 'apps/euler/2d/'
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
        
    gallery.create('gallery_2d.html')
    return gallery

def make_all():
    gallery_1d = make_1d()
    gallery_2d = make_2d()

    # make gallery of everything:
    gallery_all = Gallery(title="Gallery of all applications")
    gallery_all.sections = gallery_1d.sections + gallery_2d.sections 
    gallery_all.create('gallery_all.html')

if __name__ == "__main__":
    make_all()
