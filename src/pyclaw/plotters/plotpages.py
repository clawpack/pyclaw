
"""
Module plotpages

Utilities for taking a set of plot files and creating a set of html and/or
latex/pdf pages displaying the plots.
"""


import os, time, string, glob

# Clawpack logo... not used on plot pages currently.
clawdir = os.getenv('CLAW')
if clawdir is not None:
    logo = os.path.join(clawdir,'doc/images/clawlogo.jpg')
    if not os.path.isfile(logo):
        logo = None


#===========================
class PlotPagesData(object):
#===========================
    
    def __init__(self):
        self.plotdir = 'plots'
        self.overwrite = True
        self.verbose = True

        self.latex = True                # make latex files for figures
        self.latex_fname = 'plots'       # name of latex file to create
        self.latex_title = 'Plots'       # title on top of latex file
        self.latex_itemsperpage = 'all'  # number of items on each page
        self.latex_itemsperline = 2      # number of items on each line
        self.latex_framesperpage = 'all' # number of frames on each page
        self.latex_framesperline = 2     # number of frames on each line
        self.latex_figsperline = 'all'   # number of figures on each line
        self.latex_makepdf = False       # run pdflatex on latex file
        self.latex_preplots = None       # latex to for top of page before plots

        self.html = True                # make html files for figures
        self.html_index_fname = '_PlotIndex.html'   # name of html index file
        self.html_index_title = 'Plot Index'    # title on top of index file
        self.html_homelink = None       # link to here from top of index file
        self.html_itemsperline = 2      # number of items on each line
        self.html_preplots = None       # html to for top of page before plots
        self.html_movie = True          # make html with java script for movie
        self.html_eagle = False         # use EagleClaw titles on html pages?

        self.gif_movie = False          # make animated gif movie of frames

        self.timeframes_framenos = 'all'
        self.timeframes_frametimes = {}
        self.timeframes_fignos = 'all'
        self.timeframes_fignames = {}
    
        self.pageitem_list = []


    def new_pageitem(self):
        """
        Create a new PageItem to be printed on this page
        """
        pageitem = PageItem()
        self.pageitem_list.append(pageitem)
        return pageitem

    def make_html(self):
        plots2html(self)
        path_to_html_index = os.path.join(os.path.abspath(self.plotdir), \
                                   self.html_index_fname)
        print_html_pointers(path_to_html_index)

    def make_latex(self):
        plots2latex(self)

    def make_pages(self):
        if self.latex:
            self.make_latex()
        if self.html:
            self.make_html()

    def make_timeframes_latex(self):
        timeframes2latex(self)

    def make_timeframes_html(self):
        timeframes2html(self)
        path_to_html_index = os.path.join(os.path.abspath(self.plotdir), \
                                   self.html_index_fname)
        print_html_pointers(path_to_html_index)


#=======================
class PageItem(object):
#=======================
    
    def __init__(self):
        self.fname = ''  # full path to png or other figure file
        self.html_index_entry = 'Untitled figure'  # Name for link from
                                                  # html index page
        self.html_preitem = None   # any html to be inserted in file
                                   # just before this item.
        self.latex_preitem = None  # any latex to be inserted in file
                                   # just before this item.
    
#=======================
class HtmlIndex(object):
#=======================

    def __init__(self, fname='_Index.html', title="Index"):
        self.fname = fname
        self.file = open(fname, 'w')
        self.file.write('<html><meta http-equiv="expires" content="0">')
        self.file.write('\n<title>Index</title>')
        self.file.write('\n<body><center><h1>%s</h1></center>\n' \
                   % title)

    def add(self,text = '', link = None):
        if link:
            self.file.write("""
                <p>
                <a href="%s">%s</a>
                """ % (link,text))
        else:
            self.file.write("""
                <p>
                %s
                """ % text)

    def close(self):
        self.file.write("\n</body></html>")
        self.file.close()
        path_to_html_index = os.path.join(os.getcwd(), \
                                   self.fname)
        print_html_pointers(path_to_html_index)
                    

#======================================================================
def plots2html(plot_pages_data):
#======================================================================
    """
    take a sequence of figure files and produce an html file to display them.
    """

    print '\n-----------------------------------\n'
    print '\nCreating html pages...\n'
    startdir = os.getcwd()
    ppd = plot_pages_data
    numitems = len(ppd.pageitem_list)   # number of page items (separate plots)

    if numitems == 0:
        print '*** Warning: 0 plots to put in html file'
        return 
        
    ppd =plot_pages_data
    try:
        cd_with_mkdir(ppd.plotdir, ppd.overwrite, ppd.verbose)
    except:
        print "*** Error, aborting plots2html"
        raise


    creationtime = current_time()
    
    
    for pageitem in ppd.pageitem_list:
        splitname = os.path.splitext(pageitem.fname)
        pageitem.hname = splitname[0] + '.html'
        pageitem.ext = splitname[1]

    
    # Create the index page:
    #-----------------------
    
    html = open(ppd.html_index_fname,'w')
    
    if ppd.html_eagle:
        html.write("""
          <html><meta http-equiv="expires" content="0">
          <title>EagleClaw Plot Index</title>
          <head>
          <link type="text/css" rel="stylesheet"
                href="http://localhost:50005/eagleclaw/eagleclaw.css">
          </head>
          <eagle1>EagleClaw -- Plot Index</eagle1>
          <eagle2>Easy Access Graphical Laboratory for Exploring Conservation
          Laws</eagle2>
          <p>
          <center><eagle3>
          <a href="../eaglemenu.html">Main Menu for this run-directory
          </a></eagle3> </center><p>
        """)


    else:
        html.write('<html><meta http-equiv="expires" content="0">')
        html.write('\n<title>%s</title>' % ppd.html_index_title)
        html.write('\n<body><center><h1>%s</h1></center>\n' \
                   % ppd.html_index_title)
	homelink = getattr(ppd,'html_homelink',None)
        if homelink:
	    html.write('<center><a href="%s">Back to %s</a></center>\n' \
	               % (homelink, homelink))

    html.write('<p>\n')
    html.write('<center>Plots created: %s &nbsp;&nbsp; ' % creationtime )
    html.write('</center><p>\n')

    html.write('<p>\n<table border=0 cellpadding=5 cellspacing=5>\n')



    html.write("""<p>\n<tr><td><b>All figures:</b></td>
          <td><a href="allfigures.html">html<a> &nbsp;&nbsp;&nbsp;  </td>""")
    if ppd.latex_makepdf:
        html.write('  <td><a href="%s.pdf">%s.pdf</a></td>' \
               % (ppd.latex_fname,ppd.latex_fname))
    html.write('</tr>\n')

    html.write('<p>\n<tr><td><b>Individual Figures:</b></td> </tr>\n')
    for pageitem in ppd.pageitem_list:
        html.write("""
           <td>%s</td>
           <td><a href="%s">html</a></td>
           <td><a href="%s">%s</a></td><tr>
           """ % (pageitem.html_index_entry, \
                  pageitem.hname,\
                  pageitem.fname, pageitem.fname))
    html.write('</table>\n')
    html.write('</body></html>')

    #----------------------------------------------------------------------
    
    # allfigures.html
    #-------------------
    html = open('allfigures.html', 'w')
    html.write("""
          <html><meta http-equiv="expires" content="0">
          <title>Plots</title>
          <p>
          <h1>All figures</h1>
          <p>
          <h3><a href=%s>Return to Plot Index</a> </h3>
          <p>
          <h3>Click on a figure to enlarge</h3>
          <p>
        """ % ppd.html_index_fname)

    for pageitem in ppd.pageitem_list:
        html.write('  <a href="%s"><img src="%s" width=400></a>\n' \
                % (pageitem.hname, pageitem.fname))
    
    html.write('\n<p><h3><a href=%s>Return to Plot Index</a> </h3>' \
                % ppd.html_index_fname)
    html.write('\n</center></body></html>\n')
    html.close()
    
    
    # individual html files for each figure
    #--------------------------------------

    for j in range(len(ppd.pageitem_list)):
        pageitem = ppd.pageitem_list[j]
        html = open(pageitem.hname, 'w')
        html.write("""
              <html><meta http-equiv="expires" content="0">
              <title>%s</title>
              <p>
              <h1>%s</h1>
              <p>

              <p>
            """ % (pageitem.html_index_entry,pageitem.html_index_entry))

        html.write("""
              <p><img src="%s" ><p>  
              <h3><a href=%s>Return to Plot Index</a> 
            """ % (pageitem.fname,ppd.html_index_fname))
        if j>0:
            html.write("&nbsp; ... &nbsp;  <a href=%s>Previous Figure</a> "\
                   % ppd.pageitem_list[j-1].hname)
        if j<len(ppd.pageitem_list)-2:
            html.write("&nbsp; ... &nbsp;  <a href=%s>Next Figure</a> "\
                   % ppd.pageitem_list[j+1].hname)
        html.write("\n</h3>")
    
    html.write('\n</center></body></html>\n')
    html.close()
    
    os.chdir(startdir)
    # end of plots2html

    
#======================================================================
def print_html_pointers(path_to_html_index):
#======================================================================
    #PlotPath = os.getcwd()
    #if PlotPath[0] != '/':
        #PlotPath = '/' + PlotPath
    #PlotPath.replace('\\','/') # for windows
    print "\n--------------------------------------------------------"
    print "\nPoint your browser to:"
    print "    file://%s" % path_to_html_index
    clawdir = os.getenv('CLAW')
    if clawdir in path_to_html_index:
        path_to_html_index = path_to_html_index.replace(clawdir,'')
        print "\nOr, if you have the Clawpack server running, point your browser to:"
        print "    http://localhost:50005%s"  % path_to_html_index


#======================================================================
def timeframes2html(plot_pages_data):
#======================================================================
    """
    Replaced by plotclaw_driver below, which also handles gauges
    ... RJL, 1/1/10

    take a sequence of figure files in format frame000NfigJ.png for
    N in framenos and J in fignos, and produce an html page for each
    along with an index page and a page showing all frames together for
    each figno and all figs together for each frameno.
      plot_pages_data.timeframes_framenos  is list of frames to use,
      plot_pages_data.timeframes_frametimes  is dictionary of time for each frame
      plot_pages_data.timeframes_fignos  is list of figs to use,
      plot_pages_data.timeframes_fignames  is dictionary of fig names for index.

    """


    print '\n-----------------------------------\n'
    print '\nCreating html pages for timestepping figures...\n'

    startdir = os.getcwd()
    ppd =plot_pages_data
        
    try:
        cd_with_mkdir(ppd.plotdir, ppd.overwrite, ppd.verbose)
    except:
        print "*** Error, aborting timeframes2html"
        raise

    creationtime = current_time()
    ppd = massage_frames_data(ppd)

    framenos = ppd.timeframes_framenos
    frametimes = ppd.timeframes_frametimes
    fignos = ppd.timeframes_fignos
    fignames = ppd.timeframes_fignames
    pngfile = ppd._pngfile
    htmlfile = ppd._htmlfile
    frametimef = ppd._frametimef
    allfigsfile = ppd._allfigsfile
            
    numframes = len(framenos)
    numfigs = len(fignos)
    

    eagle = getattr(ppd,'html_eagle',False)

    
    
    # Create the index page:
    #-----------------------
    
    html = open(ppd.html_index_fname,'w')
    
    if eagle:
        html.write("""
          <html><meta http-equiv="expires" content="0">
          <title>EagleClaw Plot Index</title>
          <head>
          <link type="text/css" rel="stylesheet"
                href="http://localhost:50005/eagleclaw/eagleclaw.css">
          </head>
          <eagle1>EagleClaw -- Plot Index</eagle1>
          <eagle2>Easy Access Graphical Laboratory for Exploring Conservation
          Laws</eagle2>
          <p>
          <center><eagle3>
          <a href="../eaglemenu.html">Main Menu for this run-directory
          </a></eagle3> </center><p>
        """)


    else:
        html.write('<html><meta http-equiv="expires" content="0">')
        html.write('\n<title>%s</title>\n<body>' % ppd.html_index_title)
        html.write('\n<center><h1>%s</h1></center>\n' \
                   % ppd.html_index_title)
	homelink = getattr(ppd,'html_homelink',None)
        if homelink:
	    html.write('<center><a href="%s">Back to %s</a></center>\n' \
	               % (homelink, homelink))

    html.write('<p>\n')
    html.write('<center>Plots created: %s &nbsp;&nbsp; ' % creationtime )
    html.write('</center><p>\n')

    html.write('<p>\n<table border=0 cellpadding=5 cellspacing=5>\n')


    if ppd.latex_makepdf:
        html.write('<p><tr><td><b>pdf file:</b></td>')
        html.write('\n   <td><a href="%s.pdf">%s.pdf</a></td>' \
               % (ppd.latex_fname,ppd.latex_fname))
        html.write('</tr>\n')

    if ppd.html_movie:
        html.write('<p><tr><td><b>js Movies:</b></td>')
        for figno in fignos:
            html.write('\n   <td><a href="moviefig%s.html">%s</a></td>' \
                           % (figno,fignames[figno]))
        html.write('</tr>\n')
    if ppd.gif_movie:
        html.write('<p><tr><td><b>gif Movies:</b></td>')
        for ifig in range(len(fignos)):
            html.write('\n   <td><a href="moviefig%s.gif">%s</a></td>' \
                           % (fignos[ifig],fignames[fignos[ifig]]))
        html.write('</tr>\n')
    html.write('<p>\n<tr><td><b>All Frames:</b></td> ')
    for ifig in range(len(fignos)):
        html.write('\n   <td><a href="allframesfig%s.html">%s</a></td>' \
                       % (fignos[ifig],fignames[fignos[ifig]]))
    html.write('</tr>\n')
    html.write('<p>\n<tr><td><b>Individual Frames:</b></td> </tr>\n')

    for frameno in framenos:

        html.write('\n <tr><td>Frame %s, t = %s:</td>' \
                    % (frameno,frametimef[frameno]))
        for figno in fignos:
            figname = fignames[figno]
            html.write('\n   <td><a href="%s">%s</a></td>' \
                       % (htmlfile[frameno,figno],figname))
        if numfigs > 1:
            html.write('\n<td><a href="%s">All figures</a></td>' \
                       % allfigsfile[frameno])
        html.write('</tr>\n')
    html.write('</table>\n')
    html.write('</body></html>')

    #----------------------------------------------------------------------
    
    # allframesfigJ.html
    #-------------------
    for figno in fignos:
        html = open('allframesfig%s.html' % figno, 'w')
        html.write('<html><meta http-equiv="expires" content="0">')
        html.write('<title>Plots</title>')
        html.write('<body>\n<center><h1>All Frames -- %s</h1>\n' \
                   % fignames[figno])
        html.write('<p>\n')
        html.write('\n<p><h3><a href=%s>Plot Index</a></h3>\n' \
                   % (ppd.html_index_fname))
        html.write('<p>\n')
        html.write('<h3>Click on a figure to enlarge</h3>\n')
        html.write('<p>\n')
    
        for frameno in framenos:
            html.write('  <a href="%s"><img src="%s" width=400></a>\n' \
                % (htmlfile[frameno,figno], pngfile[frameno,figno]))
    
        html.write('\n</center></body></html>\n')
        html.close()
    
    
    # allfigsframeN.html
    #-------------------
    if numfigs > 1:
        for iframe in range(numframes):
            frameno = framenos[iframe]
            html = open(allfigsfile[frameno], 'w')
            html.write('<html><meta http-equiv="expires" content="0">')
            html.write('<title>Plots</title>')
            html.write('<body>\n<center><h3>All Figures -- Frame %s' \
                 % framenos[iframe])
            html.write('&nbsp; at time t = %s' % frametimef[frameno])
            html.write('<p>\n')

            # Write link commands to previous and next frame:

            html.write('<p> <a href="%s">' % allfigsfile[framenos[0]])
            html.write('&#060; &#060;</a> &nbsp; &nbsp;\n')
            if iframe==0:
                html.write('&#060; &nbsp; &nbsp; ')
                html.write('\n<a href="_PlotIndex.html">Index</a>  ')
                if numframes > 1:
                    html.write('&nbsp; &nbsp; <a href="%s"> &#062; </a> ' \
                        % allfigsfile[framenos[1]])

            elif iframe==numframes-1:
                if numframes > 1:
                    html.write('\n<a href="%s"> &#060; </a>  &nbsp; &nbsp; ' \
                        % allfigsfile[framenos[iframe-1]])
                html.write('\n<a href="_PlotIndex.html">Index</a>  ')
                html.write(' &nbsp; &nbsp; &#062; ')

            else:
                html.write('\n<a href="%s"> &#060; </a>  &nbsp; &nbsp; ' \
                        % allfigsfile[framenos[iframe-1]])
                html.write('\n<a href="_PlotIndex.html">Index</a>  ')
                html.write('\n&nbsp; &nbsp; <a href="%s"> &#062; </a>  &nbsp; &nbsp; ' \
                        % allfigsfile[framenos[iframe+1]])

            html.write('&nbsp; &nbsp; \n<a href="%s"> ' \
                      % allfigsfile[framenos[numframes-1]])
            html.write('&#062; &#062;</a>  \n') 

            html.write('</h3><p>\n')
            html.write('<h3>Click on a figure to enlarge</h3>\n')
            html.write('<p>\n')
    
            for figno in fignos:
                html.write('  <a href="%s"><img src="%s" width=400></a>\n' \
                        % (htmlfile[frameno,figno], pngfile[frameno,figno]))

            # list of all frames at bottom:

            html.write('\n<p><b>Other frames:</b></a> &nbsp;&nbsp;')
            for frameno2 in framenos:
                if frameno2 == frameno:
                    html.write('\n<font color=red>%i</font>&nbsp;&nbsp;' \
                               % frameno)
                else:
                    html.write('\n<a href="%s">%i</a>  &nbsp; &nbsp; ' \
                           % (allfigsfile[frameno2],frameno2))
    
            html.write('\n</center></body></html>\n')
            html.close()
    
    
    # frameNfigJ.html  -- individual files for each frame/fig combo
    #----------------
    
    for iframe in range(numframes):
        frameno = framenos[iframe]
        for figno in fignos:
            html = open(htmlfile[frameno,figno],'w')
            html.write('<html><meta http-equiv="expires" content="0">\n')
            html.write('<title>Plots</title>')
            html.write('<body><center>\n')
            html.write('\n<h3>Frame %i ' % frameno)
            if numfigs > 1:
                html.write(' &nbsp;---&nbsp; %s' % fignames[figno] )
            html.write('&nbsp;&nbsp; at time t = %s</h3>' % frametimef[frameno])
        
            # Write link commands to previous and next frame:

            html.write('<p> <a href="%s">' % htmlfile[framenos[0],figno])
            html.write('&#060; &#060;</a> &nbsp; &nbsp;\n')
            if iframe==0:
                html.write('&#060; &nbsp; &nbsp; ')
                html.write('\n<a href="_PlotIndex.html">Index</a>  ')
                if numframes > 1:
                    html.write('&nbsp; &nbsp; <a href="%s"> &#062; </a> ' \
                        % htmlfile[framenos[1],figno])

            elif iframe==numframes-1:
                if numframes > 1:
                    html.write('\n<a href="%s"> &#060; </a>  &nbsp; &nbsp; ' \
                        % htmlfile[framenos[iframe-1],figno])
                html.write('\n<a href="_PlotIndex.html">Index</a>  ')
                html.write(' &nbsp; &nbsp; &#062; ')

            else:
                html.write('\n<a href="%s"> &#060; </a>  &nbsp; &nbsp; ' \
                        % htmlfile[framenos[iframe-1],figno])
                html.write('\n<a href="_PlotIndex.html">Index</a>  ')
                html.write('\n &nbsp; &nbsp;<a href="%s"> &#062; </a>  &nbsp; &nbsp; ' \
                        % htmlfile[framenos[iframe+1],figno])

            html.write('&nbsp; &nbsp; \n<a href="%s"> ' \
                      % htmlfile[framenos[numframes-1],figno])
            html.write('&#062; &#062;</a>  \n') 
        
            # image:
            html.write('\n\n <p><img src="%s"><p>  \n ' \
                        % pngfile[frameno,figno])

            html.write('\n\nImage source: &nbsp; %s'  \
                   % os.path.join(os.getcwd(),pngfile[frameno,figno]))

            # list of all figures at bottom of page:

            if numfigs > 1:
                html.write('\n<p><b>Other figures at this time:</b> &nbsp;&nbsp;')
                for figno2 in fignos:
                    if figno2 == figno:
                        html.write('\n<font color=red>%s</font>&nbsp;&nbsp;' \
                               % fignames[figno])
                    else:
                        html.write('\n<a href="%s">%s</a>  &nbsp; &nbsp; ' \
                           % (htmlfile[frameno,figno2],fignames[figno2]))
                html.write('\n<a href="%s"> All Figures </a>' \
                     % allfigsfile[frameno])

            # list of all frames at bottom of page:

            html.write('\n<p><b>Other frames:</b></a> &nbsp;&nbsp;')
            for frameno2 in framenos:
                if frameno2 == frameno:
                    html.write('\n<font color=red>%i</font>&nbsp;&nbsp;' \
                               % frameno)
                else:
                    html.write('\n<a href="%s">%i</a>  &nbsp; &nbsp; ' \
                           % (htmlfile[frameno2,figno],frameno2))
            html.write('\n<a href="allframesfig%s.html">  All Frames </a>' \
                     % figno)
        
            html.write('\n<p><h3><a href=%s>Plot Index</a></h3>' \
                      % (ppd.html_index_fname))
            if eagle:
                html.write("""<p><h3><a href="../eaglemenu.html">Main Menu for
                this run-directory</a></h3>  """)
            html.write('</center></body></html>')
            html.close()
    
    # moviefigJ.html
    #-------------------
    if ppd.html_movie:
        for figno in fignos:
            html = open('moviefig%s.html' % figno, 'w')
            text = htmlmovie(plot_pages_data,pngfile,framenos,figno)
            html.write(text)
            html.close()
    
    

    os.chdir(startdir)
    # end of timeframes2html
    

#=====================================
def htmlmovie(plot_pages_data,pngfile,framenos,figno):
#=====================================
    """
    Input:
     pngfile: a dictionary indexed by (frameno,figno) with value the
              corresponding png file for this figure.
     framenos: a list of frame numbers to include in movie
     figno: integer with the figure number for this movie.

    Returns:
     text for an html file that incorporates javascript to loop through the 
          plots one after another.  
    
    New 6/7/10: The html page also has buttons for controlling the movie.

    The parameter iterval below is the time interval between loading
    successive images and is in milliseconds.

    The img_width and img_height parameters do not seem to have any effect.
    """


    text = """
           <html>
           <head>
           <script language="Javascript">
           <!---
           var num_images = %s; """ % len(framenos)

    text += """
           var img_width = 800;
           var img_height = 600;
           var interval = 300;
           var images = new Array();


        function preload_images()
        {
            t = document.getElementById("progress");
            """

    i = 0
    for frameno in framenos:
        i = i+1
        text += """
        t.innerHTML = "Preloading image ";
        images[%s] = new Image(img_width, img_height);
        images[%s].src = "%s";
        """ % (i,i,pngfile[frameno,figno])
    text += """
        t.innerHTML = "";
        }

        function tick()
        {
          frame += 1;
          if (frame > num_images+1)
              frame = 1;

          document.movie.src = images[frame].src;
          tt = setTimeout("tick()", interval);
        }

        function startup()
        {
          preload_images();
          frame = 1;
          document.movie.src = images[frame].src;
        }
        function rewind()
        {
          frame = 1;
          document.movie.src = images[frame].src;
        }
        function start()
        {
          tt = setTimeout("tick()", interval);
        }
        function pause()
        {
          clearTimeout(tt);
        }
        function restart()
        {
          tt = setTimeout("tick()", interval);
        }
        function slower()
        {
          interval = interval / 0.7;
        }
        function faster()
        {
          interval = interval * 0.7;
        }

        // --->
        </script>
        </head>
        <body onLoad="startup();">

        <form>
        &nbsp;&nbsp;&nbsp;
        <input type="button" value="Start movie" onClick="start()">
        <input type="button" value="Pause" onClick="pause()">
        &nbsp;&nbsp;&nbsp;
        <input type="button" value="Rewind" onClick="rewind()">
        &nbsp;&nbsp;&nbsp;
        <input type="button" value="Slower" onClick="slower()">
        <input type="button" value="Faster" onClick="faster()">
        &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
        <a href="%s">Plot Index</a>
        </form>

        <p><div ID="progress"></div></p>
          <img src="%s" name="movie"/>

        </body>
        </html>
        """ % (plot_pages_data.html_index_fname,pngfile[framenos[0],figno])

    return text
    # end of htmlmovie


#======================================================================
def plots2latex(plot_pages_data):
#======================================================================
    """
    Take a list of figure files and produce latex file to display them.
    So far only works with time frames, not with gauges or other plots.
    """

    print '\n-----------------------------------\n'
    print '\nCreating latex file...\n'
    startdir = os.getcwd()
    ppd = plot_pages_data
    plotdir = ppd.plotdir
    numitems = len(ppd.pageitem_list)   # number of page items (separate plots)

    if numitems == 0:
        print '*** Warning: 0 plots to put in latex file'
        print 'No latex file generated'
        return 
        

    ppd =plot_pages_data

    try:
        cd_with_mkdir(ppd.plotdir, ppd.overwrite, ppd.verbose)
    except:
        print "*** Error, aborting plots2latex"
        raise


    creationtime = current_time()
    
    latexfile = open(ppd.latex_fname + '.tex', 'w')
    
    # latex header
    #-------------

    latexfile.write(r"""
        \documentclass[11pt]{article}
        \usepackage{graphicx}
        \setlength{\textwidth}{7.5in}
        \setlength{\oddsidemargin}{-0.5in}
        \setlength{\evensidemargin}{-0.5in}
        \setlength{\textheight}{9.2in}
        \setlength{\voffset}{-1in}
        \setlength{\headsep}{5pt}
        \begin{document}
        \begin{center}{\Large\bf %s}\vskip 5pt
        """ % ppd.latex_title)

    latexfile.write(r"""
        \bf Plots created {\tt %s} in directory: \vskip 5pt
        \verb+%s+
        \end{center}
        \vskip 5pt
        """ % (creationtime, startdir))

    # latex layout
    #-------------

    itemsperline = ppd.latex_itemsperline
    if itemsperline == 'all': itemsperline = numitems
    itemsperpage = ppd.latex_itemsperpage
    if itemsperpage == 'all': itemsperpage = numitems

    # width each plot must be:
    fwidth = 0.95/itemsperline

    # latex for each item:
    #---------------------

    itemlinecnt = 0
    itempagecnt = 0
    for pageitem in ppd.pageitem_list:
        if itempagecnt >= itemsperpage:
            latexfile.write('\\newpage \n')                       
            itempagecnt = 0
            itemlinecnt = 0
        elif itemlinecnt >= itemsperline:
            latexfile.write('\\vskip 10pt \n')                       
            itemlinecnt = 0
        itemlinecnt += 1
        itempagecnt += 1
        if pageitem.latex_preitem:
            latexfile.write(pageitem.latex_preitem)
        latexfile.write('\\includegraphics[width=%s\\textwidth]{%s}\n' \
                            % (fwidth,pageitem.fname))
        #latexfile.write('\\vskip 10pt \n')                       
    latexfile.write('\\end{document}\n')                         
    latexfile.close()
    print "\nLatex file created:  " 
    print "  %s/%s.tex" % (plotdir, ppd.latex_fname)
    print "\nUse pdflatex to create pdf file" 

    if ppd.latex_makepdf:
        try:
            os.system('pdflatex %s' % ppd.latex_fname)
            print "\nSuccessfully created pdf file:  %s/%s.pdf" \
                   % (plotdir, ppd.latex_fname)
        except:
            print  '*** pdflatex command failed'

    os.chdir(startdir)
    # end of plots2latex
    

#======================================================================
def cd_with_mkdir(newdir, overwrite=False, verbose=True):
#======================================================================

    newdir = os.path.abspath(newdir)
    if os.path.isfile(newdir):
        print "*** Error in cd_with_mkdir: directory specified is a file"
        raise
    elif (os.path.isdir(newdir) & overwrite):
        if verbose:
            print "Directory '%s' " % newdir
            print "    already exists, files may be overwritten "
    elif (os.path.isdir(newdir) & (not overwrite)):
        print "*** Error in cd_with_mkdir"
        print "Directory already exists:\n  ",newdir
        print "Remove directory with \n '  rm -r %s' " % newdir
        print "  and try again, or set overwrite=True "
        raise
    else:
        try:
            os.mkdir(newdir)
            if verbose:
                print "Created directory:\n   ", newdir
        except:
            print "*** Error in cd_with_mkdir"
            print "Cannot make directory: \n  ",newdir
            raise
    try:
        os.chdir(newdir)
    except:
        print "*** Error in cd_with_mkdir"
        print "Cannot change directory to \n  ",newdir


#======================================================================
def cd_plotdir(plot_pages_data):
#======================================================================

    verbose = False
    ppd = plot_pages_data
    if os.path.isfile(ppd.plotdir):
        print "*** Error in cd_plotdir: plotdir specified is a file"
        raise
    elif (os.path.isdir(ppd.plotdir) & ppd.overwrite):
        if verbose:
            print "Directory '%s' " % ppd.plotdir
            print "    already exists, files may be overwritten "
    elif (os.path.isdir(ppd.plotdir) & (not ppd.overwrite)):
        print "Directory '%s'" % ppd.plotdir
        print "  already exists"
        print "Remove directory with \n '  rm -r %s' " % ppd.plotdir
        print "  and try again, or set overwrite=True "
        print "*** Error in cd_plotdir"
        raise
    else:
        try:
            os.mkdir(ppd.plotdir)
        except:
            print "Cannot make directory ",ppd.plotdir
            print  "*** Error in cd_plotdir"
            raise
    try:
        os.chdir(ppd.plotdir)
    except:
        print "*** Error trying to cd to ",ppd.plotdir


#=====================================
def massage_frames_data(plot_pages_data):
#=====================================
    ppd = plot_pages_data

    try:
        framenos = ppd.timeframes_framenos
        frametimes = ppd.timeframes_frametimes
        fignos = ppd.timeframes_fignos
        fignames = ppd.timeframes_fignames
        prefix = getattr(ppd, 'timeframes_prefix', 'frame')
    except:
        print '*** Error: timeframes not set properly'
        return

    startdir = os.getcwd()
        
        
    if framenos == 'all' or fignos == 'all':
        # need to determine which figures exist
        files = glob.glob('%s*.png' % prefix)
        np = len(prefix)
        if framenos == 'all':
            framenos = set()
            for file in files:
                frameno = int(file[np:(np+4)])
                framenos.add(frameno)
            framenos = list(framenos)
            framenos.sort()
        if fignos == 'all':
            fignos = set()
            for file in files:
                figno = int(os.path.splitext(file)[0][(np+7):])
                fignos.add(figno)
            fignos = list(fignos)
            fignos.sort()

    for figno in fignos:
        if not fignames.has_key(figno):
            fignames[figno] = 'Solution'
            
    numframes = len(framenos)
    numfigs = len(fignos)

    if len(framenos) == 0:
        print '*** Warning: 0 frames to print'
    if len(fignos) == 0:
        print '*** Warning: 0 figures to print each frame'

    pngfile = {}
    htmlfile = {}
    frametimef = {}
    allfigsfile = {}
    #print '    Making png and html files for %i frames:' % numframes, framenos
    for frameno in framenos:
        framef = string.zfill(frameno,4)
        try:
            ftime = frametimes[frameno]
        except:
            ftime = '?'

        if ftime == '?':
            ftimef = ftime
        elif ((ftime == 0) | ((ftime > 0.001) & (ftime < 1000))):
            ftimef = '%9.5f' % ftime
        else:
            ftimef = '%12.5e' % ftime
        frametimef[frameno] = ftimef
        framef = string.zfill(frameno,4)
        for figno in fignos:
            pngfile[frameno,figno] = '%s%sfig%s.png'  % (prefix,framef,figno)
            htmlfile[frameno,figno] = '%s%sfig%s.html' % (prefix,framef,figno)
        allfigsfile[frameno] = 'allfigs%s%s.html' % (prefix,framef)

    ppd.timeframes_framenos = framenos
    ppd.timeframes_fignos = fignos
    ppd.timeframes_fignames = fignames
    ppd._pngfile = pngfile
    ppd._htmlfile = htmlfile
    ppd._frametimef = frametimef
    ppd._allfigsfile = allfigsfile
    return ppd




#======================================================================
def timeframes2latex(plot_pages_data):
#======================================================================
    """
    take a sequence of figure files in format frame000NfigJ.png for
    N in framenos and J in fignos, and produce a latex file containing
    them all.
      plot_pages_data.timeframes_framenos  is list of frames to use,
      plot_pages_data.timeframes_frametimes  is dictionary of time for each frame
      plot_pages_data.timeframes_fignos  is list of figs to use,
      plot_pages_data.timeframes_fignames  is dictionary of fig names for index.
      plot_pages_data.timeframes_prefix  is the string indicating how the 
                             files are named  ('frame' by default).
    """

    print '\n-----------------------------------\n'
    print 'Creating latex file...'

    startdir = os.getcwd()

    ppd =plot_pages_data
    try:
        cd_with_mkdir(ppd.plotdir, ppd.overwrite, ppd.verbose)
    except:
        print "*** Error, aborting timeframes2latex"
        raise

    creationtime = current_time()
    ppd = massage_frames_data(ppd)

    plotdir = ppd.plotdir

    framenos = ppd.timeframes_framenos
    frametimes = ppd.timeframes_frametimes
    fignos = ppd.timeframes_fignos
    fignames = ppd.timeframes_fignames
    pngfile = ppd._pngfile
            
    numframes = len(framenos)
    numfigs = len(fignos)
    
    latexfile = open(ppd.latex_fname + '.tex', 'w')
    
    # latex header
    #-------------

    latexfile.write(r"""
        \documentclass[11pt]{article}
        \usepackage{graphicx}
        \setlength{\textwidth}{7.5in}
        \setlength{\oddsidemargin}{-0.5in}
        \setlength{\evensidemargin}{-0.5in}
        \setlength{\textheight}{9.2in}
        \setlength{\voffset}{-1in}
        \setlength{\headsep}{5pt}
        \begin{document}
        \begin{center}{\Large\bf %s}\vskip 5pt
        """ % ppd.latex_title)

    latexfile.write(r"""
        \bf Plots created {\tt %s} in directory: \vskip 5pt
        \verb+%s+
        \end{center}
        \vskip 5pt
        """ % (creationtime, startdir))

    # latex layout
    #-------------

    # determine how many plots should appear on each page and line:
    framesperpage = ppd.latex_framesperpage 
    if framesperpage == 'all':
        framesperpage = len(framenos)
    framesperline = ppd.latex_framesperline 
    if framesperline == 'all':
        framesperline = len(framenos)
    figsperline = ppd.latex_figsperline      
    if figsperline == 'all':
        figsperline = len(fignos)
    if (figsperline < len(fignos)) & (framesperline > 1):
        print '*** Incompatible layout: resetting framesperline to 1'
        framesperline = 1
    totalperline = framesperline * figsperline      
    if totalperline < 1:
        print '*** Warning: 0 figures per line requested in latex file'
        print 'No latex file generated due to format error'
        return 

    # width each plot must be:
    fwidth = 0.95/totalperline

    framecnt = 0
    for frameno in framenos:
        #latexfile.write('\\centerline{\Large Frame %s at time = %s' \
        #       % (frameno frametime[frameno])
        if framecnt >= framesperpage:
            latexfile.write('\\newpage \n')                       
            framecnt = 0
        elif framecnt >= framesperline:
            latexfile.write('\\vskip 10pt \n')                       
            framecnt = 0
        framecnt += 1
        figcnt = 0
        for figno in fignos:
            if figcnt >= figsperline:
                latexfile.write('\\vskip 10pt \n')                       
                figcnt = 0
            figcnt += 1
            latexfile.write('\\includegraphics[width=%s\\textwidth]{%s}\n' \
                            % (fwidth,pngfile[frameno,figno]))
        #latexfile.write('\\vskip 10pt \n')                       
    latexfile.write('\\end{document}\n')                         
    latexfile.close()

    print "\nLatex file created:  " 
    print "  %s/%s.tex" % (plotdir, ppd.latex_fname)
    print "\nUse pdflatex to create pdf file" 
    if ppd.latex & ppd.latex_makepdf:
        try:
            os.system('pdflatex %s' % ppd.latex_fname)
        except:
            print  '*** pdflatex command failed'
        print "\nSuccessfully created pdf file:  %s/%s.pdf" \
                % (plotdir, ppd.latex_fname)

    os.chdir(startdir)
    # end of timeframes2latex

    

#============================
def test(makeplots = True):
#============================
    try:
        from pylab import linspace,clf,plot,title,savefig,mod
    except:
        print '*** Error: could not import pylab'
        return

    ppd = PlotPagesData()

    ppd.plotdir = 'plots'

    ppd.html = True

    ppd.latex = True 
    ppd.latex_itemsperline = 2
    ppd.latex_itemsperpage = 4
    ppd.latex_makepdf = False

    # create test figures:
    x = linspace(0,1,201)
    for n in range(6):
        fname = 'plot%s.png' % n
        fname_savefig = os.path.join(ppd.plotdir, fname)
        if makeplots:
            clf()
            y = x**n
            plot(x,y)
            title('$f(x) = x^%s$' % n)
            savefig(fname_savefig)
        pid = ppd.new_pageitem()
        pid.fname = fname
        pid.html_index_entry = "Plot of x^%s" % n
        if mod(n,2) == 0:
            pid.latex_preitem = r"""
              \vskip 5pt \noindent{\large\bf Plot of $x^%s$}\vskip 2pt""" % n
    
    ppd.make_pages()


#============================
def clawtest():
#============================
    html_index = HtmlIndex(fname='vary_mx_index.html', \
          title='Results from running vary_mx.py')
    html_index.add(text = 'Experiments varying mx')

    for mx in [50, 100]:
        ppd = PlotPagesData()
    
        outdir = 'output.mx%s' % mx
        ppd.plotdir = outdir
        ppd.overwrite = True
    
        ppd.html = True
        ppd.html_index_title = 'Clawpack Plots with mx = %s' % mx
    
        ppd.latex = True 
        ppd.latex_makepdf = False
    
        ppd.timeframes_framenos = 'all'
        ppd.timeframes_frametimes = {}
        ppd.timeframes_fignos = 'all'
        ppd.timeframes_fignames = {}
    
        ppd.make_timeframes_html()
        ppd.make_timeframes_latex()

        # update global index:
        mx_text = 'mx = %s' % mx
        mx_index = os.path.join(outdir, ppd.html_index_fname)
        html_index.add(text = mx_text, link = mx_index)

    html_index.close()

#-----------------------------
def current_time(addtz=False):
#-----------------------------
    # determine current time and reformat:
    time1 = time.asctime()
    year = time1[-5:]
    day = time1[:-14]
    hour = time1[-13:-5]
    current_time = day + year + ' at ' + hour
    if addtz:
        current_time = current_time + ' ' + time.tzname[time.daylight]
    return current_time



#======================================================================
def plotclaw2html(plotdata):
#======================================================================

    """
    Create and html index and html pages for each figure created from the
    specified plotdata.

    Assumes the following types of figures may exist:
       time frame figures of the form frame000NfigJ.png 
       gauge figures of the form gauge000NfigJ.png 
       other each_run type figures of the form figJ.png
       other figures can be specified in a dictionary plotdata.otherfigs

      plotdata.timeframes_framenos  is list of frames to use,
      plotdata.timeframes_frametimes  is dictionary of time for each frame
      plotdata.timeframes_fignos  is list of figs to use,
      plotdata.timeframes_fignames  is dictionary of fig names for index.
      plotdata.gauges_gaugenos  is list of gauges to use,
      plotdata.gauges_fignos  is list of figs to use,
      plotdata.gauges_fignames  is dictionary of fig names for index.
      plotdata.eachrun_fignos  is list of figs to use,
      plotdata.eachrun_fignames  is dictionary of fig names for index.

    """


    print '\n-----------------------------------\n'
    print '\nCreating html pages for figures...\n'

    startdir = os.getcwd()
        
    try:
        cd_with_mkdir(plotdata.plotdir, plotdata.overwrite, plotdata.verbose)
    except:
        print "*** Error, aborting timeframes2html"
        raise

    creationtime = current_time()
    plotdata = massage_frames_data(plotdata)
    plotdata = massage_gauges_data(plotdata)

    framenos = plotdata.timeframes_framenos
    frametimes = plotdata.timeframes_frametimes
    fignos = plotdata.timeframes_fignos
    fignames = plotdata.timeframes_fignames
    pngfile = plotdata._pngfile
    htmlfile = plotdata._htmlfile
    frametimef = plotdata._frametimef
    allfigsfile = plotdata._allfigsfile
    gauge_pngfile = plotdata._gauge_pngfile
    gauge_htmlfile = plotdata._gauge_htmlfile
    gauge_allfigsfile = plotdata._gauge_allfigsfile
            
    numframes = len(framenos)
    numfigs = len(fignos)
    

    eagle = getattr(plotdata,'html_eagle',False)

    
    
    # Create the index page:
    #-----------------------
    
    html = open(plotdata.html_index_fname,'w')
    
    if eagle:
        html.write("""
          <html><meta http-equiv="expires" content="0">
          <title>EagleClaw Plot Index</title>
          <head>
          <link type="text/css" rel="stylesheet"
                href="http://localhost:50005/eagleclaw/eagleclaw.css">
          </head>
          <eagle1>EagleClaw -- Plot Index</eagle1>
          <eagle2>Easy Access Graphical Laboratory for Exploring Conservation
          Laws</eagle2>
          <p>
          <center><eagle3>
          <a href="../eaglemenu.html">Main Menu for this run-directory
          </a></eagle3> </center><p>
        """)


    else:
        html.write('<html><meta http-equiv="expires" content="0">')
        html.write('\n<title>%s</title>' % plotdata.html_index_title)
        html.write('\n<body><center><h1>%s</h1></center>\n' \
                   % plotdata.html_index_title)
	homelink = getattr(plotdata,'html_homelink',None)
        if homelink:
	    html.write('<center><a href="%s">Back to %s</a></center>\n' \
	               % (homelink, homelink))

    html.write('<p>\n')
    html.write('<center>Plots created: %s &nbsp;&nbsp; ' % creationtime )
    html.write('</center><p>\n')

    html.write('<p>\n<b>Go to:</b>\n')
    #html.write('<p>\n&nbsp;&nbsp; <a href="#timeframes">Time Frames</a>\n')
    gaugenos = plotdata.gauges_gaugenos
    numgauges = len(gaugenos)
    if (numgauges>0) & (len(plotdata.gauges_fignos)>0):
        html.write('&nbsp;&nbsp; <a href="#gauges">Gauges</a>\n')
    html.write('&nbsp;&nbsp; <a href="#eachrun">Other plots</a>\n')

    html.write('<p>\n<a name="timeframes"><h3>Time frames:</h3></a>\n')
    html.write('<p>\n<table border=0 cellpadding=5 cellspacing=5>\n')


    if plotdata.latex_makepdf:
        html.write('<p><tr><td><b>pdf file:</b></td>')
        html.write('\n   <td><a href="%s.pdf">%s.pdf</a></td>' \
               % (plotdata.latex_fname,plotdata.latex_fname))
        html.write('</tr>\n')

    if plotdata.html_movie:
        html.write('<p><tr><td><b>js Movies:</b></td>')
        for figno in fignos:
            html.write('\n   <td><a href="moviefig%s.html">%s</a></td>' \
                           % (figno,fignames[figno]))
        html.write('</tr>\n')
    if plotdata.gif_movie:
        html.write('<p><tr><td><b>gif Movies:</b></td>')
        for ifig in range(len(fignos)):
            html.write('\n   <td><a href="moviefig%s.gif">%s</a></td>' \
                           % (fignos[ifig],fignames[fignos[ifig]]))
        html.write('</tr>\n')
    html.write('<p>\n<tr><td><b>All Frames:</b></td> ')
    for ifig in range(len(fignos)):
        html.write('\n   <td><a href="allframesfig%s.html">%s</a></td>' \
                       % (fignos[ifig],fignames[fignos[ifig]]))
    html.write('</tr>\n')
    html.write('<p>\n<tr><td><b>Individual Frames:</b></td> </tr>\n')

    for frameno in framenos:

        html.write('\n <tr><td>Frame %s, t = %s:</td>' \
                    % (frameno,frametimef[frameno]))
        for figno in fignos:
            figname = fignames[figno]
            html.write('\n   <td><a href="%s">%s</a></td>' \
                       % (htmlfile[frameno,figno],figname))
        if numfigs > 1:
            html.write('\n<td><a href="%s">All figures</a></td>' \
                       % allfigsfile[frameno])
        html.write('</tr>\n')
    html.write('</table>\n')

    # Gauges:
    #----------------
    gaugenos = plotdata.gauges_gaugenos
    numgauges = len(gaugenos)
    fignos = plotdata.gauges_fignos
    fignames = plotdata.gauges_fignames
    if (numgauges>0) & (len(fignos)>0):
        html.write('<p>\n<a name="gauges"><h3>Gauges:</h3></a>\n')
        html.write('<p>\n<table border=0 cellpadding=5 cellspacing=5>\n')
        html.write('<p>\n<tr><td><b>All Gauges:</b></td> ')
        for ifig in range(len(fignos)):
            html.write('\n   <td><a href="allgaugesfig%s.html">%s</a></td>' \
                           % (fignos[ifig],fignames[fignos[ifig]]))
        html.write('</tr>\n')
        html.write('<p>\n<tr><td><b>Individual Gauges:</b></td> </tr>\n')
    
        for gaugeno in gaugenos:
    
            html.write('\n <tr><td>Gauge %s:</td>' \
                        % (gaugeno))
            for figno in fignos:
                figname = fignames[figno]
                html.write('\n   <td><a href="%s">%s</a></td>' \
                           % (gauge_htmlfile[gaugeno,figno],figname))
            if numfigs > 1:
                html.write('\n<td><a href="%s">All figures</a></td>' \
                           % gauge_allfigsfile[gaugeno])
            html.write('</tr>\n')
        html.write('</table>\n')
    #else:
    #   html.write('<p>None\n')   

    html.write('<p>\n<a name="eachrun"><h3>Other plots:</h3></a>\n')
    if len(plotdata.otherfigure_dict)==0:
        html.write('<p>None\n')  
    else:
        html.write('<p><ul>\n')  
        for name in plotdata.otherfigure_dict.iterkeys():
            otherfigure = plotdata.otherfigure_dict[name]
            fname = otherfigure.fname
            extension = os.path.splitext(fname)[1]
            if extension not in ['.png','.jpg']:
                print "*** Error: unrecognized extension in ",fname
                print "*** Use .png or .jpg"
            else:
                makefig = otherfigure.makefig
                if makefig:
                    if type(makefig)==str:
                        try:
                            exec(makefig)
                        except:
                            print "*** Problem executing makefig "
                            print "    for otherfigure ",name
                    else:
                        try:
                            makefig(plotdata)
                        except:
                            print "*** Problem executing makefig function"
                            print "    for otherfigure ",name
                    try:
                        from pylab import savefig
                        savefig(fname)
                    except:
                        print "*** Problem importing pylab or executing savefig"
                html.write('<p><li><a href="%s">%s</a>\n' %(fname,name))  
        html.write('<p></ul>\n')  

    
    html.write('</body></html>')

    # end of index
    #----------------------------------------------------------------------

    fignos = plotdata.timeframes_fignos
    fignames = plotdata.timeframes_fignames
    
    # allframesfigJ.html
    #-------------------
    for figno in fignos:
        html = open('allframesfig%s.html' % figno, 'w')
        html.write('<html><meta http-equiv="expires" content="0">')
        html.write('<title>Plots</title>')
        html.write('<body>\n<center><h1>All Frames -- %s</h1>\n' \
                   % fignames[figno])
        html.write('<p>\n')
        html.write('\n<p><h3><a href=%s>Plot Index</a></h3>\n' \
                   % (plotdata.html_index_fname))
        html.write('<p>\n')
        html.write('<h3>Click on a figure to enlarge</h3>\n')
        html.write('<p>\n')
    
        for frameno in framenos:
            html.write('  <a href="%s"><img src="%s" width=400></a>\n' \
                % (htmlfile[frameno,figno], pngfile[frameno,figno]))
    
        html.write('\n</center></body></html>\n')
        html.close()
    
    
    # allfigsframeN.html
    #-------------------
    if numfigs > 1:
        for iframe in range(numframes):
            frameno = framenos[iframe]
            html = open(allfigsfile[frameno], 'w')
            html.write('<html><meta http-equiv="expires" content="0">')
            html.write('<title>Plots</title>')
            html.write('<body>\n<center><h3>All Figures -- Frame %s' \
                 % framenos[iframe])
            html.write('&nbsp; at time t = %s' % frametimef[frameno])
            html.write('<p>\n')

            # Write link commands to previous and next frame:

            html.write('<p> <a href="%s">' % allfigsfile[framenos[0]])
            html.write('&#060; &#060;</a> &nbsp; &nbsp;\n')
            if iframe==0:
                html.write('&#060; &nbsp; &nbsp; ')
                html.write('\n<a href="_PlotIndex.html">Index</a>  ')
                if numframes > 1:
                    html.write('&nbsp; &nbsp; <a href="%s"> &#062; </a> ' \
                        % allfigsfile[framenos[1]])

            elif iframe==numframes-1:
                if numframes > 1:
                    html.write('\n<a href="%s"> &#060; </a>  &nbsp; &nbsp; ' \
                        % allfigsfile[framenos[iframe-1]])
                html.write('\n<a href="_PlotIndex.html">Index</a>  ')
                html.write(' &nbsp; &nbsp; &#062; ')

            else:
                html.write('\n<a href="%s"> &#060; </a>  &nbsp; &nbsp; ' \
                        % allfigsfile[framenos[iframe-1]])
                html.write('\n<a href="_PlotIndex.html">Index</a>  ')
                html.write('\n&nbsp; &nbsp; <a href="%s"> &#062; </a>  &nbsp; &nbsp; ' \
                        % allfigsfile[framenos[iframe+1]])

            html.write('&nbsp; &nbsp; \n<a href="%s"> ' \
                      % allfigsfile[framenos[numframes-1]])
            html.write('&#062; &#062;</a>  \n') 

            html.write('</h3><p>\n')
            html.write('<h3>Click on a figure to enlarge</h3>\n')
            html.write('<p>\n')
    
            for figno in fignos:
                html.write('  <a href="%s"><img src="%s" width=400></a>\n' \
                        % (htmlfile[frameno,figno], pngfile[frameno,figno]))

            # list of all frames at bottom:

            html.write('\n<p><b>Other frames:</b></a> &nbsp;&nbsp;')
            for frameno2 in framenos:
                if frameno2 == frameno:
                    html.write('\n<font color=red>%i</font>&nbsp;&nbsp;' \
                               % frameno)
                else:
                    html.write('\n<a href="%s">%i</a>  &nbsp; &nbsp; ' \
                           % (allfigsfile[frameno2],frameno2))
    
            html.write('\n</center></body></html>\n')
            html.close()
    
    
    # frameNfigJ.html  -- individual files for each frame/fig combo
    #----------------
    
    for iframe in range(numframes):
        frameno = framenos[iframe]
        for figno in fignos:
            html = open(htmlfile[frameno,figno],'w')
            html.write('<html><meta http-equiv="expires" content="0">\n')
            html.write('<title>Plots</title>')
            html.write('<body><center>\n')
            html.write('\n<h3>Frame %i ' % frameno)
            if numfigs > 1:
                html.write(' &nbsp;---&nbsp; %s' % fignames[figno] )
            html.write('&nbsp;&nbsp; at time t = %s</h3>' % frametimef[frameno])
        
            # Write link commands to previous and next frame:

            html.write('<p> <a href="%s">' % htmlfile[framenos[0],figno])
            html.write('&#060; &#060;</a> &nbsp; &nbsp;\n')
            if iframe==0:
                html.write('&#060; &nbsp; &nbsp; ')
                html.write('\n<a href="_PlotIndex.html">Index</a>  ')
                if numframes > 1:
                    html.write('&nbsp; &nbsp; <a href="%s"> &#062; </a> ' \
                        % htmlfile[framenos[1],figno])

            elif iframe==numframes-1:
                if numframes > 1:
                    html.write('\n<a href="%s"> &#060; </a>  &nbsp; &nbsp; ' \
                        % htmlfile[framenos[iframe-1],figno])
                html.write('\n<a href="_PlotIndex.html">Index</a>  ')
                html.write(' &nbsp; &nbsp; &#062; ')

            else:
                html.write('\n<a href="%s"> &#060; </a>  &nbsp; &nbsp; ' \
                        % htmlfile[framenos[iframe-1],figno])
                html.write('\n<a href="_PlotIndex.html">Index</a>  ')
                html.write('\n &nbsp; &nbsp;<a href="%s"> &#062; </a>  &nbsp; &nbsp; ' \
                        % htmlfile[framenos[iframe+1],figno])

            html.write('&nbsp; &nbsp; \n<a href="%s"> ' \
                      % htmlfile[framenos[numframes-1],figno])
            html.write('&#062; &#062;</a>  \n') 
        
            # image:
            html.write('\n\n <p><img src="%s"><p>  \n ' \
                        % pngfile[frameno,figno])

            html.write('\n\nImage source: &nbsp; %s'  \
                   % os.path.join(os.getcwd(),pngfile[frameno,figno]))

            # list of all figures at bottom of page:

            if numfigs > 1:
                html.write('\n<p><b>Other figures at this time:</b> &nbsp;&nbsp;')
                for figno2 in fignos:
                    if figno2 == figno:
                        html.write('\n<font color=red>%s</font>&nbsp;&nbsp;' \
                               % fignames[figno])
                    else:
                        html.write('\n<a href="%s">%s</a>  &nbsp; &nbsp; ' \
                           % (htmlfile[frameno,figno2],fignames[figno2]))
                html.write('\n<a href="%s"> All Figures </a>' \
                     % allfigsfile[frameno])

            # list of all frames at bottom of page:

            html.write('\n<p><b>Other frames:</b></a> &nbsp;&nbsp;')
            for frameno2 in framenos:
                if frameno2 == frameno:
                    html.write('\n<font color=red>%i</font>&nbsp;&nbsp;' \
                               % frameno)
                else:
                    html.write('\n<a href="%s">%i</a>  &nbsp; &nbsp; ' \
                           % (htmlfile[frameno2,figno],frameno2))
            html.write('\n<a href="allframesfig%s.html">  All Frames </a>' \
                     % figno)
        
            html.write('\n<p><h3><a href=%s>Plot Index</a></h3>' \
                      % (plotdata.html_index_fname))
            if eagle:
                html.write("""<p><h3><a href="../eaglemenu.html">Main Menu for
                this run-directory</a></h3>  """)
            html.write('</center></body></html>')
            html.close()
    
    # moviefigJ.html
    #-------------------
    if plotdata.html_movie:
        for figno in fignos:
            html = open('moviefig%s.html' % figno, 'w')
            text = htmlmovie(plotdata,pngfile,framenos,figno)
            html.write(text)
            html.close()
    
    

    #----------------------------------------------------------------------
    fignos = plotdata.gauges_fignos
    fignames = plotdata.gauges_fignames
    
    # allgaugesfigJ.html
    #-------------------
    for figno in fignos:
        html = open('allgaugesfig%s.html' % figno, 'w')
        html.write('<html><meta http-equiv="expires" content="0">')
        html.write('<title>Plots</title>')
        html.write('<body>\n<center><h1>All Gauges -- %s</h1>\n' \
                   % fignames[figno])
        html.write('<p>\n')
        html.write('\n<p><h3><a href=%s>Plot Index</a></h3>\n' \
                   % (plotdata.html_index_fname))
        html.write('<p>\n')
        html.write('<h3>Click on a figure to enlarge</h3>\n')
        html.write('<p>\n')
    
        for gaugeno in gaugenos:
            html.write('  <a href="%s"><img src="%s" width=400></a>\n' \
                % (gauge_htmlfile[gaugeno,figno], gauge_pngfile[gaugeno,figno]))
    
        html.write('\n</center></body></html>\n')
        html.close()
    
    
    # allfigsgaugeN.html
    #-------------------
    if numfigs > 1:
        for igauge in range(numgauges):
            gaugeno = gaugenos[igauge]
            html = open(gauge_allfigsfile[gaugeno], 'w')
            html.write('<html><meta http-equiv="expires" content="0">')
            html.write('<title>Plots</title>')
            html.write('<body>\n<center><h3>All Figures -- Gauge %s' \
                 % gaugenos[igauge])
            html.write('<p>\n')

            # Write link commands to previous and next gauge:

            html.write('<p> <a href="%s">' % gauge_allfigsfile[gaugenos[0]])
            html.write('&#060; &#060;</a> &nbsp; &nbsp;\n')
            if igauge==0:
                html.write('&#060; &nbsp; &nbsp; ')
                html.write('\n<a href="_PlotIndex.html">Index</a>  ')
                if numgauges > 1:
                    html.write('&nbsp; &nbsp; <a href="%s"> &#062; </a> ' \
                        % gauge_allfigsfile[gaugenos[1]])

            elif igauge==numgauges-1:
                if numgauges > 1:
                    html.write('\n<a href="%s"> &#060; </a>  &nbsp; &nbsp; ' \
                        % gauge_allfigsfile[gaugenos[igauge-1]])
                html.write('\n<a href="_PlotIndex.html">Index</a>  ')
                html.write(' &nbsp; &nbsp; &#062; ')

            else:
                html.write('\n<a href="%s"> &#060; </a>  &nbsp; &nbsp; ' \
                        % gauge_allfigsfile[gaugenos[igauge-1]])
                html.write('\n<a href="_PlotIndex.html">Index</a>  ')
                html.write('\n&nbsp; &nbsp; <a href="%s"> &#062; </a>  &nbsp; &nbsp; ' \
                        % gauge_allfigsfile[gaugenos[igauge+1]])

            html.write('&nbsp; &nbsp; \n<a href="%s"> ' \
                      % gauge_allfigsfile[gaugenos[numgauges-1]])
            html.write('&#062; &#062;</a>  \n') 

            html.write('</h3><p>\n')
            html.write('<h3>Click on a figure to enlarge</h3>\n')
            html.write('<p>\n')
    
            for figno in fignos:
                html.write('  <a href="%s"><img src="%s" width=400></a>\n' \
                        % (gauge_htmlfile[gaugeno,figno], gauge_pngfile[gaugeno,figno]))

            # list of all gauges at bottom:

            html.write('\n<p><b>Other gauges:</b></a> &nbsp;&nbsp;')
            for gaugeno2 in gaugenos:
                if gaugeno2 == gaugeno:
                    html.write('\n<font color=red>%i</font>&nbsp;&nbsp;' \
                               % gaugeno)
                else:
                    html.write('\n<a href="%s">%i</a>  &nbsp; &nbsp; ' \
                           % (gauge_allfigsfile[gaugeno2],gaugeno2))
    
            html.write('\n</center></body></html>\n')
            html.close()
    
    
    # gaugeNfigJ.html  -- individual files for each gauge/fig combo
    #----------------
    
    for igauge in range(numgauges):
        gaugeno = gaugenos[igauge]
        for figno in fignos:
            html = open(gauge_htmlfile[gaugeno,figno],'w')
            html.write('<html><meta http-equiv="expires" content="0">\n')
            html.write('<title>Plots</title>')
            html.write('<body><center>\n')
            html.write('\n<h3>Gauge %i ' % gaugeno)
            if numfigs > 1:
                html.write(' &nbsp;---&nbsp; %s' % fignames[figno] )
        
            # Write link commands to previous and next gauge:

            html.write('<p> <a href="%s">' % gauge_htmlfile[gaugenos[0],figno])
            html.write('&#060; &#060;</a> &nbsp; &nbsp;\n')
            if igauge==0:
                html.write('&#060; &nbsp; &nbsp; ')
                html.write('\n<a href="_PlotIndex.html">Index</a>  ')
                if numgauges > 1:
                    html.write('&nbsp; &nbsp; <a href="%s"> &#062; </a> ' \
                        % gauge_htmlfile[gaugenos[1],figno])

            elif igauge==numgauges-1:
                if numgauges > 1:
                    html.write('\n<a href="%s"> &#060; </a>  &nbsp; &nbsp; ' \
                        % gauge_htmlfile[gaugenos[igauge-1],figno])
                html.write('\n<a href="_PlotIndex.html">Index</a>  ')
                html.write(' &nbsp; &nbsp; &#062; ')

            else:
                html.write('\n<a href="%s"> &#060; </a>  &nbsp; &nbsp; ' \
                        % gauge_htmlfile[gaugenos[igauge-1],figno])
                html.write('\n<a href="_PlotIndex.html">Index</a>  ')
                html.write('\n &nbsp; &nbsp;<a href="%s"> &#062; </a>  &nbsp; &nbsp; ' \
                        % gauge_htmlfile[gaugenos[igauge+1],figno])

            html.write('&nbsp; &nbsp; \n<a href="%s"> ' \
                      % gauge_htmlfile[gaugenos[numgauges-1],figno])
            html.write('&#062; &#062;</a>  \n') 
        
            # image:
            html.write('\n\n <p><img src="%s"><p>  \n ' \
                        % gauge_pngfile[gaugeno,figno])

            html.write('\n\nImage source: &nbsp; %s'  \
                   % os.path.join(os.getcwd(),gauge_pngfile[gaugeno,figno]))

            # list of all figures at bottom of page:

            if numfigs > 1:
                html.write('\n<p><b>Other figures at this time:</b> &nbsp;&nbsp;')
                for figno2 in fignos:
                    if figno2 == figno:
                        html.write('\n<font color=red>%s</font>&nbsp;&nbsp;' \
                               % fignames[figno])
                    else:
                        html.write('\n<a href="%s">%s</a>  &nbsp; &nbsp; ' \
                           % (gauge_htmlfile[gaugeno,figno2],fignames[figno2]))
                html.write('\n<a href="%s"> All Figures </a>' \
                     % gauge_allfigsfile[gaugeno])

            # list of all gauges at bottom of page:

            html.write('\n<p><b>Other gauges:</b></a> &nbsp;&nbsp;')
            for gaugeno2 in gaugenos:
                if gaugeno2 == gaugeno:
                    html.write('\n<font color=red>%i</font>&nbsp;&nbsp;' \
                               % gaugeno)
                else:
                    html.write('\n<a href="%s">%i</a>  &nbsp; &nbsp; ' \
                           % (gauge_htmlfile[gaugeno2,figno],gaugeno2))
            html.write('\n<a href="allgaugesfig%s.html">  All Gauges </a>' \
                     % figno)
        
            html.write('\n<p><h3><a href=%s>Plot Index</a></h3>' \
                      % (plotdata.html_index_fname))
            if eagle:
                html.write("""<p><h3><a href="../eaglemenu.html">Main Menu for
                this run-directory</a></h3>  """)
            html.write('</center></body></html>')
            html.close()

    os.chdir(startdir)
    # end of plotclaw2html
    

#=====================================
def massage_gauges_data(plot_pages_data):
#=====================================
    ppd = plot_pages_data

    try:
        gaugenos = ppd.gauges_gaugenos
        fignos = ppd.gauges_fignos
        fignames = ppd.gauges_fignames
        prefix = getattr(ppd, 'gauges_prefix', 'gauge')
    except:
        print '*** Error: gauges not set properly'
        return

    startdir = os.getcwd()
        

    for figno in fignos:
        if not fignames.has_key(figno):
            fignames[figno] = 'Solution'
            
    numgauges = len(gaugenos)
    numfigs = len(fignos)

    #if len(gaugenos) == 0:
    #    print '*** Warning: 0 gauges to print'
    #if len(fignos) == 0:
    #    print '*** Warning: 0 figures to print each gauge'

    pngfile = {}
    htmlfile = {}
    allfigsfile = {}
    for gaugeno in gaugenos:
        gaugef = string.zfill(gaugeno,4)
        for figno in fignos:
            pngfile[gaugeno,figno] = '%s%sfig%s.png'  % (prefix,gaugef,figno)
            htmlfile[gaugeno,figno] = '%s%sfig%s.html' % (prefix,gaugef,figno)
        allfigsfile[gaugeno] = 'allfigs%s%s.html' % (prefix,gaugef)

    ppd.gauges_gaugenos = gaugenos
    ppd.gauges_fignos = fignos
    ppd.gauges_fignames = fignames
    ppd._gauge_pngfile = pngfile
    ppd._gauge_htmlfile = htmlfile
    ppd._gauge_allfigsfile = allfigsfile
    return ppd


#============================================
def plotclaw_driver(plotdata, verbose=False):
#============================================
    """
    The ClawPlotData object plotdata will be initialized by a call to
    function setplot unless plotdata.setplot=False.  

    If plotdata.setplot=True then it is assumed that the current directory
    contains a module setplot.py that defines this function.

    If plotdata.setplot is a string then it is assumed this is the name of
    a module to import that contains the function setplot.

    If plotdata.setplot is a function then this function will be used.
    """

    import glob, sys, os
    import numpy as np
    from pyclaw.plotters.data import ClawPlotData
    from pyclaw.plotters import frametools, gaugetools, plotpages

    datadir = os.getcwd()  # assume data files in this directory

    if not sys.modules.has_key('matplotlib'):
        print '*** Error: matplotlib not found, no plots will be done'
        return plotdata
        
    if not isinstance(plotdata,ClawPlotData):
        print '*** Error, plotdata must be an object of type ClawPlotData'
        return plotdata

    plotdata._mode = 'printframes'

    plotdata = frametools.call_setplot(plotdata.setplot, plotdata)

    try:
        plotdata.rundir = os.path.abspath(plotdata.rundir)
        plotdata.outdir = os.path.abspath(plotdata.outdir)
        plotdata.plotdir = os.path.abspath(plotdata.plotdir)

        framenos = plotdata.print_framenos # frames to plot
        gaugenos = plotdata.print_gaugenos # gauges to plot
        fignos = plotdata.print_fignos     # figures to plot at each frame
        fignames = {}                      # names to use in html files

        rundir = plotdata.rundir       # directory containing *.data files
        outdir = plotdata.outdir       # directory containing fort.* files
        plotdir = plotdata.plotdir     # where to put png and html files
        overwrite = plotdata.overwrite # ok to overwrite?
        msgfile = plotdata.msgfile     # where to write error messages
    except:
        print '*** Error in printframes: plotdata missing attribute'
        print '  *** plotdata = ',plotdata
        return plotdata

    if fignos == 'all':
        fignos = plotdata._fignos
        #for (figname,plotfigure) in plotdata.plotfigure_dict.iteritems():
        #    fignos.append(plotfigure.figno)


    # filter out the fignos that will be empty, i.e.  plotfigure._show=False.
    plotdata = frametools.set_show(plotdata)
    fignos_to_show = []
    for figname in plotdata._fignames:
        figno = plotdata.plotfigure_dict[figname].figno
        if (figno in fignos) and plotdata.plotfigure_dict[figname]._show:
            fignos_to_show.append(figno)
    fignos = fignos_to_show
        
    # figure out what type each figure is:
    fignos_each_frame = []
    fignos_each_gauge = []
    fignos_each_run = []
    for figno in fignos:
        figname = plotdata._figname_from_num[figno]
        if plotdata.plotfigure_dict[figname].type == 'each_frame':
            fignos_each_frame.append(figno)
        if plotdata.plotfigure_dict[figname].type == 'each_gauge':
            fignos_each_gauge.append(figno)
        if plotdata.plotfigure_dict[figname].type == 'each_run':
            fignos_each_run.append(figno)
        

    rootdir = os.getcwd()

    # annoying fix needed when EPD is used for plotting under cygwin:
    if rootdir[0:9] == 'C:\cygwin' and outdir[0:9] != 'C:\cygwin':
        outdir = 'C:\cygwin' + outdir
        plotdata.outdir = outdir
    if rootdir[0:9] == 'C:\cygwin' and rundir[0:9] != 'C:\cygwin':
        rundir = 'C:\cygwin' + rundir
        plotdata.rundir = rundir
    if rootdir[0:9] == 'C:\cygwin' and plotdir[0:9] != 'C:\cygwin':
        plotdir = 'C:\cygwin' + plotdir
        plotdata.plotdir = plotdir

    try:
        os.chdir(rundir)
    except:
        print '*** Error: cannot move to run directory ',rundir
        print 'rootdir = ',rootdir
        return plotdata


    if msgfile != '':
        sys.stdout = open(msgfile, 'w')
        sys.stderr = sys.stdout


    try:
        plotpages.cd_plotdir(plotdata)
    except:
        print "*** Error, aborting plotframes"
        return plotdata


    framefiles = glob.glob(os.path.join(plotdir,'frame*.png')) + \
                    glob.glob(os.path.join(plotdir,'frame*.html'))
    if overwrite:
        # remove any old versions:
        for file in framefiles:
            os.remove(file)
    else:
        if len(framefiles) > 1:
            print "*** Remove frame*.png and frame*.html and try again,"
            print "  or use overwrite=True in call to printframes"
            return plotdata


    try:
        os.chdir(outdir)
    except:
        print '*** Error plotclaw_driver: cannot move to outdir = ',outdir
        return plotdata


    fortfile = {}
    pngfile = {}
    frametimes = {}

    for file in glob.glob('fort.q*'):
        frameno = int(file[7:10])
        fortfile[frameno] = file
        for figno in fignos_each_frame:
            pngfile[frameno,figno] = 'frame' + file[-4:] + 'fig%s.png' % figno
    
    if len(fortfile) == 0:
        print '*** No fort.q files found in directory ', os.getcwd()
        return plotdata
    
    # Discard frames that are not from latest run, based on
    # file modification time:
    framenos = frametools.only_most_recent(framenos, plotdata.outdir)

    numframes = len(framenos)

    print "Will plot %i frames numbered:" % numframes, framenos
    print 'Will make %i figure(s) for each frame, numbered: ' \
          % len(fignos_each_frame), fignos_each_frame

    #fignames = {}
    #for figname in plotdata._fignames:
        #figno = plotdata.plotfigure_dict[figname].figno
        #fignames[figno] = figname
    # use new attribute:
    fignames = plotdata._figname_from_num

    # Only grab times by loading in time
    for frameno in framenos:
        frametimes[frameno] = plotdata.gettime(frameno, plotdata.outdir)
    # for frameno in framenos:
    #     frametimes[frameno] = plotdata.getframe(frameno, plotdata.outdir).t

    plotdata.timeframes_framenos = framenos
    plotdata.timeframes_frametimes = frametimes
    plotdata.timeframes_fignos = fignos_each_frame
    plotdata.timeframes_fignames = fignames

    # Gauges:
    # -------
    gaugenos = plotdata.print_gaugenos
    if gaugenos == 'all':
        # Read gauge numbers from setgauges.data if it exists:
        setgauges = gaugetools.read_setgauges(datadir)
        gaugenos = setgauges.gaugenos

    plotdata.gauges_gaugenos = gaugenos
    plotdata.gauges_fignos = fignos_each_gauge
    plotdata.gauges_fignames = fignames


    # Make html files for time frame figures:
    # ---------------------------------------

    os.chdir(plotdir)

    if plotdata.html:
        #plotpages.timeframes2html(plotdata)
        plotpages.plotclaw2html(plotdata)
        pass
    
    # Make png files for all frames and gauges:
    # -----------------------------------------

    if not plotdata.printfigs:
        print "Using previously printed figure files"
    else:
        print "Now making png files for all figures..."

        for frameno in framenos:
            frametools.plotframe(frameno, plotdata, verbose)
            print 'Frame %i at time t = %s' % (frameno, frametimes[frameno])

        for gaugeno in gaugenos:
            gaugetools.plotgauge(gaugeno, plotdata, verbose)
            print 'Gauge %i ' % gaugeno


    if plotdata.latex:
        plotpages.timeframes2latex(plotdata)
    

    # Movie:
    #-------
    
    if plotdata.gif_movie:
        print 'Making gif movies.  This may take some time....'
        for figno in fignos_each_frame:
            try:
                os.system('convert -delay 20 frame*fig%s.png moviefig%s.gif' \
                   % (figno,figno))
                print '    Created moviefig%s.gif' % figno
            except:
                print '*** Error creating moviefig%s.gif' % figno
    
    os.chdir(rootdir)

    # print out pointers to html index page:
    path_to_html_index = os.path.join(os.path.abspath(plotdata.plotdir), \
                               plotdata.html_index_fname)
    plotpages.print_html_pointers(path_to_html_index)

    # reset stdout for future print statements
    sys.stdout = sys.__stdout__

    return plotdata
    # end of printframes

