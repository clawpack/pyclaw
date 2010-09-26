#!/usr/bin/env python 

# convert CLAWPACK source code, data files, READMEs, etc. into html.
# for example,
#   clawcode2html filename.f
# generates
#   filename.f.html
# If filename.f.html already exists, you will be prompted before
# overwriting unless the -f or --force flag is used, e.g.
#   clawcode2html --force filename.f
#
# Based on mathcode2html.py from 
#    http://www.amath.washington.edu/~rjl/mathcode2html
# with some modifications for CLAWPACK.
#
# The environment variable CLAW must be properly set to the root
# of the claw directory before using this script.
# 
# Most of code is put into <pre> environment.
# Any comments enclosed by begin_html and end_html are taken out 
# of the <pre> environment and indented properly (using a table) 
# to match the surrounding code.

# By default, the html comments are offset in blue font, except if the 
# input file has extension .txt or no extension, in which case it's 
# left as black.  The default_color can be changed below.
# Also, the begin_html line may contain [color: red] for example,
# to determine the color of this comment.

# You may modify the recognized extensions and comment characters below.

# Allows many latex commands in comments if jsMath is used on webpage.
# In this case make sure the variable jsMathScript is properly set below.
# If you don't want to include a call to the jsMath script on the html
# page, invoke mathcode2html with the --nojsmath option.


#---------------------------------------------------------------------------
# Copyright (C) 2008   Randall J. LeVeque 
#
# Distributed as part of Clawpack, under the BSD license
# See www.clawpack.org
#---------------------------------------------------------------------------


import sys,os,glob
import string,re
import time
try:
    from optparse import OptionParser
except:
    print 'You must use a more recent version of Python'
    sys.exit(1)

# parse command line arguments:
parser = OptionParser()
parser.add_option("-f", "--force",
                  action="store_true", dest="forced", default=False,
                  help="force action even if it overwrites html file")
parser.add_option("--nojsmath",
                  action="store_false", dest="nojsMath", default=False,
                  help="don't include call to jsMath script in html")
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")
parser.add_option("--noheader",
                  action="store_false", dest="header", default=True,
                  help="suppress printing header at top of html page")
parser.add_option("--dropext",
                  action="store_true", dest="dropext", default=False,
                  help="drop the file extension before adding .html for output")
parser.add_option("--eagle",
                  action="store_true", dest="eagle", default=False,
                  help="load the eagleclaw.css style")

(options, infiles) = parser.parse_args()

forced = options.forced
nojsMath = options.nojsMath
verbose = options.verbose
header = options.header
dropext = options.dropext 
eagle = options.eagle 



# CLAWPACK environment variables:

try:
    clawdir = os.getenv('CLAW')  # assumes environment variable set properly
except:
    print "CLAW environment variable not set"
    print "You need to run setenv.py before setup.py"
    sys.exit(1)

if clawdir == None:
    print "CLAW environment variable not set"
    print "You need to run setenv.py before setup.py"
    sys.exit(1)


# Create addresses for links on html pages.
# change the lines below if you want to point to a webpage home instead
# of the local file system  (e.g.  clawaddr = "http://www.mywebpage/claw")

clawaddr = 'http://localhost:50005'
#clawaddr = 'http://kingkong.amath.washington.edu/claw4'


# Set comment characters for different programming languages:
# Augment or modify as desired.
commentchar = { '.f'  : ['!', '#'], \
                '.f95'  : ['!','#'], \
                '.m'  : '%', \
                '.py' : '#', \
                '.data' : ['#','=:'], \
                '.txt': None, \
                'Makefile': '#', \
                ''    : None}

firstfort = ['c','*','C','!']   # valid fortran .f comment char's in col. 1
firstfort95 = ['!']   # valid fortran .f95 comment char's in col. 1

leadingindent = ''    # additional indentation for webpage, if desired

default_color = 'blue'


try:
    infiles.remove('clawcode2html.py')  # can't apply this code to itself!
except:
    pass


for infilename in infiles:

    # check if this is a code file of a recognized language.
    ext = os.path.splitext(infilename)[1]
    if infilename == 'Makefile':
        ext = 'Makefile'  # special case
    
    if not commentchar.has_key(ext):
        print "  "
        print "  Warning: Unrecognized extension, will proceed"
        print "           with no replacement of comment characters"
        commentchar[ext] = None 
    
    
    # open input and output files:
    #-----------------------------
    
    try:
        ifile = open(infilename,'r')
    except:
        print "File not found:", infilename
        sys.exit(1)

    # Search for [use: ...] statements in file:
    all_lines = ifile.read()
    regexp = re.compile(r"\[use:[^\]]*jsMath")
    result = regexp.search(all_lines)

    # use jsMath if [use:jsMath] found and option --nojsmath was
    # not used in call:
    usejsMath = (result is not None) & (not nojsMath)

    regexp = re.compile(r"\[use:(?P<cssfile>[^\]]*).css\]")
    result = regexp.search(all_lines)
    if result is not None:
        cssfile = result.group('cssfile').strip() + '.css'
    else:
        cssfile = None

    
    ifile.seek(0)  # return to start of file
    lines = ifile.readlines()
    
    if dropext:
        infileroot = os.path.splitext(infilename)[0]
        outfilename = infileroot + '.html'
    else:
        outfilename = infilename + '.html'

    if (glob.glob(outfilename) != []) & (not forced):
        sys.stdout.write('  OK to overwrite %s?  '  %  outfilename)
        answer = raw_input()
        if answer not in ['y','Y','yes']:
            print '  Aborting!'
            sys.exit(1)
    
    ofile = open(outfilename,'w')
    
    
    # start creating html file:
    #--------------------------
    
    if verbose:
        print '  Converting ', infilename, ' to ', outfilename

    
    ofile.write("""<!-- **** DO NOT EDIT THIS FILE *** 
         This file was generated automatically from file
              %s
         using $CLAW/doc/clawcode2html.py         -->""" % infilename)

    ofile.write('\n\n<html>\n<title> %s </title>\n\n'  % outfilename)
    if eagle:
        ofile.write("""
          <head>
          <link type="text/css" rel="stylesheet"
          href="%s/eagleclaw/eagleclaw.css">
          </head> 
          <body>
        """ % clawaddr)
    elif cssfile:
        ofile.write("""
          <head>
          <link type="text/css" rel="stylesheet"
          href="%s/%s">
          </head> 
          <body>
        """ % (clawaddr, cssfile))
    else:
        ofile.write(
           """
           <BODY TEXT="#000000" BGCOLOR="#FFFFFF" LINK="#0000FF" 
                 VLINK="#5500DD" ALINK="#FF0000">
           <font FACE="HELVETICA,ARIAL">
           """)

    
    # determine time and reformat:
    time1 = time.asctime()
    year = time1[-5:]
    day = time1[:-14]
    hour = time1[-13:-5]
    creationtime = day + year + ' at ' + hour

    # put full file name with path into a comment for future reference:
    fullinfilename = os.path.join(os.getcwd(),infilename)
    ofile.write('\n<!-- Created from the file %s -->\n'  % fullinfilename)
    ofile.write('<!-- Date: %s -->\n\n'  % creationtime)
    
    
    if header:
        #ofile.write('<table bgcolor="#DDEEEE"> <tr> <td>\n')
        ofile.write('<table bgcolor="#FFEE99"> <tr> <td>\n')
#       ofile.write('&nbsp;<font size=6> %s </font> </td>\n' % outfilename)
    
#       ofile.write("""<td>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
#           <a href="http://www.clawpack.org"><img
#           src="%s/doc/images/clawlogo.jpg"
#           width=100 alt="CLAWPACK"></a>""" % clawaddr)
        ofile.write("""<tr> <td><a href="http://www.clawpack.org"><img
            src="%s/doc/images/clawlogo.jpg"
            width=100 alt="CLAWPACK"></a>
            &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
            &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;""" % clawaddr)
        ofile.write("""<font size=6> %s </font>
            </td><td>&nbsp;&nbsp;&nbsp;&nbsp; """ % outfilename)

        ofile.write('&nbsp;&nbsp; </td> </tr> <tr> <td>\n')
    
        ofile.write('&nbsp;Source file: &nbsp;&nbsp;<a href="%s">%s </a>\n' \
                     % (infilename,infilename))
        ofile.write('</td></tr><tr><td>\n')
        ofile.write('&nbsp;Directory:  &nbsp; %s \n' %  os.getcwd())
        ofile.write('</td></tr><tr><td>\n')
        ofile.write('&nbsp;Converted:  &nbsp; %s \n' % creationtime)
        ofile.write('&nbsp; using <a href="%s/doc/clawcode2html.html">clawcode2html</a>\n'\
                    % clawaddr)
        ofile.write('</td></tr><tr><td>\n')
        ofile.write('<font color="#BB3300"> &nbsp;This documentation file will \n')
        ofile.write('not reflect any later changes in the source file. </font>\n')
        ofile.write('</td></tr></table></font>\n')
        ofile.write('<p>\n')
    
    
    if usejsMath:
        
        # Set jsMathScript to the file or URL of the java script load.js.

        # proper location for using webserver started by executing
        # "python startserver.py" in the $CLAW directory:
        jsMathScript = "%s/doc/load.js"  % clawaddr

        # script location for posting on kingkong server:
        # jsMathScript = "http://kingkong.amath.washington.edu/claw/doc/load.js"
    
        # script location for posting on rjl's webpage:
        # jsMathScript = "http://www.amath.washington.edu/~rjl/jsMath/easy/load.js"
    
        ofile.write("    <!-- jsMath stuff:  -->\n")
        ofile.write(" <SCRIPT SRC= %s> </SCRIPT>\n"  % jsMathScript)

    
        ofile.write("""
            <!-- Warning message if script doesn't work:  -->
            $\phantom{******** If you see this on the webpage then the 
                               browser could not locate *********}$<br>
            $\phantom{******** the jsMath file load.js *********}$<p>
            \n""")
        
        # define any latex macros that you want to use in jsMath: 
        
        ofile.write("""
          <!- latex macros: -->
          $\\newcommand{\\vector}[1]{\\left[\\begin{array}{c} #1 \\end{array}\\right]}$ 
          $\\newenvironment{matrix}{\\left[\\begin{array}{cccccccccc}} {\\end{array}\\right]}$ 
          $\\newcommand{\\A}{{\\cal A}}$
          $\\newcommand{\\W}{{\\cal W}}$
          \n""")
        #
        # end of jsMath stuff
    

    # start writing input file to output file...
    
    # code is in pre-formatted environment:
    ofile.write('<pre> \n')
    
    insidehtml = 0;   # set to 1 when we're processing html comments
    lineno = 0;       # line number counter for error message
    
    for line in lines:
    
        lineno += 1
    
        if string.count(line,"begin_html"):
            regexp = re.compile(r"\[color:(?P<color>[^\]]*)\]")
            result = regexp.search(line)
	    if result:
	        font_color = result.group('color')
	    else:
	        font_color = default_color
    
            if insidehtml:
                print '  Error at line ', lineno, '\n'
                print '  Unexpected begin_html  when already in html mode\n'
                print '  Missing end_html? \n'
                sys.exit(1)
    
            # switch out of pre-formatted mode and create table to indent
            ofile.write('</pre>\n<table><tr><td>\n')
    
            # The first column of the table is spaces for indentation 
            # to match surrounding source code.
            # Count how many spaces there are before the begin_html
            #  or to the first comment character preceeding that string:
            regexp = re.compile(r"(?P<spaces>[ ]*)(#|%|begin)")
            result = regexp.search(line)
            if result:
                numindent = len(leadingindent) + len(result.group('spaces'))
                numindent = numindent
                ofile.write('<pre><tt>')
                for i in range(numindent):
                    ofile.write(' ')
                ofile.write('</tt></pre>')
            else:
                print '  Strange error in clawcode2html - should not be here'
                print '     at line number ', lineno
              
    
            ofile.write('\n</td>\n')
    
            # the next column of the table has the comment itself:
            ofile.write('<td width=900>\n')
            if ext not in ('.txt', ''):
                # use colored font for comments except for text files.
                ofile.write('<font size="-1" color=%s>\n'  % font_color)
            insidehtml = 1;
    
        elif string.count(line,"end_html"):
            # switch back to pre-formatted environment
            ofile.write('</font></td></tr></table>\n')
            ofile.write('<pre> \n')
            insidehtml = 0;
    
        else:
            if insidehtml:
    
                # replace blank line in html comment by new paragraph <p>:
                blankline = (string.split(line) == [])
                if not blankline:
                   firstchar = string.split(line)[0][0]
                   if ((ext in ['.f','.f95']) & (firstchar not in firstfort)): 
                       if firstchar not in commentchar[ext]:
                           print '  Error... in line ', lineno,'\n' \
                                 '  In html but not in a comment.'\
                                 '   Forgotten "end_html" ?'
                           sys.exit(1)
    
                # fortran comments may have comment symbol in column 1
                # strip lines inside an html comment of leading symbol:
                if ext == '.f':
                    if line[0] in firstfort:
                       line = ' ' + line[1:]
                if ext == '.f95':
                    if line[0] in firstfort95:
                       line = ' ' + line[1:]
    
                # replace any comment character by ' ' 
                if commentchar[ext]:
                    for char in commentchar[ext]:
                        line = string.replace(line,char,' ')  
    
    
                blankline = (string.split(line) == [])
                if blankline:
                   line = ('<p>\n')


                # Allow wiki formatting of links:
                # -------------------------------

		# Replace [name: placemark] by <a name="placemark">
                # (to jump to a different spot in the same html file)
                regexp = re.compile(r"\[name:[ ]*(?P<placemark>[^ ^\]]*)\]")
                result = regexp.search(line)
                while result:
                    placemark = result.group('placemark')
                    oldpat = result.group()
		    newpat = '<a name="%s">' % placemark
                    line = line.replace(oldpat,newpat)
                    result = regexp.search(line)

                # Replace links of the form [code: target]
                # by html links to both target and target.html.
		# Also allows [code: target#placemark] with links to target 
		# and target.html#placemark.

                regexp = re.compile(r"\[code:[ ]*(?P<target>[^ ^\]^#]*)([#]?)" + \
                                    r"(?P<placemark>[^\]]*)\]")
                result = regexp.search(line)
                while result:
                    targetname = result.group('target')
                    placemark = result.group('placemark')
                    oldpat = result.group()
		    if placemark:
                        newpat = '<a href="%s">%s</a>' \
			             % (targetname,targetname) + \
                                 '&nbsp;<a href="%s.html#%s">[.html]</a>' \
			             % (targetname,placemark) 
		    else:
                        oldpat = result.group()
                        newpat = '<a href="%s">%s</a>' \
			             % (targetname,targetname) + \
                                 '&nbsp;<a href="%s.html">[.html]</a>' \
			             % targetname 
                    line = line.replace(oldpat,newpat)
                    result = regexp.search(line)


                # replace links of the form [link: target text]
                # by an html link from text to the target page.
                regexp = re.compile(r"\[link:[ ]?(?P<target>[^ ^\]]*)" + \
                                    r"([ ]*)(?P<text>[^\]]*)\]")
                result = regexp.search(line)
                while result:
                    targetname = result.group('target')
                    text = result.group('text')
                    oldpat = result.group()
                    if text=='': 
                        text = targetname
                    newpat = '<a href="' + targetname + '">' + \
                                text + '</a>  '
                    line = line.replace(oldpat,newpat)
                    result = regexp.search(line)


                # replace links of the form [http:etc text]
                # by an html link from text to the http page.
                regexp = re.compile(r"\[http:(?P<target>[^ ^\]]*)" + \
                                    r"([ ]*)(?P<text>[^\]]*)\]")
                result = regexp.search(line)
                while result:
                    targetname = result.group('target')
                    text = result.group('text')
                    oldpat = result.group()
                    if text=='': 
                        text = targetname[2:]
                    newpat = '<a href="http:' + targetname + '">' + \
                                text + '</a>  '
                    line = line.replace(oldpat,newpat)
                    result = regexp.search(line)


                # replace links of the form [www.etc text]
                # by an html link from text to the http://www.etc page.
                regexp = re.compile(r"\[www.(?P<target>[^ ^\]]*)" + \
                                    r"([ ]*)(?P<text>[^\]]*)\]")
                result = regexp.search(line)
                while result:
                    targetname = result.group('target')
                    text = result.group('text')
                    oldpat = result.group()
                    if text=='': 
                        text = 'www.' + targetname
                    newpat = '<a href="http://www.' + targetname + '">' + \
                                text + '</a>  '
                    line = line.replace(oldpat,newpat)
                    result = regexp.search(line)


		# special things for CLAWPACK:

                # replace links of the form [clawcode:clawpack/1d/lib/step1.f]
                # for example by links relative to clawaddr, 
		# along with a link to the .html version.
                regexp = re.compile(r"\[clawcode:[ ]*(?P<target>[^ ^\]^#]*)([#]?)" + \
                                    r"(?P<placemark>[^\]]*)\]")
                result = regexp.search(line)
                while result:
                    targetname = result.group('target').lstrip()
                    placemark = result.group('placemark')
                    oldpat = result.group()
		    if placemark:
                        newpat = '<a href="%s/%s">claw/%s</a>' \
			             % (clawaddr,targetname,targetname) + \
                                 '&nbsp;<a href="%s/%s.html#%s">[.html]</a>' \
			             % (clawaddr,targetname,placemark) 
		    else:
                        newpat = '<a href="%s/%s">%s</a>' \
			             % (clawaddr,targetname,targetname) + \
                                 '&nbsp;<a href="%s/%s.html">[.html]</a>' \
			             % (clawaddr,targetname) 
                    line = line.replace(oldpat,newpat)
                    result = regexp.search(line)


                # replace links of the form [claw:clawpack/1d/lib]
                # for example by links relative to clawaddr,
		# with no .html version.
                regexp = re.compile(r"\[claw:[ ]?(?P<target>[^ ^\]]*)" + \
                                    r"([ ]*)(?P<text>[^\]]*)\]")
                result = regexp.search(line)
                while result:
                    targetname = result.group('target')
                    text = result.group('text')
                    oldpat = result.group()
                    if text=='': 
                        text = '$CLAW/' + targetname
                    newpat = '<a href="' + clawaddr + '/' + targetname + \
                               '">' + text + '</a>  ' 
                    line = line.replace(oldpat,newpat)
                    result = regexp.search(line)

		# place text surrounded by triple braces with 
		# pre environment with background color:
		newpat = '<pre class="clawcode">'
                line = line.replace('{{{',newpat)
                line = line.replace('}}}','</pre>')

          
	    else:
                # not insidehtml - make regular comments default_color.
                # Determine if this line contains a comment and if so,
                # what column the comment starts in:

		startcomment = 1000
	        if (ext == '.f') & (line[0] in firstfort):
	            startcomment = 0
	        elif (ext == '.f95') & (line[0] in firstfort95):
	            startcomment = 0
		else:
		    if commentchar[ext]:
                        for c in commentchar[ext]:
		            commentcol = string.find(line,c)
		            if (commentcol>-1)&(commentcol<startcomment):
                                startcomment = commentcol

		if startcomment<1000:
		    line = line[0:startcomment] + \
		           '<font color="%s">'  % default_color + \
		           line[startcomment:-1] + '</font>\n'


            # output the (possibly modified) line to the output file:
            ofile.write('%s' % leadingindent+line)
    
    # Done with all lines.  Add closing stuff at bottom of html file:
    ofile.write('</pre></html>\n')
    
    ifile.close()
    ofile.close()
    
