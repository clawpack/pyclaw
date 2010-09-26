
# Makefile for making htmls from programs, data, and .txt files

# Make list of all files that should be converted to .html:
HTML = \
  Makefile.html \
  $(subst .f,.f.html,$(wildcard *.f)) \
  $(subst .f95,.f95.html,$(wildcard *.f95)) \
  $(subst .py,.py.html,$(wildcard *.py)) \
  $(subst .data,.data.html,$(wildcard *.data)) \
  $(subst .txt,.html,$(wildcard *.txt))


# Rules to make html files:
CC2HTML = python $(CLAW)/doc/clawcode2html.py --force  
# e.g. qinit.f --> qinit.f.html
%.f.html : %.f ; $(CC2HTML) $<
%.f95.html : %.f95 ; $(CC2HTML) $<
%.m.html : %.m ; $(CC2HTML) $<
%.py.html : %.py ; $(CC2HTML) $<
%.data.html : %.data ; $(CC2HTML) $<

# drop .txt extension, e.g. README.txt --> README.html
%.html : %.txt ; $(CC2HTML) --dropext --noheader  $<
Makefile.html: Makefile ; $(CC2HTML) $<

htmls: $(HTML);

clean: 
	-rm -f $(HTML)
