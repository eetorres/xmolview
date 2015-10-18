#========================================================================
# Copyrigth 2006-2015 by Edmanuel Torres
# eetorres@gmail.com
#========================================================================

include makeinclude

DIRS = src bin

all:   makeinclude
#	if [ -f bin/$(PRGS) ] ; then
#	  rm bin/$(PRGS)
#	fi
	@for dir in $(DIRS); do\
	    echo "==== Compiling in $$dir ====";\
	    (cd $$dir; make || break);\
	done
	
clean: makeinclude
	@for dir in $(DIRS); do\
	    echo "==== Cleaning in $$dir ====";\
	    (cd $$dir; make clean || break);\
	done

install:
	#rm -f ~/bin/$(PRGS)
	cp bin/$(PRGS) ~/bin

run:
	./bin/$(PRGS)
