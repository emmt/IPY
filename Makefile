# where are the sources? (automatically filled in by configure script)
srcdir=.

# these values filled in by "yorick -batch make.i" or configure script
Y_MAKEDIR=
Y_EXE=
Y_EXE_PKGS=
Y_EXE_HOME=
Y_EXE_SITE=
Y_HOME_PKG=

# ----------------------------------------------------- optimization flags

# options for make command line, e.g.-   make COPT=-g TGT=exe
COPT=$(COPT_DEFAULT)
TGT=$(DEFAULT_TGT)

# ------------------------------------------------ macros for this package

PKG_NAME=ipy
PKG_I=${srcdir}/ipy.i

OBJS=ipy.o

# change to give the executable a name other than yorick
PKG_EXENAME=yorick

# PKG_DEPLIBS=-Lsomedir -lsomelib   for dependencies of this package
PKG_DEPLIBS=
# set compiler (or rarely loader) flags specific to this package
PKG_CFLAGS=
PKG_LDFLAGS=

# list of additional package names you want in PKG_EXENAME
# (typically $(Y_EXE_PKGS) should be first here)
EXTRA_PKGS=$(Y_EXE_PKGS)

# list of additional files for clean
PKG_CLEAN=

# autoload file for this package, if any
PKG_I_START=${srcdir}/ipy-autoload.i
# non-pkg.i include files for this package, if any
PKG_I_EXTRA=${srcdir}/ipy-tests.i ${srcdir}/ipy-cost.i ${srcdir}/linop.i ${srcdir}/rgl.i

RELEASE_FILES = AUTHORS LICENSE.md Makefile NEWS README.md TODO \
	configure linop.i ipy-cost.i ipy.i ipy.c \
	ipy-autoload.i ipy-tests.i
RELEASE_NAME = $(PKG_NAME)-$(RELEASE_VERSION).tar.bz2

# -------------------------------- standard targets and rules (in Makepkg)

# set macros Makepkg uses in target and dependency names
# DLL_TARGETS, LIB_TARGETS, EXE_TARGETS
# are any additional targets (defined below) prerequisite to
# the plugin library, archive library, and executable, respectively
PKG_I_DEPS=$(PKG_I)
Y_DISTMAKE=distmake

ifeq (,$(strip $(Y_MAKEDIR)))
$(info *** WARNING: Y_MAKEDIR not defined, you may run 'yorick -batch make.i' first)
else
include $(Y_MAKEDIR)/Make.cfg
include $(Y_MAKEDIR)/Makepkg
include $(Y_MAKEDIR)/Make$(TGT)
endif

# override macros Makepkg sets for rules and other macros
# see comments in Y_HOME/Makepkg for a list of possibilities

# if this package built with mpy: 1. be sure mpy appears in EXTRA_PKGS,
# 2. set TGT=exe, and 3. uncomment following two lines
# Y_MAIN_O=$(Y_LIBEXE)/mpymain.o
# include $(Y_MAKEDIR)/Makempy

# configure script for this package may produce make macros:
# include output-makefile-from-package-configure

# reduce chance of yorick-1.5 corrupting this Makefile
MAKE_TEMPLATE = protect-against-1.5

# ------------------------------------- targets and rules for this package

# Dummy default target in case Y_MAKEDIR was not defined:
dummy-default:
	@echo >&2 "*** ERROR: Y_MAKEDIR not defined, aborting..."; false

%.o: ${srcdir}/%.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<


# simple example:
#myfunc.o: myapi.h
# more complex example (also consider using PKG_CFLAGS above):
#myfunc.o: myapi.h myfunc.c
#	$(CC) $(CPPFLAGS) $(CFLAGS) -DMY_SWITCH -o $@ -c myfunc.c

start: ${srcdir}/ipy.i ${srcdir}/ipy-cost.i
	for x in $^; do \
	  y=`basename "$$x"`; \
	  sed -E "/^ *(func|extern) +[A-Za-z][_0-9A-Za-z]*/!d;s/^ *(func|extern) +([A-Za-z][_0-9A-Za-z]*).*/autoload, \"$$y\", \2;/" < "$$x"; \
	done

tag:
	@version=`sed -e '/^ *IPY_VERSION *= *"[0-9.]*" *; *$$/!d;s/^ *IPY_VERSION *= *"\([0-9.]*\)" *; *$$/\1/' <ipy.i`; \
	echo "To tag this version, type:"; \
	echo "   git tag -m 'Release $$version of IPY' v$$version"; \
	echo "   git push"; \
	echo "   git push --tags"

release: $(RELEASE_NAME)

$(RELEASE_NAME):
	@if test "x$(RELEASE_VERSION)" = "x"; then \
	  echo >&2 "set package version: make RELEASE_VERSION=... release"; \
	else \
          dir=`basename "$(RELEASE_NAME)" .tar.bz2`; \
	  if test "x$$dir" = "x" -o "x$$dir" = "x."; then \
	    echo >&2 "bad directory name for archive"; \
	  elif test -d "$$dir"; then \
	    echo >&2 "directory $$dir already exists"; \
	  else \
	    mkdir -p "$$dir"; \
	    for src in $(RELEASE_FILES); do \
	      dst="$$dir/$$src"; \
	      if test "$$src" = "Makefile"; then \
	        sed <"$$src" >"$$dst" -e 's/^\( *Y_\(MAKEDIR\|EXE\(\|_PKGS\|_HOME\|_SITE\)\|HOME_PKG\) *=\).*/\1/'; \
	        touch -r "$$src" "$$dst"; \
	      else \
	        cp -p "$$src" "$$dst"; \
	      fi; \
	    done; \
	    rm -f "$$dir"/*~ "$$dir"/*/*~; \
	    echo "$(RELEASE_VERSION)" > "$$dir/VERSION"; \
	    tar jcf "$(RELEASE_NAME)" "$$dir"; \
	    rm -rf "$$dir"; \
	    echo "$(RELEASE_NAME) created"; \
	  fi; \
	fi;

.PHONY: clean release tag

# -------------------------------------------------------- end of Makefile
