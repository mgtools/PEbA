# Generated automatically from Makefile.in by configure.
################################################################
# Process this file with top-level configure script to produce Makefile
#
# Copyright 2000 Clark Cooper
#
#  This file is part of EXPAT.
#
#  EXPAT is free software; you can redistribute it and/or modify it
#  under the terms of the License (based on the MIT/X license) contained
#  in the file COPYING that comes with this distribution.
#
# EXPAT IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
# CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN EXPAT.
#


SHELL = /bin/sh

srcdir = .
top_srcdir = ..
prefix = /user/julie/expat1.95.2
exec_prefix = ${prefix}

bindir = ${exec_prefix}/bin
sbindir = ${exec_prefix}/sbin
libexecdir = ${exec_prefix}/libexec
datadir = ${prefix}/share
sysconfdir = ${prefix}/etc
sharedstatedir = ${prefix}/com
localstatedir = ${prefix}/var
libdir = ${exec_prefix}/lib
infodir = ${prefix}/info
mandir = ${prefix}/man
includedir = ${prefix}/include
oldincludedir = /usr/include

subdir = lib

top_builddir = ..

INSTALL = ../conftools/install-sh -c
INSTALL_PROGRAM = ${INSTALL}
INSTALL_SCRIPT = ${INSTALL_PROGRAM}
INSTALL_DATA = ${INSTALL} -m 644

host_alias = alphaev56-dec-osf4.0d
host_triplet = alphaev56-dec-osf4.0d
AS = 
CC = cc
DLLTOOL = 
LIBTOOL = $(SHELL) $(top_builddir)/libtool
LN_S = ln -s
OBJDUMP = 
PACKAGE = expat
RANLIB = ranlib
VERSION = 1.95.2

LIBRARY = libexpat.la
SOURCES = xmlparse.c xmltok.c xmlrole.c
OBJECTS = $(SOURCES:.c=.o)
LTOBJECTS = $(SOURCES:.c=.lo)

TEMPLATES = xmltok_impl.c xmltok_ns.c
APIHEADER = expat.h
HEADERS = ascii.h iasciitab.h utf8tab.h xmltok.h asciitab.h latin1tab.h \
          nametab.h xmldef.h xmlrole.h xmltok_impl.h

mkinstalldirs = $(SHELL) $(top_srcdir)/conftools/mkinstalldirs
CONFIG_HEADER = ../config.h
CONFIG_CLEAN_FILES = 

INCLUDES = -I$(srcdir) -I. -I..
DEFS = -DHAVE_CONFIG_H -DPACKAGE='"$(PACKAGE)"' -DVERSION='"$(PACKAGE)_$(VERSION)"'

CPPFLAGS = 
LDFLAGS = 
LIBS = 
CFLAGS = -g

LIBREVISION = 0
LIBCURRENT  = 1
LIBAGE      = 1

COMPILE = $(CC) $(DEFS) $(INCLUDES) $(CPPFLAGS) $(CFLAGS)
LTCOMPILE = $(LIBTOOL) --mode=compile $(CC) $(DEFS) $(INCLUDES) $(CPPFLAGS) $(CFLAGS)
CCLD = $(CC)
LINK = $(LIBTOOL) --mode=link  $(CCLD) -version-info $(LIBCURRENT):$(LIBREVISION):$(LIBAGE) $(CFLAGS) $(LDFLAGS) -o $@
DIST_COMMON =  Makefile.in


DISTFILES = $(DIST_COMMON) $(SOURCES) $(TEMPLATES) $(APIHEADER) $(HEADERS) 

TAR = gtar
GZIP_ENV = --best

all: $(LIBRARY)

.SUFFIXES: .c .lo .o
.PHONY: all clean distclean maintainer-clean

.c.o:
	$(COMPILE) -c $<

.c.lo:
	$(LTCOMPILE) -c $<

Makefile: $(srcdir)/Makefile.in  $(top_builddir)/config.status
	cd $(top_builddir) \
	  && CONFIG_FILES=$(subdir)/$@ CONFIG_HEADERS= $(SHELL) ./config.status

$(top_builddir)/config.status: $(top_builddir)/configure
	cd $(top_builddir) && $(MAKE) config.status

$(top_builddir)/config.h: $(top_builddir)/config.h.in
	cd $(top_builddir) && $(MAKE) config.h

clean:
	rm -f $(LIBRARY) *.o *.lo *~
	rm -rf .libs _libs

distclean: clean
	rm -f Makefile

maintainer-clean: distclean

check:	$(SUBDIRS)
	@echo
	@echo This package does not yet have a regression test.
	@echo

install: $(LIBRARY) $(APIHEADER)
	$(mkinstalldirs) $(libdir) $(includedir)
	$(LIBTOOL) --mode=install $(INSTALL) $(LIBRARY) $(libdir)/$(LIBRARY)
	$(INSTALL_DATA) $(APIHEADER) $(includedir)

uninstall:
	$(LIBTOOL) --mode=uninstall rm -f $(libdir)/$(LIBRARY);
	rm -f $(libdir)/$(APIHEADER)

$(LIBRARY): $(LTOBJECTS)
	$(LINK) -rpath $(libdir) $(LDFLAGS) $(LTOBJECTS)

xmlparse.o \
xmlparse.lo: xmlparse.c expat.h xmlrole.h xmltok.h $(top_builddir)/config.h

xmlrole.o \
xmlrole.lo: xmlrole.c ascii.h xmlrole.h $(top_builddir)/config.h

xmltok.o \
xmltok.lo: xmltok.c xmltok_impl.c xmltok_ns.c \
          ascii.h asciitab.h iasciitab.h latin1tab.h nametab.h utf8tab.h \
          xmltok.h xmltok_impl.h $(top_builddir)/config.h
