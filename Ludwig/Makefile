# Ludwig Makefile
#

INSTALL = ./
INSTALL_PROGRAM = $(INSTALL)

#include Makefile.conf

DIRS = src

# the sets of directories to do various things in
BUILDDIRS = $(DIRS:%=build-%)
INSTALLDIRS = $(DIRS:%=install-%)
CLEANDIRS = $(DIRS:%=clean-%)

all: $(BUILDDIRS)
$(DIRS): $(BUILDDIRS)
$(BUILDDIRS):
	$(MAKE) -C $(@:build-%=%)

# the utils need the libraries in dev built first
build-utils: build-dev

install: $(INSTALLDIRS) $(COPYEXE)
$(INSTALLDIRS):
	$(MAKE) -C $(@:install-%=%) mpi
$(COPYEXE):
	cp src/Ludwig.exe $(INSTALL)

clean: $(CLEANDIRS)
$(CLEANDIRS): 
	$(MAKE) -C $(@:clean-%=%) clean


.PHONY: subdirs $(DIRS)
.PHONY: subdirs $(BUILDDIRS)
.PHONY: subdirs $(INSTALLDIRS)
.PHONY: subdirs $(CLEANDIRS)
.PHONY: all install clean 