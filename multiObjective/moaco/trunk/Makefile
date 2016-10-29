include Makefile.in
PWD	:= $(shell /bin/pwd)
RNODE = iridiacluster
#---------------------------------------------------------------------
PROBLEM_OBJ :=$(foreach ITER, $(PROBLEM_OBJ),$(join $(PROBLEM_DIR)/,$(ITER)))
LIBMISC :=$(foreach ITER, $(LIBMISC), $(join $(LIBMISC_DIR)/, $(ITER)))

PARETO_OBJS:=$(foreach PARETO_OBJ, $(PARETO_OBJS), $(join $(PARETO_DIR)/, $(PARETO_OBJ)))
PARETO_HDRS:=$(foreach PARETO_HDR, $(PARETO_HDRS), $(join $(PARETO_DIR)/, $(PARETO_HDR)))

OBJS := $(MAIN_OBJS) $(PARETO_OBJS)
HDRS := $(MAIN_HDRS) $(PARETO_HDRS)

EXE =
ECHO = @echo "$(1)"

march=native

override LINKERFLAGS := $(OFLAGS) $(DEBUGFLAGS) \
              $(LIBPROBLEM) $(LIBMISC) $(LINKERFLAGS) -lm # Order is important !
INCLUDE_DIRS := -I$(PROBLEM_DIR) -I$(LIBMISC_DIR) -I$(PARETO_DIR) -I.
override CFLAGS += -std=gnu99 $(OFLAGS) \
		-D DEBUG=$(DEBUG) $(DEBUGFLAGS)

ifneq ($(DEBUG),0)
DEBUGFLAGS = -g3 -DMALLOC_CHECK_=3
else
DEBUGFLAGS = -DNDEBUG
endif

ifndef BINDIR
BINDIR = ./
endif

SVN_REV := $(if $(shell which svnversion 1> /dev/null 2>&1 && test `svnversion -n .` != exported  && echo 1),$(shell svnversion -n . | tee svn_version),$(shell cat svn_version 2> /dev/null))


guess-march := $(strip $(shell $(CC) $(OFLAGS) -march=$(march) -x c -S -\#\#\# - < /dev/null 2>&1  | \
                grep -m 1 -o -e "march=[^'\"]*" | \
                sed 's,march=,,'))

ifndef march
  $(error please define an architecture, e.g., 'make march=pentium')
endif

override CFLAGS += -march=$(march) -DMARCH='"$(guess-march)"'


MAINEXE=$(BINDIR)/$(MAIN)$(EXE)

.PHONY : all cleanall cleanexe clean misc problem rsync dist main

main : misc problem $(MAINEXE)

$(MAINEXE) : $(LIBMISC) $(LIBPROBLEM) $(OBJS) $(HDRS)
	@mkdir -p $(BINDIR)
	$(call ECHO,--> Building $(MAIN) version $(VERSION) <---)
	$(CC) -o $(MAINEXE) $(OBJS) $(LINKERFLAGS)

$(OBJS) : $(HDRS) $(LIBMISC) $(LIBPROBLEM)

%.o : %.c
	$(CC)  $(CFLAGS) $(INCLUDE_DIRS) -D VERSION='"$(VERSION)"' -c -o $@ $<

misc $(LIBMISC):
	$(MAKE) march=$(march) DEBUG=$(DEBUG) OFLAGS='$(OFLAGS)' $(MAKE_OPTS) -C $(MISC_MAKE_DIR)

problem $(PROBLEM_OBJ) :
	$(MAKE) march=$(march) DEBUG=$(DEBUG) OFLAGS='$(OFLAGS)' $(MAKE_OPTS) -C $(PROBLEM_MAKE_DIR)

cleanlibmisc :
	$(MAKE) clean -C $(MISC_MAKE_DIR)

cleanproblem :
	$(MAKE) clean -C $(PROBLEM_MAKE_DIR)

cleanall : cleanlibmisc cleanproblem clean cleanexe

cleanexe :
	$(call ECHO,---> Removing $(MAIN) <---)
	@$(RM) $(MAINEXE)

clean :
	$(call ECHO,---> Removing object files <---)
	@$(RM) $(BINDIR)/core.*
	@$(RM) $(OBJS)
	@$(RM) *~

all : MAKE_OPTS := all
all : misc problem
	$(MAKE) cleanexe clean
	$(MAKE) $(MAINEXE)

rsync : main
ifdef RNODE
	cd $(BASE_DIR) && \
	rsync -rlp -CIzc -L --delete --copy-unsafe-links --exclude=.svn --progress --relative \
	moaco/trunk moaco/launch libpareto libmisc/current $(PROBLEM) \
	$(RNODE):~/src/
else
	@echo "ERROR: You must specify a remote node (e.g., RNODE=b21)"
	@exit 1
endif

DIST_SRC_FILES =  Makefile Makefile.in svn_version README LICENSE PISA_LICENSE.txt *.c *.h
DIST_SRC := moaco-$(VERSION)-src

dist : DEBUG=0
dist : DEBUGFLAGS= -DNDEBUG
dist : all
	@(mkdir -p $(BASE_DIR)/$(DIST_SRC)/moaco/trunk			\
	&& rsync -rlpC --exclude=.svn --exclude=todo --exclude=TODO \
	$(DIST_SRC_FILES) $(BASE_DIR)/$(DIST_SRC)/moaco/trunk/ \
        && cd $(BASE_DIR) 						\
	&& rsync -rlLpC --exclude=todo --exclude=TODO --exclude=.svn --exclude='svn-commit.tmp' --relative \
	libpareto libmisc/current btsp ./$(DIST_SRC)/     	 	\
	&& tar cf - $(DIST_SRC) | gzip -f9 > $(DIST_SRC).tar.gz 	\
	&& rm -rf ./$(DIST_SRC)/ 					\
	&& echo "$(DIST_SRC).tar.gz created." && cd $(PWD) )

test : TMPDIR=$(shell mktemp -d)
test : dist
	(cd $(TMPDIR) \
	&& tar -xaf $(shell readlink -f $(BASE_DIR))/$(DIST_SRC).tar.gz \
	&& make -C $(DIST_SRC)/moaco/trunk && echo "It builds!") \
	|| echo "FAILED!"
