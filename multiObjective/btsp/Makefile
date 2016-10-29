include Makefile.in
include Makefile.lib

#--------------------------------------------------------------------------
OBJS :=$(foreach OBJ, $(OBJS), $(join $(SRC_DIR), $(OBJ)) )
HDRS :=$(foreach OBJ, $(HDRS), $(join $(SRC_DIR), $(OBJ)) )
LIBMISC :=$(foreach ITER, $(LIBMISC), $(join $(LIBMISC_DIR)/, $(ITER)) )

PARETO_HDRS:=$(foreach PARETO_HDR, $(PARETO_HDRS), $(join $(PARETO_DIR)/, $(PARETO_HDR)))

HDRS := $(HDRS) $(PARETO_HDRS)

TAGS = etags  --totals=yes
march=native

LINKERFLAGS  = $(OFLAGS) $(LIBMISC)
INCLUDE_DIRS := -I$(LIBMISC_DIR) -I$(PARETO_DIR) -I.
override CFLAGS += $(OFLAGS) \
		  -D DEBUG=$(DEBUG) -D PROBLEM_VERSION=\"$(PROBLEM_VERSION)\"\
		  $(DEBUGFLAGS) $(INCLUDE_DIRS) $(PROBLEMFLAGS)

ifneq ($(DEBUG), 0)
DEBUGFLAGS     += -g3 -DMALLOC_CHECK_=3
endif

SVN_REV := $(if $(shell which svnversion 1> /dev/null 2>&1 && test `svnversion -n .` != exported  && echo 1),$(shell svnversion -n . | tee svn_version),$(shell cat svn_version 2> /dev/null))

guess-march := $(strip $(shell $(CC) $(OFLAGS) -march=$(march) -x c -S -\#\#\# - < /dev/null 2>&1  | \
                grep -m 1 -o -e "march=[^'\"]*" | \
                sed 's,march=,,'))

ifndef march
  $(error please define an architecture, e.g., 'make march=pentium')
endif

override CFLAGS += -march=$(march) -DMARCH='"$(guess-march)"'

ECHO = @echo "$(1)"



.PHONY : all cleanall cleanexe clean problem  misc tags

problem : misc $(LIBPROBLEM)

tags:
	$(TAGS) $(MISC_MAKE_DIR)/*.h $(OBJS:.o=.c) $(HDRS)

$(LIBPROBLEM) : $(OBJS) $(HDRS)
	$(call ECHO,--> Creating problem object \"$(LIBPROBLEM)\" ($(PROBLEM_VERSION)) <---)
	ar rcs $(LIBPROBLEM) $(OBJS)
#	armerge $(LIBPROBLEM) ./libsolution.a $(LIBMISC)
#	rm -f ./libsolution.a
#	ld -r -o $(LIBPROBLEM) $(OBJS) $(LINKERFLAGS)


$(OBJS) : $(HDRS) $(LIBMISC)

misc $(LIBMISC):
	$(MAKE) $(LIB_OPTS) DEBUG=$(DEBUG) -C $(MISC_MAKE_DIR)


%.o : %.c
	$(CC)  $(CFLAGS)  -c -o $@ $<

cleanall : clean cleanexe

cleanexe :
	$(call ECHO,---> Removing $(LIBPROBLEM) <---)
	@$(RM) $(LIBPROBLEM)

clean :
	$(call ECHO,---> Removing object files <---)
	@$(RM) $(OBJS)
	@$(RM) *~

all :  LIB_OPTS := all
all :  cleanall problem
