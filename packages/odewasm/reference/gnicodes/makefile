CC := g++
CFLAGS := -Wall -Wstrict-overflow=2 -g
LDFLAGS := -lm

.PHONY : list shipit

## dependency-finding scheme from
## http://www.gnu.org/manual/make-3.77/html_mono/make.html#SEC42
## (see also include statement at end of this makefile)
%.d: %.c
	@$(SHELL) -ec '$(CC) -MM $(CPPFLAGS) $< \
		| sed '\''s/\($*\)\.o[ :]*/\1.o $@ : /g'\'' > $@; \
		[ -s $@ ] || rm -f $@'

sources := dr_irk2.c gni_irk2.c kepler.c solfix.c ratio-gun-eom.c \
        harmonic.c pendulum.c outer-solar.c

extra := makefile licence.txt README .htaccess

all: dr_irk2

dr_irk2 : $(sources:.c=.o)

# show all the dependencies, by reading the .d files
list :
	@for file in $(sources:.c=.d) ; do <$$file sed 's/.*://' ; done | \
	  words_per_line | sort -u ; echo $(extra) | words_per_line

# show all the source filed and their dependencies:
shipit :
	@put $(extra) \
	  $$( for file in $(sources:.c=.d) ; do \
            <$$file sed 's/.*://' ; done | \
	      words_per_line | sort -u )

include $(sources:.c=.d)