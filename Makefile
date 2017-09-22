#!/usr/bin/make

target	= kinema

src_dir		= src
include_dir	= include
bin_dir		= bin
build_dir	= $(src_dir)/build

CXX	= g++

root_config	= root-config
root_include	= $(shell $(root_config) --cflags)
root_libs	= $(shell $(root_config) --evelibs)

cflags	= -g -O3 -Wall
#dflags	= -DDEBUG

flags	= $(cflags) -I. -I$(include_dir) $(root_include) $(dflags)
libs	= $(root_libs) -lz

srcs	= $(wildcard $(src_dir)/*.cc)
deps	= $(srcs:$(src_dir)/%.cc=$(build_dir)/%.d)
objs	= $(srcs:$(src_dir)/%.cc=$(build_dir)/%.o)

echo	= echo -e

#===========================================================================#
.PHONY: all clean distclean

all: $(bin_dir)/$(target)

-include $(deps)

$(bin_dir)/$(target): $(objs)
	@ $(echo) "\e[35;1m=== Linking $@ \e[m"
	@ mkdir -p $(bin_dir)
	$(CXX) -o $@ $(objs) $(libs)

clean:
	@ $(echo) '\e[32;1m=== Cleaning \e[m'
	@ rm -rfv $(build_dir)
	@ rm -fv *.d *.o *~ \#*

distclean:
	@ $(echo) "\e[32;1m=== Cleaning \e[m"
	@ rm -rfv $(bin_dir)/$(target) $(build_dir)
	@ rm -fv *.d *.o *~ \#*

$(build_dir)/%.o: $(src_dir)/%.cc
	@ $(echo) "\e[32;1m=== Compiling $@ \e[m"
	@ mkdir -p $(build_dir)
	$(CXX) $(flags) -o $@ -MMD -c $<
