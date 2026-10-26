# The rescue plugin needs librescueranger, see "Optional Compilation of the rescue Plugin" in INSTALL. When the library is not
# found, the plugin is left out with a notice, and the rest of bcftools builds as usual.
#
# By default the library is the one pkg-config finds. To build against a particular installation, give its prefix (headers in
# <prefix>/include/rescueranger, the library in <prefix>/lib or <prefix>/lib64):
#
#   make RESCUERANGER_ROOT=/opt/rescueranger
#
# The directory of the library is put in the RPATH of the plugin, so that it is found when bcftools loads the plugin.
ifneq "$(RESCUERANGER_ROOT)" ""
RR_LIBDIR := $(patsubst %/,%,$(dir $(firstword $(wildcard $(RESCUERANGER_ROOT)/lib/librescueranger.* $(RESCUERANGER_ROOT)/lib64/librescueranger.*))))
ifeq "$(RR_LIBDIR)" ""
$(error RESCUERANGER_ROOT=$(RESCUERANGER_ROOT): expected librescueranger in $(RESCUERANGER_ROOT)/lib or $(RESCUERANGER_ROOT)/lib64)
endif
RR_CPPFLAGS := -I$(RESCUERANGER_ROOT)/include
RR_LIBS := -L$(RR_LIBDIR) -lrescueranger -Wl,-rpath,$(RR_LIBDIR)
else ifneq "$(shell pkg-config --exists rescueranger 2>/dev/null && echo yes)" ""
RR_LIBDIR := $(shell pkg-config --variable=libdir rescueranger)
RR_CPPFLAGS := $(shell pkg-config --cflags rescueranger)
RR_LIBS := $(shell pkg-config --libs rescueranger) -Wl,-rpath,$(RR_LIBDIR)
endif

ifeq "$(RR_LIBS)" ""
PLUGINS := $(filter-out plugins/rescue$(PLUGIN_EXT),$(PLUGINS))

# Said when the plugins are built, so that a missing plugin does not go unnoticed
plugins: rescue-not-built
.PHONY: rescue-not-built
rescue-not-built:
	@echo "rescue plugin not built: librescueranger not found (see INSTALL)"
else
plugins/rescue$(PLUGIN_EXT): plugins/rescue.c version.h version.c
	$(CC) $(PLUGIN_FLAGS) $(CFLAGS) $(ALL_CPPFLAGS) $(EXTRA_CPPFLAGS) $(RR_CPPFLAGS) $(LDFLAGS) -o $@ version.c $< $(PLUGIN_LIBS) $(RR_LIBS) $(LIBS)
endif
