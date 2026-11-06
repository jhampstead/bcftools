# The rescue plugin needs librescue, see "Optional Compilation of the rescue Plugin" in INSTALL. When the library is not
# found, the plugin is left out with a notice, and the rest of bcftools builds as usual.
#
# By default the library is the one pkg-config finds. To build against a particular installation, give its prefix (headers in
# <prefix>/include/rescue, the library in <prefix>/lib or <prefix>/lib64):
#
#   make RESCUE_ROOT=/opt/rescue
#
# The directory of the library is put in the RPATH of the plugin, so that it is found when bcftools loads the plugin.
ifneq "$(RESCUE_ROOT)" ""
RR_LIBDIR := $(patsubst %/,%,$(dir $(firstword $(wildcard $(RESCUE_ROOT)/lib/librescue.* $(RESCUE_ROOT)/lib64/librescue.*))))
ifeq "$(RR_LIBDIR)" ""
$(error RESCUE_ROOT=$(RESCUE_ROOT): expected librescue in $(RESCUE_ROOT)/lib or $(RESCUE_ROOT)/lib64)
endif
RR_CPPFLAGS := -I$(RESCUE_ROOT)/include
RR_LIBS := -L$(RR_LIBDIR) -lrescue -Wl,-rpath,$(RR_LIBDIR)
else ifneq "$(shell pkg-config --exists rescue 2>/dev/null && echo yes)" ""
RR_LIBDIR := $(shell pkg-config --variable=libdir rescue)
RR_CPPFLAGS := $(shell pkg-config --cflags rescue)
RR_LIBS := $(shell pkg-config --libs rescue) -Wl,-rpath,$(RR_LIBDIR)
endif

ifeq "$(RR_LIBS)" ""
PLUGINS := $(filter-out plugins/rescue$(PLUGIN_EXT),$(PLUGINS))

# Said when the plugins are built, so that a missing plugin does not go unnoticed
plugins: rescue-not-built
.PHONY: rescue-not-built
rescue-not-built:
	@echo "rescue plugin not built: librescue not found (see INSTALL)"
else
plugins/rescue$(PLUGIN_EXT): plugins/rescue.c version.h version.c
	$(CC) $(PLUGIN_FLAGS) $(CFLAGS) $(ALL_CPPFLAGS) $(EXTRA_CPPFLAGS) $(RR_CPPFLAGS) $(LDFLAGS) -o $@ version.c $< $(PLUGIN_LIBS) $(RR_LIBS) $(LIBS)
endif
