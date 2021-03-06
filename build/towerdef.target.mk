# This file is generated by gyp; do not edit.

TOOLSET := target
TARGET := towerdef
DEFS_Debug := \
	'-DNODE_GYP_MODULE_NAME=towerdef' \
	'-DUSING_UV_SHARED=1' \
	'-DUSING_V8_SHARED=1' \
	'-DV8_DEPRECATION_WARNINGS=1' \
	'-D_LARGEFILE_SOURCE' \
	'-D_FILE_OFFSET_BITS=64' \
	'-DBUILDING_NODE_EXTENSION' \
	'-DDEBUG' \
	'-D_DEBUG' \
	'-DV8_ENABLE_CHECKS'

# Flags passed to all source files.
CFLAGS_Debug := \
	-fPIC \
	-pthread \
	-Wall \
	-Wextra \
	-Wno-unused-parameter \
	-m64 \
	-g \
	-O0

# Flags passed to only C files.
CFLAGS_C_Debug :=

# Flags passed to only C++ files.
CFLAGS_CC_Debug := \
	-fno-rtti \
	-fno-exceptions \
	-std=gnu++1y \
	-fexceptions

INCS_Debug := \
	-I/home/uaio/.electron-gyp/.node-gyp/iojs-3.0.8/include/node \
	-I/home/uaio/.electron-gyp/.node-gyp/iojs-3.0.8/src \
	-I/home/uaio/.electron-gyp/.node-gyp/iojs-3.0.8/deps/openssl/config \
	-I/home/uaio/.electron-gyp/.node-gyp/iojs-3.0.8/deps/openssl/openssl/include \
	-I/home/uaio/.electron-gyp/.node-gyp/iojs-3.0.8/deps/uv/include \
	-I/home/uaio/.electron-gyp/.node-gyp/iojs-3.0.8/deps/zlib \
	-I/home/uaio/.electron-gyp/.node-gyp/iojs-3.0.8/deps/v8/include \
	-I$(srcdir)/node_modules/nan

DEFS_Release := \
	'-DNODE_GYP_MODULE_NAME=towerdef' \
	'-DUSING_UV_SHARED=1' \
	'-DUSING_V8_SHARED=1' \
	'-DV8_DEPRECATION_WARNINGS=1' \
	'-D_LARGEFILE_SOURCE' \
	'-D_FILE_OFFSET_BITS=64' \
	'-DBUILDING_NODE_EXTENSION'

# Flags passed to all source files.
CFLAGS_Release := \
	-fPIC \
	-pthread \
	-Wall \
	-Wextra \
	-Wno-unused-parameter \
	-m64 \
	-O3 \
	-fno-omit-frame-pointer

# Flags passed to only C files.
CFLAGS_C_Release :=

# Flags passed to only C++ files.
CFLAGS_CC_Release := \
	-fno-rtti \
	-fno-exceptions \
	-std=gnu++1y \
	-fexceptions

INCS_Release := \
	-I/home/uaio/.electron-gyp/.node-gyp/iojs-3.0.8/include/node \
	-I/home/uaio/.electron-gyp/.node-gyp/iojs-3.0.8/src \
	-I/home/uaio/.electron-gyp/.node-gyp/iojs-3.0.8/deps/openssl/config \
	-I/home/uaio/.electron-gyp/.node-gyp/iojs-3.0.8/deps/openssl/openssl/include \
	-I/home/uaio/.electron-gyp/.node-gyp/iojs-3.0.8/deps/uv/include \
	-I/home/uaio/.electron-gyp/.node-gyp/iojs-3.0.8/deps/zlib \
	-I/home/uaio/.electron-gyp/.node-gyp/iojs-3.0.8/deps/v8/include \
	-I$(srcdir)/node_modules/nan

OBJS := \
	$(obj).target/$(TARGET)/cpp/ETP/Point.o \
	$(obj).target/$(TARGET)/cpp/index.o \
	$(obj).target/$(TARGET)/cpp/clipper/clipper.o \
	$(obj).target/$(TARGET)/cpp/core.o \
	$(obj).target/$(TARGET)/cpp/utilz.o \
	$(obj).target/$(TARGET)/cpp/RVO2/Agent.o \
	$(obj).target/$(TARGET)/cpp/RVO2/Definitions.o \
	$(obj).target/$(TARGET)/cpp/RVO2/KdTree.o \
	$(obj).target/$(TARGET)/cpp/RVO2/Obstacle.o \
	$(obj).target/$(TARGET)/cpp/RVO2/RVOSimulator.o \
	$(obj).target/$(TARGET)/cpp/RVO2/Vector2.o \
	$(obj).target/$(TARGET)/cpp/poly2tri/poly2tri.o \
	$(obj).target/$(TARGET)/cpp/poly2tri/common/shapes.o \
	$(obj).target/$(TARGET)/cpp/poly2tri/common/utils.o \
	$(obj).target/$(TARGET)/cpp/poly2tri/sweep/advancing_front.o \
	$(obj).target/$(TARGET)/cpp/poly2tri/sweep/cdt.o \
	$(obj).target/$(TARGET)/cpp/poly2tri/sweep/sweep_context.o \
	$(obj).target/$(TARGET)/cpp/poly2tri/sweep/sweep.o

# Add to the list of files we specially track dependencies for.
all_deps += $(OBJS)

# CFLAGS et al overrides must be target-local.
# See "Target-specific Variable Values" in the GNU Make manual.
$(OBJS): TOOLSET := $(TOOLSET)
$(OBJS): GYP_CFLAGS := $(DEFS_$(BUILDTYPE)) $(INCS_$(BUILDTYPE))  $(CFLAGS_$(BUILDTYPE)) $(CFLAGS_C_$(BUILDTYPE))
$(OBJS): GYP_CXXFLAGS := $(DEFS_$(BUILDTYPE)) $(INCS_$(BUILDTYPE))  $(CFLAGS_$(BUILDTYPE)) $(CFLAGS_CC_$(BUILDTYPE))

# Suffix rules, putting all outputs into $(obj).

$(obj).$(TOOLSET)/$(TARGET)/%.o: $(srcdir)/%.cc FORCE_DO_CMD
	@$(call do_cmd,cxx,1)

$(obj).$(TOOLSET)/$(TARGET)/%.o: $(srcdir)/%.cpp FORCE_DO_CMD
	@$(call do_cmd,cxx,1)

# Try building from generated source, too.

$(obj).$(TOOLSET)/$(TARGET)/%.o: $(obj).$(TOOLSET)/%.cc FORCE_DO_CMD
	@$(call do_cmd,cxx,1)

$(obj).$(TOOLSET)/$(TARGET)/%.o: $(obj).$(TOOLSET)/%.cpp FORCE_DO_CMD
	@$(call do_cmd,cxx,1)

$(obj).$(TOOLSET)/$(TARGET)/%.o: $(obj)/%.cc FORCE_DO_CMD
	@$(call do_cmd,cxx,1)

$(obj).$(TOOLSET)/$(TARGET)/%.o: $(obj)/%.cpp FORCE_DO_CMD
	@$(call do_cmd,cxx,1)

# End of this set of suffix rules
### Rules for final target.
LDFLAGS_Debug := \
	-pthread \
	-rdynamic \
	-m64

LDFLAGS_Release := \
	-pthread \
	-rdynamic \
	-m64

LIBS :=

$(obj).target/towerdef.node: GYP_LDFLAGS := $(LDFLAGS_$(BUILDTYPE))
$(obj).target/towerdef.node: LIBS := $(LIBS)
$(obj).target/towerdef.node: TOOLSET := $(TOOLSET)
$(obj).target/towerdef.node: $(OBJS) FORCE_DO_CMD
	$(call do_cmd,solink_module)

all_deps += $(obj).target/towerdef.node
# Add target alias
.PHONY: towerdef
towerdef: $(builddir)/towerdef.node

# Copy this to the executable output path.
$(builddir)/towerdef.node: TOOLSET := $(TOOLSET)
$(builddir)/towerdef.node: $(obj).target/towerdef.node FORCE_DO_CMD
	$(call do_cmd,copy)

all_deps += $(builddir)/towerdef.node
# Short alias for building this executable.
.PHONY: towerdef.node
towerdef.node: $(obj).target/towerdef.node $(builddir)/towerdef.node

# Add executable to "all" target.
.PHONY: all
all: $(builddir)/towerdef.node

