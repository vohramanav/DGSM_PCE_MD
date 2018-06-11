# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/vohram/DGSM_PCE_MD/Codes/Kinetics/TChem/lib/cvode-2.7.0

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/vohram/DGSM_PCE_MD/Codes/Kinetics/TChem/lib/build

# Include any dependencies generated for this target.
include src/nvec_ser/CMakeFiles/sundials_nvecserial_static.dir/depend.make

# Include the progress variables for this target.
include src/nvec_ser/CMakeFiles/sundials_nvecserial_static.dir/progress.make

# Include the compile flags for this target's objects.
include src/nvec_ser/CMakeFiles/sundials_nvecserial_static.dir/flags.make

src/nvec_ser/CMakeFiles/sundials_nvecserial_static.dir/nvector_serial.o: src/nvec_ser/CMakeFiles/sundials_nvecserial_static.dir/flags.make
src/nvec_ser/CMakeFiles/sundials_nvecserial_static.dir/nvector_serial.o: /home/vohram/DGSM_PCE_MD/Codes/Kinetics/TChem/lib/cvode-2.7.0/src/nvec_ser/nvector_serial.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/vohram/DGSM_PCE_MD/Codes/Kinetics/TChem/lib/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object src/nvec_ser/CMakeFiles/sundials_nvecserial_static.dir/nvector_serial.o"
	cd /home/vohram/DGSM_PCE_MD/Codes/Kinetics/TChem/lib/build/src/nvec_ser && /usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/sundials_nvecserial_static.dir/nvector_serial.o   -c /home/vohram/DGSM_PCE_MD/Codes/Kinetics/TChem/lib/cvode-2.7.0/src/nvec_ser/nvector_serial.c

src/nvec_ser/CMakeFiles/sundials_nvecserial_static.dir/nvector_serial.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/sundials_nvecserial_static.dir/nvector_serial.i"
	cd /home/vohram/DGSM_PCE_MD/Codes/Kinetics/TChem/lib/build/src/nvec_ser && /usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/vohram/DGSM_PCE_MD/Codes/Kinetics/TChem/lib/cvode-2.7.0/src/nvec_ser/nvector_serial.c > CMakeFiles/sundials_nvecserial_static.dir/nvector_serial.i

src/nvec_ser/CMakeFiles/sundials_nvecserial_static.dir/nvector_serial.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/sundials_nvecserial_static.dir/nvector_serial.s"
	cd /home/vohram/DGSM_PCE_MD/Codes/Kinetics/TChem/lib/build/src/nvec_ser && /usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/vohram/DGSM_PCE_MD/Codes/Kinetics/TChem/lib/cvode-2.7.0/src/nvec_ser/nvector_serial.c -o CMakeFiles/sundials_nvecserial_static.dir/nvector_serial.s

src/nvec_ser/CMakeFiles/sundials_nvecserial_static.dir/nvector_serial.o.requires:

.PHONY : src/nvec_ser/CMakeFiles/sundials_nvecserial_static.dir/nvector_serial.o.requires

src/nvec_ser/CMakeFiles/sundials_nvecserial_static.dir/nvector_serial.o.provides: src/nvec_ser/CMakeFiles/sundials_nvecserial_static.dir/nvector_serial.o.requires
	$(MAKE) -f src/nvec_ser/CMakeFiles/sundials_nvecserial_static.dir/build.make src/nvec_ser/CMakeFiles/sundials_nvecserial_static.dir/nvector_serial.o.provides.build
.PHONY : src/nvec_ser/CMakeFiles/sundials_nvecserial_static.dir/nvector_serial.o.provides

src/nvec_ser/CMakeFiles/sundials_nvecserial_static.dir/nvector_serial.o.provides.build: src/nvec_ser/CMakeFiles/sundials_nvecserial_static.dir/nvector_serial.o


src/nvec_ser/CMakeFiles/sundials_nvecserial_static.dir/__/sundials/sundials_math.o: src/nvec_ser/CMakeFiles/sundials_nvecserial_static.dir/flags.make
src/nvec_ser/CMakeFiles/sundials_nvecserial_static.dir/__/sundials/sundials_math.o: /home/vohram/DGSM_PCE_MD/Codes/Kinetics/TChem/lib/cvode-2.7.0/src/sundials/sundials_math.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/vohram/DGSM_PCE_MD/Codes/Kinetics/TChem/lib/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object src/nvec_ser/CMakeFiles/sundials_nvecserial_static.dir/__/sundials/sundials_math.o"
	cd /home/vohram/DGSM_PCE_MD/Codes/Kinetics/TChem/lib/build/src/nvec_ser && /usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/sundials_nvecserial_static.dir/__/sundials/sundials_math.o   -c /home/vohram/DGSM_PCE_MD/Codes/Kinetics/TChem/lib/cvode-2.7.0/src/sundials/sundials_math.c

src/nvec_ser/CMakeFiles/sundials_nvecserial_static.dir/__/sundials/sundials_math.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/sundials_nvecserial_static.dir/__/sundials/sundials_math.i"
	cd /home/vohram/DGSM_PCE_MD/Codes/Kinetics/TChem/lib/build/src/nvec_ser && /usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/vohram/DGSM_PCE_MD/Codes/Kinetics/TChem/lib/cvode-2.7.0/src/sundials/sundials_math.c > CMakeFiles/sundials_nvecserial_static.dir/__/sundials/sundials_math.i

src/nvec_ser/CMakeFiles/sundials_nvecserial_static.dir/__/sundials/sundials_math.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/sundials_nvecserial_static.dir/__/sundials/sundials_math.s"
	cd /home/vohram/DGSM_PCE_MD/Codes/Kinetics/TChem/lib/build/src/nvec_ser && /usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/vohram/DGSM_PCE_MD/Codes/Kinetics/TChem/lib/cvode-2.7.0/src/sundials/sundials_math.c -o CMakeFiles/sundials_nvecserial_static.dir/__/sundials/sundials_math.s

src/nvec_ser/CMakeFiles/sundials_nvecserial_static.dir/__/sundials/sundials_math.o.requires:

.PHONY : src/nvec_ser/CMakeFiles/sundials_nvecserial_static.dir/__/sundials/sundials_math.o.requires

src/nvec_ser/CMakeFiles/sundials_nvecserial_static.dir/__/sundials/sundials_math.o.provides: src/nvec_ser/CMakeFiles/sundials_nvecserial_static.dir/__/sundials/sundials_math.o.requires
	$(MAKE) -f src/nvec_ser/CMakeFiles/sundials_nvecserial_static.dir/build.make src/nvec_ser/CMakeFiles/sundials_nvecserial_static.dir/__/sundials/sundials_math.o.provides.build
.PHONY : src/nvec_ser/CMakeFiles/sundials_nvecserial_static.dir/__/sundials/sundials_math.o.provides

src/nvec_ser/CMakeFiles/sundials_nvecserial_static.dir/__/sundials/sundials_math.o.provides.build: src/nvec_ser/CMakeFiles/sundials_nvecserial_static.dir/__/sundials/sundials_math.o


# Object files for target sundials_nvecserial_static
sundials_nvecserial_static_OBJECTS = \
"CMakeFiles/sundials_nvecserial_static.dir/nvector_serial.o" \
"CMakeFiles/sundials_nvecserial_static.dir/__/sundials/sundials_math.o"

# External object files for target sundials_nvecserial_static
sundials_nvecserial_static_EXTERNAL_OBJECTS =

src/nvec_ser/libsundials_nvecserial.a: src/nvec_ser/CMakeFiles/sundials_nvecserial_static.dir/nvector_serial.o
src/nvec_ser/libsundials_nvecserial.a: src/nvec_ser/CMakeFiles/sundials_nvecserial_static.dir/__/sundials/sundials_math.o
src/nvec_ser/libsundials_nvecserial.a: src/nvec_ser/CMakeFiles/sundials_nvecserial_static.dir/build.make
src/nvec_ser/libsundials_nvecserial.a: src/nvec_ser/CMakeFiles/sundials_nvecserial_static.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/vohram/DGSM_PCE_MD/Codes/Kinetics/TChem/lib/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking C static library libsundials_nvecserial.a"
	cd /home/vohram/DGSM_PCE_MD/Codes/Kinetics/TChem/lib/build/src/nvec_ser && $(CMAKE_COMMAND) -P CMakeFiles/sundials_nvecserial_static.dir/cmake_clean_target.cmake
	cd /home/vohram/DGSM_PCE_MD/Codes/Kinetics/TChem/lib/build/src/nvec_ser && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/sundials_nvecserial_static.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/nvec_ser/CMakeFiles/sundials_nvecserial_static.dir/build: src/nvec_ser/libsundials_nvecserial.a

.PHONY : src/nvec_ser/CMakeFiles/sundials_nvecserial_static.dir/build

src/nvec_ser/CMakeFiles/sundials_nvecserial_static.dir/requires: src/nvec_ser/CMakeFiles/sundials_nvecserial_static.dir/nvector_serial.o.requires
src/nvec_ser/CMakeFiles/sundials_nvecserial_static.dir/requires: src/nvec_ser/CMakeFiles/sundials_nvecserial_static.dir/__/sundials/sundials_math.o.requires

.PHONY : src/nvec_ser/CMakeFiles/sundials_nvecserial_static.dir/requires

src/nvec_ser/CMakeFiles/sundials_nvecserial_static.dir/clean:
	cd /home/vohram/DGSM_PCE_MD/Codes/Kinetics/TChem/lib/build/src/nvec_ser && $(CMAKE_COMMAND) -P CMakeFiles/sundials_nvecserial_static.dir/cmake_clean.cmake
.PHONY : src/nvec_ser/CMakeFiles/sundials_nvecserial_static.dir/clean

src/nvec_ser/CMakeFiles/sundials_nvecserial_static.dir/depend:
	cd /home/vohram/DGSM_PCE_MD/Codes/Kinetics/TChem/lib/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/vohram/DGSM_PCE_MD/Codes/Kinetics/TChem/lib/cvode-2.7.0 /home/vohram/DGSM_PCE_MD/Codes/Kinetics/TChem/lib/cvode-2.7.0/src/nvec_ser /home/vohram/DGSM_PCE_MD/Codes/Kinetics/TChem/lib/build /home/vohram/DGSM_PCE_MD/Codes/Kinetics/TChem/lib/build/src/nvec_ser /home/vohram/DGSM_PCE_MD/Codes/Kinetics/TChem/lib/build/src/nvec_ser/CMakeFiles/sundials_nvecserial_static.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/nvec_ser/CMakeFiles/sundials_nvecserial_static.dir/depend

