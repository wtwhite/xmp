# xmp

The XMP program from the 2011 Bioinformatics paper "Faster exact maximum parsimony search with XMP", by myself and Barbara Holland. Or at least I hope so -- this was resurrected from an old laptop after the original site http://www.massey.ac.nz/~wtwhite/xmp became unavailable.

I have dumped this here in the state I found it in, untouched since 2010, with just a few tweaks to persuade the serial (single-CPU) version to build and run on both Windows and Linux. The executable produced is named fastdnamp -- I don't recall now why it was renamed to (or from?) XMP, but if it was because some existing program already has that name, please let me know.

There is also a parallel (multi-CPU) version that should in theory build in any environment that has an MPI implementation (there are [free implementations](https://en.wikipedia.org/wiki/Message_Passing_Interface#Official_implementations) for commodity PCs, and it once even worked nicely on a BlueGene/L supercomputer), but I haven't tried getting that to build myself since uploading here. That is probably going to be a nightmare, but if anyone makes progress on this, please let me know! Pull requests happily accepted.

Why does `Release` appear in different places for the Windows and Linux build instructions given below? See [here](https://stackoverflow.com/a/24470998/47984).

## Windows serial version build instructions

You will need [CMake](https://cmake.org/) and a C compiler. I succeeded in using MSVC 2017, but other compilers (e.g., clang or MinGW) should work fine. If your choice of compiler lacks its own concept of Release/Debug project configurations, then the Linux build instructions may be more appropriate.

1. Clone this repo and the [useful_libs repo](https://github.com/wtwhite/useful_libs) as sibling directories:

    ```
    git clone git@github.com:wtwhite/xmp.git
    git clone git@github.com:wtwhite/useful_libs.git
    ```

2. Build useful_libs:

    Open a Command Prompt window that sets up the compiler environment variables. E.g., for MSVC 2017, there is a shortcut named "x64 Native Tools Command Prompt for VS 2017" that opens such a window. From this window, the compiler can be run with `cl`, etc.
    
    ```
    cd useful_libs\source
    mkdir build
    cd build
    cmake ..
    cmake --build . --config Release
    ```

    This should create `..\..\lib\useful.lib`.

3. Build fastdnamp:

    ```
    cd ..\..\..\xmp\src
    mkdir build
    cd build
    cmake ..
    cmake --build . --config Release
    ```

    This should create `fastdnamp.exe` in the `Release` subdirectory. Run this to see a(n incomplete) list of options.

Change `--config Release` to `--config Debug` in the above to build the slower debug version instead.

## Linux serial version build instructions

You will need [CMake](https://cmake.org/), `make` and a C compiler.

1. Clone this repo and the [useful_libs repo](https://github.com/wtwhite/useful_libs) as sibling directories:

    ```
    git clone git@github.com:wtwhite/xmp.git
    git clone git@github.com:wtwhite/useful_libs.git
    ```

2. Build useful_libs:

    ```
    cd useful_libs/source
    mkdir build
    cd build
    cmake -DCMAKE_BUILD_TYPE=Release ..
    cmake --build .
    ```

    This should create `../../lib/libuseful.a`.

3. Build fastdnamp:

    ```
    cd ..\..\..\xmp\src
    mkdir build
    cd build
    cmake -DCMAKE_BUILD_TYPE=Release
    cmake --build .
    ```

    This should create the executable `fastdnamp` in the `Release` subdirectory. Run this to see a(n incomplete) list of options.

Change `CMAKE_BUILD_TYPE=Release` to `CMAKE_BUILD_TYPE=Debug` in the above to build the slower debug version instead.

## Known issues

* There is a crash bug (assertion failure) in the TBR method for finding initial upper bounds (UBs), which is enabled by default. Disable it by specifying `--tbr=N` on the command line. This is just one of several ways to obtain an initial UB; if the remaining options are not giving what you feel is a good enough bound, you can also specify an initial bound with `b <nn>`. If this turns out to be lower than the true MP score, the program will exit without finding any solutions.
* Not all command-line options are listed when the program is run without options.
* I have not attempted to build the parallel version since uploading this -- that will probably be tricky.
* For some reason, CMake 3.11.0-rc3 on Windows detects my 64-bit compiler, OS and computer as 32-bit. If you know how to fix/override this, please let me know. CMake 3.5.1 on Linux correctly detects a 64-bit compiler.
* The SSE2-optimized version only builds with MSVC on Windows. Porting this to other x86 environments is compiler-specific and probably a fair amount of work.
