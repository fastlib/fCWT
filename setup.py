#!/usr/bin/env python

"""
setup.py file for SWIG
"""

import numpy
import os
from pathlib import Path
import platform
from setuptools import Extension, setup, find_packages
from setuptools.command.build_ext import build_ext
import shutil
import sysconfig

ROOT = Path(__file__).parent.resolve()
LIBS = ROOT / "libs"
PACKAGE_DIR = ROOT / "src" / "fcwt"

# Platform detection
PLATFORM = sysconfig.get_platform()
MACHINE = platform.machine().lower()

IS_MACOS = PLATFORM.startswith("macosx") or "darwin" in PLATFORM
IS_LINUX = PLATFORM.startswith("linux")
IS_WINDOWS = PLATFORM.startswith("win")
IS_X86 = MACHINE in ("x86_64", "amd64")

if IS_MACOS:
    os.environ.setdefault("ARCHFLAGS", f"-arch {MACHINE}")

class BuildExt(build_ext):
    def _copy_runtime_libs(self):
        PACKAGE_DIR.mkdir(parents=True, exist_ok=True)

        if IS_WINDOWS:
            shutil.copy2(LIBS / "fftw3f.dll", PACKAGE_DIR / "fftw3f.dll")

        elif IS_LINUX:
            # The bundled ELF libraries advertise their upstream SONAMEs.
            shutil.copy2(LIBS / "libfftw3fl.so", PACKAGE_DIR / "libfftw3f.so.3")
            shutil.copy2(LIBS / "libfftw3f_ompl.so", PACKAGE_DIR / "libfftw3f_omp.so.3")

    def build_extension(self, ext):
        super().build_extension(ext)
        self._copy_runtime_libs()

    def copy_extensions_to_source(self):
        super().copy_extensions_to_source()
        self._copy_runtime_libs()


# Obtain the numpy include directory.  This logic works across numpy versions.
try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()

include_dirs = ["src/fcwt", "src", "libs", numpy.get_include()]
library_dirs = ["libs"]
runtime_library_dirs = []
libraries = []
link_args = []
files2 = [
    "omp.h",
    "fftw3.h",
    "fftw3f.dll",
    "fftw3f.lib",
    "libfftw3fmac.a",
    "libfftw3f_ompmac.a",
    "libfftw3fl.so",
    "libfftw3f_ompl.so",
    "libomp.a",
]
files = ["fcwt.h", "fcwt.cpp"]

files = files + files2

plat = sysconfig.get_platform()
machine = platform.machine().lower()

if IS_MACOS:
    libraries = ["fftw3fmac", "fftw3f_ompmac"]
    comp_args = ["-O3", "-Xpreprocessor", "-fopenmp"]

    # Only Intel can use AVX.
    if IS_X86:
        comp_args.append("-mavx")

    omp_prefix = Path("/opt/homebrew/opt/libomp")
    if not omp_prefix.exists():
        omp_prefix = Path("/usr/local/opt/libomp")

    include_dirs.append(str(omp_prefix / "include"))
    library_dirs.append(str(omp_prefix / "lib"))
    omp_library = str(omp_prefix / "lib" / "libomp.dylib")
    link_args = [
        "-L" + str(omp_prefix / "lib"),
        "-Wl,-needed_library," + omp_library,
        "-Wl,-rpath," + str(omp_prefix / "lib"),
    ]

if IS_LINUX:
    libraries = ["fftw3fl", "fftw3f_ompl", "gomp"]
    comp_args = ["-mavx", "-O3", "-fopenmp"]
    runtime_library_dirs = ["$ORIGIN"]
    link_args = ["-fopenmp", "-Wl,-rpath,$ORIGIN"]


if IS_WINDOWS:
    libraries = ["fftw3f"]
    comp_args = ["/arch:AVX", "/O2", "/openmp"]


setup(
    ext_modules=[
        Extension(
            "fcwt._fcwt",
            sources=["src/fcwt/fcwt.cpp", "src/fcwt/fcwt_wrap.cxx"],
            include_dirs=include_dirs,
            library_dirs=library_dirs,
            runtime_library_dirs=runtime_library_dirs,
            libraries=libraries,
            extra_compile_args=comp_args,
            extra_link_args=link_args,
        )
    ],
    cmdclass={"build_ext": BuildExt},
    packages=find_packages(where="src"),
    package_dir={"fcwt": "src/fcwt"},
    package_data={"fcwt": files},
)

# swig -c++ -python fcwt-swig.i
