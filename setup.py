# Copyright 2012 Patrick Varilly, Stefano Angioletti-Uberti
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

# setup.py
"""DNACC: Calculation of DNA-coated colloid interactions

DNACC is a Python module to calculate the interaction details for a set of
DNA-coated colloids.  These include the free energy of binding and of
entropic repulsion and the probability of any two tethers (or groups of
tethers) to form a bond.

Currently supported tether and colloid configurations include: short dsDNA
tethers (treated as rigid rods) grafted on parallel plates (in explicit and
in mean-field representation), as well as on spheres.  An approximate
sphere-sphere potential can be obtained using the Derjaguin approximation.

The module is extensible, so it should be relatively straigthforward to add
other kinds of tethers grafted on colloids with different geometries.

If you use this module, you *must* cite the following paper::

  Varilly, Angioletti-Uberti, Mognetti and Frenkel, "A general theory
  of DNA-mediated and other valence-limited interactions", submitted (2012)
"""

# A good portion of this setup.py file is modelled after SciPy's setup.py
DOCLINES = __doc__.split("\n")

import subprocess
import sys
import os

CLASSIFIERS = """\
Development Status :: 4 - Beta
Intended Audience :: Science/Research
License :: OSI Approved :: GNU General Public License v3 (GPLv3)
Operating System :: POSIX
Operating System :: Unix
Operating System :: MacOS
Operating System :: Microsoft :: Windows
Programming Language :: C
Programming Language :: Python
Topic :: Software Development
Topic :: Scientific/Engineering

"""

MAJOR = 1
MINOR = 1
MICRO = 0
ISRELEASED = False
VERSION = '%d.%d.%d' % (MAJOR, MINOR, MICRO)


# Return the git revision as a string
def git_version():
    def _minimal_ext_cmd(cmd):
        # construct minimal environment
        env = {}
        for k in ['SYSTEMROOT', 'PATH']:
            v = os.environ.get(k)
            if v is not None:
                env[k] = v
        # LANGUAGE is used on win32
        env['LANGUAGE'] = 'C'
        env['LANG'] = 'C'
        env['LC_ALL'] = 'C'
        out = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                               env=env).communicate()[0]
        return out

    try:
        out = _minimal_ext_cmd(['git', 'rev-parse', 'HEAD'])
        GIT_REVISION = out.strip().decode('ascii')
    except Exception:
        GIT_REVISION = "Unknown"
    return GIT_REVISION


def write_version_py(filename='dnacc/version.py'):
    cnt = """# THIS FILE IS GENERATED FROM DNACC SETUP.PY
short_version = '%(version)s'
version = '%(version)s'
full_version = '%(full_version)s'
git_revision = '%(git_revision)s'
release = %(isrelease)s

if not release:
    version = full_version
"""
    FULLVERSION = VERSION
    if os.path.exists('.git'):
        GIT_REVISION = git_version()
    elif os.path.exists('dnacc/version.py'):
        # must be a source distribution, use existing version file
        from dnacc.version import git_revision as GIT_REVISION
    else:
        GIT_REVISION = "Unknown"

    if not ISRELEASED:
        FULLVERSION += '.dev-' + GIT_REVISION[:7]

    with open(filename, 'w') as f:
        f.write(cnt % {'version': VERSION,
                       'full_version': FULLVERSION,
                       'git_revision': GIT_REVISION,
                       'isrelease': str(ISRELEASED)})

try:
    from setuptools import setup, Extension
except ImportError:
    from distutils.core import setup, Extension

# On OS X, the automatic downloading and installation of numpy and scipy
# is very problematic.  Avoid it
if sys.platform == 'darwin':
    try:
        import numpy
        import scipy
    except ImportError:
        print("""

FATAL ERROR: NumPy and SciPy both need to be installed to install DNACC.
Ordinarily, they would be installed automatically, but this seems to be
very problematic in OS X.  Fortunately, there are easy-to-install binaries
available for download.

Please read the README.md file for installation instructions.

""")
        exit()

# We need NumPy installed to find out the NumPy include dirs!
# Computer Says No.
try:
    import numpy
    numpy_includes = [numpy.get_include()]
except ImportError:
    print("""

WARNING: numpy isn't installed!!!

It should be downloaded automatically, but from there on, the
installation will fail.  Simply run the installation again once
numpy has been installed.

""")
    numpy_includes = []

# Rewrite version file every time
write_version_py()

setup(
    name="dnacc",
    version=VERSION,
    author='Patrick Varilly, Stefano Angioletti-Uberti',
    author_email='pv271@cam.ac.uk, sa601@cam.ac.uk',
    description=DOCLINES[0],
    long_description="\n".join(DOCLINES[2:]),
    url="https://github.com/patvarilly/dnacc",
    download_url="https://github.com/patvarilly/dnacc/downloads",
    license='GPLv3',
    classifiers=[_ for _ in CLASSIFIERS.split('\n') if _],
    platforms=["Linux", "Unix", "Mac OS-X", "Windows", "Solaris"],
    install_requires=[
        'numpy>=1.3.0',
        'scipy>=0.7.0'],
    packages=['dnacc'],
    ext_modules=[
        Extension("dnacc.generic", ['dnacc/generic.c'],
                  include_dirs=numpy_includes)],
    scripts=['simple_dnacc'])
