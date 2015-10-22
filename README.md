# Python module for calculating DNA-coated colloid interactions and binding details
Written by Patrick Varilly and Stefano Angioletti-Uberti
Last updated: 17 Oct 2015

This program accompanies the papers:

> 1) P. Varilly, S. Angioletti-Uberti, B.M. Mognetti, and D. Frenkel. “A general theory 
> of DNA-mediated and other valence-limited colloidal interactions”. 
> The Journal of Chemical Physics, 137:094108–094122, 2012.
> arXiv:1205.6921 (2012)
>
> 2) S. Angioletti-Uberti, P. Varilly, B.M. Mognetti, A.V. Tkachenko, and D. Frenkel. 
> “Communication: A simple analytical formula for the free energy of ligand-receptor-mediated interactions”. 
> The Journal of Chemical Physics, 138:021102–021106, 2013.
> arXiv:1211.1873 (2013).

We kindly ask to cite those papers in any publication where this code is used.

You can download the Python package [here](http://github.com/downloads/patvarilly/DNACC/dnacc-1.0.1.tar.gz).

## 1. Compiling

For compiling on OS X or Windows, see below.

It should be nearly trivial to build and install this module, since I'm
using Python's built-in setup tools.  From the main directory, just run

```
python setup.py install --user
```

That will build and install the module in the folder
`~/.local/lib/pythonX.X/site-packages/dnacc`, where X.X is your version
of python.

If you want to run a script like `simple_dnacc` without installing
the whole module, make sure that the inner loops in C get built by
issuing the following command:

```
python setup.py build_ext --inplace
```

(If you're seeing the following warning:

> Could not import extension code to speed up inner loops
> -- using much slower pure Python version instead.

then issuing the above command should solve the problem)


### 1.1. Compiling in OS X

In Mac OS X, the usual automatic downloading of dependencies during a
`python setup.py install` isn't able to successfully install NumPy and SciPy
(in my machine, it's the lack of a Fortran compiler, but the SciPy docs
point to other potential problems).  So you have to download and install
NumPy and SciPy manually.

The good news: binaries are easily available
The bad news: they only work with the version of Python from
  http://www.python.org, not the version that ships with OS X!

So you have to download and install three things:

1. Python 2.7.2 from http://www.python.org/download
2. Latest version of NumPy from http://numpy.org
3. Latest version of SciPy from http://scipy.org

This is far less painful than it sounds.  In my own case, the disk images
that I ended up downloading were:

1. python-2.7.2-macosx10.6.dmg
2. numpy-1.6.1-py2.7-python.org-macosx10.6.dmg
3. scipy-0.10.0-py2.7-python.org-macosx10.6.dmg

CAREFUL: for some of these packages, the "link to the latest version" that
SourceForge suggests may be incorrect!  Do look at the full list of downloads
available and pick the one that is most appropriate to your own setup

Once you've installed python2.7, numpy and scipy, you can run the command

```
python setup.py install --user
```

The DNACC package will be installed in
`~/Library/Python/2.7/lib/python/site-packages`

#### Extra Note on GCC on Mac OS X

Thanks to Stefano Angioletti-Uberti for pointing this out.

You may or may not have problems with getting python setup.py to find gcc on
your Mac. Typically, if you installed gcc using the Xcode package provided
by Apple, its name is not simply `gcc-4.2` but something longer like
`i686-apple-darwin11-llvm-gcc-4.2`. This is a link to the gcc compiler, whose
full path name is (in my system, but I just used the automatic installation
done by Xcode)

```
path-to-compiler = "/Developer/usr/llvm-gcc-4.2/bin/llvm-gcc-4.2"
```

You can create a softlink that corrects this problem as follows:

```
sudo ln -s path-to-compiler /usr/bin/gcc-4.2
sudo ln -s path-to-compiler /usr/local/bin/gcc-4.2
```

### 1.2 Compiling on Windows

Thanks to Bortolo M. Mognetti for these instructions.

With a DOS prompt go into the dnacc directory and launch `setup.py install`

Be sure that NumPy and Scipy have been installed. A C++ compiler is also
required. A Visual C++ 2008 Express edition is available for free [here](www.microsoft.com/visualstudio/en-us/products/2008-editions/express).

## 2. Test programs

All the calculations that I have ever run on DNA-coated colloids are
included in the `examples` folder.  Generally, these have a `.py` file
that you should run.  They produce lots of data which can then be plotted.
For example, in `competing_linkages`, type

```
python competing_linkages.py
for f in *.gp; do gnuplot $f; done
```

MORE DETAILED DESCRIPTIONS TO COME

## 3. General usage

For detailed instructions, take a look at [the HTML documentation](http://github.com/downloads/patvarilly/DNACC/dnacc-1.0.1-docs.tar.gz).

These are generated using Sphinx as follows:

```
cd docs
make clean
make html
firefox _build/html/index.html
```


    
