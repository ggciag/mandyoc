How to install
==============

MANDYOC installation is very simple and it consists of installing both `PETSc`_
and MANDYOC.

.. warning::
	The following installation steps work for both Linux and macOS machines
	**only** and no tests were made to install MANDYOC on Windows machines yet.

Dependencies
------------

* PETSc_
* gcc
* make
* git

*Optional* dependencies:

* `gfortran`

PETSc Installation
------------------

MANDYOC requires the `PETSc`_ library to run.
The first step is to **download** the latest release of PETSc from `PETSc website`_
or **clone** the repository into your machine::

	git clone -b release https://gitlab.com/petsc/petsc.git $HOME/petsc

By default, we will download/clone in ``~/petsc``.

Second, **configure the PETSc build** and set up the installation directory.
By default, we will install PETSc in ``~/petsc``.

.. code-block:: bash

	cd $HOME/petsc
	./configure \
	  PETSC_DIR=$HOME/petsc \
	  PETSC_ARCH=arch-label-optimized \
	  --with-debugging=0 \
	  --with-cc=gcc \
	  --with-cxx=g++ \
	  --download-fblaslapack \
	  --download-mpich \
	  --download-hdf5 \
	  --download-superlu_dist \
	  --download-metis \
	  --download-parmetis \
	  --download-mumps \
	  --download-scalapack \
	  --download-cmake \
	  COPTFLAGS='-O3 -march=native -mtune=native' \
	  CXXOPTFLAGS='-O3 -march=native -mtune=native'


.. note::

	If using ``gfortran`` optional dependency add the options
	``--with-fc=gfortran`` and ``FOPTFLAGS='-O3 -march=native -mtune=native'``
	to the PETSc build configuration above.

.. note::

	If you are build a development version of MANDYOC you can build
	a **debug version** of PETSc by setting ``--with-debugging=1`` and removing
	the ``COPTFLAGS``, ``CXXOPTFLAGS`` (and ``FOPTFLAGS``) flags.
	In this case, you may set ``PETSC_ARCH=arch-label-debug``.

**Check** the installation with:

.. code-block::

	make all check

Or follow the instructions that pop up on the terminal.

For further information about the PETSc library, check the `PETSc website`_.

MANDYOC Installation
--------------------

To install the MANDYOC in your machine,  you need to **clone or download  the latest release** of the
code from the `MANDYOC repository page`_.

To clone the repository, navigate to the directory you wish to install MANDYOC and type:

.. code-block:: bash

   git clone https://github.com/ggciag/mandyoc

Next, **build and install** MANDYOC by running::

	make all

.. note::

	To print MANDYOC runtime options, run mandyoc with `-flags` command line
	argument.



.. _PETSc: https://www.mcs.anl.gov/petsc/
.. _PETSc website: https://www.mcs.anl.gov/petsc/download/index.html
.. _PETSc repository: https://bitbucket.org/petsc/petsc/src/maint/
.. _MANDYOC repository page: https://github.com/ggciag/mandyoc