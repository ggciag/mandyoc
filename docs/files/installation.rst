How to install
==============

*Mandyoc* installation is very simple and it consists of installing both `PETSc`_
and *Mandyoc*.

.. warning::
	The following installation steps work for both Linux and macOS machines
	**only** and no tests were made to install *Mandyoc* on Windows machines yet.

Dependencies
------------

* PETSc_
* gcc
* make
* git

**Optional** dependencies:

* gfortran
* pytest
* numpy
* os

PETSc Installation
------------------

*Mandyoc* requires the `PETSc`_ library to run.
The first step is to **download** the latest release of PETSc from `PETSc website`_
or **clone** the repository into your machine.

Choose the path to your PETSC installation and clone the repository::

	cd /path/to/petsc
	git clone -b release https://gitlab.com/petsc/petsc.git $HOME/petsc

Second, **configure the PETSc build** and set up the installation directory.

.. code-block:: bash

	cd path/to/petsc
	./configure \
	  PETSC_DIR=/path/to/petsc \
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

	If you are build a development version of *Mandyoc* you can build
	a **debug version** of PETSc by setting ``--with-debugging=1`` and removing
	the ``COPTFLAGS``, ``CXXOPTFLAGS`` (and ``FOPTFLAGS``) flags.
	In this case, you may set ``PETSC_ARCH=arch-label-debug``.

**Check** the installation with:

.. code-block::

	make all check

Or follow the instructions that pop up on the terminal.

For further information about the PETSc library, check the `PETSc website`_.

Finally, *add a symlinks* of `mpirun` to `~/.local/bin`

.. code-block::

	ln -s /path/to/pets/arch-label-optimized/bin/mpirun ~/.local/bin/mpirun


*Mandyoc* Installation
----------------------

To install the *Mandyoc* in your machine, you need to **clone or download the latest release** of the code from the `Mandyoc repository page`_.
To clone the repository, navigate to the directory you wish to install *Mandyoc* and type:

.. code-block:: bash

   git clone https://github.com/ggciag/mandyoc

Before to install Mandyoc, you mast *set an env variable* which indicates the path to PETSc installation folder:

.. code-block:: bash

	export PETSC_DIR=/path/to/petsc

Next, *build and install Mandyoc* by running::

	make all

.. note::

	To print *Mandyoc* runtime options, run mandyoc with `-flags` command line
	argument.

**Check** Mandyoc installation with:

.. code-block::

	make test

Examples
--------

The benchmarks and other experiments are located in the `examples <https://github.com/ggciag/mandyoc/tree/main/examples>`_ folder of the Mandyoc repository.

Inside each example folder, you find a ``README.md`` file with detailed explanation and instrutions on how to run the experiment.
First, you need to run the python script file named ``generate_input_files.py`` to generate the :ref:`input files<inputfiles>` needed by Mandyoc.
Then, you may execute `mandyoc` directly from a terminal command or update the bash script ``run.sh`` accordingly to your setup and execute it to run the experiment.


.. _PETSc: https://www.mcs.anl.gov/petsc/
.. _PETSc website: https://www.mcs.anl.gov/petsc/download/index.html
.. _PETSc repository: https://bitbucket.org/petsc/petsc/src/maint/
.. _Mandyoc repository page: https://github.com/ggciag/mandyoc
