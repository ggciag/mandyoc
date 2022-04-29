How to install
==============

*Mandyoc* installation is very simple and it consists of installing both `PETSc`_
and *Mandyoc*.

.. warning::
	The following installation steps work for both Linux and macOS machines
	**only** and no tests were made to install *Mandyoc* on Windows machines yet.

.. _Dependencies:

Dependencies
------------

To build *Mandyoc*, the following requirements are needed:

* PETSc_ (currently tested on version v3.15.5)
* gcc
* make
* git (recommended, but not strictly needed)

If you do not already have a PETSc installation, you may need a Fortran compiler.

Additionally, the following additional software is needed to run the examples
that come with *Mandyoc*:

* Python 3.5+

With required python packages:

* numpy
* matplotlib
* pandas
* jupyterlab

To run the tests, some additional python packages are required:

* pytest
* numpy

To build the documentation, further python packages are necessary:

* sphinx
* sphinx_rtd_theme

PETSc Installation
------------------

*Mandyoc* requires the PETSc library to run.
Check out the `PETSc installation`_ for details.

.. note::

	The following steps its a example installation with minimum requirements
	to run *Mandyoc*.

**Requirements:**

* ``CMake``

The first step is to **download** PETSc (v3.15.5) release from `PETSc website`_
or **clone** the repository into your machine.
*Mandyoc* might work with latest release of PETSc, but this is not guaranteed
since new versions might introduce breaking changes.

Clone the repository to your desired location::

	git clone -b v3.15.5 https://gitlab.com/petsc/petsc.git

Second, **configure the PETSc build** and set up the installation directory.

.. code-block:: bash

	cd path/to/petsc
	./configure \
	  PETSC_DIR=/path/to/petsc \
	  PETSC_ARCH=arch-label-optimized \
	  --with-debugging=0 \
	  --with-cc=gcc \
	  --with-cxx=g++ \
	  --with-fc=gfortran \
	  --download-fblaslapack \
	  --download-mpich \
	  --download-mumps \
	  --download-scalapack \
	  --download-parmetis \
	  --download-metis

.. note::
	This example installation uses the ``gfortran`` compiler.
	You may use another Fortran compiler changing the ``--with-fc`` flag.
	It is also possible to configure PETSc without Fortran by using ``--with-fc=0``
	and changing ``--download-fblaslapack`` to ``--download-f2cblaslapack``.

.. note::

	By default, *Mandyoc* uses direct solvers (LU and Cholesky) provided by `MUMPS`_.
	This requires additional external packages. Refer to `PETSc documentation`_
	for further information.

.. note::

	If you want to build a development version of *Mandyoc*
	its recommended to build a **debug version** of PETSc
	by setting ``--with-debugging=1``.
	In this case, you may set ``PETSC_ARCH=arch-label-debug``.

.. note::

	If you prefer *openmpi*, you need to swith ``--download-mpich`` to ``--download-openmpi``.

**Check** the installation with:

.. code-block::

	make all check

Or follow the instructions that pop up on the terminal.

For further information about the PETSc library, check out the `PETSc website`_.

Finally, add a symlink of `mpirun` to `~/.local/bin`:

.. code-block::

	ln -s /path/to/petsc/arch-label-optimized/bin/mpirun ~/.local/bin/mpirun

.. note::

    Make sure the directory ``~/.local/bin`` exists, otherwise the above
    command will fail.
	You can create it by running ``mkdir -p ~/.local/bin``.


*Mandyoc* Installation
----------------------

To install the *Mandyoc* in your machine, you need to **clone or download the latest release** of the code from the `Mandyoc repository`_.
To clone the repository, navigate to the directory you wish to install *Mandyoc* and type:

.. code-block:: bash

   git clone https://github.com/ggciag/mandyoc
   cd mandyoc

Before to install Mandyoc, you mast *set an environment variable* which indicates the path to PETSc installation folder:

.. code-block:: bash

	export PETSC_DIR=/path/to/petsc

*Build Mandyoc* by running:

.. code-block::

	make all

Next, *install Mandyoc* with:

.. code-block::

	make install

By default, it will be installed in ``~/.local/bin``.

.. note::

	Make sure the directory ``~/.local/bin`` exists, otherwise the above command will fail.
	You can change the installation location setting ``INSTALL_PATH`` variable by running:

	.. code-block::

		make INSTALL_PATH=/path/to/install/mandyoc install

.. note::

	To print *Mandyoc* runtime options, run mandyoc with ``-flags`` command line
	argument.

**Check** Mandyoc installation with:

.. code-block::

	make test

.. note::

	You need python and some python packages to run the last commmand succesfully.
	Check out requirements in `Dependencies`_ section.

Docker Container
----------------

We provide a `Docker container`_ image for *Mandyoc*.
Docker is an implementation of container virtualization.
Citing their documentation "it is a lightweight, standalone, executable package of software
that includes everything needed to run an application:
code, runtime, system tools, system libraries and settings".

Visit the `Dockerhub Mandyoc repository`_ to find out more on how to use the container to run *Mandyoc*.

.. note::

	To use the *Mandyoc* docker image, it is required to install the Docker Engine.
	Find out more on `Install Docker Engine`_ page.

Examples
--------

The benchmarks and other experiments are located in the `examples <https://github.com/ggciag/mandyoc/tree/main/examples>`_ folder of the Mandyoc repository.

Inside each example folder, you find a Jupyter notebook with detailed explanation and instructions on how to run the experiment.



.. _PETSc: https://petsc.org/release/
.. _PETSc installation: https://petsc.org/release/install/
.. _PETSc website: https://petsc.org/release/download/
.. _PETSc documentation: https://petsc.org/main/docs/manualpages/Mat/MATSOLVERMUMPS.html
.. _Mandyoc repository: https://github.com/ggciag/mandyoc
.. _MUMPS: http://mumps.enseeiht.fr/
.. _Docker container: https://www.docker.com/resources/what-container
.. _Dockerhub Mandyoc repository: https://hub.docker.com/r/ggciag/mandyoc
.. _Install Docker Engine: https://docs.docker.com/engine/install/
