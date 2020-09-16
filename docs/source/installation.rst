.. _install:

**************
Installation
**************


.. _install-main:

Main installation method
===========================

1. **Fork** 🍴 the ```pyEPR top-level repository```_ on
   GitHub. (`How to fork a GitHub repo?`_). Share some love by
   **staring** :star: `pyEPR`_.
2. **Clone** 👇 your forked repository locally. (`How to clone
   a GitHub repo?`_). Setup the ``pyEPR`` python code by following
   `Installation and Python Setup`_.
3. **Tutorials** Learn how to use using the `jupyter notebook
   tutorials`_
4. **Stay up to date** Enjoy and make sure to git add the master remote
   branch
   ``git remote add MASTER_MINEV git://github.com/zlatko-minev/pyEPR.git``
   `(help?)`_.
5. **Cite ``pyEPR``** `arXiv:1902.10355`_ and enjoy!  🎂

.. _``pyEPR top-level repository``: https://github.com/zlatko-minev/pyEPR
.. _How to fork a GitHub repo?: https://help.github.com/en/articles/fork-a-repo
.. _pyEPR: https://github.com/zlatko-minev/pyEPR/
.. _How to clone a GitHub repo?: https://help.github.com/en/articles/cloning-a-repository
.. _Installation and Python Setup: #installation-of-pyepr
.. _jupyter notebook tutorials: https://github.com/zlatko-minev/pyEPR/tree/master/_tutorial_notebooks
.. _(help?): https://stackoverflow.com/questions/11266478/git-add-remote-branch
.. _`arXiv:1902.10355`: https://arxiv.org/abs/1902.10355

.. _install-via_pip:

Installing locally with via pip
===============================

In the future, ``pyEPR`` can be installed using the Python package manager `pip <http://www.pip-installer.org/>`_.


However, for the moment, we recommend a local develper instalation, which allows for fast upgrades. We are still in active development.
Perform the steps in the :ref:`install-main` section.
What you could do, once you have the local clone git, is to install pyEPR locally. Navigate to the local root folder of the repo.

First, in bash, upgrade python ``pip``

.. code-block:: bash

    python -m pip install -U pip

Now we can locally install the pyEPR module.

.. code-block:: bash

    python -m pip install -r requirements.txt -e .


.. _install-via_conda:

Installing via conda
====================

For Python 3.6+, installation via conda will be supported in a future version.


In the meantime, if you are using conda, you can locally install from the cloned repo.
Perform the steps in the :ref:`install-main` section.
Now, you can use

.. code-block:: bash

    python -m pip install -e .