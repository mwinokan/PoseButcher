.. PoseButcher documentation master file, created by
   sphinx-quickstart on Wed Mar 20 14:31:51 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

=========================
PoseButcher Documentation
=========================

Installation
============

PoseButcher has several dependencies, which can all be managed via a pip installation:

::

   pip install --upgrade posebutcher

However, open3d and PyGAMer can occasionally pose problems. The latter is only necessary to generate protein surfaces, which can be solved on a different machine, or omitted from calculations. In these cases prefer to install posebutcher on its own:

::

   pip install --upgrade posebutcher --no-dependencies

Getting started
===============

See the example notebook_ on github.

.. _notebook: https://github.com/mwinokan/PoseButcher/blob/main/example.ipynb

.. toctree::
   :maxdepth: 2
   :caption: Contents:
