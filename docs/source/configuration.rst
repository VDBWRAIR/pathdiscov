=============
Configuration
=============

At this time the configuration files are a bit confusing. In a later release of the
pipeline this will be fixed. But for now, here is a brief overview of how the three 
configuration files interact:

* When you first :doc:`install <install>`, you copy ``pathdiscov/files/config.yaml.base``
  to ``pathdiscov/files/config.yaml`` and edit it to suit your needs
* When you run ``python setup.py install`` the ``config.yaml`` and ``sample.param.base``
  files are copied from the location where you edited it to the virtualenv's 
  installation directory::

    pathdiscov/lib/python*/site-packages/pathdiscov/files/

* When you run :doc:`scripts/pathdiscov_cli` the virtualenv's copy of 
  ``config.yaml`` is used to replace the values inside of the virtualenv's copy 
  of ``sample.param.base`` and that is used to create the :ref:`paramtxt` 
  file used by the pipeline.

.. _config-yaml-base:

`pathdiscov/files/config.yaml.base <../../../pathdiscov/files/config.yaml.base>`_
=================================================================================================


This file allows you to set the defaults that define where your databases are located
under the databases directory.
GENOMEDIR will be replaced with ~/databases(at some point that will be configurable)

When the pipeline is installed

.. code-block:: bash

   $> python setup.py install

pathdiscov/files/config.yaml.base is copied to pathdiscov/files/config.yaml
and GENOMEDIR is replaced with ~/databases

.. _sample-param-base:

`pathdiscov/files/sample.param.base <../../../pathdiscov/files/sample.param.base>`_
===================================================================================================


Then when the pipeline is run with pathdiscov_cli pathdiscov/files/config.yaml is used to modify
`pathdiscov/files/sample.param.base <../../../pathdiscov/files/sample.param.base>`_ which is copied to the python installation under
pathdiscov/lib/python*/site-packages/pathdiscov/files/sample.param

This file is then in turn used to create the :ref:`paramtxt` file that run_standard_stable4.pl uses

The file is comprised of sections defining parameters for each :doc:`pipeline stage <stages/index>`
Each stage's documentation contains the information about the configuration available

.. _paramtxt:

param.txt
=========

The :ref:`paramtxt` file is created each time you run the :doc:`scripts/pathdiscov_cli`

It is placed inside of the directory you specify with the ``--outdir`` parameter inside
of the input directory.

The :ref:`paramtxt` has the same format as the :ref:`sample-param-base` file.
