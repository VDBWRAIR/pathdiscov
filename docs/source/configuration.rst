=============
Configuration
=============

At this time the configuration files are a bit confusing. In a later release of the
pipeline this will be fixed. But for now, here is a brief overview of how the three 
configuration files interact:

* When you first :doc:`install <install>`, you copy ``usamriidPathDiscov/files/config.yaml.base``
  to ``usamriidPathDiscov/files/config.yaml`` and edit it to suit your needs
* When you run ``python setup.py install`` the ``config.yaml`` and ``sample.param.base``
  files are copied from the location where you edited it to the virtualenv's 
  installation directory::

    usamriidPathDiscov/lib/python*/site-packages/usamriidPathDiscov/files/

* When you run :doc:`scripts/usamriidPathDiscov_cli` the virtualenv's copy of 
  ``config.yaml`` is used to replace the values inside of the virtualenv's copy 
  of ``sample.param.base`` and that is used to create the :ref:`paramtxt` 
  file used by the pipeline.

.. _config-yaml-base:

`usamriidPathDiscov/files/config.yaml.base <../../../usamriidPathDiscov/files/config.yaml.base>`_
=================================================================================================


This file allows you to set the defaults that define where your databases are located
under the databases directory.
GENOMEDIR will be replaced with ~/databases(at some point that will be configurable)

When the pipeline is installed

.. code-block:: bash

   $> python setup.py install

usamriidPathDiscov/files/config.yaml.base is copied to usamriidPathDiscov/files/config.yaml
and GENOMEDIR is replaced with ~/databases

.. _sample-param-base:

`usamriidPathDiscov/files/sample.param.base <../../../usamriidPathDiscov/files/sample.param.base>`_
===================================================================================================


Then when the pipeline is run with usamriidPathDiscov_cli usamriidPathDiscov/files/config.yaml is used to modify
`usamriidPathDiscov/files/sample.param.base <../../../usamriidPathDiscov/files/sample.param.base>`_ which is copied to the python installation under
usamriidPathDiscov/lib/python*/site-packages/usamriidPathDiscov/files/sample.param

This file is then in turn used to create the param.txt file that run_standard_stable4.pl uses

The file is comprised of sections defining parameters for each :doc:`pipeline stage <stages/index>`
Each stage's documentation contains the information about the configuration available

.. _paramtxt:

param.txt
=========

The param.txt file is created each time you run the :doc:`scripts/usamriidPathDiscov_cli`

It is placed inside of the directory you specify with the ``--outdir`` parameter inside
of the input directory.

The param.txt has the same format as the :ref:`sample-param-base` file.
