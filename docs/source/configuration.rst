=============
Configuration
=============

.. _config-yaml-base:

`usamriidPathDiscov/files/config.yaml.base <../../../usamriidPathDiscov/files/config.yaml.base>`_
=================================================================================================


This file allows you to set the defaults that define where your databases are located
under the databases directory.
GENOMEDIR will be replaced with ~/databases(at some point that will be configurable)

When the pipeline is installed

.. code-block:: bash

   #> python setup.py install

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
