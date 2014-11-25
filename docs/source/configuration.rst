Configuration
=============

usamriidPathDiscov/files/config.yml.base
----------------------------------------

This file allows you to set the defaults that define where your databases are located
under the databases directory.
GENOMEDIR will be replaced with ~/databases(at some point that will be configurable)

When the pipeline is installed

.. code-block:: bash

   #> python setup.py install

usamriidPathDiscov/files/config.yaml.base is copied to usamriidPathDiscov/files/config.yaml
and GENOMEDIR is replaced with ~/databases
Then when the pipeline is run with usamriidPathDiscov_cli that config.yaml is used to modify
usamriidPathDiscov/files/sample.param.base which is copied to the python installation under
usamriidPathDiscov/lib/python*/site-packages/usamriidPathDiscov/files/sample.param
