=========
Uninstall
=========

The pipeline comes with an uninstaller that you can use to either completely remove the pipeline and its dependencies(except system packages) or it will just uninstall the pipeline from the virtualenv

Partial Uninstall
=================

This is useful if you just want to uninstall the pipeline from the virtualenv. This is nice since you will not have to wait for some of the bigger packages to recompile.

This is especially useful for testing

.. code-block:: bash

    ./uninstall.sh

Full Uninstall
==============

This completely removes all dependencies and the pipeline. Basically, ``rm -rf`` specific directories related to installation

.. code-block:: bash

    ./uninstall.sh -full

**Notes**: 

1. This does not remove the system packages your administrator had to install during the :ref:`install-system-packages`
2. Neither of these uninstall commands will remove the established directory that was cloned during installation so you will still have the usamriidPathDiscov directory. You will have to do ``rm -rf`` on this directory if you want to remove absolutely everything associated with this pipeline.
