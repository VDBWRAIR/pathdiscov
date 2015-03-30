Testing/Devlopment
==================

Tests are composed of the full functional tests and the smaller unit tests

Full Functional
---------------

Functional tests are all located under the tests directory and can be run
as follows:

.. code-block:: bash

    nosetests -vx tests

Unit Tests
----------

Unit tests are all located under the usamriidPathDiscov/tests directory
and can be run as follows:

.. code-block:: bash

    nosetests -vx usamriidPathDiscov/tests

Setup
-----

#. Install pipeline as normal
#. Follow `<rikkcdna/Readme.rst>`_ to setup test databases
#. Run nosetests

    The following will run all test scripts that start with ``test``

    .. code-block:: bash

        $> nosetests -vx
