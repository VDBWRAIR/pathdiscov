Testing/Development
===================

Tests are composed of the full functional tests and the smaller unit tests

Full Functional
---------------

Functional tests are all located under the tests directory.
Functional tests are tests that run the entire pipeline for you
and compare the output with what is expected from start to finish.

The functional tests should use the smaller rikkcdna database as well
as the testData input files so they run in a matter of a few minutes
each.

They can be run as follows:

.. code-block:: bash

    nosetests -vx tests

Unit Tests
----------

Unit tests are all located under the usamriidPathDiscov/tests directory.
Unit tests are smaller in size and test the smaller logic units of the
pipeline such as individual python functions or the individual pipeline
stages.

Unit tests should run quickly and complete in a matter of seconds.

They can be run as follows:

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
