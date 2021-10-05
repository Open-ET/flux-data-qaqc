
Automated testing with `pytest <https://docs.pytest.org/en/6.2.x/contents.html#>`__
==========================================================================================

Software tests are automatically run each time a change to ``flux-data-qaqc`` is made on the master branch in GitHub using this `GitHub Actions workflow <https://github.com/Open-ET/flux-data-qaqc/actions/workflows/fluxdataqaqc_tests.yml>`__.  Automated tests help spot potential bugs early so that they be identified and corrected efficiently resulting in an improved user experience.  

Running tests manually
^^^^^^^^^^^^^^^^^^^^^^

``pytest`` is required to run software tests that are provided. You can install ``pytest`` with PIP:

.. code-block:: bash

    pip install pytest

The tests utilize the example flux input data and ``flux-data-qaqc`` configuration files that are provided with the software whether installed from PyPI or GitHub. These files can be found `here <https://github.com/Open-ET/flux-data-qaqc/tree/master/examples>`__.

To run the tests, navivgate to the root directory of the source code (from the command line or shell) and run pytest:

.. code-block:: bash

    pytest
    
This will print out basic test results, usage of ``pytest`` plugins and command line options can be used for getting more information out of the tests.

