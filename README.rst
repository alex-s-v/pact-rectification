========
Overview
========

.. start-badges

.. list-table::
    :stub-columns: 1

    * - docs
      - |docs|
    * - tests
      - | |travis| |requires|
        | |codecov|
    * - package
      - | |version| |wheel| |supported-versions| |supported-implementations|
        | |commits-since|
.. |docs| image:: https://readthedocs.org/projects/pact-rectification/badge/?style=flat
    :target: https://readthedocs.org/projects/pact-rectification
    :alt: Documentation Status

.. |travis| image:: https://api.travis-ci.org/alex-s-v/pact-rectification.svg?branch=master
    :alt: Travis-CI Build Status
    :target: https://travis-ci.org/alex-s-v/pact-rectification

.. |requires| image:: https://requires.io/github/alex-s-v/pact-rectification/requirements.svg?branch=master
    :alt: Requirements Status
    :target: https://requires.io/github/alex-s-v/pact-rectification/requirements/?branch=master

.. |codecov| image:: https://codecov.io/github/alex-s-v/pact-rectification/coverage.svg?branch=master
    :alt: Coverage Status
    :target: https://codecov.io/github/alex-s-v/pact-rectification

.. |version| image:: https://img.shields.io/pypi/v/rectification.svg
    :alt: PyPI Package latest release
    :target: https://pypi.org/project/rectification

.. |wheel| image:: https://img.shields.io/pypi/wheel/rectification.svg
    :alt: PyPI Wheel
    :target: https://pypi.org/project/rectification

.. |supported-versions| image:: https://img.shields.io/pypi/pyversions/rectification.svg
    :alt: Supported versions
    :target: https://pypi.org/project/rectification

.. |supported-implementations| image:: https://img.shields.io/pypi/implementation/rectification.svg
    :alt: Supported implementations
    :target: https://pypi.org/project/rectification

.. |commits-since| image:: https://img.shields.io/github/commits-since/alex-s-v/pact-rectification/v0.1.0.svg
    :alt: Commits since latest release
    :target: https://github.com/alex-s-v/pact-rectification/compare/v0.1.0...master



.. end-badges

Simple rectification simulation software.

* Free software: MIT license

Installation
============

::

    pip install rectification

You can also install the in-development version with::

    pip install https://github.com/alex-s-v/pact-rectification/archive/master.zip


Documentation
=============


https://pact-rectification.readthedocs.io/


Development
===========

To run the all tests run::

    tox

Note, to combine the coverage data from all the tox environments run:

.. list-table::
    :widths: 10 90
    :stub-columns: 1

    - - Windows
      - ::

            set PYTEST_ADDOPTS=--cov-append
            tox

    - - Other
      - ::

            PYTEST_ADDOPTS=--cov-append tox
