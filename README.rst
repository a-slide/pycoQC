
pycoQC
======

**Python 3 package for Jupyter Notebook, computing metrics and
generating simple QC plots from Oxford Nanopore technologies (ONT)
Albacore basecaller **

Installation
------------

.. code:: python

    %%bash
    get 

Dependencies
~~~~~~~~~~~~




Usage
-----

.. code:: python

    import pylab as pl # Namespace containing nupy + matplotlib
    %pylab inline
    pl.rcParams['figure.figsize'] = 20, 7
    pl.rcParams['font.family'] = 'sans-serif'
    pl.rcParams['font.sans-serif'] = ['DejaVu Sans']
    pl.style.use('ggplot')


.. parsed-literal::

    Populating the interactive namespace from numpy and matplotlib



.. code:: python

    # Import pycoQC main class
    from pycoQC import pycoQC





.. parsed-literal::

    <unbound method pycoQC.channels_activity>


