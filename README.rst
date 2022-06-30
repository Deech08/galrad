Milky Way Radiation Field and Ionization Photon Flux Tool
---------------------------------------------------------

.. image:: http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat
    :target: http://www.astropy.org
    :alt: Powered by Astropy Badge

This package is meant to provide an easy way to estimate the Ionizing photon flux 
around the CGM environment of Milky Way, including contributions from the 
Magellanic Clouds. In the future this will include an updated radiation field 
field model, but for now provides a simple photon flux calculation based on 
models from Barger et al. (2013) ,Bland-Hawthorn et al. (2019), 
Antwi-Danso et al. (2020). The radiation field SED available for now is from 
Fox et al. (2005), which is assumed to be the radiaiton field shape for the Milky 
Way as well as for the Magellanic Clouds. 

An extragalactic background field can also be included. 

For use with Cloudy, you will need to point your cloudy input file to the radiation 
field file (found in galrad/data/MW.0kpc.dat) as well as the normalization (photon flux provided by this package). 
We recommend including an extragalactic background independently through built in 
Cloudy framework. 


You can install the package using the following command:

    pip install -q git+https://github.com/Deech08/galrad.git

For a quickstart example see the following `Google Colab Notebook <https://colab.research.google.com/drive/1PXw585xXJJIY836WkVHHtCHXUVRAJRI3?usp=sharing>`_





License
-------

This project is Copyright (c) Dhanesh Krishnarao (DK) and licensed under
the terms of the BSD 3-Clause license. This package is based upon
the `Astropy package template <https://github.com/astropy/package-template>`_
which is licensed under the BSD 3-clause license. See the licenses folder for
more information.


Contributing
------------

We love contributions! galrad is open source,
built on open source, and we'd love to have you hang out in our community.

**Imposter syndrome disclaimer**: We want your help. No, really.

There may be a little voice inside your head that is telling you that you're not
ready to be an open source contributor; that your skills aren't nearly good
enough to contribute. What could you possibly offer a project like this one?

We assure you - the little voice in your head is wrong. If you can write code at
all, you can contribute code to open source. Contributing to open source
projects is a fantastic way to advance one's coding skills. Writing perfect code
isn't the measure of a good developer (that would disqualify all of us!); it's
trying to create something, making mistakes, and learning from those
mistakes. That's how we all improve, and we are happy to help others learn.

Being an open source contributor doesn't just mean writing code, either. You can
help out by writing documentation, tests, or even giving feedback about the
project (and yes - that includes giving feedback about the contribution
process). Some of these contributions may be the most valuable to the project as
a whole, because you're coming to the project with fresh eyes, so you can see
the errors and assumptions that seasoned contributors have glossed over.

Note: This disclaimer was originally written by
`Adrienne Lowe <https://github.com/adriennefriend>`_ for a
`PyCon talk <https://www.youtube.com/watch?v=6Uj746j9Heo>`_, and was adapted by
galrad based on its use in the README file for the
`MetPy project <https://github.com/Unidata/MetPy>`_.
