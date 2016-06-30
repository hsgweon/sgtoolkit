#!/usr/bin/env python

import os
from distutils.core import setup

__version__ = os.environ.get("VERSION", "1.0.1")

setup(name ='sgtoolkit',
      version = __version__,
      description ='Command line utilities for illumina sequences data analysis',
      author='Hyun Soon Gweon',
      author_email='soonio@gmail.com',
      url='https://github.com/hsgweon/sgtoolkit',
      download_url="http://github.com/hsgweon/sgtoolkit/blob/master/dist/sgtoolkit-%s.tar.gz" % __version__,
      packages = ["sgtoolkit"],
      scripts = ['sgtoolkit/sgtk_getreadpairslist.py', 
                 'sgtoolkit/sgtk_prepseqs.py',
                 'sgtoolkit/sgtk_processseqs.py',
                 'sgtoolkit/sgtk_uc2otutable.py',
                 'sgtoolkit/sgtk_getsamplelistfromfasta.py',
                 'sgtoolkit/sgtk_subsampler.py']
     )
