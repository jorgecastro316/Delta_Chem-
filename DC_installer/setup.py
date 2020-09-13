#!/usr/bin/env python

# setup.py.in.distutils

from distutils.core import setup
import platform


if platform.system() == 'Linux':
    doc_dir = '/usr/local/share/doc/deltachem'
else:
    try:
        from win32com.shell import shellcon, shell
        homedir = shell.SHGetFolderPath(0, shellcon.CSIDL_APPDATA, 0, 0)
        appdir = 'deltachem'
        doc_dir = os.path.join(homedir, appdir)
    except:
        pass

long_desc = \
"""
"""

setup(name='deltachem',
      version='1.0',
      author='Jorge L Castro and Thibault Terencio',
      author_email='jcastromda316@gmail.com',
      maintainer='Jorge L Castro',
      maintainer_email='jcastromda316@gmail.com',
      #url='',
      #description='',
      #long_description=long_desc,
      #download_url='',
      #classifiers=[''],
      #platforms=[''],
      #license='',
      install_requires=['numpy >= 1.18.1','matplotlib >= 3.1.3','scipy >= 1.4.1, <= 1.4.99','prettytable >= 0.7.2',
                  'pyfiglet >= 0.8.post0','sympy >= 0.7.3, <=  0.7.99'],
      python_requires='>=3.6',
      packages=['deltachem'],
      package_data={'':['./*.dat','./*.txt'],'deltachem':['test_materials/*.xyz']},
      
     )
