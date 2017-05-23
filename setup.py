#!/bin/env/python

if __name__ == "__main__":
    from distutils.core import setup
    from glob import glob

    setup(name='mcdiff',
        version='0.001',
        description='Program to find diffusion parameters with Monte Carlo routine.',
        #author='An Ghysels',
        #author_email='An.Ghysels@UGent.be',
        #url='http://molmod.ugent.be/code/',
        package_dir = {'mcdiff': 'lib'},
        packages = ['mcdiff','mcdiff.tools','mcdiff.permeability'],
        scripts=glob("scripts/*"),
        classifiers=[
            #'Development Status :: 3 - Alpha',
            #'Environment :: Console',
            #'Intended Audience :: Science/Research',
            #'License :: OSI Approved :: GNU General Public License (GPL)',
            #'Operating System :: POSIX :: Linux',
            'Programming Language :: Python',
        ],
    )

