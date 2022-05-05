from setuptools import setup


setup(
    name='ligtools',
    version='0.1.0',    
    description='LIGO tools package',
    url='',
    author='Ligo Scientific Collaboration (LSC) and Duang',
    packages=['ligotools'],
    install_requires=['matplotlib',
                      'numpy',   
                      'scipy',
                      'h5py',
                      'json',
                      'IPython'                  
                      ],

    classifiers=[
        'Development Status :: 1 - Planning',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',  
        'Operating System :: POSIX :: Linux',        
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],
)
