from distutils.core import setup

setup(
        name='PYACDCPF',
        version='1.0',
        author='Roni Irnawan',
        author_email='roni.irnawan@gmail.com',
        description='Runs a sequential ac/dc power flow.',
        url='https://github.com/rwl/PYACDCPF',
        packages=['pyacdcpf'],
        package_data={
              'pyacdcpf': ['Cases/PowerflowAC/*.py','Cases/PowerflowDC/*.py','*.html','*.ipynb']
                },
        classifiers=[
            'Development Status :: 1 - Production/Stable',
            'Environment :: Console',
            'Intended Audience :: Developers',
            'Intended Audience :: Education',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: BSD License',
            'Natural Language :: English',
            'Operating System :: OS Independent',
            'Programming Language :: Python',
            'Programming Language :: Python :: 3.3',
            'Programming Language :: Python :: 3.4',
            'Programming Language :: Python :: 3.5',
            'Topic :: Scientific/Engineering',
        ],
      
      )