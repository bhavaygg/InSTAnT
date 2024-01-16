
from distutils.core import setup
setup(
  name = 'sc-instant',
  packages=['InSTAnT'],
    package_dir={
        'InSTAnT': 'InSTAnT',
    },
  version = '1.0.2',      
  license='MIT',        
  description = 'InSTAnT is a toolkit to identify gene pairs which are d-colocalized from single molecule measurement data.',   
  author = 'Anurendra Kumar, Bhavay Aggarwal',                   
  author_email = 'anu.ankesh@gmail.com, bhavayaggarwal07@gmail.com',      
  url = 'https://github.com/anurendra, https://github.com/Chokerino',   
  keywords = ['InSTAnT', 'Python'],   
  install_requires=[            
          'numpy',
          'scipy',
          'pandas',
          'pickle',
          'scikit-learn',
      ],
  classifiers=[
    'Development Status :: 3 - Alpha',      
    'Intended Audience :: Developers',      
    'Topic :: Software Development :: Build Tools',
    'License :: OSI Approved :: MIT License',   
    'Programming Language :: Python :: 3',      
    'Programming Language :: Python :: 3.4',
    'Programming Language :: Python :: 3.5',
    'Programming Language :: Python :: 3.6',
  ],
  long_description=open('README.md').read(),
  long_description_content_type='text/markdown',
)
