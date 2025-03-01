from distutils.core import setup
setup(
  name = 'bio-IGLoo',
  packages = ['IGLoo', 'IGLoo.scripts'],#, 'IGLoo.materials', 'IGLoo.materials.personalized_ref', 'IGLoo.materials.gene_annotations', 'IGLoo.materials.gene_annotations.GRCh38'],
  version = '0.1.0',
  license='MIT',
  description = 'The toolkits to recover the IG region from LCL dataset.',
  author = 'Mao-Jan Lin',
  author_email = 'mj.maojanlin@gmail.com',
  url = 'https://github.com/maojanlin/IGLoo',
  download_url = 'https://github.com/maojanlin/IGLoo/tarball/master',
  keywords = ['Immunoglobulin', 'IG', 'assembly', 'Lymphoblastoid Cell Line', 'LCL'],
  install_requires=[
          'numpy',
          'pysam',
          'pandas',
          'matplotlib',
          'seaborn',
          'scikit-learn',
          'scipy'
      ],
  include_package_data=True,
  package_data={'IGLoo': ['materials/*', 'materials/gene_annotations/*', 'materials/gene_annotations/GRCh38/*', 'materials/personalized_ref/*', 'scripts/*.sh']},
  zip_safe = False,
  classifiers=[
    'Development Status :: 3 - Alpha',
    'Intended Audience :: Developers',
    'Topic :: Software Development :: Build Tools',
    'License :: OSI Approved :: MIT License',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.6',
    'Programming Language :: Python :: 3.7',
    'Programming Language :: Python :: 3.8',
    'Programming Language :: Python :: 3.9',
  ],
  entry_points={"console_scripts": ["IGLoo = IGLoo.IGLoo:main"],},
)



