from setuptools import setup

setup(name='domain_classifier',
      version='0.0.1',
      description='Use PFAM domains to classify DNA or proteins to taxonomic domain',
      author='Asaf Peer',
      author_email='asaf.peer@jax.org',
      license='MIT',
      packages=['domain_classifier'],
      scripts=[
        'bin/build_database.py',
        'bin/predict_domain.py',
        'bin/filter_kraken_db.py'], 
      install_requires=[
        'biopython',
        'sqlite3',
        ],
      zip_safe=False)
