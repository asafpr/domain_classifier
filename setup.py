from setuptools import setup
import os.path

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(name='domain_classifier',
      version='0.0.4',
      description='Use PFAM domains to classify DNA or proteins to taxonomic domain',
      author='Asaf Peer',
      author_email='asaf.peer@jax.org',
      license='MIT',
      packages=['domain_classifier'],
      scripts=[
        'bin/build_domains_DB.py',
        'bin/predict_domain.py',
        'bin/filter_kraken_db.py'], 
      install_requires=[
        'biopython',
        'apsw',
        ],
      zip_safe=False,
      url='https://github.com/asafpr/domain_classifier',
      long_description=read("README.md"),
      long_description_content_type="text/markdown",)
