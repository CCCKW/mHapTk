import re
from setuptools import setup,find_packages
version = '1.3'

# README as the PyPI description: drop the R Markdown header, and point images at GitHub
# since relative paths do not resolve on PyPI
with open('README.md', encoding='utf-8') as fh:
    long_description = fh.read()
long_description = re.sub(r'\A---\n.*?\n---\s*\*\*\*\s*', '', long_description, flags=re.S)
long_description = long_description.replace(
    'src="images/', 'src="https://raw.githubusercontent.com/CCCKW/mHapTk/%s/images/' % version)

setup(name='mHapTk',
      version=version,
      description='mHapTk is a tool kit for analysis of DNA methylation haplotypes. It has 5 sub-commands: tanghulu, stat, genomeWide, R2 and MHBDiscovery',
      long_description=long_description,
      long_description_content_type='text/markdown',
      author='ckw',
      author_email='1353595807@qq.com',
      packages=find_packages(),
      entry_points={
             'console_scripts': ['mhaptk = src.mhaptk:main']
        }
      )
