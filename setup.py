from setuptools import setup,find_packages
version = '1.3'
setup(name='mHapTk',
      version='1.3',
      description='mHapTk is a tool kit for analysis of DNA methylation haplotypes. It has 5 sub-commands: tanghulu, stat, genomeWide, R2 and MHBDiscovery',
      author='ckw',
      author_email='1353595807@qq.com',
      packages=find_packages(),
      entry_points={
             'console_scripts': ['mhaptk = src.mhaptk:main']
        }
      )

