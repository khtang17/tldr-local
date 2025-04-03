from setuptools import setup, find_packages

setup(
    name="tldr-local",
    version="0.8.0",
    packages=[
        'tldr_local', 
        'tldr_local.bbfilter',  
        'tldr_local.cluster',  
        'tldr_local.strain',    
        'tldr_local.helpers',   
        'tldr_local.fine_tranche', 
    ],
    install_requires=[  
       "tqdm",
       "rdkit-pypi",
    ],
    entry_points={  # Define the console scripts here
        'console_scripts': [
            'bb_filter = tldr_local.bbfilter.main:main',  
            'cluster = tldr_local.cluster.main:main',
            'strain = tldr_local.strain.main:main',
            'fine_tranche = tldr_local.fine_tranche.main:main'
        ],
    },
    include_package_data=True, 
    author="Your Name",
    author_email="your.email@example.com",
    description="A short description of your package",
    long_description=open('README.md').read(),
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/tldr-local",
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],

)
