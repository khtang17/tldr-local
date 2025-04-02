from setuptools import setup, find_packages

setup(
    name="tldr-local",
    version="0.8.0",
    packages=[
        'tldr_local',  # Add the root package
        'tldr_local.bbfilter',  # Add the submodule 'bbfilter'
        'tldr_local.cluster',   # Add the submodule 'cluster'
        'tldr_local.strain',    # Add the submodule 'strain'
        'tldr_local.helpers',   # Add the submodule 'helpers'
        'tldr_local.fine_tranche',  # Add the submodule 'fine_tranche'
    ],
    install_requires=[  # Add any dependencies your package needs
        # Example: "requests",
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
