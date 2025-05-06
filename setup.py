from setuptools import setup, find_packages

setup(
    name="cna-inferer",
    version="0.1.0",
    description="Infer copy number alterations (CNAs) from single-cell RNA-seq data",
    author="Chenming Pu",
    author_email="your.email@example.com",
    url="https://github.com/chenming-pu/cna-inferer",
    packages=find_packages(exclude=["tests", "notebooks", "data"]),
    install_requires=[
        "scanpy>=1.9.0",
        "pandas>=1.3.0",
        "numpy>=1.21.0",
        "pyensembl>=1.1.3",
        "gtfparse>=1.2.1",
        "polars>=0.16.0",
        "ruptures>=1.1.6",
        "joblib>=1.1.0"
    ],
    python_requires='>=3.8',
    entry_points={
        'console_scripts': [
            # once you add a cna_inferer/cli.py with a `main()` function:
            'cna-infer=cna_inferer.cli:main'
        ]
    },
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
)
