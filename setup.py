from setuptools import setup, find_packages

setup(
    name="psosp",
    version="1.1.0",  
    packages=find_packages(),
    author="Mujie Zhang",  
    author_email="3224413991@qq.com",  
    description="PSOSP (Prophage SOS-dependency Predictor)",
    long_description=open('README.md', encoding='utf-8').read(),
    long_description_content_type="text/markdown",
    url="https://github.com/mujiezhang/PSOSP",  
    license="MIT",  
    entry_points={
        'console_scripts': [
            'psosp = psosp:main',
        ],
    },
    
    package_data={
        'psosp': ['files/19-motifs-meme.txt', 'files/uniprot_swiss_prot_LexA.dmnd'],
    },
    include_package_data=True,
) 
