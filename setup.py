from setuptools import setup

setup(
    name="psosp",
    version="1.1.0",  
    py_modules=['psosp'], 
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
    # include_package_data=True tells setuptools to use MANIFEST.in
    include_package_data=True,
    package_data={
        '': ['files/*', 'test/*']
    },
) 
