from setuptools import setup, find_packages

setup(
    name="psosp",
    version="1.1.2",
    # Use find_packages() to correctly identify the 'psosp' directory as a package
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
            # Point to the module, which will be executed via -m
            'psosp = psosp.__main__:main',
        ],
    },
    # This is crucial for including non-Python files specified in MANIFEST.in
    include_package_data=True,
    # package_data is an alternative, but MANIFEST.in is generally more robust
    # We can leave this commented out or remove it if MANIFEST.in is correct
    # package_data={
    #     'psosp': ['files/*', 'test/data/*']
    # },
) 