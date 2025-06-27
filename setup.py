from setuptools import setup, find_packages

setup(
    name="psosp",
    version="1.1.0",  # 确保这个版本号和您希望发布的版本一致
    packages=find_packages(),
    author="Mujie Zhang",  # 您的名字
    author_email="3224413991@qq.com",  # 您的邮箱
    description="PSOSP (Prophage SOS-dependency Predictor)",
    long_description=open('README.md', encoding='utf-8').read(),
    long_description_content_type="text/markdown",
    url="https://github.com/mujiezhang/PSOSP",  # 您的项目GitHub链接
    license="MIT",  # 建议选择一个开源许可证，如MIT
    entry_points={
        'console_scripts': [
            'psosp = psosp:main',
        ],
    },
    # 这一步很关键，用于包含非代码文件
    package_data={
        'psosp': ['files/19-motifs-meme.txt', 'files/uniprot_swiss_prot_LexA.dmnd'],
    },
    include_package_data=True,
) 