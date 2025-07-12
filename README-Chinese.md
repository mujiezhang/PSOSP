# PSOSP: 原噬菌体SOS依赖性预测器
PSOSP (<b>P</b>rophage <b>SOS</b> dependency <b>P</b>redictor) 通过分析LexA蛋白与靶DNA的结合异质性指数（_HI_），预测原噬菌体的诱导模式，将原噬菌体分类为SOS依赖型（SdPs）和SOS非依赖型（SiPs）。

![PSOSP-LOGO](https://github.com/user-attachments/assets/3f02c2e6-5596-4293-970f-45e60139e6a9)

![GitHub Release](https://img.shields.io/github/v/release/mujiezhang/psosp?style=flat)
[![Conda Downloads](https://img.shields.io/conda/dn/bioconda/psosp?style=flat&label=Bioconda%20downloads&color=%2397CA00)](https://anaconda.org/bioconda/psosp)
[![GitHub License](https://img.shields.io/github/license/mujiezhang/psosp?style=flat)](https://anaconda.org/bioconda/psosp)
[![conda release](https://anaconda.org/bioconda/psosp/badges/version.svg)](https://anaconda.org/bioconda/psosp)

## 目录
<!-- TOC -->
- [PSOSP: 原噬菌体SOS依赖性预测器](#psosp-原噬菌体sos依赖性预测器)
- [Webserver](#webserver)
- [背景](#背景)
  - [基本原理](#基本原理)
  - [工作流程](#工作流程)
  - [实验验证](#实验验证)
  - [对噬菌体分离的意义](#对噬菌体分离的意义)
- [软件使用](#软件使用)
  - [输入要求](#输入要求)
  - [依赖](#依赖)
  - [**安装**](#安装)
  - [**输入文件**](#输入文件)
  - [**如何运行**](#如何运行)
  - [**输出结果**](#输出结果)
- [引用](#引用)

<!-- /TOC -->

## Webserver
[![Total Predictions](https://img.shields.io/badge/dynamic/json?color=brightgreen&label=Online%20predictions&query=$.total_predictions&url=https%3A%2F%2Fvee-lab.sjtu.edu.cn%2FPSOSP%2Fbadge.php)](https://vee-lab.sjtu.edu.cn/PSOSP/)

我们提供了一个在线平台（PSOSP）用于快速预测噬菌体诱导模式：**https://vee-lab.sjtu.edu.cn/PSOSP**

## 背景

### 基本原理
温和噬菌体作为前噬菌体整合到细菌宿主基因组中。在正常条件下，LexA蛋白与前噬菌体中的SOS box结合，抑制噬菌体相关基因的表达，维持溶原状态。当受到外部刺激（如暴露于DNA损伤剂）时，RecA蛋白被激活，导致LexA自切割并从SOS box上解离。这解除了对前噬菌体的抑制，触发温和噬菌体进入裂解周期，从而促进其复制。

![psosp-theory](https://github.com/user-attachments/assets/654a77e1-dbb6-44bb-9719-0fe4fca7519c)

### 工作流程
- 鉴定LexA和经典SOS box（CBS）：扫描宿主基因组以鉴定LexA蛋白和位于lexA基因上游的经典SOS box（CSBs）
  
- 异质性指数（_HI_）计算：识别细菌基因组中的PSBs (Potential SOS Boxes)，计算每个PSB的异质性指数（_HI_），并通过Mean Shift聚类结果建立分类阈值（_HI<sub>c1</sub>_ 和 _HI<sub>c2</sub>_）
  
- 在前噬菌体中扫描PSB：在前噬菌体启动子区域扫描PSBs并确定最小 _HI_（_HI<sub>min</sub>_）
  
- 前噬菌体分类：通过比较 _HI<sub>min</sub>_ 与阈值评估LexA与前噬菌体启动子区域的结合能力
   - _HI<sub>min</sub>_ ≤ _HI<sub>c1</sub>_ → **SdP**（SOS依赖型前噬菌体）
   - _HI<sub>min</sub>_ ≥ _HI<sub>c2</sub>_ → **SiP**（SOS非依赖型前噬菌体）
   - _HI<sub>c1</sub>_ < _HI<sub>min</sub>_ < _HI<sub>c2</sub>_ → **SuP**（SOS不确定型前噬菌体）

![PSOSP_workflow](https://github.com/user-attachments/assets/c3a0334f-ce4f-4533-960a-ae8c19b71514)

### 实验验证
我们使用14个经过实验验证的噬菌体（囊括10个病毒科，其中2个属于*Peduoviridae*、3个属于*Inoviridae*，另外9个属于不同的新科）验证PSOSP的准确性，其宿主覆盖7个细菌属（*Salmonella、Escherichia、Vibrio、Pseudomonas、Serratia_J、Hafnia*和*Shewanella*），囊括3个细菌目（Enterobacterales、Enterobacterales_A和Pseudomonadales）。值得注意的是，PSOSP对这些噬菌体的所有预测与实验证据完全一致，证明了该工具的通用性和可靠性。

![experiment_validation](https://github.com/user-attachments/assets/f39bc3c6-a18b-4bf4-9459-d14176d76289)

### 对噬菌体分离的意义
我们建议未来的噬菌体分离工作可以首先使用PSOSP来确定噬菌体类型
  - 对于<strong>SdPs</strong>，传统的SOS诱导剂（如MMC、UV）仍然适用。
  - 对于<strong>SiPs</strong>，应考虑SOS非依赖型诱导剂，如DPO、C4-HSL、EDTA和pyocyanin，或物理因素，如不同的盐度、温度和pH。

  
## 软件使用
### 输入要求
**对于宿主**：
- 宿主类别：**PSOSP主要适用于Gammaproteobacteria**，包括*Vibrio cholerae, Pseudomonas aeruginosa, Yersinia pestis, Escherichia coli, Salmonella enterica, Shigella spp*和 *Klebsiella spp*在内的多类重要致病菌，可以通过PSOSP在线网站的[**Statistics**](https://vee-lab.sjtu.edu.cn/PSOSP/statistics.html)页面查看适用的细菌属。
- 基因组质量：**建议使用完整度高于90%的宿主基因组**，因为低质量的基因组可能会丢失LexA蛋白降低预测准确度。可以使用[**CheckM2**](https://github.com/chklovski/CheckM2)评估宿主基因组的完整度。
- 多contig基因组：如果宿主基因组由多个contig组成，请确保输入的宿主基因组文件包含所有contigs。

**对于前噬菌体**：

- 基因组质量：PSOSP使用CheckV进行质量评估。**对于CheckV评估的完整度>90%的噬菌体，PSOSP预测结果比较可信；若完整度低于90%，则预测准确度降低**。如果确定输入的噬菌体基因组是完整的，可以忽略输出文件中的CheckV结果（有些情况下CheckV评估的完整度不够准确）。
- 多个输入：输入的噬菌体基因组文件可以包含多个前噬菌体序列。
- 宿主关联：输入的噬菌体必须是整合在对应宿主基因组内或者能侵染对应宿主的噬菌体。预测不匹配的噬菌体-宿主之间的调控关系没有意义。

### 依赖
* PSOSP是一个Python脚本，依赖于：
```
DIAMOND=2.1.8
MEME=5.5.5
Python=3.12
scikit-learn=1.6.1
Prodigal=2.6.3
biopython=1.85
checkv=1.0.3
```

### 安装
(1) conda (**推荐，最简单的安装方式**)

  - 安装conda并添加官方镜像 (如已安装, 请跳过)
    ```
    wget -c https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh
    source ~/.bashrc
    conda config --add channels bioconda
    conda config --add channels conda-forge
    ```
  - 安装 PSOSP
    ```
    conda create -n psosp psosp
    conda activate psosp
    ```
    测试安装是否成功:`psosp test`
    使用方法: ```psosp -h```

  - 如果需要CheckV结果，请下载CheckV数据库
    ```
    checkv download_database ./
    ```

(2) git (**请先安装上述依赖**)
  ```
  git clone https://github.com/mujiezhang/PSOSP.git
  cd PSOSP
  pip install -e .
  ```
  测试安装是否成功:`psosp test`


### 输入文件
PSOSP需要两个文件作为输入：
* ```-hf```：宿主基因组（fasta格式）
* ```-vf```：原噬菌体基因组（fasta格式）

其他参数：
* ```-wd```：保存结果文件的目录
* ```-faa```：宿主蛋白质序列（fasta格式，可选）
* ```-db```：checkv参考数据库路径（可选）

### 如何运行
* 通过Conda安装
```
psosp -hf /path/to/host-genome.fasta -vf /path/to/virus-genome.fasta -wd output_dir -db /path/to/checkv-db
```

用 [**github**](https://github.com/mujiezhang/PSOSP/tree/main/psosp/test/) 或者 [**zenodo**](https://zenodo.org/records/15795217/files/test.zip) 中的示例数据测试:
```
psosp -hf test/data/host_wp2.fna -vf test/data/virus_wp2-phage-sp1-sp2-sp3.fna -wd test/test-result -db /path/to/checkv-db
```


### 输出结果
在此示例中，PSOSP的分析结果将写入`test/test-result`目录，其结构如下：
```
test/test-result
├── virus_wp2-phage-sp1-sp2-sp3_checkv
├── host_wp2.fna.faa_lexa_blast.tsv
├── host_wp2.fna_whole_genome_HI.tsv
├── host_wp2_prodigal.faa
╰── virus_wp2-phage-sp1-sp2-sp3_prediction.tsv
```

1. `virus_wp2-phage-sp1-sp2-sp3_checkv`: 此目录包含[**CheckV**](https://bitbucket.org/berkeleylab/checkv/src/master/)的结果
2. `host_wp2.fna.faa_lexa_blast.tsv`: LexA蛋白的blast结果
3. `host_wp2.fna_whole_genome_HI.tsv`: 使用MeanShift方法对宿主基因组中所有潜在SOS box的HI聚类结果
4. `host_wp2_prodigal.faa`: 由prodigal产生的蛋白质序列
5. `virus_wp2-phage-sp1-sp2-sp3_prediction.tsv`: **PSOSP的预测结果**

**`virus_wp2-phage-sp1-sp2-sp3_prediction.tsv`**文件内容:

|                **host**                |                **virus**                |                **prediction_result**                |                **prediction_quality**                |                **completeness**                |                **contamination**                |                **viral-HI(min)**                |                **box-seq**                |                **box-seq_start_pos**                |                **box-seq_strand**                |                **confidence_window_lower**                |                **confidence_window_upper**                |                **blast_status**                |                **fimo_status**                |
|:--------------------------------------:|:----------------------------------------:|:------------------------------------------------:|:------------------------------------------------:|:--------------------------------------------:|:--------------------------------------------:|:--------------------------------------------:|:----------------------------------------:|:------------------------------------------------:|:--------------------------------------------:|:----------------------------------------------------:|:----------------------------------------------------:|:--------------------------------------------:|:------------------------------------------:|
|               host_wp2.fna              |                    sp1                   |          SiP (SOS-independent Prophage)          |                      High                         |                    100.0                     |                     0.0                      |                   14.2017                    |         CACTGTATTATTATACCACA              |                     29652                      |                      -                       |                      11.8522                       |                      13.3752                       |                   Blast_OK                   |                   Fimo_OK                  |
|               host_wp2.fna              |                    sp2                   |           SdP (SOS-dependent Prophage)           |                      High                         |                    100.0                     |                     0.0                      |                   11.7239                    |         TTATGTATGTATATTCAGCA              |                     17757                      |                      -                       |                      11.8522                       |                      13.3752                       |                   Blast_OK                   |                   Fimo_OK                  |
|               host_wp2.fna              |                    sp3                   |          SiP (SOS-independent Prophage)          |                      High                         |                    99.62                     |                     0.0                      |                   15.4715                    |         CACTGTATAAAAAAACATAC              |                     1225                       |                      -                       |                      11.8522                       |                      13.3752                       |                   Blast_OK                   |                   Fimo_OK                  |

- `host`: 输入的宿主文件名
- `virus`: 输入FASTA文件中噬菌体名称
- `prediction_result`: PSOSP预测的诱导模式. **SiPs**：SOS非依赖型前噬菌体；**SuPs**：SOS不确定型前噬菌体；**SdPs**：SOS依赖型前噬菌体
- `prediction quality`: 完整度在90%-100%之间或预测为SdPs的噬菌体为**`High`**；完整度在50%-90%之间为**`Medium`**；完整度低于50%为**`Low`**
- `completeness`: 噬菌体完整度（CheckV）
- `contamination`: 噬菌体污染度（CheckV）
- `viral-HI(min)`: 噬菌体基因组中的最小 _HI_
- `box-seq`: 最小 _HI_ 对应的SOS box的序列
- `box-seq_start_pos`: 最小 _HI_ 对应的SOS box的序列在噬菌体上的起始位置
- `box-seq_strand`: +（正链）或-（负链）
- `confidence_window_lower`: 阈值 _HI<sub>c1</sub>_（_HI<sub>min</sub>_ ≤ _HI<sub>c1</sub>_ → `SdP`）
- `confidence_window_upper`: 阈值 _HI<sub>c2</sub>_（_HI<sub>min</sub>_ ≥ _HI<sub>c2</sub>_ → `SiP`）
- `blast_status`: 'Blast_OK'（找到LexA同源物）或"-"（无同源物）
- `fimo_status`: 'Fimo_OK'（在LexA上游检测到SOS box）或"-"（未检测到SOS box）

## 引用
......
