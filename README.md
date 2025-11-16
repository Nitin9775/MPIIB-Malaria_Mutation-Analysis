content = """
# Kelch13 C580Y Mutation Analysis Pipeline (Plasmodium falciparum)
### Author: Nitin Sharma (BTech Biotechnology, NIT Allahabad)
📧 nitin.20230035@mnnit.ac.in

## 🧬 Project Overview
This repository contains a complete, fully automated pipeline for identifying the Kelch13 C580Y drug-resistance mutation in Plasmodium falciparum using real sequencing data.

## 💡 Scientific Background
Artemisinin resistance in Plasmodium falciparum is strongly associated with mutations in the Kelch13 (K13) gene.  
The most important mutation is **C580Y (Cysteine → Tyrosine)**.

| Feature | Value |
|--------|-------|
| Contig | NC_004331.3 |
| Gene | Kelch13 (PfK13) |
| Codon 580 genomic positions | 1725260–1725262 |
| Wild-type codon | ACA → Threonine (T) |
| Mutant codon | TAT/TAC → Tyrosine (Y) |

## 📁 Repository Structure
    
├── pipeline.sh     
├── analyze.py       
├── environment.yml       
├── README.md       
│       
├── data/          
│ └── README.md       
│          
├── example/           
│ ├── example.vcf          
│ └── C580Y_report.txt         
│           
└── docs/            
└── workflow_diagram.png          


## ⚙️ Installation
git clone https://github.com/Nitin9775/MPIIB-Project.git

cd MPIIB-Project
conda env create -f environment.yml
conda activate malaria-k13


## 🚀 Running the Pipeline
Place FASTQ files into `data/`.

Run:
bash pipeline.sh

## 🔬 Mutation Detection
python3 analyze.py round_1.vcf Pfalciparum_3D7.fasta

makefile
Produces:
C580Y_report.txt


## 🧰 Tools Used
- minimap2
- samtools
- bcftools
- sra-tools
- medaka
- Python 3

## 👨‍🎓 Author
Nitin Sharma  
📧 nitin.20230035@mnnit.ac.in  
GitHub: https://github.com/Nitin9775

## 📄 License
MIT License © 2025 Nitin Sharma
"""

with open("/mnt/data/README.md", "w") as f:
    f.write(content)

"/mnt/data/README.md created"
