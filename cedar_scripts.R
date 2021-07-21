##### mkref script #####
#!/bin/bash
#SBATCH --mail-user=cara.e.ellis@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --time=3:00:00
#SBATCH --account=def-ubcdrg
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G

#script for creating a reference transcriptome
echo "Running on hostname `hostname`"
echo "Starting run at: `date`"

export PATH=/home/caraee/cellranger-4.0.0:$PATH

wget ftp://ftp.ensembl.org/pub/release-101/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
gunzip Mus_musculus.GRCm38.dna.primary_assembly.fa.gz

wget ftp://ftp.ensembl.org/pub/release-101/gtf/mus_musculus/Mus_musculus.GRCm38.101.gtf.gz
gunzip Mus_musculus.GRCm38.101.gtf.gz

cellranger mkgtf Mus_musculus.GRCm38.101.gtf Mus_musculus.GRCm38.101.filtered.gtf \
--attribute=gene_biotype:protein_coding \
--attribute=gene_biotype:lincRNA \
--attribute=gene_biotype:antisense \
--attribute=gene_biotype:IG_LV_gene \
--attribute=gene_biotype:IG_V_gene \
--attribute=gene_biotype:IG_V_pseudogene \
--attribute=gene_biotype:IG_D_gene \
--attribute=gene_biotype:IG_J_gene \
--attribute=gene_biotype:IG_J_pseudogene \
--attribute=gene_biotype:IG_C_gene \
--attribute=gene_biotype:IG_C_pseudogene \
--attribute=gene_biotype:TR_V_gene \
--attribute=gene_biotype:TR_V_pseudogene \
--attribute=gene_biotype:TR_D_gene \
--attribute=gene_biotype:TR_J_gene \
--attribute=gene_biotype:TR_J_pseudogene \
--attribute=gene_biotype:TR_C_gene

cellranger mkref --genome=mm10 \
--fasta=Mus_musculus.GRCm38.dna.primary_assembly.fa \
--genes=Mus_musculus.GRCm38.101.filtered.gtf \
--ref-version=GRCm38.p6 \
--nthreads=6

echo " Job finished with exit code $? at: `date`"

##### trimgalore script #####
#!/bin/bash
#SBATCH --mail-user=cara.e.ellis@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --time=8:00:00
#SBATCH --account=def-ubcdrg
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G

#script for trimming  Fluidigm fastq files with cutadapt

module load fastqc/0.11.9
module load python/3.8.2

virtualenv --no-download $SLURM_TMPDIR/env
source $SLURM_TMPDIR/env/bin/activate
pip install --no-index --upgrade pip
pip install --no-index --upgrade cutadapt

echo "Running on hostname `hostname`"
echo "Starting run at: `date`"

cd /home/caraee/scratch/CE282_scRNAseq/SingleCell/Kieffer/200824_NS500668_0851_AHVJF5BGXF/SC200825fastqs/outs/fastq_path/SC67

for i in *.fastq.gz; do
cutadapt --no-indels -e 0.1 -g ^AAGCAGTGGTATCAACGCAGAGTACATGGG \
-j 4 -o /home/caraee/scratch/CE282_scRNAseq/trimmed/CE282_621666/ ${i}

done

echo " Job finished with exit code $? at: `date`"

cd /home/caraee/scratch/CE282_scRNAseq/SingleCell/Kieffer/200825_NS500668_0852_AHVWJKBGXF/SC200825fastqs/outs/fastq_path/SC67

for i in *.fastq.gz; do
cutadapt --no-indels -e 0.1 -g ^AAGCAGTGGTATCAACGCAGAGTACATGGG \
-j 4 -o /home/caraee/scratch/CE282_scRNAseq/trimmed/CE282_621666/ ${i}

done



##### count script #####
#!/bin/bash
#SBATCH --mail-user=cara.e.ellis@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --time=24:00:00
#SBATCH --account=def-ubcdrg
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G

#script for creating a reference transcriptome
echo "Running on hostname `hostname`"
echo "Starting run at: `date`"

export PATH=/home/caraee/cellranger-4.0.0:$PATH

cellranger count --id=run_count_CE282_621665 \
--fastqs=/home/caraee/scratch/CE282_scRNAseq/SingleCell/Kieffer/200824_NS500668_0851_AHVJF5BGXF/SC200825fastqs.trimmed/SC67/ \
--sample=CE282_621665 \
--transcriptome=/scratch/caraee/CE282_scRNAseq/SingleCell/Kieffer/mm10

count --id=SC67__CE282_621665 --transcriptome=/brcwork/sequence/10x_data/refdata-cellranger-mm10-2.1.0/ 
  --fastqs=SC200825fastqs.trimmed/SC67 --sample=CE282_621665 --jobmode=sge --maxjobs=100


echo " Job finished with exit code $? at: `date`"