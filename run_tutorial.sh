main(){
    readonly READEMPTION=reademption
    readonly READEMPTION_ANALYSIS_FOLDER=READemption_analysis

    create_project
    store_environment_variable
    download_fasta
    modify_fasta_header
    download_annotation
    download_reads
    align
    #convert_bam_to_sam
    coverage
    gene_quanti
    deseq
    viz_align
    viz_gene_quanti
    viz_deseq
}
create_project(){
    ${READEMPTION} create --project_path ${READEMPTION_ANALYSIS_FOLDER}
}

store_environment_variable(){
    FTP_SOURCE=ftp://ftp.ncbi.nih.gov/genomes/archive/old_refseq/Bacteria/Salmonella_enterica_serovar_Typhimurium_SL1344_uid86645/
}

download_fasta(){
    wget -O READemption_analysis/input/reference_sequences/NC_016810.fa $FTP_SOURCE/NC_016810.fna
    wget -O READemption_analysis/input/reference_sequences/NC_017718.fa $FTP_SOURCE/NC_017718.fna
    wget -O READemption_analysis/input/reference_sequences/NC_017719.fa $FTP_SOURCE/NC_017719.fna
    wget -O READemption_analysis/input/reference_sequences/NC_017720.fa $FTP_SOURCE/NC_017720.fna
}

modify_fasta_header(){
    sed -i "s/>/>NC_016810.1 /" READemption_analysis/input/reference_sequences/NC_016810.fa
    sed -i "s/>/>NC_017718.1 /" READemption_analysis/input/reference_sequences/NC_017718.fa
    sed -i "s/>/>NC_017719.1 /" READemption_analysis/input/reference_sequences/NC_017719.fa
    sed -i "s/>/>NC_017720.1 /" READemption_analysis/input/reference_sequences/NC_017720.fa
}

download_annotation(){
    wget -P READemption_analysis/input/annotations ${FTP_SOURCE}/*gff
}

download_reads(){
    wget -P READemption_analysis/input/reads http://reademptiondata.imib-zinf.net/InSPI2_R1.fa.bz2
    wget -P READemption_analysis/input/reads http://reademptiondata.imib-zinf.net/InSPI2_R2.fa.bz2
    wget -P READemption_analysis/input/reads http://reademptiondata.imib-zinf.net/LSP_R1.fa.bz2
    wget -P READemption_analysis/input/reads http://reademptiondata.imib-zinf.net/LSP_R2.fa.bz2
}

align(){
    ${READEMPTION} align -p 4 --poly_a_clipping -f ${READEMPTION_ANALYSIS_FOLDER}
}

convert_bam_to_sam(){
    samtools view -h -o InSPI2_R1_alignments_final.sam READemption_analysis/output/align/alignments/InSPI2_R1_alignments_final.bam
    samtools view -h -o InSPI2_R2_alignments_final.sam READemption_analysis/output/align/alignments/InSPI2_R1_alignments_final.bam
    samtools view -h -o  LSP_R1_alignments_final.sam READemption_analysis/output/align/alignments/LSP_R1_alignments_final.bam
    samtools view -h -o  LSP_R2_alignments_final.sam READemption_analysis/output/align/alignments/LSP_R2_alignments_final.bam
 	     
}

coverage(){
    ${READEMPTION} coverage -p 4 -f ${READEMPTION_ANALYSIS_FOLDER}
}

gene_quanti(){
    ${READEMPTION} gene_quanti -p 4 --features CDS,tRNA,rRNA -f ${READEMPTION_ANALYSIS_FOLDER}
}

deseq(){
    ${READEMPTION} deseq \
   -l InSPI2_R1.fa.bz2,InSPI2_R2.fa.bz2,LSP_R1.fa.bz2,LSP_R2.fa.bz2 \
   -c InSPI2,InSPI2,LSP,LSP -f ${READEMPTION_ANALYSIS_FOLDER}
}

viz_align(){
    ${READEMPTION} viz_align -f ${READEMPTION_ANALYSIS_FOLDER}
}

viz_gene_quanti(){
    ${READEMPTION} viz_gene_quanti -f ${READEMPTION_ANALYSIS_FOLDER}
}

viz_deseq(){
    ${READEMPTION} viz_deseq -f ${READEMPTION_ANALYSIS_FOLDER}
}

main
