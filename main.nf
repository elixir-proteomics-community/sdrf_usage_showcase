#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.sdrf = ""
params.fasta = ""
params.mzml_dir = ""
params.out_dir = ""
params.maxed_missed_cleavages = 2
params.max_charge = 4

process create_xtandem_xmls {
    publishDir "${params.out_dir}", mode:'copy'

    input:
    path(sdrf)
    path(fasta)

    output:
    path "*.tandem_input.xml"
    path "*_taxonomy.xml"

    """
    mkdir -p xtandem_confs/
    # Create xtandem input files. All pathes are relative to the current directory
    # as the files will be linked into the next process's working directory
    python -m sdrf_convert ${sdrf} tandem -r ${params.mzml_dir} -f ${fasta} -x ./ -o ./ -m ${params.maxed_missed_cleavages} -c ${params.max_charge}
    """
}

process xtandem {
    publishDir "${params.out_dir}/idents", mode:'copy'

    input:
    path input_conf

    output:
    path "*.mzML.xml"

    """
    xtandem ${input_conf}
    """
}

process peptide_shaker {
    publishDir "${params.out_dir}/idents", mode:'copy'

    input:
    path(sdrf)
    path(fasta)
    path(idents)

    output:
    path("./shaker.psdb")

    """
    # Create peptide shaker identifications files
    idparamscli_params=(\$(python -m sdrf_convert ${sdrf} identificationparameterscli))
    sh -c 'peptide-shaker eu.isas.peptideshaker.cmd.IdentificationParametersCLI "\${idparamscli_params}" -out ./id_params.par'

    # Move the links into one directory
    mkdir -p ./idents
    mv ${idents} ./idents

    # Run peptide shaker
    peptide-shaker eu.isas.peptideshaker.cmd.PeptideShakerCLI -reference sdrf-hupo-showcase -fasta_file ${fasta} -identification_files ./idents -spectrum_files ${params.mzml_dir} -id_params ./id_params.par -out "./shaker.psdb"
    """
}

process flashlfq {
    publishDir "${params.out_dir}/quant", mode:'copy'

    input:
    path(sdrf)
    path(peptide_shaker_psdb)

    output:
    path("./quant/*")

    """
    # Read tolerance from SDRF
    ppm=\$(python -m sdrf_convert ${sdrf} flashlfq)
    # Convert peptide shaker results to tsv: https://github.com/smith-chem-wisc/FlashLFQ/issues/125
    peptide-shaker eu.isas.peptideshaker.cmd.ReportCLI -in ${peptide_shaker_psdb} -out ${params.out_dir}/idents -reports 11
    # Run FlashLFQ
    sh -c 'FlashLFQ --idt "/home/myfolder/msms.txt" --rep "./quant" "\${ppm}"'
    """
}

workflow {
    // Parse input parameters
    sdrf = Channel.fromPath(params.sdrf)
    fasta = Channel.fromPath(params.fasta)
    out_dir = Channel.fromPath(params.out_dir)

    // Run the workflow
    create_xtandem_xmls(sdrf, fasta)
    xtandem(create_xtandem_xmls.out[0].flatten())  
    peptide_shaker(sdrf, fasta, xtandem.out.collect())
    flashlfq(sdrf, peptide_shaker.out)
}
