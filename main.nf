#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.sdrf = ""
params.fasta = ""
params.mzml_dir = ""
params.out_dir = ""
params.max_missed_cleavages = 2
params.max_charge = 4
params.max_threads = 4

process create_xtandem_xmls {
    publishDir "${params.out_dir}/xtandem_params", mode:'copy'

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
    python -m sdrf_convert ${sdrf} tandem -r ${params.mzml_dir} -f ${fasta} -x ./ -o ./ -m ${params.max_missed_cleavages} -c ${params.max_charge} -t ${params.max_threads}
    """
}

process xtandem {
    cpus 4
    publishDir "${params.out_dir}/idents", mode:'copy'

    input:
    path input_conf

    output:
    path "*.t.xml"

    """
    xtandem ${input_conf}
    """
}

process peptide_shaker {
    cpus 4
    publishDir "${params.out_dir}/idents", mode:'copy'

    input:
    path(sdrf)
    path(fasta)
    path(idents)

    output:
    path("./shaker.psdb")

    shell:
    '''
    # Create peptide shaker identifications files
    idparamscli_params="$(python -m sdrf_convert !{sdrf} identificationparameterscli)"
    sh -c "peptide-shaker eu.isas.peptideshaker.cmd.IdentificationParametersCLI ${idparamscli_params} -fasta_target_decoy 1 -fasta_decoy_tag 'DECOY_' -fasta_decoy_type 0 -fasta_decoy_file_tag '-target_decoy' -out ./id_params.par"

    # Move the links into one directory
    mkdir -p ./idents
    mv !{idents} ./idents
    
    # Run peptide shaker
    peptide-shaker eu.isas.peptideshaker.cmd.PeptideShakerCLI -reference sdrf-hupo-showcase -fasta_file !{fasta} -identification_files ./idents -spectrum_files !{params.mzml_dir} -id_params ./id_params.par -out ./shaker.psdb -threads !{params.max_threads}
    '''
}

process flashlfq {
    cpus 4
    publishDir "${params.out_dir}/quant", mode:'copy'

    input:
    path(sdrf)
    path(peptide_shaker_psdb)

    output:
    path("./QuantifiedPeaks.tsv")
    path("./QuantifiedPeptides.tsv")
    path("./QuantifiedProteins.tsv")

    shell:
    '''
    # Convert peptide shaker results to tsv: https://github.com/smith-chem-wisc/FlashLFQ/issues/125
    peptide-shaker eu.isas.peptideshaker.cmd.ReportCLI -in !{peptide_shaker_psdb} -out_reports !{params.out_dir}/idents -reports 11

    # Read tolerance from SDRF
    ppm="$(python -m sdrf_convert !{sdrf} flashlfq)"

    # Run FlashLFQ
    FlashLFQ --idt "!{params.out_dir}/idents/sdrf-hupo-showcase_Extended_PSM_Report.txt" --rep "!{params.mzml_dir}" --out ./ ${ppm} --thr !{params.max_threads}
    '''
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
