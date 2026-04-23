
// Get allele details. Mainly exists because if I edited reformat_vars then the pipeline would run all the models again
process GET_ALLELES {
    label 'single_cpu_job'
    
    input:
    path(vcf)
    
    output:
    path("${params.analysis_name}.alleles.tsv"), emit: alleles

    publishDir "$params.resdir/int_files/"

    script:
    """
    bcftools query -f "%ID\t%POS\t%REF\t%ALT\n" ${vcf} > ${params.analysis_name}.alleles.tsv
    """
}

process NOWBIALL_DETAILS {

    label 'small_job'
    
    input:
    path(alleles)
    path(recoding)
    path(rconf)
    
    output:
    path("${params.analysis_name}.multi_now_bi.new_allele_info.tsv"), emit: nowbi_dat

    publishDir "$params.resdir/int_files/"
    
    script:
    """
    class_nowbi.R ${rconf} ${params.analysis_name} ${task.cpus}
    """
}


process GET_PLOTDATA {
    label 'medium_job_increasing'
    
    input:
    path(model_results)
    path(vardetails)
    path(normscores)
    path(vartypes)
    path(alleles)
    path(nowbi_dat)
    path(recodegt)
    path(rconf)

    output:
    path("${params.analysis_name}.snps.tsv.gz"), emit: snpsfile
    path("${params.analysis_name}.non_snps.tsv.gz"), emit: otherfile
    path("${params.analysis_name}.model_warnings.tsv.gz")

    publishDir "$params.resdir/final_stats/", pattern: '*.{tsv.gz}'
    publishDir "$params.resdir/final_stats/", pattern: '*.model_warnings.tsv.gz'


    scratch true

    script:
    """
    get_plotdata.nf.R ${rconf} ${params.analysis_name} ${task.cpus}
    """
}

process GET_PLOTDATA_MULTI {
    label 'small_job'
    
    input:
    path(model_results)
    path(vardetails)
    path(alleles)
    path(normscores)
    path(rconf)

    output:
    path("${params.analysis_name}.multiallelic.tsv.gz"), emit: multifile
    path("${params.analysis_name}.multiallelic.model_warnings.tsv.gz")

    publishDir "$params.resdir/final_stats/"

    //scratch true

    script:
    """
    get_plotdata.multi.nf.R ${rconf} ${params.analysis_name} ${task.cpus}
    """
}


process SIGNIF_FRAGS {
    label 'medium_job_increasing'
    
    input:
    tuple val(chrname), path(filelist)
    path(snpsfile)
    path(otherfile)
    path(normfactors)
    path(rconf)
    
    output:
    path("*signif_frag_scores.txt.gz")
    tuple val(chrname), path("*signif_frag_chunks.txt.gz"), emit: signif_chunks
    path("*signif_frag_positions.txt.gz")

    publishDir "$params.resdir/signif_frags/", pattern: '*.txt.gz'

    script:
    """
    if [ -f chunk.${chrname}_1.tsv.gz ]; then
        gzip -cd ${snpsfile} | gawk -v chr="${chrname}" 'NR == 1 {print; next} \$2 == chr' > ${chrname}.snps.tsv
        get_signif_frags.R ${rconf} ${chrname} ${normfactors} ${params.analysis_name} ${task.cpus}
    else
        touch ${chrname}.${params.analysis_name}.signif_frag_scores.txt.gz
        touch ${chrname}.${params.analysis_name}.signif_frag_chunks.txt.gz
        touch ${chrname}.${params.analysis_name}.signif_frag_positions.txt.gz
    fi
    """
}

process SIGNIF_FRAGS_MULTI {
    label 'medium_job_increasing'
    
    input:
    tuple val(chrname), path(vartable)
    path(multifile)
    path(multichunks)
    path(normfactors)
    path(rconf)

    output:
    path("${chrname}.${params.analysis_name}.multi.signif_frag_scores.txt.gz")
    path("${chrname}.${params.analysis_name}.multi.signif_frag_chunks.txt.gz")
    path("${chrname}.${params.analysis_name}.multi.signif_frag_positions.txt.gz")

    publishDir "$params.resdir/signif_frags_multi/", pattern: '*.txt.gz'

    script:
    """
    get_signif_frags.multi.R ${rconf} ${chrname} ${normfactors} ${params.analysis_name} ${task.cpus}
    """
}

process VAR_PAIRS {
    label 'medium_job_increasing'

    input:
    tuple val(chrname), path(signif_chunks)
    path(snpsfile)
    path(otherfile)
    path(rconf)

    output:
    path("${params.analysis_name}.${chrname}.variant_pairs.tsv.gz")

    publishDir "$params.resdir/variant_pairs/"

    script:
    """
    var_pairs.R ${rconf} ${chrname} ${params.analysis_name} ${task.cpus}
    """
}