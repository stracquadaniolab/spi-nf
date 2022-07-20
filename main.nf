// enabling nextflow DSL v2
nextflow.enable.dsl=2

// define some parameters shared by all processes
timestampName = (new Date().format("yyyy-MM-dd")) + "-${params.name}"
params.resultsDir = params.outdir + "/" + timestampName

motd = """
--------------------------------------------------------------------------
spi-nf ($workflow.manifest.version)
--------------------------------------------------------------------------
Session ID   : $workflow.sessionId
--------------------------------------------------------------------------
Parameters
--------------------------------------------------------------------------
Genome structure            : $params.genome.structure
Genome reconstructions      : $params.genome.reconstruction
Lambda                      : $params.model.lambda
Nu                          : $params.model.nu
B                           : $params.model.b
Max configs                 : $params.model.configs
Results dir                 : $params.resultsDir
Reweighting                 : $params.model.reweight.run
Reweighting trajectories    : $params.model.reweight.trajectories
Reweighting lambda          : $params.model.reweight.lambda [$params.model.reweight.lambda_radius]
Reweighting nu              : $params.model.reweight.nu [$params.model.reweight.nu_radius]
Reweigthing b               : $params.model.reweight.b  [$params.model.reweight.b_radius]
--------------------------------------------------------------------------
Environment information
--------------------------------------------------------------------------
Container                   : $workflow.container
Config files                : $workflow.configFiles
Project dir                 : $workflow.projectDir
Work dir                    : $workflow.workDir
Launch dir                  : $workflow.launchDir
Command line                : $workflow.commandLine
Repository                  : $workflow.repository
CommitID                    : $workflow.commitId
Revision                    : $workflow.revision
--------------------------------------------------------------------------
"""

log.info motd

process MODEL_SAMPLING {

    input:
    tuple path(genome_structure), val(lam), val(nu), val(b), val(max_config)

    output:
    tuple val(lam), val(nu), val(b), val(max_config), path('profile.txt'), path('trajectories.txt'), path('run-info.txt')

    script:
    """
    spi-sample.py -i ${genome_structure} --lb ${lam} --nu ${nu} --b ${b} -m ${max_config} -o profile.txt -t trajectories.txt --run-info run-info.txt
    """

    stub:
    """
    touch profile.txt
    touch trajectories.txt
    touch run-info.txt
    """

}

process MODEL_REWEIGHTING {

     input:
     tuple path(genome_structure), path(trajectories), val(lam), val(nu), val(b)

     output:
     path('profile.txt')
     
     script:
     """
     spi-reweight.py -i ${genome_structure} -t ${trajectories} --lb ${lam} --nu ${nu} --b ${b} -o profile.txt --rl ${params.model.reweight.lambda_radius} --rnu ${params.model.reweight.nu_radius} --rb ${params.model.reweight.b_radius}
     """

     stub:
     """
     touch profile.txt
     """

 }

process MODEL_EVALUATION {

    publishDir "${params.resultsDir}", mode: 'copy'

    input:
    tuple path(genome_structure), path(genome_reconstruction)
    path('profiles.txt')

    output:
    path('loglik.txt')

    script:
    """
    spi-eval.py -i ${genome_structure} -t ${genome_reconstruction} -p profiles.txt -o loglik.txt
    """

    stub:
    """
    touch loglik.txt
    """

}

process REWEIGHTING_EVALUATION {

    publishDir "${params.resultsDir}", mode: 'copy'

    input:
    tuple path(genome_structure), path(genome_reconstruction)
    path('profiles.txt')

    output:
    path('reweight-loglik.txt')

    script:
    """
    spi-eval.py -i ${genome_structure} -t ${genome_reconstruction} -p profiles.txt -o reweight-loglik.txt
    """

    stub:
    """
    touch reweight-loglik.txt
    """

}


process COMPRESS_PROFILES {
    publishDir "${params.resultsDir}", mode: 'copy'

    input:
    path(profiles_file)

    output:
    path "profiles.gz"

    script:
    """
    cat ${profiles_file} | gzip > profiles.gz
    """

    stub:
    """
    touch profiles.gz
    """
}

process TRAJECTORIES_EVALUATION {
    publishDir "${params.resultsDir}", mode: 'copy'

    input:
    path(trajectories_file)

    output:
    path "trajectories-stats.txt"

    script:
    """
    spi-trajectories-stats.py -i ${trajectories_file} -o trajectories-stats.txt
    """

    stub:
    """
    touch trajectories-stats.txt
    """
}

process COMPRESS_TRAJECTORIES {
    publishDir "${params.resultsDir}", mode: 'copy'

    input:
    path(trajectories_file)

    output:
    path "trajectories.gz"

    script:
    """
    cat ${trajectories_file} | gzip > trajectories.gz
    """

    stub:
    """
    touch trajectories.gz
    """
}

process DECOMPRESS_TRAJECTORIES {
    input:
    path(trajectories_file)

    output:
    path "trajectories.txt"

    script:
    """
    gunzip -c ${trajectories_file} > trajectories.txt
    """

    stub:
    """
    touch trajectories.txt
    """
}

process COMPRESS_REWEIGHTED_PROFILES {
    publishDir "${params.resultsDir}", mode: 'copy'

    input:
    path(trajectories_file)

    output:
    path "reweight-profiles.gz"

    script:
    """
    cat ${trajectories_file} | gzip > reweight-profiles.gz
    """

    stub:
    """
    touch reweight-profiles.gz
    """

}

process TELEMETRY {
    publishDir "${params.resultsDir}", mode: 'copy'

    output:
    path('run.info.txt')

    script:
    """
    echo '${motd}' > run.info.txt
    """

    stub:
    """
    touch run.info.txt
    """
}


workflow SPI{

    // parameters channel
    genome_structure_ch = channel.fromPath(params.genome.structure)
    genome_reconstruction_ch = channel.fromPath(params.genome.reconstruction)
    genome_ch = genome_structure_ch.combine(genome_reconstruction_ch)

    if (params.model.reweight.run != "only"){
        // model parameters 
        model_lam_ch = channel.from(params.model.lambda)
        model_nu_ch = channel.from(params.model.nu)
        model_b_ch = channel.from(params.model.b)
        model_config_ch = channel.from(params.model.configs)
        
        // build parameters grid
        sampling_params_ch = genome_structure_ch.combine(model_lam_ch)
                                                .combine(model_nu_ch)
                                                .combine(model_b_ch)
                                                .combine(model_config_ch)

        // sample genomes
        MODEL_SAMPLING(sampling_params_ch)
        

        
        // ensemble channel
        ensemble_ch = MODEL_SAMPLING.out.multiMap{ lam, nu, b, max_config, profile, trajectories, run_info -> 
                                                    profiles: file(profile) 
                                                    trajectories: file(trajectories)
                                                    run_info: file(run_info)
                                            }
        ensemble_profiles_ch = ensemble_ch.profiles.collectFile(name: "profiles.txt")
        ensemble_trajectories_ch = ensemble_ch.trajectories.collectFile(name: "trajectories.txt")
        ensemble_run_info_ch = ensemble_ch.run_info.collectFile(name: "sampling-info.txt", storeDir: params.resultsDir)

        // evaluate parameters likelihood
        MODEL_EVALUATION(genome_ch, ensemble_profiles_ch)
        TRAJECTORIES_EVALUATION(ensemble_trajectories_ch)

        // publish files
        COMPRESS_PROFILES(ensemble_profiles_ch)
        COMPRESS_TRAJECTORIES(ensemble_trajectories_ch)

    }



    if (params.model.reweight.run != "no"){
        
        // load precomputed trajectories
        if (params.model.reweight.run == "only"){
            ensemble_trajectories_ch = DECOMPRESS_TRAJECTORIES(file(params.model.reweight.trajectories))
        }

        // histogram reweighting
        model_rlam_ch = channel.from(params.model.reweight.lambda)
        model_rnu_ch = channel.from(params.model.reweight.nu)
        model_rb_ch = channel.from(params.model.reweight.b)

        reweight_params_ch = genome_structure_ch.combine(ensemble_trajectories_ch)
                                                .combine(model_rlam_ch)
                                                .combine(model_rnu_ch)
                                                .combine(model_rb_ch)

        MODEL_REWEIGHTING(reweight_params_ch)
        reweighting_profiles_ch = MODEL_REWEIGHTING.out.collectFile(name: "reweigthed-trajectories.txt")
        REWEIGHTING_EVALUATION(genome_ch, reweighting_profiles_ch)
        COMPRESS_REWEIGHTED_PROFILES(reweighting_profiles_ch)
        
    }

    // saving run information
    TELEMETRY()
}

workflow {
    SPI()
}

