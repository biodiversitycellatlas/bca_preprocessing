process MTX_TO_H5AD {
    publishDir "${params.outdir}/anndata/${meta.id}/${meta.datatype}", mode: 'copy'
    tag "${meta.id} | ${meta.mapping_method} | ${meta.datatype}"
    label 'process_low'

    container 'oras://community.wave.seqera.io/library/scanpy:1.12--45f1dccaf83880df'
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(mtx), path(barcodes), path(features)

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    path "versions.yml",             emit: versions

    script:
    """
    #!/usr/bin/env python3
    import scanpy as sc
    import pandas as pd
    import anndata as ad
    from scipy.io import mmread

    # Read the sparse matrix
    X = mmread("${mtx}").tocsr()

    # Read Barcodes
    obs_names = pd.read_csv("${barcodes}", header=None)[0].astype(str).values
    obs = pd.DataFrame(index=obs_names)

    # Read Features
    var_names = pd.read_csv("${features}", header=None, sep='\\t')[0].astype(str).values
    var = pd.DataFrame(index=var_names)

    # Create AnnData object
    if X.shape[0] == len(var_names) and X.shape[1] == len(obs_names):
        adata = ad.AnnData(X=X.T, obs=obs, var=var)
    else:
        adata = ad.AnnData(X=X, obs=obs, var=var)

    # Clean up names and metadata
    adata.var_names_make_unique()
    adata.obs['sample_id'] = "${meta.id}"

    # Write compressed H5AD
    adata.write_h5ad("${meta.id}_${meta.datatype}.h5ad", compression="gzip")

    with open("versions.yml", "w") as f:
        f.write(f'"${task.process}":\\n')
        f.write(f'    python: {__import__("platform").python_version()}\\n')
        f.write(f'    scanpy: {sc.__version__}\\n')
    """
}
