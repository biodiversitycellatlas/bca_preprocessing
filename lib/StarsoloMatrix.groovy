/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    STARsolo count matrix
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Which of STARsolo's raw matrices the pipeline analyses. STAR always writes the
    unique-gene counts to raw/matrix.mtx; every other --soloMultiMappers type it is given
    (EM, Uniform, PropUnique, Rescue) adds raw/UniqueAndMult-<type>.mtx, the unique counts
    plus the multi-gene reads distributed by that rule. The first such type listed in
    star_soloMultiMappers is the matrix analysed downstream; with "Unique" alone it is
    matrix.mtx.

    STAR writes the multi-gene matrices to raw/ only, so every matrix directory this
    pipeline writes (filtered_secondderiv/filtered, filtered_<file>) holds the analysed
    counts as matrix.mtx again, and only raw/ needs the name resolved here.

    Called from workflow code only: a bare 'params' in a process script is hashed by
    object identity and would break -resume (see MappingRateCheck).

    See "STARsolo multi-gene reads" in docs/CONFIGURATION_PARAMETERS.md.
----------------------------------------------------------------------------------------
*/

class StarsoloMatrix {

    static final String UNIQUE_FILE = 'matrix.mtx'

    /** star_soloMultiMappers as STARSOLO_ALIGN resolves it: the parameter, else the protocol's default. */
    static String setting(Map params) {
        def value = params.star_soloMultiMappers ?: params.seqtech_parameters?.get(params.protocol)?.star_soloMultiMappers
        return (value ?: 'Unique').toString()
    }

    /** The multi-gene type whose matrix is analysed, e.g. 'EM', or null for unique counts only. */
    static String multiGeneType(Map params) {
        return setting(params).tokenize(' ').find { type -> type != 'Unique' }
    }

    /** File name of the analysed matrix inside STARsolo's raw/ directory. */
    static String rawFileName(Map params) {
        def type = multiGeneType(params)
        return type ? "UniqueAndMult-${type}.mtx" : UNIQUE_FILE
    }

    /**
    * Whether the UMI-per-cell curve has to be derived from the analysed matrix. STAR's own
    * UMIperCellSorted.txt counts unique-gene UMIs only, so with a multi-gene matrix a cutoff
    * read off it would be applied to totals on a different basis.
    */
    static boolean curveFromMatrix(Map params) {
        return multiGeneType(params) != null && params.cellfilter_method in ['second_derivative', 'manual_cutoff']
    }
}
