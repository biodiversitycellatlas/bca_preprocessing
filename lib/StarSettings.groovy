/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    STAR settings resolved from several sources
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    The strand that is counted, and STAR's extra arguments.

    Both take the values they need rather than 'params': a bare 'params' in a process
    script is hashed by object identity and would break -resume (see MappingRateCheck).
----------------------------------------------------------------------------------------
*/

class StarSettings {

    /**
    * The strand counted for a sample (STAR --soloStrand; alevin-fry and the antisense
    * metrics follow it): star_soloStrand when it is set. Otherwise 'Reverse' for a sample
    * read from an Ultima CRAM, whose reads are antisense to the gene, and 'Forward' for
    * every other sample.
    */
    static String strand(Object soloStrand, Object cramType, Object inputFormat) {
        if (soloStrand) {
            return soloStrand.toString()
        }
        return (cramType == 'ultima' && inputFormat == 'cram') ? 'Reverse' : 'Forward'
    }

    /**
    * STAR's extra arguments: the protocol's (conf/seqtech_parameters.config), then the
    * custom ones (star_extraargs). STAR refuses an option given twice, so a custom option
    * replaces the protocol's option of the same name, at its place, with all its values.
    */
    static String mergeArgs(Object protocolArgs, Object customArgs) {
        def merged = new LinkedHashMap<String, List<String>>()
        [protocolArgs, customArgs].each { args ->
            String option = null
            (args ?: '').toString().tokenize().each { token ->
                if (token.startsWith('--')) {
                    option = token
                    merged[option] = []
                } else if (option != null) {
                    merged[option] << token
                } else {
                    throw new IllegalArgumentException(
                        "star_extraargs must start with a STAR option ('--name value ...'), but starts with '${token}'.")
                }
            }
        }
        return merged.collect { option, values -> ([option] + values).join(' ') }.join(' ')
    }
}
