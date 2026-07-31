# Extract literal container image references from a Nextflow module file.
#
# Handles both forms used in this pipeline:
#   container "oras://community.wave.seqera.io/library/fastqc:0.12.1--104d..."
#   container "${ task.ext.use_gpu ? 'us.gcr.io/.../cellbender:0.3.2' :
#       workflow.containerEngine in ['singularity', 'apptainer'] ? 'https://...' :
#       'community.wave.seqera.io/library/cellbender_webcolors:156d...' }"
#
# Commented-out directives (// container "...") are ignored, because the regex
# anchors `container` to the start of the line. Fully dynamic directives such as
# `container { demuxafy_sif }` yield no output -- the caller reports those as
# skipped, since there is no literal reference to validate.

BEGIN { collecting = 0; quotes = 0; block = "" }

# Start of a container directive: `container` at the start of a line, followed
# by a quote or an opening brace.
/^[ \t]*container[ \t]*["'{]/ && collecting == 0 {
    collecting = 1; quotes = 0; block = ""
}

collecting {
    block = block $0 "\n"
    # gsub returns the number of substitutions; replacing " with " is a no-op
    # that lets us count double quotes on this line.
    quotes += gsub(/"/, "\"")

    # An even number of double quotes means the directive's string literals are
    # balanced, so the (possibly multi-line) directive is complete.
    if (quotes % 2 == 0) {
        emit(block)
        collecting = 0; block = ""
    }
}

# Print every quoted literal in the block that looks like an image reference.
function emit(text,   rest, lit) {
    rest = text
    # Literals containing $ are Groovy interpolation, not a usable reference.
    while (match(rest, /"[^"$]+"|'[^'$]+'/)) {
        lit = substr(rest, RSTART + 1, RLENGTH - 2)
        rest = substr(rest, RSTART + RLENGTH)
        # A registry reference always has a path separator; this filters out
        # bare words such as the 'singularity' / 'apptainer' engine names.
        if (lit ~ /\//) print lit
    }
}
