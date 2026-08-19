/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Input-size dependent resource requests
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Groovy classes under lib/ are added to the pipeline classpath automatically, so
    every module can call these without an include.

    A `withLabel:` tier has to be sized for the largest dataset the label will ever
    see, and every smaller task then reserves that same amount and leaves most of it
    idle. Where a process's memory provably tracks the size of its own input, a
    closure over that size tracks with it instead.

    The coefficients are not hard-coded here. They live in `params.dynamic_memory`,
    measured from previous runs by `bin/resource_efficiency.py --dynamic`, so
    retuning never means editing a module. With no entry for a process, the cap is
    returned unchanged and the process behaves exactly as its label always did.

    The same allocation also has to be divided up *inside* a task: `bamSortRam()`
    turns it into STAR's `--limitBAMsortRAM`, so that arithmetic lives here under
    unit test rather than inside the module.

    See docs/RESOURCE_TUNING.md.
----------------------------------------------------------------------------------------
*/

import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths

class BcaResources {

    /**
     * The file a staged input actually stands for.
     *
     * Inside a process -- in a directive exactly as much as in the script block --
     * Nextflow does not hand over the path it read the input from. It hands over a
     * `TaskPath`: a Path whose name is the name the file will be staged under, and
     * whose file system provider is not the one that owns the file. Every
     * `java.nio.file.Files` call on such a path raises ProviderMismatchException, so
     * a size read straight off it comes back as 0 and the request quietly collapses
     * onto the label -- the exact silent failure this class exists to avoid.
     *
     * `toRealPath()` is the way back to the file itself. On an ordinary Path it is a
     * harmless no-op, and anything that refuses it is returned untouched.
     */
    static Path realPath(Path path) {
        try {
            return path.toRealPath()
        }
        catch (Exception ignored) {
            return path
        }
    }

    /**
     * Bytes of one staged input, whatever shape it arrives in.
     *
     * Deliberately avoids Groovy's `size()` extension on Path: this runs while the
     * directive is being evaluated, and a missing or unreadable input must yield 0
     * rather than throw, or the whole run dies before a single task is submitted.
     * (`size()` would also read a directory as 0, which is what a genome index is.)
     */
    static long fileSize(Object entry) {
        try {
            if (entry == null) {
                return 0L
            }
            Path path
            if (entry instanceof Path) {
                path = realPath((Path) entry)
            }
            else if (entry instanceof File) {
                path = ((File) entry).toPath()
            }
            else {
                path = Paths.get(entry.toString())
            }
            return Files.isDirectory(path) ? directorySize(path) : Files.size(path)
        }
        catch (Exception ignored) {
            return 0L
        }
    }

    /** Recursive byte count, so a staged matrix directory measures its contents. */
    static long directorySize(Path root) {
        long total = 0L
        try {
            Files.walk(root).withCloseable { stream ->
                for (Path path : stream) {
                    try {
                        if (!Files.isDirectory(path)) {
                            total += Files.size(path)
                        }
                    }
                    catch (Exception ignored) {
                    }
                }
            }
        }
        catch (Exception ignored) {
        }
        return total
    }

    /** Total bytes across a single input or a collection of them. */
    static long totalBytes(Object files) {
        long total = 0L
        if (files instanceof Collection) {
            for (Object entry : (Collection) files) {
                total += fileSize(entry)
            }
        }
        else if (files instanceof Object[]) {
            for (Object entry : (Object[]) files) {
                total += fileSize(entry)
            }
        }
        else {
            total = fileSize(files)
        }
        return total
    }

    /**
     * A memory request scaled from the size of a process's own input.
     *
     * The model is anchored on a measured point rather than an abstract
     * coefficient, because an anchor is something a reader can sanity-check:
     * at `ref_gb` of input the process needed `mem_gb`, and from there memory grows
     * as size**exponent.
     *
     * @param spec    params.dynamic_memory entry: exponent, ref_gb, mem_gb and
     *                optionally floor_gb and cap_gb. Null or incomplete means
     *                "not measured".
     * @param files   The input this process's memory scales with.
     * @param attempt task.attempt, so the retry ladder still escalates.
     * @param labelGb The process's own label tier, in GB. Used verbatim when there
     *                is no spec, which is what makes an unmeasured process behave
     *                exactly as before, and as the ceiling unless the spec raises it
     *                with `cap_gb` -- which it must whenever a process was already
     *                close to exhausting its tier.
     * @param residentFiles Optional. Input the process must hold in memory
     *                regardless of how much data it processes -- a genome index,
     *                typically. Its size is not scaled; it raises the floor, so a
     *                small input paired with a large reference can never be handed
     *                less memory than the reference alone needs. Without this,
     *                scaling STARSOLO_ALIGN on read volume alone would hand a small
     *                sample a request too small to load the index at all.
     * @return A memory string such as "24 GB", which Nextflow parses.
     */
    static String scaledMemory(Map spec, Object files, int attempt, Number labelGb,
                               Object residentFiles = null, Number defaultFloorGb = 1) {
        int gb = Math.max(1, labelGb as int)

        if (spec && spec.exponent != null && spec.ref_gb && spec.mem_gb) {
            int cap = Math.max(1, (spec.cap_gb ?: labelGb) as int)
            long bytes = totalBytes(files)
            // A zero here means the input could not be measured. Falling through to
            // the label is the safe reading: it is the behaviour it already had.
            if (bytes > 0L) {
                double inputGb = bytes / (1024.0d * 1024.0d * 1024.0d)
                double estimate = (spec.mem_gb as double) *
                        Math.pow(inputGb / (spec.ref_gb as double), spec.exponent as double)

                int floor = Math.max(1, (spec.floor_gb ?: defaultFloorGb) as int)
                if (residentFiles != null) {
                    long residentBytes = totalBytes(residentFiles)
                    if (residentBytes > 0L) {
                        double factor = (spec.resident_factor ?: 1.3d) as double
                        double residentGb = (residentBytes / (1024.0d * 1024.0d * 1024.0d)) * factor
                        floor = Math.max(floor, (int) Math.ceil(residentGb))
                    }
                }

                gb = Math.min(Math.max((int) Math.ceil(estimate), floor), cap)
            }
        }

        return "${gb * Math.max(1, attempt)} GB"
    }

    /** Bytes from a MemoryUnit (`task.memory`), a plain number, or a numeric string. */
    static long asBytes(Object value) {
        if (value == null) {
            return 0L
        }
        if (value instanceof Number) {
            return ((Number) value).longValue()
        }
        try {
            // MemoryUnit, reached by duck typing rather than an import, so this
            // class stays compilable without Nextflow on the classpath.
            return value.toBytes() as long
        }
        catch (Exception ignored) {
        }
        try {
            return Long.parseLong(value.toString().trim())
        }
        catch (Exception ignored) {
            return 0L
        }
    }

    /** What STAR needs on top of the index and the sort buffer: alignment buffers,
     *  the solo structures, and the usual container overhead. */
    private static final long BAMSORT_SLACK_BYTES = 4L * 1024 * 1024 * 1024

    /** Asked for instead of a negative budget when the index alone fills the
     *  allocation, so that STAR fails with the figure it needs rather than on
     *  a number this file made up. */
    private static final long BAMSORT_FLOOR_BYTES = 1024L * 1024 * 1024

    /**
     * STAR's `--limitBAMsortRAM`, in bytes.
     *
     * Left at 0, STAR sets this ceiling from the size of the genome index, which has
     * nothing to do with the memory the job was actually given: on a generous
     * allocation it sorts in far less RAM than it could, and on a tight one it
     * reserves memory the scheduler never granted and the task is killed. The budget
     * that is actually correct is what remains of the allocation once the index --
     * which STAR keeps resident for the whole run -- and its own working memory are
     * accounted for.
     *
     * @param override      params.star_limitBAMsortRAM. Anything above zero pins the
     *                      value and is used verbatim, which is how a user answers a
     *                      STAR error that names the exact figure it wants.
     * @param taskMemory    task.memory, the allocation to divide up. Absent, there is
     *                      nothing to derive from and STAR's own behaviour is left in
     *                      place.
     * @param residentFiles The genome index.
     * @return `[bytes: <long>, note: <String>]` -- the value for STAR, and a
     *         one-line account of where it came from for the task log, so an OOM can
     *         be read back without re-deriving the arithmetic by hand.
     */
    static Map bamSortRam(Object override, Object taskMemory, Object residentFiles) {
        long pinned = asBytes(override)
        if (pinned > 0L) {
            return [bytes: pinned, note: "${pinned} pinned by params.star_limitBAMsortRAM".toString()]
        }

        long alloc = asBytes(taskMemory)
        if (alloc <= 0L) {
            return [bytes: 0L,
                    note: "0: no memory directive to derive from, leaving STAR its own default".toString()]
        }

        long indexBytes = totalBytes(residentFiles)
        long bytes
        String note
        if (indexBytes <= 0L) {
            // The index could not be measured. Guessing high here would get the job
            // OOM-killed, so take a conservative fraction of the allocation instead.
            bytes = (alloc * 60L).intdiv(100L)
            note = "${bytes}: genome index size unknown, using 60% of ${alloc}"
        }
        else {
            bytes = alloc - indexBytes - BAMSORT_SLACK_BYTES
            note = "${bytes} = allocation ${alloc} - genome index ${indexBytes} - slack ${BAMSORT_SLACK_BYTES}"
        }

        if (bytes < BAMSORT_FLOOR_BYTES) {
            // A genome index larger than the allocation leaves nothing to sort in.
            // Ask for a modest buffer rather than a negative one and let STAR report
            // the real shortfall -- its error names the exact figure needed.
            bytes = BAMSORT_FLOOR_BYTES
            note = "${bytes}: floored at 1 GB, the genome index leaves no room in this allocation"
        }

        return [bytes: bytes, note: note.toString()]
    }
}
