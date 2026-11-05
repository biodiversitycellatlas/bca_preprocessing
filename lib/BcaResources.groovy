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

    Time follows the same model for the steps that stream their reads
    (`params.dynamic_time`, `scaledTime()`), and `tieredCpus()` steps a thread count
    up above an input size (`params.dynamic_cpus`).

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

    /**
     * Total bytes across a single input or a collection of them.
     *
     * Recurses, because a multi-file input (a `path` that received a collected list)
     * reaches the directive as a list of its own, so `[star_logs, cell_stats]` is a list
     * of lists. Read as one entry, such a list measures as 0.
     */
    static long totalBytes(Object files) {
        long total = 0L
        if (files instanceof Collection) {
            for (Object entry : (Collection) files) {
                total += totalBytes(entry)
            }
        }
        else if (files instanceof Object[]) {
            for (Object entry : (Object[]) files) {
                total += totalBytes(entry)
            }
        }
        else {
            total = fileSize(files)
        }
        return total
    }

    /**
     * The files directly inside a staged directory whose name ends with `suffix`.
     *
     * For an input of which only some files are held in memory: a Kraken2 database
     * loads its `*.k2d` files, while a database built in place may also hold the
     * `library/` and `taxonomy/` it was built from, many times larger. Anything that
     * cannot be listed yields an empty list, which measures as 0.
     */
    static List<Path> filesEndingWith(Object dir, String suffix) {
        List<Path> found = []
        try {
            Path path = dir instanceof Path ? realPath((Path) dir) : Paths.get(dir.toString())
            Files.list(path).withCloseable { stream ->
                for (Path entry : stream) {
                    if (entry.fileName.toString().endsWith(suffix) && !Files.isDirectory(entry)) {
                        found << entry
                    }
                }
            }
        }
        catch (Exception ignored) {
        }
        return found
    }

    /**
     * The files at `relPaths` under a staged directory, those that exist.
     *
     * For a directory of which a step reads only named files: STARsolo's raw/ holds
     * matrix.mtx next to every UniqueAndMult-*.mtx, and only the analysed one is loaded;
     * a Solo.out directory holds every feature's matrices. A path that is absent or a
     * directory is left out, and anything that cannot be resolved yields an empty list,
     * which measures as 0.
     */
    static List<Path> filesIn(Object dir, List relPaths) {
        List<Path> found = []
        try {
            Path path = dir instanceof Path ? realPath((Path) dir) : Paths.get(dir.toString())
            for (Object rel : relPaths) {
                if (!rel) {
                    continue
                }
                Path entry = path.resolve(rel.toString())
                if (Files.isRegularFile(entry)) {
                    found << entry
                }
            }
        }
        catch (Exception ignored) {
        }
        return found
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
            Double estimate = powerLaw(spec, 'mem_gb', totalBytes(files))
            // Null means the input could not be measured. Falling through to the
            // label is the safe reading: it is the behaviour it already had.
            if (estimate != null) {
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

    /**
     * The anchored power law shared by the scaled requests: at `ref_gb` of input the
     * process needed `spec[valueKey]`, and from there it grows as size**exponent.
     * Null when the input measured 0 bytes, which means it could not be measured.
     */
    private static Double powerLaw(Map spec, String valueKey, long bytes) {
        if (bytes <= 0L) {
            return null
        }
        double inputGb = bytes / (1024.0d * 1024.0d * 1024.0d)
        return (spec[valueKey] as double) *
                Math.pow(inputGb / (spec.ref_gb as double), spec.exponent as double)
    }

    /**
     * A time limit scaled from the size of a process's own input.
     *
     * The streaming FASTQ and BAM steps take time in proportion to the reads they
     * pass over, so a fixed limit sized for an ordinary sample kills a deep one
     * halfway through. Same anchored model as scaledMemory(): at `ref_gb` of input
     * the process took `time_h` hours.
     *
     * Unlike memory, the floor defaults to the label rather than to a token minimum:
     * an entry without `floor_h` can lengthen a limit but never shorten one, so a
     * step that was finishing in time keeps finishing in time.
     *
     * @param spec       params.dynamic_time entry: exponent, ref_gb, time_h and
     *                   optionally floor_h and cap_h. Null or incomplete means
     *                   "not measured".
     * @param files      The input this process's run time scales with.
     * @param attempt    task.attempt, so the retry ladder still escalates.
     * @param labelHours The process's current limit, in hours (0.75 for 45 min). Used
     *                   verbatim when there is no spec, and as the ceiling unless the
     *                   spec raises it with `cap_h`.
     * @return A duration string in whole minutes, such as "45m", which Nextflow parses.
     */
    static String scaledTime(Map spec, Object files, int attempt, Number labelHours) {
        double label = Math.max(1.0d / 60.0d, labelHours as double)
        double hours = label

        if (spec && spec.exponent != null && spec.ref_gb && spec.time_h) {
            Double estimate = powerLaw(spec, 'time_h', totalBytes(files))
            if (estimate != null) {
                double floor = (spec.floor_h ?: label) as double
                double cap = (spec.cap_h ?: label) as double
                hours = Math.min(Math.max(estimate, floor), cap)
            }
        }

        long minutes = Math.max(1L, (long) Math.ceil(hours * 60.0d - 1e-9d))
        return "${minutes * Math.max(1, attempt)}m"
    }

    /**
     * A CPU count that steps up above an input size.
     *
     * For a multi-threaded tool whose wall time is worth buying down on a large input
     * but not a small one. The spec is `[above_gb: X, cpus: N]`; with no threshold set,
     * or an input that cannot be measured, the request is `baseCpus`, as before.
     *
     * @return The CPU count, already multiplied by task.attempt.
     */
    static int tieredCpus(Map spec, Object files, int attempt, int baseCpus) {
        int cpus = Math.max(1, baseCpus)
        // instanceof rather than != null: a ConfigObject answers a missing key with an
        // empty ConfigObject, not null.
        if (spec && spec.above_gb instanceof Number && spec.cpus instanceof Number) {
            long bytes = totalBytes(files)
            double threshold = (spec.above_gb as double) * 1024.0d * 1024.0d * 1024.0d
            if (bytes > 0L && bytes > threshold) {
                cpus = Math.max(1, spec.cpus as int)
            }
        }
        return cpus * Math.max(1, attempt)
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
