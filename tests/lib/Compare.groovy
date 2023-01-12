import groovy.io.FileType
import java.nio.file.Files
import java.nio.file.Paths
import java.util.zip.GZIPInputStream

public class Compare {
    /**
     * Checks that the contents of the two dirs are identical, using a recursive
     * traversal making md5 comparisons on each file pair. Any files missing or different
     * in either directory will cause this function to return false.
     * In order to allow the caller to specify complex behaviour for comparing different
     * types of files, a 'compare' closure is taken. This closure will be passed the
     * expected and actual file for each file pair and should return a `Result` indicating
     * the outcome of the comparison.
     */
    public static def dirsConform(File expectedDir, File actualDir, Closure compare) {
        if (!expectedDir.exists()) {
            System.err.println "Directory `${expectedDir}` doesn't exist"
            return false
        }
        if (!actualDir.exists()) {
            System.err.println "Directory `${actualDir}` doesn't exist"
            return false
        }

        // Get the set of every relative file path in either directory.
        def relativePaths = new HashSet<String>()
        relativePaths.addAll(relativePathsIn(expectedDir))
        relativePaths.addAll(relativePathsIn(actualDir))

        def missingFiles = []
        def unexpectedFiles = []
        def incorrectFiles = []

        // Compare each pair of files which share a relative path in their respective dir.
        relativePaths.each { path ->
            def expectedFile = new File(expectedDir, path)
            def actualFile = new File(actualDir, path)

            switch(compare(expectedFile, actualFile)) {
                case Result.MISSING:
                    missingFiles << path
                    break
                case Result.UNEXPECTED:
                    unexpectedFiles << path
                    break
                case Result.INCORRECT:
                    incorrectFiles << path
                    break
                default:
                    break
            }
        }

        // Output differences if there are any.
        if (missingFiles.size() > 0 || unexpectedFiles.size() > 0 || incorrectFiles.size() > 0) {
            reportFailedFiles(
                missingFiles, "${missingFiles.size()} files were expected but not produced"
            )
            reportFailedFiles(
                unexpectedFiles, "${unexpectedFiles.size()} files were produced that were not expected"
            )
            reportFailedFiles(
                incorrectFiles, "${incorrectFiles.size()} files were produced but their contents were not as expected"
            )
            return false
        }
        return true
    }

    /**
     * The outcome of comparing an expected file to what was actually produced.
     */
    public enum Result {
        ACCEPT, // The comparison found no issues.
        MISSING, // The expected file exists but the actual file doesn't.
        UNEXPECTED, // The actual file exists but the expected file doesn't.
        INCORRECT // The actual file's contents don't match the expected file.
    }

    /**
     * Returns the uncompressed contents of a gzipped file, as a string.
     */
    public static def readGzip(file) {
        def inputStream = new FileInputStream(file)
        def gzipStream = new GZIPInputStream(inputStream)
        return gzipStream.text
    }

    /**
     * Returns the uncompressed contents of a gzipped file, as a stream of lines.
     */
    public static def streamGzipLines(file) {
        def inputStream = new FileInputStream(file)
        def gzipStream = new GZIPInputStream(inputStream)
        def reader = new InputStreamReader(gzipStream)
        def bufferedReader = new BufferedReader(reader)
        return bufferedReader.lines()
    }
    
    /**
     * Returns the set of every relative path for which there is a file in the given directory.
     */
    private static def relativePathsIn(dir) {
        def set = new HashSet<String>()
        def dirPath = dir.toPath()
        dir.eachFileRecurse (FileType.FILES) { file ->
            def path = dirPath.relativize(file.toPath())
            set.add("$path")
        }
        return set
    }

    /**
     * Prints the paths of any files that failed their comparison with the expected files.
     */
    private static def reportFailedFiles(paths, explanation) {
        if (paths.size() > 0) {
            System.err.println "${explanation}:"
            paths.each { path ->
                System.err.println "$path"
            }
        }
        System.err.println ""
    }
}
