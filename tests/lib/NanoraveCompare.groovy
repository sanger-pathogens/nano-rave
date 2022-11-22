class NanoraveCompare {
    /**
     * A comparison function that takes into account the file extension. It:
     * - Ignores .log files.
     * - Only checks the names of .html, .json, .png, and .vcf.gz.tbi files.
     * - Compares only the data lines of .vcf.gz files.
     * - Compares the uncompressed contents of other .gz files.
     * - Compares the contents of everything else.
     * This function is designed to be passed as a closure to Compare.dirsConform.
     */
    public static def compare(expected, actual) {
        def pathStr = "$expected"
        // Early return to ignore these files completely.
        if (pathStr.endsWith(".log")) {
            return Compare.Result.ACCEPT
        }
        // Check both files exist, or return the right kind of error.
        if (!actual.exists()) {
            return Compare.Result.MISSING
        }
        if (!expected.exists()) {
            return Compare.Result.UNEXPECTED
        }
        // Early return to ignore the contents of these files.
        if (
            pathStr.endsWith(".html") ||
            pathStr.endsWith(".json") ||
            pathStr.endsWith(".png") ||
            pathStr.endsWith(".vcf.gz.tbi")
        ) {
            return Compare.Result.ACCEPT
        }
        // Compare the contents of the files.
        if (pathStr.endsWith(".vcf.gz")) {
            def expectedLines = Compare.streamGzipLines(expected)
            def actualLines = Compare.streamGzipLines(actual)
            if (!Vcf.compare(expectedLines, actualLines)) {
                return Compare.Result.INCORRECT
            }
        } else if (pathStr.endsWith(".gz")) {
            if (Compare.readGzip(expected).md5() != Compare.readGzip(actual).md5()) {
                return Compare.Result.INCORRECT
            }
        } else if (expected.text.md5() != actual.text.md5()) {
            return Compare.Result.INCORRECT
        }
        // Accept them if they got through all that.
        return Compare.Result.ACCEPT
    }
}
