public class Vcf {
    /**
     * Given two streams of lines in VCF format, returns true if they are close enough to be considered
     * a match for testing purposes. This means that the meta and header lines are ignored and
     * certain data fields have specialised comparison behaviour.
     */
    public static def compare(expected, actual) {
        def expectedLines = expected.filter( line -> !line.startsWith("#") ).iterator()
        def actualLines = actual.filter( line -> !line.startsWith("#") ).iterator()
        while (expectedLines.hasNext() && actualLines.hasNext()) {
            def expectedLine = expectedLines.next()
            def actualLine = actualLines.next()
            if (!dataLinesMatch(expectedLine, actualLine)) {
                return false
            }
        }
        // Return false if one of the streams finished before the other.
        if (expectedLines.hasNext() || actualLines.hasNext()) {
            return false
        }
        return true
    }

    /**
     * Decides whether two data lines in VCF format are effectively a match. This means
     * allowing the QUAL field to vary within a certain range.
     */
    private static def dataLinesMatch(expectedLine, actualLine) {
        def expectedFields = expectedLine.split("\t")
        def actualFields = actualLine.split("\t")
        // Return false if they have different numbers of fields, which is an invariant
        // for the `vcfDataFieldsMatch` method and clearly means they don't match anyway.
        if (expectedFields.length != actualFields.length) {
            return false
        }
        // Compare each field individually.
        def linesMatch = true
        for (def i = 0; i < expectedFields.length; i++) {
            linesMatch = linesMatch && dataFieldsMatch(expectedFields[i], actualFields[i], i)
        }
        return linesMatch
    }

    /**
     * Decides whether two data fields match, using their index to select specialised behaviour
     * for some of the fields.
     */
    private static def dataFieldsMatch(expectedField, actualField, index) {
        switch(index) {
            case 5:
                // QUAL field: if they don't match, attempt to parse them into floats
                // and allow them to be within a certain range of each other.
                return expectedField == actualField || compareFloats(expectedField, actualField)
            case 7:
                // INFO field: ignore this field as it contains arbitrary key-value pairs so
                // it's hard to reason strongly about it.
                return true
            default:
                // Do a simple comparison if there is no special behaviour for this index.
                return expectedField == actualField
        }
    }

    /**
     * Returns true if the two strings passed in can be parsed into floating point numbers
     * with a maximum of maxDiff difference between them.
     */
    private static def compareFloats(expectedField, actualField, maxDiff=0.05) {
        def expectedFloat
        def actualFloat
        try {
            expectedFloat = Float.parseFloat(expectedField)
            actualFloat = Float.parseFloat(actualField)
        } catch (NullPointerException | NumberFormatException ex) {
            // Fail if the fields cannot be parsed into floats.
            return false
        }
        // Check the floating-point values are within range of each other.
        def diff = Math.abs(expectedFloat - actualFloat)
        return diff <= maxDiff
    }
}
