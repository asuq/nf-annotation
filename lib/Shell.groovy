/** Quote data as one shell argument; never interpret original sample identifiers. */
class Shell {
    static String quote(Object value) {
        if (value == null) { throw new IllegalArgumentException('Missing shell argument') }
        return "'" + value.toString().replace("'", "'\"'\"'") + "'"
    }
}
