# Reading and writing the TSVs the flair R scripts share, with the same suffix
# driven compression the Python side uses: a path ending .gz or .bz2 is
# compressed, anything else is plain.
#
# Reading needs nothing added.  R's file() and data.table::fread both detect
# compression from the file itself, so read.table and fread already accept a
# compressed path.
#
# Writing does need it.  write.table given a file name opens a plain file
# whatever the name ends in, so a table written to x.tsv.gz would be an
# uncompressed file wearing a .gz name.  The connection has to be opened here
# instead.

# what --compress appends; gzip is the one flair writes, though the writer below
# keeps whatever the given path already ends in
tsvCompressSuffix <- ".gz"

tsvCompressedPath <- function(path, compress) {
    # the output name a run writes, given whether --compress was asked for
    if (isTRUE(compress)) {
        paste0(path, tsvCompressSuffix)
    } else {
        path
    }
}

tsvWriteConnection <- function(path, name = path) {
    # name is the file the table is written as, which states the compression; path is
    # what is opened, a temporary file whose own name states nothing
    if (grepl("\\.gz$", name)) {
        gzfile(path, "w")
    } else if (grepl("\\.bz2$", name)) {
        bzfile(path, "w")
    } else if (grepl("\\.bgz$", name)) {
        # gzfile would write ordinary gzip under a name that states the block
        # format, which is worse than rejecting
        stop(sprintf("cannot write bgzf: %s", name))
    } else {
        file(path, "w")
    }
}

tsvWriteTable <- function(x, path, ...) {
    # write.table with the tab separated, unquoted conventions flair uses, through
    # a connection that compresses when the name states it.  The table is written to
    # a temporary file beside it and renamed when it is complete, as the Python side
    # writes its outputs, so a failed run leaves no partial table
    tmp <- paste0(path, ".tmp.", Sys.getpid())
    con <- tsvWriteConnection(tmp, path)
    written <- FALSE
    on.exit({
        close(con)
        if (written) file.rename(tmp, path) else unlink(tmp)
    })
    write.table(x, file = con, sep = "\t", quote = FALSE, ...)
    written <- TRUE
}

tsvFread <- function(path, ...) {
    # data.table's fread needs the R.utils package to open a compressed file
    # itself, so a compressed one is decompressed through a pipe instead, the way
    # the Python side does it
    if (grepl("\\.(gz|bgz)$", path)) {
        data.table::fread(cmd = sprintf("gzip -dc %s", shQuote(path)), ...)
    } else if (grepl("\\.bz2$", path)) {
        data.table::fread(cmd = sprintf("bzip2 -dc %s", shQuote(path)), ...)
    } else {
        data.table::fread(path, ...)
    }
}
