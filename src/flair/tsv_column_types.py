"""Column converters shared by the flair TSV formats.

A converter is the (parse, format) pair a typeMap entry takes, so a column whose
text holds several values is a list on the row and is joined again on write.
Callers then never split or join a column themselves, and the separator is
stated once, here and in the format that uses it.
"""

def separated_list_type(sep):
    "a column holding values joined by sep; the row carries the list of them"
    return (lambda value: [] if value == '' else value.split(sep),
            lambda values: sep.join(values))


# comma separated, the usual case: read names, isoform ids, filter tags
commaListType = separated_list_type(',')

# '; ' separated, used where the values themselves contain commas
semicolonSpaceListType = separated_list_type('; ')
