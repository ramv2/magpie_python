import sys

def parse_composition_file(composition_file_str, file_type='csv',
                           verbose=False):
    """
    :param composition_file_str: Name of the file containing composition
    information of the materials.

    :param file_type: Composition file type, default is csv

    :param verbose: Flag to be used mainly for debugging purposes

    :return entries: A list of dictionaries containing <Element name,
    fraction> as <key,value> pairs.
    """
    entries = []
    split_str = ","
    if (file_type == "ssv"):
        split_str = " "
    elif (file_type == "tsv"):
        split_str = "\t"


    try:
        composition_file = open(composition_file_str,'r')
    except IOError:
        print "Input composition file {} doesn't exist!!! Please make sure " \
              "you specify the correct file name".format(composition_file_str)
        sys.exit(1)
    else:
        if (verbose):
            tmp_str = "{} opened successfully. \nStarted reading " \
                      "data..\n".format(composition_file_str)
            print tmp_str
        entry = {}
        for line in composition_file.readlines():
            words = line.strip().split(split_str)
            if (len(words) % 2 != 0):
                print "Invalid line of values in file."
                print line
                sys.exit(1)

            for i in xrange(0, len(words), 2):
                if (words[i].upper() in entry):
                    entry[words[i].upper()] += float(words[i+1])
                else:
                    entry[words[i].upper()] = float(words[i+1])

            total = sum(entry.values())
            for item in entry:
                entry[item] /= total
            entries.append(entry)

    return entries

if __name__ == "__main__":
    dict = parse_composition_file("sample_composition_file.txt")
    print dict