from shutil import which


def is_tool(name):
    """Check whether `name` is on PATH and marked as executable."""
    return which(name) is not None


def check_program_install(list_names):
    flag_violate = False
    for name in list_names:
        if is_tool(name) == False:
            print(name, "is a prerequisite program, please install it before running biastools")
            flag_violate = True
    if flag_violate:
        exit(1)


def bool2str(flag):
    if flag:
        return "1"
    else:
        return "0"


def catch_assert(parser, message):
    print('\n', message, '\n')
    parser.print_usage()
    exit(1)


