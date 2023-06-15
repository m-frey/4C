import argparse
import os
import common_utils as utils


def check_header_guards(look_cmd, allerrors):
    num_wrong_header_guards = 0
    for file_path in utils.files_changed(look_cmd):
        file_ending = os.path.splitext(os.path.basename(file_path))[1]

        if file_path == "" or (file_ending != ".H" and file_ending != ".h"):
            continue

        file_name = os.path.splitext(os.path.basename(file_path))[0]
        define_name = file_name.replace("-", "_").replace(".", "_").upper() + "_H"

        template_header_if = "#ifndef " + define_name + "\n"
        template_header_define = "#define " + define_name + "\n"

        def has_wrong_header_guards():
            all_lines = utils.file_contents(file_path)
            for line_num, line in enumerate(all_lines):
                if (
                    line == template_header_if
                    and all_lines[line_num + 1] == template_header_define
                ):
                    return False
            return True

        if has_wrong_header_guards():
            allerrors.append("Wrong header guard in " + file_path)
            num_wrong_header_guards += 1

    return num_wrong_header_guards


def main():
    # build command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--diff_only",
        action="store_true",
        help="Add this tag if only the difference to HEAD should be analyzed. This flag should be used as a pre-commit hook. Otherwise all files are checked.",
    )
    parser.add_argument(
        "--out",
        type=str,
        default=None,
        help="Add this tag if the error message should be written to a file.",
    )
    args = parser.parse_args()

    # flag, whether only touched files should be checked
    diff_only = args.diff_only

    # error file (None for sys.stderr)
    errfile = args.out
    errors = 0
    allerrors = []
    try:
        if diff_only:
            look_cmd = "git diff --name-only --cached --diff-filter=MRAC"
        else:
            look_cmd = "git ls-files"
        errors += check_header_guards(look_cmd, allerrors)
    except ValueError:
        print("Something went wrong! Check the error functions in this script again!")
        errors += 1

    utils.pretty_print_error_report(
        "Wrong header guards in the following files:", allerrors, errfile
    )
    return errors


if __name__ == "__main__":
    import sys

    sys.exit(main())
