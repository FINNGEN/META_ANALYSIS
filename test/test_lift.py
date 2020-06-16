from scripts.lift import create_parser


def test_parser_check_trivial():
    """
    Trivial test
    """
    parser = create_parser()
    arg = parser.parse_args(["file"])
    assert arg.var is None


def test_parser_check_with_args():
    """
    Parse with sample arguments
    """
    parser = create_parser()
    arg = parser.parse_args(["-chr", "CHR",
                             "-pos", "POS",
                             "-ref", "REF",
                             "-alt", "ALT",
                             "-chain_file", "CHAIN_FILE",
                             "-no_clean",
                             "-tmp_path", "TMP_PATH",
                             "file"])
    assert arg.chr == 'CHR'
    assert arg.alt == 'ALT'
    assert arg.chain_file == 'CHAIN_FILE'
    assert arg.file == 'file'
    assert arg.no_clean
    assert arg.pos == 'POS'
    assert arg.ref == 'REF'
    assert arg.tmp_path == 'TMP_PATH'
    assert arg.var is None
