

def _discard(_saved, _dropped):
    with open(_saved, 'r'), open(_dropped, 'r') as saved, dropped:
        try:
            dropped_id = next(dropped)
        except StopIteration:
            warnings.warn("discard file %s was empty.")
            return None
            # just symlink
        for read in SeqIO.parse(saved, 'fastq'):
            if str(read.id) == dropped_id:
                # read gets filtered because its ID was in discard so don't yield
                dropped_id = next(dropped)
            else:
                yield read
def discard(saved, disc, outpath):
    res = discard(saved, disc)
    if res:
        SeqIO.write(res, outpath, 'fastq')
    else:
         os.symlink(os.path.abspath(saved), os.path.abspath(outpath))

def main():
    parser = argparse.ArgumentParser(
        description = 'Drop discarded reads from the other pair.'
    )
    parser.add_argument(
        '--saved',
        help='The file with the unmapped reads.'
    )
    parser.add_argument(
        '--other-discard',
        help='The discard file for the *opposite* member of the pair.
        If --saved is R1, this should be R2.discard, etc.'
    )
    parser.add_argument(
        '--out',
        help='The output path'
    )
    args = parser.parse_args()
    discard(args.saved, args.other_discard, args.out)


