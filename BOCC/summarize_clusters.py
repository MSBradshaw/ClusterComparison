from BOCC import load_clusters, summarize_clusters
import argparse


def get_args():
    parser = argparse.ArgumentParser()


    parser.add_argument('--coms',
                        dest='coms',
                        required=True,
                        help='tsv file where each row in a community, first area is com name followed by a tab ' +
                             'seporated list of community members (not all rows will be the same length)')

    parser.add_argument('--alpha',
                        dest='alpha',
                        type=float,
                        default=0.05,
                        help='threshold for significance in GO enrichment')

    parser.add_argument('--mg2',
                        dest='mg2',
                        type=str,
                        help='threshold for significance in GO enrichment')

    parser.add_argument('--out',
                        dest='out',
                        type=str,
                        required=True,
                        help='output file, should be a .tsv)')

    args = parser.parse_args()
    return args

args = get_args()
coms = load_clusters(args.coms)
df = summarize_clusters(coms, args.mg2, args.alpha)
df.to_csv(args.out, sep='\t', index=False)
