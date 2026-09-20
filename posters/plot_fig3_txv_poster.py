"""
Fig 3 (poster variant): therapeutic-enhanced screening strategies.

Identical to plot_fig3_txv.py, output directed to the poster folder.
"""
import argparse

import plot_fig3_txv


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--resfolder', default='results/v2.2.6_baseline')
    parser.add_argument('--outpath', default='figures/poster/fig3_txv.png')
    args = parser.parse_args()
    plot_fig3_txv.plot_fig3(resfolder=args.resfolder, outpath=args.outpath)
    print('Done.')
