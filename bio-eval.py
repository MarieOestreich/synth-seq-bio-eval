from utils import *
import pandas as pd
import argparse
import numpy as np
from tabulate import tabulate
from matplotlib import pyplot as plt
import os

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--real_data",
        "-rd",
        type=str,
        required=True,
        help="path to real data .csv",
    )
    parser.add_argument(
        "--synth_data",
        "-sd",
        type=str,
        required=True,
        help="path to synthetic data .csv",
    )
    parser.add_argument(
        "--description",
        "-d",
        type=str,
        default="my-bio-eval",
        required=True,
        help="keyword to use as file prefix",
    )
    parser.add_argument(
        "--verbose",
        "-v",
        type=bool,
        default=True,
        help="whether or not to create an HTML report with the results",
    )
        
    args = parser.parse_args()
    return args

def main():
    args = parse_arguments()

    if not os.path.isdir('./plots/'):
        os.mkdir('./plots/')
    tables = []
    # read in real data:
    real_data = pd.read_csv(args.real_data)
    # read in synthetic data:
    fake_data = pd.read_csv(args.synth_data)
    # create Bioeval object:
    bioeval = BioEval(real_data, fake_data)
    real, fake = bioeval.summary()
    stats = bioeval.count_stats()
    try:
        bioeval.DE_test()
        DE = ''
    except:
        DE = 'DE-test was not run, since at least 2 labels are required'

    coex_out = bioeval.coex_test()
    fig = bioeval._plot_coex_test()
    plt.tight_layout()
    plt.savefig('./plots/'+args.description+'_coexpression.jpg')
    if args.verbose:
        tables.append(real)
        tables.append(fake)
        tables.append(round(stats, 2))
        
        # print('label counts real data:')
        # print(tabulate(real, headers='keys', tablefmt='psql'))
        # print('')
        # print('label counts fake data:')
        # print(tabulate(fake, headers='keys', tablefmt='psql'))
        # print('')
        # print('general count statistics:')
        # print(tabulate(stats, headers='keys', tablefmt='psql'))

        # print('')
        # print('')
        # print('DE_test')
        if DE == '':
            tables.append(bioeval.stats1)
            tables.append(round(bioeval.stats2, 2))
            
            # print(tabulate(bioeval.stats1, headers='keys', tablefmt='psql'))
            # print(tabulate(bioeval.stats2, headers='keys', tablefmt='psql'))
            fig = bioeval.plot_false_DE(which = 'up')
            plt.savefig('./plots/'+args.description+'_bp_false_up.jpg')
            fig = bioeval.plot_false_DE(which = 'down')
            plt.savefig('./plots/'+args.description+'_bp_false_down.jpg')
        else:
            tables.append(pd.DataFrame())
            tables.append(pd.DataFrame())
            print(DE)

        # print('')
        # print('')
        # print('co-expression test')
        # print(tabulate(coex_out, headers='keys', tablefmt='psql'))
        tables.append(round(coex_out, 2))
        write_report(tables, args.description)



if __name__ == "__main__":
    main()