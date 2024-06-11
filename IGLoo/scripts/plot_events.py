import pandas as pd
import matplotlib.pyplot as plt
import argparse



def plot_vdj(fn_csv, sample_id, output_name):
    df_sample = pd.read_csv(fn_csv, sep='\t', header=None)
    plt.title(sample_id + " VDJ combinations")
    plt.pie(x=df_sample[1], labels=df_sample[0].values, autopct=lambda p: '{:.0f}'.format((p/100)*df_sample[1].sum()))
    plt.savefig(output_name, bbox_inches='tight')



def main(arguments=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-rpt', '--event_split_report', help='the report listing number of reads supporting each event.', required=True)
    parser.add_argument('-id',  '--sample_ID', help='sample ID', default="sample", required=False)
    parser.add_argument('-out', '--output_name', help='output file name')
    args = parser.parse_args(arguments)
    
    event_report = args.event_split_report
    sample_id = args.sample_ID
    output_name = args.output_name
    if output_name == None:
        output_name = event_report + '.pdf'

    plot_vdj(event_report, sample_id, output_name)



if __name__ == "__main__":
    main()
