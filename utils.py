import numpy as np
import pandas as pd
import scipy.stats
import statsmodels.stats.multitest
import matplotlib.pyplot as plt
from IPython.display import display
import collections

__all__ = ['BioEval', 'LabelError', 'write_report']

class LabelError(Exception):
    """Exception raised when attemting to run DE-gene test with only 1 class.

    Attributes:
        message -- explanation of the error
    """

    def __init__(self, message="At least 2 different labels/classes are required to run a DE-gene analysis"):
        self.message = message
        super().__init__(self.message)

class BioEval:

    def __init__(self, original_data, synth_data):
        """Biology-Based Evaluation Class.
        
        Args:
            original_data (pd.DataFrame):
                The real data used for training. Genes are rows (indeces), samples are columns.
            synth_data (pd.DataFrame):
                The generated synthetic data. Genes are rows, samples are columns.
            original_label (pd.DataFrame):
                Sample labels of the original data.
            synth_label (pd.DataFrame):
                Sample labels of the synthetic data.
        """

        self.int2str = lambda x: str(x)
        self._original_data = original_data
        self._synth_data = synth_data
        self._original_label = pd.DataFrame([[x for x in original_data.index], [self.int2str(l) for l in self._original_data.iloc[:, 0].values.tolist()]],
                                            index=["sample", "label"]).T
        self._synth_label = pd.DataFrame([[x for x in synth_data.index], [self.int2str(l) for l in self._synth_data.iloc[:, 0].values.tolist()]],
                                            index=["sample", "label"]).T

        self._original_data = self._original_data.iloc[:, 1:].T
        self._synth_data = self._synth_data.iloc[:, 1:].T
        self._coex_res = None
        self._DE_res = None
        self.stats1 = None
        self.stats2 = None
        self.true_down = None
        self.true_up = None
        self.false_down = None
        self.false_up = None

        

    def summary(self):
        labels_real = collections.Counter(self._original_label.label.values)
        labels_real = dict(sorted(labels_real.items()))
        labels_fake = collections.Counter(self._synth_label.label.values)
        labels_fake = dict(sorted(labels_fake.items()))
        df_real = pd.DataFrame()
        df_fake = pd.DataFrame()
        for cls, count in labels_real.items():
            row = pd.DataFrame([[cls, count, int(count/sum(labels_real.values())*100)]], columns=['class', 'n samples', '%'+' of total'])
            df_real = pd.concat([df_real, row])

        for cls, count in labels_fake.items():
            row = pd.DataFrame([[cls, count, int(count/sum(labels_fake.values())*100)]], columns=['class', 'n samples', '%'+' of total'])
            df_fake = pd.concat([df_fake, row])
        
        
        return df_real, df_fake
    
    def count_stats(self, genes='all'):
        stats = pd.DataFrame()
        if genes=='all':
            x1 = self._original_data.values.flatten()
            x2 = self._synth_data.values.flatten()
        else:
            x1 = self._original_data.loc[genes, :].values.flatten()
            x2 = self._synth_data.loc[genes, :].values.flatten()

        stats = pd.concat([stats, pd.DataFrame([[x1.min(), x1.mean(), np.median(x1), x1.max()]], columns=['min', 'mean', 'median', 'max'], index=['real data'])])
        stats = pd.concat([stats, pd.DataFrame([[x2.min(), x2.mean(), np.median(x2), x2.max()]], columns=['min', 'mean', 'median', 'max'], index=['synthetic data'])])
        return stats

    def plot_false_DE(self, which = 'up'):
        if which == 'up':
            d1 = self.false_up
        else:
            d1 = self.false_down
        d = {}
        colors = []
        for k in d1.keys():
            x1 = self._original_data.loc[list(d1[k]), :].values.flatten()
            x2 = self._synth_data.loc[list(d1[k]), :].values.flatten()
            d[k + '_real'] = x1
            d[k + '_synth'] = x2
            if len(colors) % 4==0:
                colors += ['white', 'white']
            else:
                colors += ['grey', 'grey']

        fig, ax = plt.subplots()
        bplot = ax.boxplot(d.values(), sym='', showfliers=False, patch_artist=True)
        for patch, color in zip(bplot['boxes'], colors):
            patch.set_facecolor(color)

        ax.set_xticklabels(d.keys())
        plt.xticks(rotation = 90)
        plt.title('Falsely ' + which + '-regulated genes')
        plt.tight_layout()
        # plt.show()
        return fig

    # def plot_DE(self, which = 'up'):
    #     if which == 'up':
    #         d1 = self.false_up
    #         d2 = self.true_up
    #     else:
    #         d1 = self.false_down
    #         d2 = self.true_down
    #     d = {}
    #     for k in d1.keys():
    #         x1 = self._original_data.loc[list(d1[k]), :].values.flatten()
    #         x2 = self._synth_data.loc[list(d1[k]), :].values.flatten()
    #         d[k + '_delta_false_' + which] = x1-x2
    #         x1 = self._original_data.loc[list(d2[k]), :].values.flatten()
    #         x2 = self._synth_data.loc[list(d2[k]), :].values.flatten()
    #         d[k + '_delta_true_' + which] = x1-x2


    #     fig, ax = plt.subplots()
    #     ax.boxplot(d.values(), sym='', showfliers=False)
    #     ax.set_xticklabels(d.keys())
    #     plt.xticks(rotation = 90)
    #     plt.show()

    def DE_stats(self, which = 'up'):
        _stats = pd.DataFrame()
        if which == 'up':
            d1 = self.false_up
            d2 = self.true_up
        else:
            d1 = self.false_down
            d2 = self.true_down

        for k in d1.keys():
            genes1 = list(d1[k])
            genes2 = list(d2[k])
            stats1 = self.count_stats(genes=genes1)
            stats2 = self.count_stats(genes=genes2)
            stats1 = round(stats1.iloc[:2, :].diff(), 0).iloc[1:2,:]
            stats2 = round(stats2.iloc[:2, :].diff(), 0).iloc[1:2,:]
            stats1.index = [k + '_false_' + which]
            stats2.index = [k + '_true_' + which]
            _stats = pd.concat([_stats, stats1, stats2])
        return _stats


    def _pearson_corr(self, mat):

        """Compute piar-wise Pearson Correlation Coefficients for the rows of a given data matrix.
        
        Args:
            mat:
                A numeric matrix.
        """

        # calc Pearson correlation:
        corrmat = np.corrcoef(mat)
        # set diagonal and one triangle to zero:
        corrmat = np.tril(corrmat, -1)
        # make pandas data frame:
        corrmat = pd.DataFrame(data=corrmat, index=mat.index, columns=mat.index)
        corrmat.index.name = None
        # reshape to be a three column data frame (gene 1, gene 2, correlation value) isntead of correlation matrix:
        corrmat = corrmat.stack().reset_index()
        corrmat.columns = ['Gene1', 'Gene2', 'Correlation']
        # Remove rows with correlations <= 0. Other than those introduced by np.tril, this also removes true zero correlations. 
        # Since we only focus on positive correlations here, this doesn't matter.
        corrmat = corrmat[corrmat.Correlation > 0]
        # Concatenate the two genes of each row in alphabetical order to allow comparison between datasets:
        corrmat["Merged"] = [x + y if x < y else y + x for x, y in zip(corrmat["Gene1"], corrmat["Gene2"])]
        return (corrmat)

    def coex_test(self):

        """ Test the preservation of gene co-expressions between the original and the synthtic data.
        Returns a dataframe of the percentage of preserved edges and the RMSD of the correlation values for a series of correlation cutoffs.
        Also plots the results.        
        """

        corrdf1 = self._pearson_corr(mat=self._original_data)
        corrdf2 = self._pearson_corr(mat=self._synth_data)
        corrdf1.index = corrdf1["Merged"]
        corrdf2.index = corrdf2["Merged"]
        res = pd.DataFrame(columns=['cutoff', 'perc_preserved', 'rmsd', 'edges_real', 'edges_fake', 'real-fake-ratio'])
        for c in range(0, 100, 1):
            c = c / 100
            # filter for cutoff:
            corrdf1_filt = corrdf1[corrdf1.Correlation >= c]
            corrdf2_filt = corrdf2[corrdf2.Correlation >= c]
            # in case one (or both) dataset now has no edges left:
            if corrdf1_filt.shape[0] == 0 or corrdf2_filt.shape[0] == 0:
#                 return (0, float('inf'))
                res_row = pd.DataFrame([[c, 0, 'NA', corrdf1_filt.shape[0], corrdf2_filt.shape[0], 'None']], columns=['cutoff', 'perc_preserved', 'rmsd', 'edges_real', 'edges_fake', 'real-fake-ratio'])
                # res = res.append(res_row)
                res = pd.concat([res, res_row])
                continue
            # merge based on edges and remove entries that contain NAs (i.e., edges that were not preserved):
            common = pd.concat([corrdf1_filt, corrdf2_filt], axis=1).dropna()
            # of preserved edges only keep weigth (correlation vlaue) from each data set:
            out = common.iloc[0:common.shape[0], [2, 6]]
            out.columns = ["Correlation1", "Correlation2"]
            # set the total edge number to be that of the bigger data frame:
            if corrdf1_filt.shape[0] > corrdf2_filt.shape[0]:
                edge_num = corrdf1_filt.shape[0]
            else:
                edge_num = corrdf2_filt.shape[0]

            # calculate percentage of preserved edges on the number of all edges:
            perc_preserved_edges = 100 * out.shape[0] / edge_num
            # calculate root mean squared deviation for the correlation values of the preserved edges (the lower the better):
            rmsd = np.sqrt(np.mean((out["Correlation1"] - out["Correlation2"]) ** 2))
            res_row = pd.DataFrame([[c, perc_preserved_edges, rmsd, corrdf1_filt.shape[0], corrdf2_filt.shape[0], round(corrdf1_filt.shape[0]/corrdf2_filt.shape[0], 2)]], columns=['cutoff', 'perc_preserved', 'rmsd', 'edges_real', 'edges_fake', 'real-fake-ratio'])
            # res = res.append(res_row)
            res = pd.concat([res, res_row])

        # res = res.dropna()
        self._coex_res = res
        # self._plot_coex_test()
        return res

    def _plot_coex_test(self):
        """
            Plot results of coex_test().
        """

        if self._coex_res is None:
            print("run method coex_test() first")
        else:
            res = self._coex_res[self._coex_res.rmsd != 'NA']
            fig, (ax1, ax2) = plt.subplots(2, sharex=True)

            ax1.plot(res['cutoff'], res['perc_preserved'], label='perc_preserved', color='orange')
            ax1.set(ylabel='preserved edges [%]')
            ax1.legend()
            ax2.plot(res['cutoff'], res['rmsd'], label='rmsd')
            ax2.set(xlabel='Pearson Correlation Cutoff', ylabel='RMSD')
            ax2.legend()
            return fig

    def _get_label_comparisons(self, sample_to_label):

        """
            Get array of all required group comparisons.
        
        Args:
            sample_to_label (pd.DataFrame):
                A two column data frame with on column giving the sample ID and the other the label.
        """

        comparisons = []
        if len(set(sample_to_label["label"])) == 1:
            raise LabelError()

        for i in range(len(set(sample_to_label["label"]))):
            condition1 = list(set(sample_to_label["label"]))[i].replace(" ", "_")
            for j in range(i + 1, len(set(sample_to_label["label"]))):
                condition2 = list(set(sample_to_label["label"]))[j].replace(" ", "_")
                if condition1 < condition2:
                    comparisons.append(str(condition1) + "_VS_" + str(condition2))
                else:
                    comparisons.append(str(condition2) + "_VS_" + str(condition1))
        return (comparisons)

    def _wt_pair(self, m1, m2, p_threshold):

        """Calculate Wilcoxon test for two vectors."""

        out = pd.DataFrame(columns=["gene", "lesser", "greater"])
        for i in range(m1.shape[1]):
            v1 = m1.iloc[:, i].to_numpy()  # add assertion that genes (columns) are in the same order in both data frames
            v2 = m2.iloc[:, i].to_numpy()

            if sum(v1) != 0 and sum(v2) != 0:
                lesserw, lesserp = scipy.stats.ranksums(v1, v2, alternative="less")
                greaterw, greaterp = scipy.stats.ranksums(v1, v2, alternative="greater")
            else:
                lesserp = None
                greaterp = None
            res = pd.DataFrame([[m1.columns[i], lesserp, greaterp]], columns=["gene", "lesser", "greater"])
            # out = out.append(res)
            out = pd.concat([out, res])
        out = out.dropna(axis=0)

        # produces only NaN ?
        # _, out["lesser"], _, _ = statsmodels.stats.multitest.multipletests(pvals = out["lesser"], method = "fdr_bh")
        # _, out["greater"], _, _ = statsmodels.stats.multitest.multipletests(pvals = out["greater"], method = "fdr_bh")

        DE = pd.DataFrame(columns=['gene', 'DE'])
        for index, row in out.iterrows():
            if (row['lesser'] <= p_threshold):
                de_row = pd.DataFrame([[row['gene'], -1]], columns=['gene', 'DE'])
            elif (row['greater'] <= p_threshold):
                de_row = pd.DataFrame([[row['gene'], 1]], columns=['gene', 'DE'])
            else:
                de_row = pd.DataFrame([[row['gene'], 0]], columns=['gene', 'DE'])
            # DE = DE.append(de_row)
            DE = pd.concat([DE, de_row])

        return DE

    def _wt(self, counts, sample_to_label, comparisons, p_threshold=0.01):

        """Calculate Wilcoxon test for all label comparisons."""
        out = dict()

        for comp in comparisons:
            cond1 = comp.split("_VS_")[0]
            cond2 = comp.split("_VS_")[1]

            samples1 = sample_to_label[sample_to_label["label"] == cond1].loc[:, "sample"].to_numpy()
            samples2 = sample_to_label[sample_to_label["label"] == cond2].loc[:, "sample"].to_numpy()

            counts1 = counts.loc[:, samples1]
            counts2 = counts.loc[counts1.index.to_numpy(), samples2]

            res = self._wt_pair(m1=counts1.T, m2=counts2.T, p_threshold=p_threshold)
            up = set(res[res['DE'] == 1].loc[:, 'gene'])
            down = set(res[res['DE'] == -1].loc[:, 'gene'])
            negative = set(res[res['DE'] == 0].loc[:, 'gene'])
            out[comp] = {"up": up,
                         "down": down,
                         "negative": negative}

        return out

    def _DE_preservation(self, original, synthetic):

        """
        Metric for how well the differentially expressed genes are preserved between the real and the synthetic data.
        A values from [0,1], 1 being the optimum.
        """

        up_original = original['up']
        RU = len(up_original) # real up
        down_original = original['down']
        RD = len(down_original) # real down
        up_synth = synthetic['up']
        SU = len(up_synth) # synth up
        down_synth = synthetic['down']
        SD = len(down_synth) # synth down

        # true up
        TU = up_synth.intersection(up_original)
        # true down
        TD = down_synth.intersection(down_original)
        # false up
        FU = up_synth.difference(up_original)
        # false down
        FD = down_synth.difference(down_original)
        
        TN = len(original['negative'].intersection(synthetic['negative']))
        # false negatives
        FN = len(synthetic['negative'].difference(original['negative']))
        # true positives (note: positive = DE)
        TP = len(TU) + len(TD)
        # false positives
        FP = len(FU) + len(FD)
        
        # True-False ratio up
        PCU = round(len(TU)/len(FU), 2)
        # True-False ratio down
        PCD = round(len(TD)/len(FD), 2)
        # Precision
        PRCSN = TP/(TP+FP)
        # Recall
        RCLL = TP/(TP+FN)
        # F1
        F1 = 2*(PRCSN*RCLL)/(PRCSN+RCLL)
        
        stats = pd.DataFrame([[len(TU), len(TD)], [len(FU), len(FD)], [SU, SD], [RU, RD],[PCU, PCD]],
                             index = ['TRUE', 'FALSE', 'FAKE TOT.', 'REAL TOT.', 'T-F-RATIO'],
                             columns = ['UP', 'DOWN']).T
        stats2 = pd.DataFrame([PRCSN, RCLL, F1], index=['Precision', 'Recall', 'F1']).T

        return F1, stats, stats2, TU, TD, FU, FD

    def DE_test(self, p_threshold=0.05):

        """Test preservation of differentially expressed genes in synthetic dataset.
        
        Args:
            p_threshold (float):
                Maximum p-value for a gene to be considered DE. Default is 0.05.
        """

        comparisons_original = self._get_label_comparisons(sample_to_label=self._original_label)
        comparisons_synth = self._get_label_comparisons(sample_to_label=self._synth_label)
        if (len(set(comparisons_original).intersection(set(comparisons_synth))) != len(set(comparisons_original))):
            print(
                "The synthetic data does not have the same set of labels as the original data. Will proceed to only compare DE genes between common labels.")

        comparisons = set(comparisons_original).intersection(set(comparisons_synth))
        wt_original = self._wt(counts=self._original_data,
                                sample_to_label=self._original_label,
                                comparisons=comparisons,
                                p_threshold=p_threshold)
        wt_synth = self._wt(counts=self._synth_data,
                             sample_to_label=self._synth_label,
                             comparisons=comparisons,
                             p_threshold=p_threshold)

        score = []
        stats1 = pd.DataFrame()
        stats2 = pd.DataFrame()
        true_up = {}
        false_up = {}
        true_down = {}
        false_down = {}
        for comp in comparisons:
            _score, _stats1, _stats2, _true_up, _true_down, _false_up, _false_down = self._DE_preservation(original=wt_original[comp], synthetic=wt_synth[comp])
            true_up[comp] = _true_up
            false_up[comp] = _false_up
            true_down[comp] = _true_down
            false_down[comp] = _false_down
            
            _stats1.index = [comp + "_"+ x for x in list(_stats1.index)]
            stats1 = pd.concat([stats1, _stats1])
            _stats2.index = [comp]
            stats2 = pd.concat([stats2, _stats2])
            score.append(_score)
        col_means = pd.DataFrame(stats2.mean(axis=0), columns=['mean']).T
        stats2 = pd.concat([stats2, col_means])
        self.stats1 = stats1
        self.stats2 = stats2
        # display(stats1)
        # display(stats2)
        self.true_down = true_down
        self.true_up = true_up
        self.false_down = false_down
        self.false_up = false_up
        self._DE_res = sum(score) / len(score)        
        return None



def write_report(tables, page_title='title'):
    # 1. Set up multiple variables to store the titles, text within the report
    page_title_text=page_title

    title = 'Bio-eval report'
    general_info = 'General Information'
    label_freqs = 'Label Frequencies'
    real_data = 'Real Data'
    synth_data = 'Synthetic Data'
    gen_count_stats = 'General Count Statistics'

    DE = 'Differential Expression'
    description_false_DE = 'Boxplot of expression values of genes that have been identified as up-/down-regulated in the synthetic data, but that are not up-/down-regulated in the real data.'
    fasle_de = 'False DE' # plots
    plot1 = '.\\plots\\' + page_title + '_bp_false_down.jpg'
    plot2 = '.\\plots\\' + page_title + '_bp_false_up.jpg'
    plot3 = '.\\plots\\' + page_title + '_coexpression.jpg'

    coexpr = 'Co-Expression' # table + plot

    # 2. Combine them together using a long f-string
    html = f'''
        <html>
            <head>
                <title>{page_title_text}</title>
            </head>
            <body>
                <h1>{title}</h1>

                <h2>{general_info}</h2>
                <h3>{label_freqs}</h3>
                <p>{real_data}</p>
                {tables[0].to_html()}
                <p>{synth_data}</p>
                {tables[1].to_html()}
                <p>{''}</p>
                <h2>{gen_count_stats}</h2>
                {tables[2].to_html()}
                <p>{''}</p>

                <h2>{DE}</h2>
                {tables[3].to_html()}
                {tables[4].to_html()}
                <p>{''}</p>
                <h2>{fasle_de}</h2>
                <p>{description_false_DE}</p>
                <img src={plot1}>
                <p>{''}</p>
                <img src={plot2}>
                <p>{''}</p>
                <h2>{coexpr}</h2>
                {tables[5].to_html()}
                <img src={plot3}>
            </body>
        </html>
        '''
    # 3. Write the html string as an HTML file
    with open('report.html', 'w') as f:
        f.write(html)